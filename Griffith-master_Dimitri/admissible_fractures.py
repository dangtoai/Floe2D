import os, sys, traceback
import argparse
import matplotlib.pyplot as plt
from math import pi

from griffith import fracture_iterators, log
from griffith.mesh import Fracture, Mesh, Broken_Mesh, NotAdmissibleFracture, ShouldBeAdmissibleFracture, Broken_Mesh_Linear_Tip
from griffith.geometry import Point, Segment
from griffith import problem_data

def run(raw_args=None):
  parser = argparse.ArgumentParser(description='Plots admissible segments, mainly to debug Fractured_Mesh class')
  subparsers = parser.add_subparsers(dest="mode")
  
  parser_BPI = subparsers.add_parser('BPI', help='For Boundary Point Iterator')
  parser_BPI.add_argument('-bs', '--boundary_step', type=float, help="Boundary step for fracture discretisation")
  
  parser_IPI = subparsers.add_parser('IPI', help='For Interior Point Iterator')
  parser_IPI.add_argument('-is', '--interior_step', type=float, help="Interior step for fracture discretisation")
  
  parser_BSI = subparsers.add_parser('BSI', help='For Boundary Segment Iterator')
  group_BSI = parser_BSI.add_mutually_exclusive_group()
  group_BSI.add_argument('-bs', '--boundary_step', type=float, help="Boundary step for fracture discretisation")
  group_BSI.add_argument('-bp', '--boundary_point', type=float, nargs=2, help="Boundary point for fracture initiation")
  parser_BSI.add_argument('-as', '--angular_step', type=float, default=pi/4, help="Angle step for fracture discretisation")
  parser_BSI.add_argument('-ls', '--lengh_step', type=float, default=10, help="Lengh step for fracture discretisation")
  
  parser_ISI = subparsers.add_parser('ISI', help='For Interior Segment Iterator')
  parser_ISI.add_argument('-is', '--interior_step', type=float, help="Interior step for fracture discretisation")
  parser_ISI.add_argument('-as', '--angular_step', type=float, default=pi/4, help="Angle step for fracture discretisation")
  parser_ISI.add_argument('-ls', '--lengh_step', type=float, default=10, help="Lengh step for fracture discretisation")
  parser_ISI.add_argument('-isa', '--interior_start_angle', type=float, default=0, help="Interior start angle for fracture discretisation")
  parser_ISI.add_argument('-ifi', '--interior-fast-init', action='store_true', help="Interior fast initialisation")
  
  parser_EFI = subparsers.add_parser('EFI', help='For Existing Fracture Iterator')
  parser_EFI.add_argument('-p', '--points', type=float, nargs='*', help="Fracture points for fracture initiation")
  parser_EFI.add_argument('-as', '--angular_step', type=float, default=pi/4, help="Angle step for fracture discretisation")
  parser_EFI.add_argument('-ls', '--lengh_step', type=float, default=10, help="Lengh step for fracture discretisation")
  
  parser_PI = subparsers.add_parser('PI', help='For Polyline Iterator')
  group_PI = parser_PI.add_mutually_exclusive_group()
  group_PI.add_argument('-bs', '--boundary_step', type=float, help="Boundary step for fracture discretisation")
  group_PI.add_argument('-bp', '--boundary_point', type=float, nargs=2, help="Boundary point for fracture initiation")
  parser_PI.add_argument('-ns', '--number_segments', type=int, default=2, help="Number of segments for the polyline")
  parser_PI.add_argument('-as', '--angular_step', type=float, default=pi/4, help="Angle step for fracture discretisation")
  parser_PI.add_argument('-ls', '--lengh_step', type=float, default=10, help="Lengh step for fracture discretisation")
  
  parser.add_argument('-m', '--mesh_file', default='mesh/square.msh', help="Name of the mesh file")
  parser.add_argument('-o', '--output_directory', default='admissible', help="Output directory")
  parser.add_argument('-p', '--plot', action='store_true', help="Plot result")
  
  args = parser.parse_args(raw_args)
  
  mesh = Mesh(args.mesh_file)
  
  if not os.path.isdir(args.output_directory):
    os.mkdir(args.output_directory)
  else:
    for filename in os.listdir(args.output_directory):
      os.remove(args.output_directory + '/' + filename)
  
  
  #########
  # Logging
  #########
  logger = log.Log(args.output_directory + '/' + 'admissible_fracture.log', level=log.INFO)
  log_queue = logger._log_queue
  
  
  #####################
  # Discretization Data
  #####################
  if args.mode == 'EFI':
    if len(args.points) < 4 or len(args.points) % 2 == 1:
      parser.error('Invalid data for points')
  
  for arg in ['boundary_step', 'boundary_point', 'interior_step', 'angular_step', 'lengh_step', 'interior_start_angle', 'interior_fast_init']:
    if not hasattr(args, arg):
      setattr(args, arg, None)
  fracture_discretization =  problem_data.Fracture_Discretization(args.angular_step, args.lengh_step, boundary_step=args.boundary_step, boundary_point=args.boundary_point, 
                                                                  interior_step=args.interior_step, interior_start_angle=args.interior_start_angle, interior_fast_init=args.interior_fast_init)
  
  
  ####################
  # Behaviour function
  ####################
  def plot_fracture(ind, fracture):
    strind = tuple_to_string(ind)
    try:
      frac_mesh = Broken_Mesh_Linear_tip(fracture, mesh)
    except NotAdmissibleFracture as e:
      log_queue.put(('WARNING', '{} : {} #NOT ADMISSIBLE\n'.format(strind, fracture)))
      fig, ax = mesh.plot()
      fracture.plot((fig, ax))
    except NotImplementedError:
      log_queue.put(('WARNING', '{} : {} #NOT IMPLEMENTED\n'.format(strind, fracture)))
    #except Exception as e:
    #  log_queue.put(('ERROR', '{} : {} #ERROR\n'.format(strind, fracture)))
    #  fig, ax = mesh.plot()
    #  fracture.plot((fig, ax))
    else:
      log_queue.put(('INFO', '{} : {}\n'.format(strind, fracture)))
      fig, ax = frac_mesh.plot()
    
    fig.savefig(args.output_directory + '/{}.svg'.format(strind))
    plt.close(fig)
  
  def compute_fracture(ind, fracture):
    strind = tuple_to_string(ind)
    #import pudb
    #pudb.set_trace()
    try:
      frac_mesh = Broken_Mesh_Linear_Tip(fracture, mesh)
    except NotAdmissibleFracture as e:
      log_queue.put(('DEBUG', '{} : {} #NOT ADMISSIBLE\n'.format(strind, fracture)))
    except ShouldBeAdmissibleFracture as e:
      log_queue.put(('WARNING', '{} : {} #SHOULD BE ADMISSIBLE\n'.format(strind, fracture)))
    #except Exception as e:
      #log_queue.put(('ERROR', '{} : {} #ERROR\n'.format(strind, fracture)))
    except NotImplementedError:
      log_queue.put(('WARNING', '{} : {} #NOT IMPLEMENTED\n'.format(strind, fracture)))
    else:
      log_queue.put(('INFO', '{} : {}\n'.format(strind, fracture)))
  
  def tuple_to_string(tple):
    result = "{}".format(tple[0])
    for i in tple[1:]:
      result += "-{}".format(i)
    return result
  
  
  ###################
  # Fracture Iterator
  ###################
  if args.plot:
    do_fracture = plot_fracture
  else:
    do_fracture = compute_fracture
    
  if args.mode == 'BPI':
    fig, ax = mesh.plot()
    for pointsegment in fracture_iterators.Admissible_Boundary_Point(fracture_discretization, mesh.boundary_mesh):
      point, segment = pointsegment
      if args.plot:
        point.plot((fig, ax), marker='+', color='red')
    if args.plot:
      fig.savefig(args.output_directory + '/boundary_points.svg')

  elif args.mode == 'IPI':
    fig, ax = mesh.plot()
    for point in fracture_iterators.Admissible_Interior_Point(fracture_discretization, mesh.boundary_mesh):
      point.plot((fig, ax), marker='+', color='red')
    if args.plot:
      fig.savefig(args.output_directory + '/interior_points.svg')
  
  elif args.mode == 'BSI':
    if args.boundary_point:
      for i, fracture in enumerate(fracture_iterators.Admissible_Fractures_From_Fixed_Boundary_Point(fracture_discretization, mesh, fracture_discretization.boundary_point)):
        do_fracture([i], fracture)
  
    else:
      assert args.boundary_step
      for i, fracture in enumerate(fracture_iterators.Admissible_Fractures_From_Boundary(fracture_discretization, mesh)):
        do_fracture([i], fracture)
  
  elif args.mode == 'ISI':
    #import pudb
    #pudb.set_trace()
    for i, fracture in enumerate(fracture_iterators.Admissible_Fractures_From_Interior(fracture_discretization, mesh)):
      do_fracture([i], fracture)
  
  elif args.mode == 'EFI':
    points = [Point(x, y) for x, y in zip(*[iter(args.points)]*2)]
    segments = [Segment(p1, p2) for p1, p2 in zip(points, points[1:])]
    old_fracture = Fracture(segments, mesh)
    for i, bigger_fracture in enumerate(fracture_iterators.Admissible_Fractures_From_Fracture(fracture_discretization, mesh, old_fracture)):
      do_fracture([i], bigger_fracture)
  
  elif args.mode == 'PI':
    def one_more_segment(iterator):
      for segment in iterator:
        ind, segment = segment
        for i, bigger_segment in enumerate(fracture_iterators.Admissible_Fractures_From_Fracture(fracture_discretization, mesh, segment)):
          try:
            new_ind = ind + (i,)
          except TypeError:
            new_ind = (ind, i)
          yield new_ind, bigger_segment
  
    def gen(nbr_segments):
      if args.boundary_step:
        pl_iterator = enumerate(fracture_iterators.Admissible_Fractures_From_Boundary(fracture_discretization, mesh))
      else:
        pl_iterator = enumerate(fracture_iterators.Admissible_Fractures_From_Fixed_Boundary_Point(fracture_discretization, mesh, fracture_discretization.boundary_point))
      for i in range(nbr_segments - 1):
        pl_iterator = one_more_segment(pl_iterator)
      return pl_iterator
  
    for fracture in gen(args.number_segments):
      ind, fracture = fracture
      do_fracture(ind, fracture)

  logger.exit()

if __name__ == '__main__':
  run()
