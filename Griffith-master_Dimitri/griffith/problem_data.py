import sys
import os
from math import exp
import plotrc as plot_options
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull
import shapely.geometry as sg 
from Func import Unit_vect
# from griffith.mesh import Mesh, Boundary_Mesh, Triangle

from .geometry import Point, Circle, rotation, dist
sys.path.append(os.path.abspath('..'))


#####################
# TIME DISCRETIZATION
#####################
class Regular_Time_Discretization:
    def __init__(self, time_step, time_min=None):
        self._time_step = time_step
        self._time = time_min if time_min else time_step
        self._raise = False

    def __iter__(self):
        return self

    def __next__(self):
        time = self._time
        self._time += self._time_step

        if self._raise:
            raise StopIteration

        if time > 1:
            time = 1
            self._raise = True
        return time


#########################
# FRACTURE DISCRETIZATION
#########################
class Fracture_Discretization:
    def __init__(self, angular_step, lengh_step, boundary_step=None, boundary_point=None, interior_step=None, interior_fast_init=True, interior_fast_step=True, interior_start_angle=None):
        self.angular_step = angular_step
        self.lengh_step = lengh_step
        self.min_lengh = lengh_step  # FIXME
        self.max_lengh = np.inf #2*lengh_step  # np.inf  # FIXME

        self.interior_step = interior_step
        self.interior_fast_step = interior_fast_step
        self.interior_fast_init = interior_fast_init
        if interior_start_angle:
            self.interior_start_angle = interior_start_angle
        else:
            self.interior_start_angle = 1.57079632679  # pi/2
        self.boundary_step = boundary_step
        if boundary_point:
            print((boundary_point))
            self.boundary_point = Point(*boundary_point)
            assert not boundary_step
        else:
            self.boundary_point = None


#######################
# BOUNDARY DISPLACEMENT
#######################
class Boundary_Displacement:
    """
    Abstract class for boundary displacement.
    Personalized boundary displacement should inherit from it, and implement __call__ function.
    """

    def __call__(self, p):
        raise NotImplementedError


class Linear_Displacement(Boundary_Displacement):
    def __init__(self, traction_coefficient=1):
        self._traction_coefficient = traction_coefficient

    def __call__(self, time, p):
        return self._traction_coefficient*time*self._func(p)


class Constant_Displacement_On_Y(Linear_Displacement):
    r"""
    | | -> |     |
    """

    def __init__(self, traction_coefficient=1, abscissa_mid=0.9):
        super().__init__(traction_coefficient)
        self._abscissa_mid = abscissa_mid

    def _func(self, p):
        if p.x > self._abscissa_mid:  # and p.y > 40 and p.y < 50:
            return np.array([1, 0.])
        # if p.x < 0.1 and p.y > 90 and p.y < 95: return np.array([-1, 0])
        return np.array([0., 0.])
        # if p.x < 0.1 and p.y > 40 and p.y < 60:
        # return np.array([1, 0])
        # return np.array([-1, 0])  # a supprimer si Neumann


class Boundary_Displacement_by_percussion(Linear_Displacement):
    """
    The boundary displacement is computed thanks to a mass-spring deformation data
    Boundary data contains -(xi,yi) the nodes of the mass-spring system and the deformation field u1(x,y),u2(x,y)
    We compute the P1 approximation on the boundary by projecting the data on the boundary
    """

    def __init__(self, traction_coefficient=1., boundary_data=None):
        #self.some_threshold = some_threshold
        super().__init__(traction_coefficient=1)
        self.boundary_data = boundary_data

    def check_Dirichlet(self, p):
        # Implement the logic to check if point p is inside the specific region
        # Return True if it is inside; otherwise, return False
        # Example: check if the x-coordinate of p is greater than some threshold
        return dist(p, Contact_region=np.array([98, 23])) < 5

    def _func(self, p):
        # Implement the logic to calculate the displacement vector at point p
        # based on a specific function for this region
        # Return the displacement vector as a numpy array
        # Example: return the displacement as [0, f(p.y)], where f is some function
        data = self.boundary_data

        # extract the data of interpolation points and data

        # Function to calculate the angle of a point with respect to the centroid
        def angle_with_centroid(point, centroid):
            x, y = point
            cx, cy = centroid
            return np.arctan2(y - cy, x - cx)

        # Function to rearrange points to form a proper polygon
        def rearrange_polygon_points(xdata, ydata):
            # Calculate the centroid
            centroid = (np.mean(xdata), np.mean(ydata))
            # Create a list of points and calculate angles with respect to the centroid
            points = list(zip(xdata, ydata))
            angles = [angle_with_centroid(point, centroid) for point in points]
            # Sort the points based on the angles
            sorted_points = [point for _, point in sorted(zip(angles, points))]
            return zip(*sorted_points)

        def rearrange_polygon_data(xdata, ydata, z1data, z2data):
            x_reordered, y_reordered = rearrange_polygon_points(xdata, ydata)
            # Rearrange z1data based on the new order of x and y
            z1_reordered = np.array([z1data[np.where((xdata == x) & (ydata == y))[0][0]] for x, y in zip(x_reordered, y_reordered)])
            # Rearrange z2data based on the new order of x and y
            z2_reordered = np.array([z2data[np.where((xdata == x) & (ydata == y))[0][0]] for x, y in zip(x_reordered, y_reordered)])

            return x_reordered, y_reordered, z1_reordered, z2_reordered

        xdata, ydata = data[:, 0], data[:, 1]

        x_reordered, y_reordered, z1_reordered, z2_reordered = rearrange_polygon_data(xdata, ydata, data[:,2], data[:,3])

        xdata = np.array(x_reordered)
        ydata = np.array(y_reordered)
        z1data = np.array(z1_reordered)
        z2data = np.array(z2_reordered)

        polygon = sg.Polygon(np.column_stack((xdata, ydata)))

        interp_function_x = LinearNDInterpolator(list(zip(xdata, ydata)), z1data)
        interp_function_y = LinearNDInterpolator(list(zip(xdata, ydata)), z2data)
        
        def project_on_polygon(point, polygon):
            if polygon.contains(sg.Point(point)):
                return point
            nearest_point = polygon.exterior.interpolate(polygon.exterior.project(sg.Point(point)))
            return np.array([nearest_point.x, nearest_point.y])

        def f_x(x_eval, y_eval):
            interpolated_value = interp_function_x(x_eval, y_eval)
            nan_mask = np.isnan(interpolated_value)
            center = np.array([np.mean(xdata), np.mean(ydata)])
            if np.any(nan_mask):
                projected_point = project_on_polygon(sg.Point(x_eval,y_eval), polygon)
                projected_point += Unit_vect(projected_point, center) /10000.
                interpolated_value = interp_function_x(projected_point[0], projected_point[1])

            return interpolated_value
            # return np.interp(x_eval, xdata, z1data) + np.interp(y_eval, ydata, z1data)

        def f_y(x_eval,y_eval):
            interpolated_value = interp_function_y(x_eval, y_eval)
            nan_mask = np.isnan(interpolated_value)
            center = np.array([np.mean(xdata), np.mean(ydata)])
            if np.any(nan_mask):
                projected_point = project_on_polygon(sg.Point(x_eval,y_eval), polygon)
                projected_point += Unit_vect(projected_point, center) /10000.
                interpolated_value = interp_function_y(projected_point[0], projected_point[1])
            return interpolated_value
            # return np.interp(x_eval, xdata, z2data) + np.interp(y_eval, ydata, z2data)

        def Dirichlet_function(x,y):
            return np.array([f_x(x,y), f_y(x,y)])
        
        # if dist(p, Point(98, 23)) < 5:
        #     return Dirichlet_function(p.x, p.y)
        return Dirichlet_function(p.x, p.y)


class Linear_Displacement_On_Y(Linear_Displacement):
    r"""
    | | -> / \
    """

    def __init__(self, traction_coefficient=1, abscissa_mid=0.5, y_min=0, c_min=1, y_max=100, c_max=0):
        super().__init__(traction_coefficient)
        self._abscissa_mid = abscissa_mid
        self.c_min, self.c_max = c_min, c_max
        self.y_min, self.y_max = y_min, y_max

    def _func(self, p):
        c_y = self.c_min + (self.c_max - self.c_min) / \
            (self.y_max - self.y_min)*p.y
        if p.x > self._abscissa_mid:
            return np.array([c_y, 0])
        return np.array([-c_y, 0])


class Picewise_Linear_Displacement_On_Y(Linear_Displacement):
    r"""
    | | -> ||
    | | -> / \
    """

    def __init__(self, traction_coefficient=1, abscissa_mid=0.5, y_min=0, y_max=100):
        super().__init__(traction_coefficient)
        self._abscissa_mid = abscissa_mid
        self._y_min = y_min
        self._y_max = y_max
        self._amplitude = y_max - y_min

    def _func(self, p):
        if p.y > self._y_max:
            return np.array([0, 0])
        if p.x > self._abscissa_mid:
            return np.array([(self._y_max - p.y)/self._amplitude, 0])
        return np.array([-(self._y_max - p.y)/self._amplitude, 0])


class Rotation_Displacement_On_Y(Boundary_Displacement):
    r"""
    #| | -> _ _ (angle = pi/2) 
    #| | -> / \ (angle = pi/4)
    This displacement is not linear.
    """

    def __init__(self, angle=1, abscissa_mid=0.5, point_left=(0, 100), point_right=(100, 100)):
        self._angle = angle
        self._abscissa_mid = abscissa_mid
        self._point_left = Point(*point_left)
        self._point_right = Point(*point_right)

    def __call__(self, time, p):
        if p.x > self._abscissa_mid:
            return (rotation(self._point_right, p, self._angle*time) - p).array
        return (rotation(self._point_left, p, -self._angle*time) - p).array


##############
# STIFF TENSOR
##############
class Identity_Tensor:
    def __init__(self):
        pass

    def tproduct(self, matrix):
        return matrix


class Lame_Tensor:
    def __init__(self, lambda_, mu):
        self._lambda = lambda_
        self._mu = mu

    def tproduct(self, matrix):
        return 2*self._mu*matrix + self._lambda*matrix.trace()*np.array(((1, 0), (0, 1)))

    @classmethod
    def _init_with_Young_Poisson(cls, E, nu):
        return cls(E*nu/((1+nu)*(1-2*nu)), 0.5*E/(1+nu))


# Classical ice tensor
lame_tensor_ice = Lame_Tensor._init_with_Young_Poisson(8.95, 0.295)


#################
# TOUGHNESS FIELD
#################
class Constant_Toughness:
    def __init__(self, k):
        self.k = k

    def __call__(self, p):
        return self.k

    def plot(self, figax):
        pass


class Japan_Toughness:
    def __init__(self, k1, k2, center, radius):
        self.k1 = k1
        self.k2 = k2
        self.circle = Circle(Point(*center), radius)

    def __call__(self, p):
        if self.circle.has_point(p):
            return self.k2
        return self.k1

    def plot(self, figax=None):
        self.circle.plot(figax, **plot_options.circle_inclusion)
        return figax


class Smooth_Japan_Toughness(Japan_Toughness):
    def __init__(self, k1, k2, center, radius):
        super().__init__(k1, k2, center, radius)
        self.sigma = 0.25

    def __call__(self, p):
        d = dist(p, self.circle.center)
        if d < self.circle.radius:
            # /(self.sigma*sqrt(2*pi))
            return self.k1 + (self.k2 - self.k1)*exp(-(self.circle.radius - d)**2/2*self.sigma**2)
        return self.k1


##############
# PROBLEM DATA
##############
class Discretization_Data:
    def __init__(self, mesh, time_discretization, fracture_discretization, tip_enrichement=False):
        self.mesh = mesh
        self.time_discretization = time_discretization
        self.fracture_discretization = fracture_discretization
        self.tip_enrichement = tip_enrichement

        if fracture_discretization.boundary_point:
            assert mesh.has_point_on_boundary(
                self.fracture_discretization.boundary_point)


class Physical_Data:
    def __init__(self, stiffness_tensor, toughness_field, boundary_displacement, initial_fracture=None):
        self.stiffness_tensor = stiffness_tensor
        self.toughness_field = toughness_field
        self.boundary_displacement = boundary_displacement
        self.initial_fracture = initial_fracture
