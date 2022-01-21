from typing import List, Tuple, Dict, Set
from copy import deepcopy
from random import choice

class UndirectedGraph:
	"""Define a generic undirected graph where vertices are numbered from 0"""
	
	def __init__(self, numVertices: int = 0, listEdges: List[Tuple[int, int]] = []):
		if numVertices < 0:
			raise ValueError("numVertices is negative")
		for v1, v2 in listEdges:
			if v1 not in range(numVertices) or v2 not in range(numVertices):
				raise ValueError(f"Edge {(v1, v2)} has invalid vertex")
			if v1 == v2:
				raise ValueError(f"{(v1, v2)} is a loop")

		self._graphDict: Dict[int, Set[int]] = { i:set() for i in range(numVertices) }	# Graph as dict of vertices
		self._graphSet: Set[Tuple[int, int]] = set()	# Graph as set of edges
		self._removedVertices: List[int] = []	# Removed vertices, for recycling numbers
		for v1, v2 in listEdges:
			self.AddEdge(v1, v2)

	@property
	def NumVertices(self):
		return len(self._graphDict)

	@property
	def NumEdges(self):
		return len(self._graphSet)
	
	@property
	def AsDict(self):
		"""A deep copy of the graph as dictionary with key as vertex, value as set of vertices as edges"""
		return deepcopy(self._graphDict)

	@property
	def AsSet(self):
		"""A deep copy of the graph as set of pairs of vertices as edges"""
		return deepcopy(self._graphSet)

	def IsEdge(self, v1: int, v2: int):
		"""Check if (v1,v2) is an edge"""
		return (v1, v2) in self._graphSet or (v2, v1) in self._graphSet

	def VertexDegree(self, v: int):
		"""Return number of neighbors of vertex v"""
		if v not in self._graphDict:
			raise ValueError(f"{v} is not a vertex in graph")
		return len(self._graphDict[v])

	def AddEdge(self, v1: int, v2: int):
		"""Add an edge from v1 to v2"""
		if self.IsEdge(v1, v2):
			raise ValueError(f"{(v1, v2)} is already an edge")
		if v1 not in self._graphDict or v2 not in self._graphDict:
			raise ValueError(f"One of the vertices not in graph")
		if v1 == v2:
			raise ValueError(f"Trying to add a loop")

		self._graphDict[v1].add(v2)
		self._graphDict[v2].add(v1)
		self._graphSet.add((min(v1, v2), max(v1, v2)))

	def AddVertex(self, associations: List[int] = []) -> int:
		"""Add a new vertex to graph, optionally listing all of its associations with existing vertices
		as new edges. Return new vertex's number"""
		for v in associations:
			if v not in self._graphDict:
				raise ValueError(f"{v} is not a vertex in graph")

		newVertex = self._removedVertices.pop(0) if len(self._removedVertices) else self.NumVertices
		self._graphDict[newVertex] = set()
		for v in associations:
			self.AddEdge(v, newVertex)
		return newVertex

	def RemoveEdge(self, v1: int, v2: int):
		"""Remove the edge from v1 to v2"""
		if not self.IsEdge(v1, v2):
			raise ValueError(f"{(v1, v2)} is not an edge in graph")
		self._graphDict[v1].remove(v2)
		self._graphDict[v2].remove(v1)
		self._graphSet.remove((min(v1, v2), max(v1, v2)))

	def RemoveVertex(self, v: int):
		"""Remove vertex v and all its associated edges from graph"""
		if v not in self._graphDict:
			raise ValueError(f"{v} is not a vertex of graph")
		for v2 in self.AsDict[v]:
			self.RemoveEdge(v, v2)
		if v < self.NumVertices - 1:
			self._removedVertices.append(v)
		del self._graphDict[v]

	def GetConnectedSubgraphs(self):
		"""Return list of sets of vertices in which each set is a connected subgraph (composante connexe en francais)"""
		
		def getIsland(v: int):
			"""Return the connected subgraph that contains vertex v as set of vertices"""
			queue = [v]
			island = {v}
			i = 0
			while i < len(queue):
				vertex = queue[i]
				for neighbor in self._graphDict[vertex]:
					if neighbor not in island:
						queue.append(neighbor)
						island.add(neighbor)
				i += 1
			return island
		
		subgraphs = list()
		vertices = set(self._graphDict.keys())
		while len(vertices) > 0:
			l = vertices.copy()
			island = getIsland(l.pop())
			vertices -= island
			subgraphs.append(island)
		return subgraphs

	def RouteInspection(self):
		"""Return a path as list of vertices that traverse all edges of graph, provided graph has total connectivity (Chinese Postman Problem)"""

		def addDummies(oddVertices: List[int]) -> UndirectedGraph:
			"""Create a copy of current graph then add dummy vertices to eliminate odd vertices. Return
			new graph"""
			graph = deepcopy(self)
			treated = set()
			for i, v1 in enumerate(oddVertices[:-1]):
				if v1 in treated:
					continue
				for v2 in oddVertices[i+1:]:
					if v2 not in treated and graph.IsEdge(v1, v2):
						graph.AddVertex([v1, v2])
						treated.add(v1)
						treated.add(v2)
			return graph

		def getTrail(graph: UndirectedGraph, v: int, trail: List[int] = None) -> List[int]:
			"""Return Eulerian trail of graph with dummy vertices if any, starting from vertex v"""
			if trail is None:
				trail = list()
			trail.append(v)
			if graph.NumEdges == 0:
				return trail
			else:
				queue = sorted(graph._graphDict[v], key = lambda x: graph.VertexDegree(x))
				for next in queue:
					graph.RemoveEdge(v, next)
					t = getTrail(graph, next, trail)
					if t:
						return t
					graph.AddEdge(v, next)
					trail.pop()
				return None

		def removeDummies(trail: List[int], dummies: Set[int]):
			"""Correct path containing dummies and also remove them from current graph"""
			trailCopy = trail.copy()
			for v in trailCopy:
				if v in dummies:
					trail.remove(v)

		def filterDuplicates(trail: List[int]):
			"""Remove unnecessary duplicated edges at the end of trails"""
			trailCopy = trail.copy()
			edgesTraveled = set()
			for i, v1 in enumerate(trailCopy[:-1]):
				if len(edgesTraveled) == len(self._graphSet):
					break
				v2 = trail[i+1]
				mini, maxi = min(v1, v2), max(v1, v2)
				if (mini, maxi) not in edgesTraveled:
					edgesTraveled.add((mini, maxi))
			else:
				return
			i += 1
			while i < len(trail):
				trail.pop(i)

		if len(self.GetConnectedSubgraphs()) > 1:
			raise Exception("Graph has multiple connected subgraphs")

		vertices = list(self._graphDict.keys())
		oddVertices = [v for v in vertices if self.VertexDegree(v) % 2]
		graph = self
		if len(oddVertices) == 0:
			# All vertices are even degree, start anywhere
			start = choice(vertices)
		elif len(oddVertices) == 2:
			# Two odd vertices, start at any of the two
			start = choice(oddVertices)
		else:
			# More than two odd vertices, add dummy vertices to eliminate all odd vertices
			start = choice(oddVertices)
			graph = addDummies(oddVertices)

		trail = getTrail(graph, start)
		if graph is not self:
			dummies = set(graph._graphDict.keys()) - set(self._graphDict.keys())
			removeDummies(trail, dummies)
		filterDuplicates(trail)
		return trail


if __name__ == '__main__':
	from itertools import combinations
	g = UndirectedGraph(4)
	for v1, v2 in combinations(range(4), 2):
		g.AddEdge(v1, v2)
	
	"""
	ABOVE IS A COMPLETE GRAPH WITH 4 VERTICES AS FOLLOWING
	0----1
	|\  /|
	| \/ |
	| /\ |
	|/  \|
	3----2
	"""

	print(g.RouteInspection())
