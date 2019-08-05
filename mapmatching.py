import math
import gc
import time

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import networkx as nx
from geopy.point import Point
from geopy.distance import great_circle as distance

def mapmatch(points, times, graph):
    """Matches a series of gps points to defined trail segments.

    Args:
      points (array-like): Point objects for each gps point.
                           Could be (lat, lon) tuples, [lon, lat] lists,
                           or geopy.point.Point() instances.
      times (array-like): Time from start of workout for each gps point.
      graph (networkx graph): Graph of trail network.

    Returns:
      A list, same length as points and times, with matched points 
      along route. The matching algorithm reduces the number of
      gps points to match, so the entries that are skipped contain
      None.

    """
    start = time.time()

    # Eliminate gps points that are not sufficiently far from previous 
    # gps point. If gps reports movement, it might just be random error.
    mapmatcher = MapMatcher()
    points_reduced, times_reduced = mapmatcher.reduced_points(points, times)
    print('Reduced number of gps points from %0.0f to %0.0f'
        % (len(times), len(times_reduced)))
    T = len(times_reduced)
    C = graph.number_of_nodes()

    # NOT NECESSARY NOW: Construct matrix of distances between segment midpoints
    #route_distances = np.zeros((C,C),dtype=float)
    #for ci in range(0,C):
    #    for cj in range(0,C):
    #        route_distances[ci,cj] = nx.shortest_path_length()

    print("Assembling emission probability matrix")
    E = mapmatcher.emission_probability(
        points_reduced,
        graph
    )

    #probs_nonzero = [np.nonzero(E[t,:])[0] for t in range(0,T)]
    #probs_nonzero = [np.nonzero(E[t,:]) for t in range(0,T)]
 
    print("Assembling transition probability matrix")
    print("%0.0f transitions, %0.0f segments" % (T-1, C))
    A = mapmatcher.transition_probability(
        points_reduced,
        times_reduced,
        graph,
        E
        #probs_nonzero
    )

    print("Performing viterbi process")
    v = Viterbi(E, A)
    node_ids_matched = v.viterbi()
    node_dict = dict(graph.nodes.data())
    print(time.time()-start)
    return ({times_reduced[t]: Point(node_dict[node_ids_matched[t]]['pos'])  \
                for t in range(0,T)},
            points_reduced)

class Segment(object):
    """One straight segment of trail.

    Longer class information....
    Longer class information....

    Attributes:
      node_start: The point where the segment originates.
      node_end: The point where the segment terminates.
    """

    def __init__(self, node_start, node_end):
        """Constructor for the Segment object
        
        Args:
          node_start (point): The point where the segment originates.
          node_end (point): The point where the segment terminates.
        """
        self.node_start = node_start
        self.node_end = node_end

    def get_length(self):
        """Returns the length of the segment using geodesic distance"""
        return distance(self.node_start, self.node_end).meters

    def get_midpoint(self):
        """Returns midpoint of segment as a Point object"""
        return Point(
            0.5*(self.node_start.latitude+self.node_end.latitude),
            0.5*(self.node_start.longitude+self.node_end.longitude)           
        )

    def find_closest_point(self, emission_point):
        """Returns closest point to emission point along segment

        TODO: Documentation
        """ 
        node_start = self.node_start
        node_end = self.node_end
        lat_emission = emission_point.latitude
        lon_emission = emission_point.longitude
        if node_start.longitude == node_end.longitude:
            # Special case: perfectly N-S, or no movement.
            lon = node_start.longitude
            lat_min = min([node_start.latitude, node_end.latitude])
            lat_max = max([node_start.latitude, node_end.latitude])
            if lat_emission <= lat_min:
                return Point(lat_min, lon)
            elif lat_emission >= lat_min:
                return Point(lat_max, lon)
            else:
                return Point(lat_emission, lon)
        else:
            # find closest point along this horizontal or sloped line
            node_start_is_first = \
                0 == np.argmin([node_start.longitude, node_end.longitude])
            if node_start_is_first:
                lon_1 = node_start.longitude
                lat_1 = node_start.latitude
                lon_2 = node_end.longitude
                lat_2 = node_end.latitude
            else:
                lon_1 = node_end.longitude
                lat_1 = node_end.latitude
                lon_2 = node_start.longitude
                lat_2 = node_start.latitude
            slope = (lat_2-lat_1)/(lon_2-lon_1)
            if slope != 0.0:
                slope_perp = -1/slope
                b = lat_1 - slope*lon_1
                b_perp = lat_emission - slope_perp*lon_emission
                lon_perp = (b-b_perp)/(slope_perp-slope)
                lat_perp = slope*lon_perp + b
            else:  # completely horizontal line
                lon_perp = lon_emission
                lat_perp = lat_1
            if lon_perp <= lon_1:
                return Point(lat_1, lon_1)
            elif lon_perp >= lon_2:
                return Point(lat_2, lon_2)
            else:
                return Point(lat_perp, lon_perp)


class Candidate(object):
    """Candidate for corrected gps coordinate on trail network.

    TODO: doc

    """
    
    def __init__(self, point, dist, prob):
        """Constructor for the Candidate object
   
        Args:
          point (point): The point on route 
            that may match the gps point.
          dist (float): The distance between the candidate 
            and the emission (m).
          prob (float): The probability that this candidate 
            produced this emission.
        """
        self.point = point
        self.dist = dist
        self.prob = prob

class Viterbi(object):
    """Implementation of the Viterbi algorithm.

    Finds the most likely path through a Hidden Markov Model.

    Attributes:
      E: Emission probability matrix.
      A: Transition probability matrix.
      initial_probs: State probabilities at the first timestep, often assumed.

    """

    def __init__(self, E, A, initial_probs=None):
        """Constructor for the Viterbi object
        
        Args:
          E (ndarray): Emission probability matrix.
          A (ndarray): Transition probability matrix.
          initial_probs (ndarray): Probabilities at the first timestep, often assumed.
        """
        self.E = E
        self.A = A
        if initial_probs is not None:
            self.initial_probs = initial_probs
        else:
            self.initial_probs = E[0,:].todense()

    def viterbi(self):
        """Perform the Viterbi process and return a list of state numbers

        TODO: 
            -doc
            -Catch spots where A does not have a valid transition (avoid log err)
        """
        E = self.E
        A = self.A
        initial_probs = self.initial_probs
        T, S = np.shape(E)
        P = np.full((T, S), -np.inf, dtype=np.float32)
        back = np.zeros((T, S), dtype=np.int)
        z_pred = np.zeros((T,), dtype=np.int)
        for s in np.nonzero(initial_probs)[1]:
            P[0,s] = 2*initial_probs[:,s]
        for t in range(1, T):
            cols_nonzero_Et = E.indices[E.indptr[t]:E.indptr[t+1]]
            #cols_nonzero_At1 = np.ediff1d(A[t-1].indptr).nonzero()[0]
            cols_nonzero_At1 = np.unique(A[t-1].indices)
            inds_nonzero_cur = np.intersect1d(
                cols_nonzero_Et,
                cols_nonzero_At1)
            inds_nonzero_prev = E.indices[E.indptr[t-1]:E.indptr[t]]
            #indices = np.argmax(
            #    [[P[t-1,s_i]+A[t-1][s_i,s_j] for s_j in inds_nonzero_cur]  \
            #      for s_i in inds_nonzero_prev], 
            #    axis=0)
            for s_j in inds_nonzero_cur:
                #cols_nonzero_Et1 = E.indices[E.indptr[t-1]:E.indptr[t]]
                #rows_nonzero_At1 = A[t-1].indices[
                #    A[t-1].indptr[s_j]:A[t-1].indptr[s_j+1]]
                #inds_nonzero_prev = np.intersect1d(
                #    cols_nonzero_Et1,
                #    rows_nonzero_At1)
                index = np.argmax([P[t-1,s_i]  \
                                  +A[t-1][s_i,s_j]  \
                                  for s_i in inds_nonzero_prev])
                back[t,s_j] = inds_nonzero_prev[index]
                P[t,s_j] = P[t-1,back[t,s_j]]  \
                          +A[t-1][back[t,s_j],s_j]  \
                          +E[t,s_j]
                #P[s_j,t] = np.max([P[s_i,t-1]  \
                #                  +A[t-1][s_i,s_j]  \
                #                  +E[t,s_j]  \
                #           for s_i in inds_nonzero_prev])
            if not np.any(P[t,:] > -np.inf):
                print("No  path forward at timestep number %0.0f" % (t,))

        # Find the state associated with the max probability at last timestep.
        z_pred[T-1] = np.argmax([P[T-1,k] for k in range(0, S)])
        
        # Find the states associated with the most likely path.
        for t in range(1,T)[::-1]:
            z_pred[t-1] = back[t,z_pred[t]]
        return z_pred

import sys
from pympler import muppy, summary
class MapMatcher(object):
    """Master class for matching gps coordinates to route coordinates

    Attributes:
      SIGMA_Z: A float representing the standard deviation of gps error.
      BETA: A float representing the tolerance for changing gps route 
        length during matching.
      DIST_ERR_TOL: Maximum distance between gps point and matched
        route point to consider.
      #pts: list of gps points
      #pts_reduced: list of gps points after eliminating points too close together

    #Returns:
    #    A dict with two ndarrays for corrected gps points and associated times.

    TODO: 
        -complete doc.
    """

    SIGMA_Z = 5.0
    BETA = 2.0
    DIST_ERR_TOL = 50.0 

    def __init__(self, *args, **kwargs): 
        """Constructor for the MapMatcher object
        
        TODO:
          * Doc
        """

        self.SIGMA_Z = kwargs.get('SIGMA_Z', 5.0)
        self.BETA = kwargs.get('BETA', 2.0)
        self.DIST_ERR_TOL = kwargs.get('DIST_ERR_TOL', 50.0)
        #self.segments = [Segment(route_points[i], route_points[i+1])
        #    for i in range(0,len(route_points)-1)]

    def reduced_points(self, points, times):
        """Eliminates gps points that are too close together

        Args:
          points (array-like): Point objects for each gps point.
                               Could be (lat, lon) tuples, [lon, lat] lists,
                               or geopy.point.Point() instances.
          times (array-like): Time from start of workout for each gps point.
        """
        if type(points[0]) in (Point, tuple):
            points_conv = points
        elif type(points[0]) in (list, np.ndarray):
            points_conv = [(point[1], point[0]) for point in points]
        else:
            # Raise error
            return

        times_reduced_list = [times[0]]
        points_reduced_list = [points_conv[0]] 
        t1 = 0 # index of 1st point
        # index of 2nd point (may not be far enough away from first point)
        t2 = 1 
        for t in range(0,len(times)-1):
            dist = distance(points_conv[t1], points_conv[t2]).meters
            if dist < 2*self.SIGMA_Z:
                t2 = t2 + 1
            else:
                # we are sufficiently far from previous point. 
                points_reduced_list.append(points_conv[t2])
                times_reduced_list.append(times[t2])
                t1 = t2
                t2 = t1 + 1
        return (
            np.array(points_reduced_list),
            np.array(times_reduced_list)
        )

    def emission_probability(self, points, graph):
        """Assembles and returns an emission probability matrix

        Args:
          points (array-like): Point objects for each gps point.
          graph (networkx graph): graph of trail network. 
        
        Returns: 
          A CSR matrix (good for sparse matrices with row slicing) 
          of size (points, nodes) containing the log-probability 
          that each gps point belongs at each node. 
        """

        T = len(points)
        C = graph.number_of_nodes()
        #E = np.zeros((T, C), dtype=float)
        #E = np.full((T, C), -np.inf, dtype=np.float32)  # logprobs
        data = []
        row_ind = []
        col_ind = []
        # find all nodes that are within the distance tolerance
        for c, posdata in graph.nodes(data=True):
            for t in range(0,T):
                dist = distance(
                    points[t], 
                    posdata['pos']).meters
                if dist < self.DIST_ERR_TOL:
                    #E[t,c] = math.exp(-0.5*(dist/self.SIGMA_Z)**2)
                    #E[t,c] = -0.5*(dist/self.SIGMA_Z)**2
                    data.append(-0.5*(dist/self.SIGMA_Z)**2)
                    row_ind.append(t)
                    col_ind.append(c)
            #for c in range(0,C):
            #    #candidate = candidates[c]
            #    dist = distance(
            #        point,
            #        candidate).meters
            #    if dist < self.DIST_ERR_TOL:
            #        E[t,c] = ((2*math.pi)**0.5*self.SIGMA_Z)**-1 \
            #            * math.exp(-0.5*(dist/self.SIGMA_Z)**2) 
            #if sum(E[t,:]) > 0:
            #    E[t,:] = E[t,:] / sum(E[t,:])
            #else:

            #if not np.any(E[t,:] > -np.inf):
            #    print("No candidates at timestep number %0.0f" % (t,))
        E = csr_matrix((data, (row_ind, col_ind)), 
                       shape=(T, C), 
                       dtype=np.float32)  # logprobs
        return E

    def transition_probability(self, points, times, graph, E):  #probs_nonzero):
        """Assembles and returns a transition probability matrix

        Args:
          points (array-like): Point objects for each gps point.
          times (array-like): Time since start of workout for each gps point.
          graph (networkx graph): Graph of trail network.
          probs_nonzero (array-like): At each time, the indices of the 
            emission probability matrix which are nonzero.
        
        Returns: 
          A ndarray of size (points-1, nodes, nodes) containing 
          the probability, for each timestep, of transitioning from one 
          candidate node to another. 

        """
        T = len(points)
        C = graph.number_of_nodes()
        #route_distances_dict = nx.shortest_path_length(graph, weight='length')
        A = []  # will be a list of CSC matrices (sparse) for each timestep
        for t in range(0,T-1):
            if t % 20 == 0:
                print(t)
            dist_between_emissions = distance(
                points[t],
                points[t+1]).meters
            data = []
            rowind = []
            colind = []
            #for i in E.getrow(t).nonzero()[1]:
            for i in E.indices[E.indptr[t]:E.indptr[t+1]]:
                # Distance I would cover if running 9 m/s 
                # (3 min/mile, as fast as I will ever move)
                dist_9 = 9.0*(times[t+1]-times[t])
                # Distance that is larger than dist between emissions
                # by a factor. 
                # (paper uses 2000 m, but I think a factor makes sense)
                dist_dt = 2.0*dist_between_emissions
                radius = min(dist_9, dist_dt)
                subgraph = nx.ego_graph(graph, i,
                    radius=radius, distance='length')
                route_distances_i = nx.shortest_path_length(
                    subgraph, source=i, weight='length')
                #route_distances_i = nx.shortest_path_length(
                #    graph,
                #    source=i,
                #    weight='length')
                #for j in E.getrow(t+1).nonzero()[1]: 
                for j in route_distances_i.keys(): 
                    dt = abs(dist_between_emissions - route_distances_i[j])
                    #dt_frac = dt/dist_between_emissions
                    # reduce transition probability if route and 
                    # gps distances are different
                    data.append(-dt/self.BETA)
                    rowind.append(i)
                    colind.append(j)
                # Normalize the transition probabilities to sum to 1.0?
            A.append(csr_matrix(
                (data,(rowind, colind)), 
                shape=(C, C), 
                dtype=np.float32))
            #if not np.any(A[t]>-np.inf):
            if A[t].getnnz() == 0:
                print("No transitions available at timestep %0.0f" % (t,))
            #gc.collect()
        return A
