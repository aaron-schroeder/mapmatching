from geopy.point import Point
from geopy.distance import great_circle as distance
import numpy as np

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
            self.initial_probs = E[0,:]

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
        P = np.zeros((S, T))
        P.fill(-99999.0) #  Using log-probabilities, so this is filling in for zero.
        back = np.zeros((S, T))
        z_pred = np.zeros((T,), dtype=np.int)

        for s in np.nonzero(initial_probs):
            P[s,0] = np.log(initial_probs[s]**2)
        for t in range(1, T):
            #inds_nonzero_cur = np.nonzero(A[t-1,np.nonzero(E[t-1,:])[0],:])[0]
            inds_nonzero_cur = np.nonzero(E[t,:])[0]
            for s_j in inds_nonzero_cur:
                inds_nonzero_prev = np.intersect1d(
                    np.nonzero(E[t-1,:])[0],
                    np.nonzero(A[t-1,:,s_j])[0]
                )
                P[s_j,t] = np.max([P[s_i,t-1]  \
                           +np.log(A[t-1,s_i,s_j])  \
                           +np.log(E[t,s_j])  \
                           for s_i in inds_nonzero_prev])
                index = np.argmax([P[s_i,t-1]  \
                           +np.log(A[t-1,s_i,s_j])  \
                           for s_i in inds_nonzero_prev])
                back[s_j,t] = np.nonzero(E[t-1,:])[0][index]
            if not np.any(P[:,t]):
                print("No  path forward at timestep number %0.0f" % (t,))

        # Find the state associated with the max probability at last timestep.
        z_pred[T-1] = np.argmax([P[k,T-1] for k in range(0, S)])
        
        # Find the states associated with the most likely path.
        for t in range(1,T)[::-1]:
            z_pred[t-1] = back[z_pred[t],t]

        return z_pred
