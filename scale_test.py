"""Matches a simulated, drifted gps trace to the database route.

Samples points along a database route, which equate to ground truth
measurements. Using these ground truth points as a base, we simulate 
gps drift effects to create a simulated gps trace. We use the
map matching package to attempt to match the simulated gps points 
to the ground truth points along the route. We summarize performance
with speed and accuracy metrics.

TODO: 
    * Determine relationship between sigma_z (radial std dev)
      and x- and y- std devs. Consider that radius cannot be negative.
    * Add route smoothing mechanism and test effect of route simplicity
      on performance (time and accuracy).
    * Predict instantaneous hiking rate based on hiker intensity and
      route elevation profile. Right now we assume a constant pace.
    * Better docs.
"""

import os
import sys
import math
import pickle
from statistics import median, mean

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import django

from simplification.cutil import simplify_coords
from geopy.distance import great_circle as distance
from geopy.point import Point

sys.path.insert(0, '/home/aaronsch/webapps/aarontakesawalk/trailzealot')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "trailzealot.settings")
django.setup()
from boulderhikes.models import Activity, Route
from mapmatching import mapmatch
from utils.hikehelper import get_bearing_lat_lon

def fractional_index(seq, val):
    """Finds fractional index where a value occurs in a sequence.

    Sequence should be monotonically increasing, and does not need to
    contain exact value.

    Example:
        >>> my_seq = [0.0, 2.0, 5.0, 7.0, 8.1]
        >>> fractional_index(my_seq, 1.5)
        0.75
        >>> fractional_index(my_seq, 3.5)
        1.5
        >>> fractional_index(my_seq, 5.0)
        2.0
    """

    r = np.where(np.diff(np.sign(seq - val)) != 0)
    idx = r + (val - seq[r]) / (seq[r + np.ones_like(r)] - seq[r])
    idx = np.append(idx, np.where(seq == val))
    idx = np.sort(idx)
    return(idx[0])

# Get route coordinate list and determine length of each segment,
# defined as the line connecting two consecutive route coordinates.
route_slug = "green-mountain-gregory-canyon-ranger-trails"
route = Route.objects.get(slug=route_slug)
latlons_rt = np.array(route.get_latlon_coords())
nseg_rt = len(latlons_rt)-1
seg_lengths_rt = [distance(
    latlons_rt[i][::-1], latlons_rt[i+1][::-1]).meters
    for i in range(0,nseg_rt)]

# Generate ground truth points along route.
# We sample at a constant time rate, and we assume a constant pace.
# Later, will use predicted pace based on hiker intensity and grades.
delta_t = 5.0  # seconds
sampling_dist = (2.0*3600/1609.34)*delta_t  # meters
route_dists_cum = np.insert(np.cumsum(seg_lengths_rt),0,0.0)
gps_dists_cum = np.arange(0, route_dists_cum[-1], sampling_dist)
n_gps = len(gps_dists_cum)
idxs = [fractional_index(route_dists_cum, gps_dist)
    for gps_dist in gps_dists_cum]
latlons_truth = np.zeros((n_gps,2))
for i,idx in enumerate(idxs):
    i_st = math.floor(idx)
    i_ed = math.ceil(idx)
    if i_st != i_ed:
        latlons_truth[i,:] = (i_ed-idx)*latlons_rt[i_st]  \
                            +(idx-i_st)*latlons_rt[i_ed]
    else:
        latlons_truth[i,:] = latlons_rt[i_st]

# Simulate gps readings at ground truth points.
# We use a Gauss-Markov model for drift. At each point, 
# drift has a random component and a memory of the previous drift value.
beta = 1.0/20  # (sec^-1)^-1 = seconds
s = beta*delta_t
sigma_z = 5.0  # meters
nu_xs = np.random.normal(0, sigma_z/2**0.5, n_gps)
nu_ys = np.random.normal(0, sigma_z/2**0.5, n_gps)
delta_xs = [nu_xs[0]]
delta_ys = [nu_ys[0]]
for i in range(1, n_gps):
    delta_xs.append(nu_xs[i] + (1-s)*delta_xs[i-1])
    delta_ys.append(nu_ys[i] + (1-s)*delta_ys[i-1])
ratio = np.sin(np.radians(latlons_truth[0,1]))
latlons_gps = np.zeros(np.shape(latlons_truth))
for i in range(0, n_gps):
    #print('%d: (%0.1f, %0.1f), (%0.1f, %0.1f)' % (i,nu_xs[i],delta_xs[i],nu_ys[i],delta_ys[i]))
    p=Point(latlons_truth[i,::-1])
    dist=(delta_xs[i]**2 + delta_ys[i]**2)**0.5
    bearing=np.degrees(np.arctan2(delta_xs[i], ratio*delta_ys[i]))
    d=distance(meters=dist)
    p_delta=d.destination(p, bearing=bearing)
    latlons_gps[i,0]=p_delta.longitude
    latlons_gps[i,1]=p_delta.latitude
times_gps = np.array([i*delta_t for i in range(0, n_gps)])

# Perform map matching procedure. Generate a graph of the trail network,
# and pass the graph to the mapmatching algorithm, along with the route
# coordinates and simulated gps trace coordinates.
graph = nx.Graph() 
for i in range(0,len(latlons_rt)):
    graph.add_node(
        i,
        pos=(latlons_rt[i][1], latlons_rt[i][0])
    )
for i in range(0,nseg_rt):
    graph.add_edge(
        i,
        i+1,
        key=i,
        length=distance(
            latlons_rt[i][::-1],
            latlons_rt[i+1][::-1]).meters
    )
points_matched, emission_points_matched = mapmatch(
    latlons_gps,
    times_gps,
    graph)

# Pickle the results
match_dict = {
    'beta': beta,
    'sigma_z': sigma_z,
    'delta_t': delta_t,
    'latlons_rt': latlons_rt,
    'latlons_truth': latlons_truth,
    'latlons_gps': latlons_gps,
    'times_gps': times_gps,
    'graph': graph,
    'points_matched': points_matched,
    'emission_points_matched': emission_points_matched
}
pickle.dump(match_dict, open('%s_%d.p' % (route_slug, delta_t), 'wb'))
