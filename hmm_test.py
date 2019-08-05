"""
TODO: 
  -Figure out issue with transition matrix. 
   Why does the match keep happening in a far-away place?
  -Add a check to ensure there are some matching candidates, 
  and do something different if there are not.
"""

import math
import sys
import os
import time

import numpy as np
from geopy.point import Point
from geopy.distance import great_circle as distance
import networkx as nx
import django
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/aaronsch/webapps/aarontakesawalk/trailzealot')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "trailzealot.settings")
django.setup()
from boulderhikes.models import Activity, Route
from mapmatching import mapmatch, Segment, Viterbi, MapMatcher
from myutils import find_closest_point, get_route_dist
from utils.hikehelper import get_bearing_lat_lon

"""
SPECIFY ALGORITHM PARAMETERS
These variables determine how the algorithm behaves. 
These are everything.
"""
# Choose a value for SIGMA_Z: the standard deviation of the gps location, meters
# Describes the uncertainty in each location measurement.
# Defines this algorithms tolerance to gps points being far away from route points.
# Also helps define minimum point-to-point distance between consecutive gps latlon points to match to route.
# A higher value implies less certainty in gps location, and more tolerance for distant route matches.
# SIGMA_Z = 1.4826*median_t(z - x) ~ 2.97 to 3.5 m, smaller for more points included.
SIGMA_Z = 5.0  
#SIGMA_Z = 3.5 

# Choose a number of points to skip from reduced-point device gps tracks. Generally set to 0, unless speeding up algorithm.
num_skip = 10 

# Choose a value for beta: the tolerance for the difference between point-to-point gps distances and distances along route.
# A smaller value is indicative as a lower tolerance for this difference.
# beta = median_t(dt)/ln(2) ~ 0.149 to 0.33. Lower if more data points? Not bad estimating! 
# In other words, a ballpark for beta is the median of the ACTUAL difference between route and point distances.
#BETA = 0.3
#BETA = 0.2
#BETA = 0.5
BETA = 2.0

# Maximum distance between gps point (which may be wrong), and route point (which theoretically is not wrong at all)
# Making this number too small will result in a) the algorithm only choosing nearby points (not flexible), or b) not finding a path at all
DIST_ERR_TOL = 50.0 # meters. The paper uses 200m, but I suspect that is a bit high for my application.

# User gps coordinates
"""
# (Desktop version)
import pickle
with open('activity_green-gregory.p', "rb") as fp:   # Unpickling
    data_dict_emission = pickle.load(fp)
activity = data_dict_emission
emission_times_unreduced = np.array(activity['times'])
emission_latlons_unreduced = np.array(activity['latlons'])
"""
# (Server version)
activity = Activity.objects.get(number=185) 
emission_times_unreduced = np.array(activity.src_data['times'][0:200])
emission_points_unreduced = np.array([Point(latlon[::-1]) 
    for latlon in activity.get_latlon_coords()[0:200]])
#emission_times_unreduced = np.array(activity.src_data['times'])
#emission_points_unreduced = np.array([Point(latlon[::-1]) 
#    for latlon in activity.get_latlon_coords()])

## Eliminate gps points that are not sufficiently far from last gps point. 
## (If gps reports movement, it might just be random error)
#mapmatcher = MapMatcher()
#emission_points_reduced, emission_times_reduced = \
#    mapmatcher.reduced_points(
#        emission_points_unreduced, 
#        emission_times_unreduced)
#
#print('Reduced number of gps points from %0.0f to %0.0f' 
#    % (len(emission_times_unreduced), len(emission_times_reduced)))

# Route coordinates
"""
# (Desktop version)
with open('route_green-mountain-gregory-canyon-ranger-trails.p', "rb") as fp:   # Unpickling
    data_dict_route = pickle.load(fp)
num_skip_route = 0 # if > 0, reduces number of route points we are matching
node_latlons = np.array(data_dict_route['latlons'][::num_skip_route+1])
"""
# (Server version)
route = Route.objects.get(slug="green-mountain-gregory-canyon-ranger-trails")
num_skip_route = 0 # if > 0, reduces number of route points we are matching
latlons = route.get_latlon_coords()[:1000:num_skip_route+1]
#latlons = route.get_latlon_coords()[::num_skip_route+1]
node_points = [Point(latlon[::-1]) for latlon in latlons]

#segments = [Segment(node_points[i], node_points[i+1])
#    for i in range(0,len(node_points)-1)]

graph = nx.Graph()  # nx.DiGraph() for directional graph
for i in range(0,len(latlons)):
    graph.add_node(
        i, 
        pos=(node_points[i].latitude, node_points[i].longitude)
    )
for i in range(0,len(latlons)-1):
    graph.add_edge(
        i,
        i+1,
        key=i,
        length=distance(
            node_points[i],
            node_points[i+1]).meters
    )

points_matched, emission_points_matched = mapmatch(
    emission_points_unreduced,
    emission_times_unreduced,
    graph)

## Skip a specified number of already-reduced gps points
#stride = num_skip + 1
#emission_times_reduced_skipped = emission_times_reduced[::stride]
#emission_points_reduced_skipped = emission_points_reduced[::stride]
##emission_cadences = np.array(activity.src_data['cadences'][::stride]) # used once I am capable of using cadence data
#print('Further reduced number of gps points: %0.0f, after skipping every %0.0f points' % (len(emission_times_reduced_skipped),num_skip))
#
## Store some info for later
#dists_gps_cum = {emission_times_reduced[0]: 0.0}
#for t in range(1,len(emission_times_reduced)):
#    dists_gps_cum[emission_times_reduced[t]] =  \
#        dists_gps_cum[emission_times_reduced[t-1]] + \
#        distance(emission_points_reduced[t-1], 
#                 emission_points_reduced[t]).meters
#
## handy variables
#nT = len(emission_times_reduced_skipped) 
#nSeg = len(segments)
#dists_cum = np.cumsum([segment.get_length() for segment in segments])
#
## Construct array of distances (along route) between route segment midpoints
#route_distance_array = np.array([
#        0.5*(segments[i].get_length() + segments[i+1].get_length()) \
#        for i in range(0,nSeg-1)
#    ])
## Construct matrix of distances between segment midpoints
#route_distances = np.zeros((nSeg,nSeg),dtype=float)
#for i in range(0,nSeg):
#    route_distances[i,i+1:] = np.cumsum(route_distance_array[i:])
#    route_distances[i+1:,i] = np.cumsum(route_distance_array[i:])
#
##candidates = np.zeros((nT,nSeg,2))
#
## Construct list of candidates corresponding to the midpoint of each segment
#candidates = [segment.get_midpoint() for segment in segments]
#
#print("Assembling emission probability matrix")
##candidates = {} 
#E = mapmatcher.emission_probability(
#    emission_points_reduced_skipped,
#    candidates
#)
#
#probs_nonzero = [np.nonzero(E[t,:])[0] for t in range(0,T)]
#
#print("Assembling transition probability matrix")
#print("%0.0f transitions, %0.0f segments" % (nT-1, nSeg))
#A = mapmatcher.transition_probability(
#    emission_points_reduced_skipped,
#    candidates,
#    probs_nonzero
#)
#
##A = {}
#A = np.zeros((nT-1,nSeg,nSeg))
#for t in range(0,nT-1):
#    if t % 20 == 0:
#        print(t)
#    dist_between_emissions = dists_gps_cum[emission_times_reduced_skipped[t+1]]  \
#                            -dists_gps_cum[emission_times_reduced_skipped[t]]
#    for i in np.nonzero(E[t,:])[0]:
#        for j in np.nonzero(E[t+1,:])[0]:
#            #A[t,i,j] = 1.0
#            dist_along_route = route_distances[i,j]
#            speed_route = dist_along_route / \
#                (emission_times_reduced_skipped[t+1] \
#               - emission_times_reduced_skipped[t])
#            dt = abs(dist_between_emissions - dist_along_route)
#            dt_frac = dt/dist_between_emissions
#            # 9 m/s is 20.1mph / 3 min/mile, which is as fast as I will ever move.
#            # paper uses 2000 m, but I am interested in 1- to 10-second data. # ADS HERE: RE-EVALUATE!!
#            if speed_route < 9.0 and dt <= 200.0:
#                # reduce transition probability if route and gps distances are different
#                p_dt = (1/BETA)*math.exp(-dt/BETA)
#                A[t, i, j] = p_dt
#        # Normalize the transition probabilities to sum to 1.0
#        if sum(A[t,i,:]) > 0.0: 
#            A[t,i,:] = A[t,i,:]/sum(A[t,i,:])
#    if not np.any(A[t]):
#        print("No transitions available at timestep %0.0f" % (t,))
#
#v = Viterbi(E, A)
#z_pred = v.viterbi()

"""
PLOT THE ROUTE PATH, THE GPS PATH, AND THE CORRECTED GPS PATH.
"""
fig1, ax = plt.subplots()
plt.plot(
    [p.longitude for p in emission_points_unreduced],
    [p.latitude for p in emission_points_unreduced],
    'b-',
    lw=1.0)
#plt.plot(
#    [p.longitude for p in emission_points_reduced],
#    [p.latitude for p in emission_points_reduced],
#    'g-',
#    lw=2.0)
plt.plot(
    [p.longitude for p in node_points],
    [p.latitude for p in node_points],
    'k-')

max_dist_err = 0.0
key_list = list(points_matched.keys())
for t in range(0, len(key_list)):
    #seg_number_selected = z_pred[t]
    #route_point_i = candidates[seg_number_selected]
    route_point_i = points_matched[key_list[t]]
    dist_err = distance(
        emission_points_matched[t],
        route_point_i).meters
    #if dist_err > max_dist_err:
    #    max_dist_err = dist_err
    plt.plot(
        [emission_points_matched[t].longitude,
            route_point_i.longitude],
        [emission_points_matched[t].latitude,
            route_point_i.latitude],
        'r-',
        lw=1.0)
    
    """
    plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'rx')
    """
ratio = (1.0/math.cos(math.radians(node_points[0].latitude)))
ratio = 1.0
ax.set_aspect(ratio)
plt.savefig('mapmatch.png',dpi=900)

import sys
sys.exit()
"""
PLOT HORIZ SPEED vs TIME, VERT SPEED vs TIME, ELEVATION vs DISTANCE
All using the route data. Final data should not use the gps points, except as comparison.
""" 
route_points = [candidates[seg_number_selected] for seg_number_selected in z_pred]
dists_gps_unreduced = [distance(
    emission_points_unreduced[t-1],
    emission_points_unreduced[t]).meters \
    for t in range(1,len(emission_times_unreduced))]
dists_gps_reduced = [distance(
    emission_points_reduced[t-1],
    emission_points_reduced[t]).meters \
    for t in range(1,len(emission_times_reduced))]
dists_gps_reduced_skipped = [distance(
    emission_points_reduced_skipped[t-1],
    emission_points_reduced_skipped[t]).meters \
    for t in range(1,nT)]
dists_matched = [distance(
    route_points[t-1],
    route_points[t]).meters \
    for t in range(1,nT)]
dists_matched_rt = [route_distances[
    z_pred[t-1],
    z_pred[t]]  \
    for t in range(1,nT)]

speed_gps_unreduced = [dists_gps_unreduced[t-1]/  \
    (emission_times_unreduced[t] - emission_times_unreduced[t-1])  \
    for t in range(1,len(emission_times_unreduced))]
speed_gps_reduced = [dists_gps_reduced[t-1]/  \
    (emission_times_reduced[t] - emission_times_reduced[t-1])  \
    for t in range(1,len(emission_times_reduced))]
speed_gps_reduced_skipped = [dists_gps_reduced_skipped[t-1]/  \
    (emission_times_reduced_skipped[t] - emission_times_reduced_skipped[t-1])  \
    for t in range(1,nT)]
speed_matched = [dists_matched[t-1]/  \
    (emission_times_reduced_skipped[t] - emission_times_reduced_skipped[t-1])  \
    for t in range(1,nT)]
speed_matched_rt = [dists_matched_rt[t-1]/  \
    (emission_times_reduced_skipped[t] - emission_times_reduced_skipped[t-1])  \
    for t in range(1,nT)]

dists_gps_unreduced_cum = np.cumsum(dists_gps_unreduced)
dists_gps_reduced_cum = np.cumsum(dists_gps_reduced)
dists_gps_reduced_skipped_cum = np.cumsum(dists_gps_reduced_skipped)
dists_matched_cum = np.cumsum(dists_matched)
dists_matched_rt_cum = np.cumsum(dists_matched_rt)

fig, axs = plt.subplots(2)
# Distance plot
axs[0].plot(
    emission_times_unreduced[1::], 
    dists_gps_unreduced_cum,
    label="Unreduced")
axs[0].plot(
    emission_times_reduced[1::],
    dists_gps_reduced_cum,
    label="Reduced")
axs[0].plot(
    emission_times_reduced_skipped[1::],
    dists_gps_reduced_skipped_cum,
    label="Skipped")
axs[0].plot(
    emission_times_reduced_skipped[1::],
    dists_matched_cum,
    label="Matched")
axs[0].plot(
    emission_times_reduced_skipped[1::],
    dists_matched_rt_cum,
    label="Matched (on-route)")
axs[0].legend(loc='upper left')
# Velocity plot
axs[1].plot(
    emission_times_unreduced[1::], 
    speed_gps_unreduced,
    label="Unreduced")
axs[1].plot(
    emission_times_reduced[1::],
    speed_gps_reduced,
    label="Reduced")
axs[1].plot(
    emission_times_reduced_skipped[1::],
    speed_gps_reduced_skipped,
    label="Skipped")
axs[1].plot(
    emission_times_reduced_skipped[1::],
    speed_matched,
    label="Matched")
axs[1].plot(
    emission_times_reduced_skipped[1::],
    speed_matched_rt,
    label="Matched (on-route)")
axs[1].set_ylim(0,5)
axs[1].legend(loc='upper left')
plt.savefig('graphs.png',dpi=900)

# STATISTICS
dts = []
ratios = []
dist_errs = []
#for t in range(0,len(emission_latlons_reduced)-1):
for t in range(0,nT-1):
    dist_between_emissions = dists_gps_cum[emission_times_reduced_skipped[t+1]] \
                           - dists_gps_cum[emission_times_reduced_skipped[t]] 
    #dist_between_emissions = get_dist_lat_lon_m(emission_latlons_reduced[t],emission_latlons_reduced[t+1])
    dist_route = route_distances[z_pred[t],z_pred[t+1]]
    #dist_route = get_route_dist(candidates,segments,dists_cum,z_pred[t],z_pred[t+1],t,t+1)
    dts.append(abs(dist_between_emissions - dist_route))
    ratios.append(abs(dist_between_emissions - dist_route)/dist_between_emissions)
    dist_errs.append(distance(
        emission_points_reduced_skipped[t],
        candidates[z_pred[t]]).meters)
print("Route differences:")
print("Maximum difference between two gps points' distance and the chosen route distance: "+str(max(dts))+" m")
print("for a ratio of: "+str(ratios[np.argmax(dts)]))
print("Maximum ratio of route distance difference to gps distance: "+str(max(ratios)))
print("Beta (approx): "+str(np.median(dts)/math.log(2)))

print("Location errors:")
print("Maximum distance between a gps point and its matching route point: "+str(max(dist_errs))+" m")
print("SIGMA_Z (approx): "+str(1.4826*np.median(dist_errs))+" m")

