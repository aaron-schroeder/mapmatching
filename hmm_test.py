from boulderhikes.models import Activity, Route
from boulderhikes.mapmatching import *
from utils.hikehelper import get_dist_lat_lon_m, get_bearing_lat_lon, get_distance_from_start_m
import numpy as np
import math

"""
SPECIFY ALGORITHM PARAMETERS
These variables determine how the algorithm behaves. 
These are everything.
"""
# Choose a value for sigma_z: the standard deviation of the gps location, meters
# Describes the uncertainty in each location measurement.
# Defines this algorithms tolerance to gps points being far away from route points.
# Also helps define minimum point-to-point distance between consecutive gps latlon points to match to route.
# A higher value implies less certainty in gps location, and more tolerance for distant route matches.
# sigma_z = 1.4826*median_t(z - x) ~ 2.97 to 3.5 m, smaller for more points included.
sigma_z = 5.0  
#sigma_z = 3.5 

# Choose a number of points to skip from reduced-point device gps tracks. Generally set to 0, unless speeding up algorithm.
num_skip = 5 

# Choose a value for beta: the tolerance for the difference between point-to-point gps distances and distances along route.
# A smaller value is indicative as a lower tolerance for this difference.
# beta = median_t(dt)/ln(2) ~ 0.149 to 0.33. Lower if more data points? Not bad estimating! 
# In other words, a ballpark for beta is the median of the ACTUAL difference between route and point distances.
#beta = 0.3
#beta = 0.2
beta = 0.5

# Maximum distance between gps point (which may be wrong), and route point (which theoretically is not wrong at all)
# Making this number too small will result in a) the algorithm only choosing nearby points (not flexible), or b) not finding a path at all
dist_err_tol = 50.0 # meters. The paper uses 200m, but I suspect that is a bit high for my application.

# ADS HERE: WANT TO REDUCE POINTS BASED ON SPACING (X * SIGMA_Z), AND THEN FURTHER REDUCE FROM THERE WITH NUM_SKIP

# DATA INPUT ---------------------------------------------------------
# User Green Mtn data

# User gps coordinates
"""
# Desktop version
import pickle
with open('activity_green-gregory.p', "rb") as fp:   # Unpickling
    data_dict_emission = pickle.load(fp)
activity = data_dict_emission
emission_times_unreduced = np.array(activity['times'])
emission_latlons_unreduced = np.array(activity['latlons'])
"""
# Server version
activity = Activity.objects.get(number=98) # I think this is ranger trail?
emission_times_unreduced = np.array(activity.src_data['times'])
emission_latlons_unreduced = np.array(activity.get_latlon_coords())

# route coordinates
"""
# Desktop version
with open('route_green-mountain-gregory-canyon-ranger-trails.p', "rb") as fp:   # Unpickling
    data_dict_route = pickle.load(fp)
num_skip_route = 0 # if > 0, reduces number of route points we are matching
node_latlons = np.array(data_dict_route['latlons'][::num_skip_route+1])
"""
# Server version
route = Route.objects.get(slug="green-mountain-gregory-canyon-ranger-trails")
num_skip_route = 0 # if > 0, reduces number of route points we are matching
node_latlons = np.array(route.get_latlon_coords()[::num_skip_route+1])

# END OF DATA INPUT ---------------------------------------------------------

# Eliminate gps points that are not sufficiently far from last gps point. 
# (If gps reports movement, it might just be random error)
print('Unreduced number of gps points: %0.0f' % len(emission_times_unreduced))
emission_times_reduced_list = [emission_times_unreduced[0]]      # initializing lists 
emission_latlons_reduced_list = [emission_latlons_unreduced[0]]  # which will be appended in loop
t1 = 0 # index of first point
t2 = 1 # index of second point (which may or may not be far enough away from first point)
#for t in range(0,len(emission_times_unreduced)-2):
for t in range(0,len(emission_times_unreduced)-1):
    point_1 = emission_latlons_unreduced[t1]
    point_2 = emission_latlons_unreduced[t2]
    emission_dist = get_dist_lat_lon_m(point_1,point_2)
    if emission_dist < 2*sigma_z:
        # increment the index of the point to check, then go through loop to again to see if we are far enough away from previous point
        t2 = t2 + 1 
        #continue
    else:
        # we are sufficiently far from previous point. Add this point to the list.
        emission_latlons_reduced_list.append(emission_latlons_unreduced[t2])
        emission_times_reduced_list.append(emission_times_unreduced[t2])
        t1 = t2     # set indices
        t2 = t1 + 1 # for next loop
    #t2 = t2 + 1 # happens no matter what. Consider moving out of if statement.
emission_times_reduced = np.array(emission_times_reduced_list)
emission_latlons_reduced = np.array(emission_latlons_reduced_list)
"""
# NOT IN USE: Ensure the last emission is included in the reduced-point scheme
emission_latlons = np.concatenate((emission_latlons_reduced,[emission_latlons[-1]]))
emission_times = np.concatenate((emission_times_reduced,[emission_times[-1]]))
"""
#print(len(emission_times))
print('Reduced number of gps points: %0.0f' % len(emission_times_reduced))

# Skip a specified number of already-reduced gps points
stride = num_skip + 1
emission_times_reduced_skipped = emission_times_reduced[::stride]
emission_latlons_reduced_skipped = emission_latlons_reduced[::stride]
#emission_cadences = np.array(activity['latlons'][::stride]) # used once I am capable of using cadence data
print('Further reduced number of gps points: %0.0f, after skipping every %0.0f points' % (len(emission_times_reduced_skipped),num_skip))

# Store some info for later
"""
#dists_gps_cum = {emission_times_unreduced[0]: 0.0}
#for t in range(0,len(emission_times_unreduced)-1):
    #dists_gps_cum[emission_times_unreduced[t+1]] = dists_gps_cum[emission_times_unreduced[t]] + \
    #    get_dist_lat_lon_m(emission_latlons_unreduced[t], emission_latlons_unreduced[t+1])
"""
dists_gps_cum = {emission_times_reduced[0]: 0.0}
for t in range(0,len(emission_times_reduced)-1):
    dists_gps_cum[emission_times_reduced[t+1]] = dists_gps_cum[emission_times_reduced[t]] + \
        get_dist_lat_lon_m(emission_latlons_reduced[t], emission_latlons_reduced[t+1])

"""
REDUCE NUMBER OF GPS POINTS TO CONSIDER
Ensure each gps point considered is spaced a minimum distance from last point
"""
# handy variables
nT = len(emission_times_reduced_skipped) 
nSeg = len(node_latlons)-1 # two nodes define one segment
dists_cum = get_distance_from_start_m(node_latlons) # point-to-point cumulative distances along route nodes. Get your mind out of the gutter.

"""
CONSTRUCT CANDIDATE OBJECTS
AND EMISSION PROBABILITY MATRIX
For each gps point, there is a dictionary of candidates that are {'latlon': [lon,lat], 'dist': dist}
"""
print("Gathering candidate info")

#candidates = np.zeros((nT,nSeg,2))
#candidates = {time:{} for time in emission_times_reduced_skipped } 
#emission_probs = {time:{} for time in emission_times_reduced_skipped}
candidates = {} 
emission_probs = {}
candidate_ti = find_closest_point(emission_latlons_reduced_skipped[0],[node_latlons[0],node_latlons[1]]) # returns a latlon pair (What is this doing?) 

for t in range(0,nT):
    emission_time = emission_times_reduced_skipped[t]
    candidates[t] = {}
    emission_probs[t] = {}
    for i in range(0,nSeg): 
       candidate_ti = find_closest_point(emission_latlons_reduced_skipped[t],[node_latlons[i],node_latlons[i+1]]) # retur`ns a latlon pair
       dist_candidate = get_dist_lat_lon_m(emission_latlons_reduced_skipped[t],candidate_ti)
       if dist_candidate < dist_err_tol: # defined at start of script 
           candidates[t][i] = {
               'latlon': candidate_ti, 
               'dist':   dist_candidate,
           }
           emission_probs[t][i] = ((2*math.pi)**0.5*sigma_z)**-1 * math.exp(-0.5*(dist_candidate/sigma_z)**2)
           #sum_emission_probs_t = sum_emission_probs_t + emission_probs_t[i]
    #candidates[emission_times[t]] = candidates_t # for the gps location at this timestep, documents all segments within the maximum distance.
    #candidates[t] = candidates_t # for the gps location at this timestep, documents all segments within the maximum distance.
   
    # Normalize the emission probabilities at the timestep so they sum to 1.0
    #sum_emission_probs_t = sum([emission_probs[t][seg_number] for seg_number in emission_probs[t].keys()])
    sum_emission_probs_t = sum(emission_probs[t].values())
    emission_probs[t] = {seg_number: emission_probs[t][seg_number]/sum_emission_probs_t for seg_number in emission_probs[t].keys() }

"""
# BRUTE FORCE NEAREST POINT SEARCH: BREAK EACH SEGMENT UP INTO PIECES AND FIND THE DISTANCE TO EACH NODE
for seg_number in range(0,nSeg):
    # Check how long the segment is
    seg_length = get_dist_lat_lon(node_latlons[seg_number],node_latlons[seg_number+1])

    # break the segment up into pieces, so the pieces are less than 5m long.
    num_pieces = int(math.ceil(seg_length/1.0))
    node_lon_arr = np.linspace(node_latlons[seg_number][0],node_latlons[seg_number+1][0],num_pieces + 1)
    node_lat_arr = np.linspace(node_latlons[seg_number][1],node_latlons[seg_number+1][1],num_pieces + 1)

    for t in range(0,nT):
        emission_time = emission_times[t]
        closest_index = np.argmin([ get_dist_lat_lon(emission_latlons[t],[node_lon_arr[i],node_lat_arr[i]]) for i in range(0,num_pieces+1) ])
        dist_candidate = get_dist_lat_lon(emission_latlons[t],[node_lon_arr[closest_index],node_lat_arr[closest_index]])
        if dist_candidate <= dist_err_tol: # close enough to route segment to consider
            candidates[emission_time][seg_number] = {
                'latlon': [ node_lon_arr[closest_index], node_lat_arr[closest_index] ],
                'dist':   dist_candidate,
            }
"""
    # ADD A CHECK TO ENSURE THERE ARE SOME MATCHING CANDIDATES, AND DO SOMETHING DIFFERENT IF THERE ARE NOT.

"""
CONSTRUCT INITIAL PROBABILITY VECTOR
For each state: prob that the initial state is that state
Then, the inital entry in P is the inital prob times the prob that (if the initial state is really that state) the observation would come true (i.e. the prob that the state would make this observation show up)
"""
print("Constructing initial probability vector")

initial_probs = np.zeros(nSeg)
P = np.zeros((nSeg,nT))      # Probability of the most likely path, at each timestep, that generates each state as its final state
back = np.zeros((nSeg,nT))   # Stores the actual most likely paths associated with probabilities in "P". 
z_pred = np.zeros(nT, dtype=np.int) # Stores the most likely path itself. The end result of our quest.
"""
initial_probs = {}
P = {0:{}}  #{} for t in range(0,nT)}      # Probability of the most likely path, at each timestep, that generates each state as its final state
back = {0: {}} #{} for t in range(0,nT)}   # Stores the actual most likely paths associated with probabilities in "P". 
z_pred = np.zeros(nT, dtype=np.int) # Stores the most likely path itself. The end result of our quest.
"""

for seg_number in candidates[0].keys():
    candidate = candidates[0][seg_number]
    dist_greatcircle = candidate['dist'] # distance from gps location to candidate route location
    initial_probs[seg_number] = ((2*math.pi)**0.5*sigma_z)**-1 * math.exp(-0.5*(dist_greatcircle/sigma_z)**2)

#for i in range(0,nSeg):
#    candidate_i = candidates[0,i]
#    dist_greatcircle = get_dist_lat_lon(emission_latlons[0],candidate_i)
#    #dist_greatcircle = dists_greatcircle[i,0]
#    if dist_greatcircle < 200: # Paper uses 200m
#        initial_probs[i] = ((2*math.pi)**0.5*sigma_z)**-1 * math.exp(-0.5*(dist_greatcircle/sigma_z)**2)

# Make the probabilities sum to 1.0. Decreases likelihood for roundoff errors that make small probabilities zero.

if sum(initial_probs) > 0.0:
    initial_probs = initial_probs/sum(initial_probs)
"""
sum_init_probs = sum(initial_probs.values())
if sum_init_probs > 0.0:
    initial_probs = {key: initial_probs[key]/sum_init_probs for key in initial_probs.keys()}
"""

# Populate the initial entries in the vitterbi arrays, which store info about the most likely paths.
for i in range(0,nSeg):
    # probability that the initial emission resulted from each state
    #   = probability of each initial state (assumed known), multiplied by probability of each state resulting in this emission?
    P[i,0] = initial_probs[i]**2
    back[i,0] = 0 # I think this is not used, just a placeholder. 
"""
for i in initial_probs.keys():
    P[0][i] = initial_probs[i]**2
    #P[i] = {emission_times[0]: initial_probs[i]**2}
    back[0][i] = 0 # I think this is not used, just a placeholder. 
"""
"""
DO EVERYTHING AT ONCE!
One loop through all time steps.
  Sub loop through all segments available at current time step.
    Sub loop through all segments available at next time step.
O(nT x nSeg x nSeg): Double nSeg -> quadruple processing time. Oops.
"""

# TIME LOOP
#A = {}
#for t in candidates:
for t in range(0,nT):
    if t % 20 == 0:
        print(t)

    if t < nT-1:
        dist_between_emissions = dists_gps_cum[emission_times_reduced_skipped[t+1]] \
                               - dists_gps_cum[emission_times_reduced_skipped[t]]
        #dist_between_emissions = get_dist_lat_lon(emission_latlons[t],emission_latlons[t+1])

        bearing_between_emissions = get_bearing_lat_lon(emission_latlons_reduced_skipped[t],
                                                        emission_latlons_reduced_skipped[t+1])

        # Initialize transition probability matrix AT THIS TIMESTEP, 
        # describing probability of going from any segment to any other:
        A_t = np.zeros((nSeg,nSeg)) 
        #A_t = {}

    candidate_list = candidates[t]

    # SEGMENT SUBLOOP TIME
    #for i in range(0,nSeg): 
    #for seg_number in range(0,nSeg): # (this approach would require catching errors where the seg_number is not in the candidate list)
    for seg_number in candidate_list:
        candidate = candidate_list[seg_number]

        # Needed? I think this is an artifact of a previous approach 
        dist_greatcircle = candidate['dist'] # distance from gps location to candidate route location
        #dist_greatcircle = get_dist_lat_lon(emission_latlons[t],candidates[t,i])

        # Emission probability matrix entry:
        # probability that a given observation (device GPS coord) came from a given state (official route segment)
        # Total matrix size = observation space (nT) x state space (nSeg)
        if t < nT-1:
            # SEGMENT SUBLOOP TIME
            #for next_seg_number in range(0,nSeg):
            for next_seg_number in candidates[t+1]:
                # ADS HERE: I changed candidates list and messed this up
                dist_route = get_route_dist(candidates,node_latlons,dists_cum,seg_number,next_seg_number,t,t+1)
                # 9 m/s is 20.1mph / 3 min/mile, which is as fast as I will ever move.
                speed_route = dist_route / (emission_times_reduced_skipped[t+1] - emission_times_reduced_skipped[t])
                if speed_route < 9.0:
                    dt = abs(dist_between_emissions - dist_route)
                    # filter for the difference between gps distance and route distance.
                    if dt <= 100: # paper uses 2000 m, but I am interested in 1- to 10-second data. # ADS HERE: RE-EVALUATE!!
                        # reduce transition probability if route and gps distances are different
                        p_dt = (1/beta)*math.exp(-dt/beta)

                        # reduce transition probability if the gps and route bearings are different
                        #latLngStart = candidate['latlon']
                        #latLngEnd = candidates[emission_times[t+1]][transition_seg_number]['latlon']
                        #bearing_between_candidates = get_bearing_lat_lon(latLngStart,latLngEnd)
                        #d_bearing = abs(bearing_between_emissions - bearing_between_candidates)
                        #if d_bearing > math.pi:
                        #    d_bearing = 2*math.pi - d_bearing
                        #p_bearing = math.exp(-d_bearing/(math.pi))
                        p_bearing = 1.0

                        A_t[seg_number,next_seg_number] = p_dt*p_bearing
                        """
                        if seg_number in A_t.keys():
                            A_t[seg_number][transition_seg_number] = p_dt*p_bearing
                        else:
                            A_t[seg_number]={transition_seg_number: p_dt*p_bearing}
                        """
        if t > 0:
            P[seg_number,t] = np.max([ P[prev_seg_number,t-1]*A_t_prev[prev_seg_number,seg_number]*emission_probs[t][seg_number] for prev_seg_number in range(0,nSeg) ])
            back[seg_number,t] = np.argmax([ P[prev_seg_number,t-1]*A_t_prev[prev_seg_number,seg_number] for prev_seg_number in range(0,nSeg) ])
            """
            #if emission_times[t] in P[seg_number].keys():
            if t in P.keys():
                if seg_number==234:
                    print(P[t].keys())
                P[t][seg_number] = np.max([ P[t-1][prev_seg_number]*A_t_prev[prev_seg_number][seg_number]*emission_probs[emission_times[t]][seg_number] for prev_seg_number in A_t_prev.keys() ])
                back[t][seg_number] = np.argmax([ P[t-1][prev_seg_number]*A_t_prev[prev_seg_number][seg_number] for prev_seg_number in A_t_prev.keys() ]) # is this a correct implementation?
            else:
                P[t]= {seg_number: np.max([ P[t-1][prev_seg_number]*A_t_prev[prev_seg_number][seg_number]*emission_probs[emission_times[t]][seg_number] for prev_seg_number in A_t_prev.keys() ])}
                back[t]= {seg_number: np.argmax([ P[t-1][prev_seg_number]*A_t_prev[prev_seg_number][seg_number] for prev_seg_number in A_t_prev.keys() ])} # is this a correct implementation?
            """              
    # ADS HERE
    # PLAYING WITH FIRE: Ensure each row of P sums to 1.0, reducing change of roundoff errors down the line.
    P[:,t] = P[:,t]/np.sum(P[:,t])
    """
    sum_P = sum(P[t].values())
    P[t] = { seg_number: P[t][seg_number]/sum_P for seg_number in P[t].keys() }
    """

    A_t_prev = A_t # store the transition matrix from this timestep for use in the next timestep.
   
    # Check if any paths forward exist from each segment at this timestep, and if so, normalize the path probabilities to sum to 1.0   
    for seg_num in range(0,nSeg):
        if np.sum(A_t_prev[seg_num])>0.0: # Check if a path forward exists from here, and if so, normalize it to sum to 1.0
            A_t_prev[seg_num] = A_t_prev[seg_num]/np.sum(A_t_prev[seg_num])
        #else:
        #    print("Zero sum A_t at time = "+str(emission_times[t])+", seg_num = "+str(seg_num)) # get rid? Sounding an alarm for nothing?
    """
    for seg_number in A_t_prev.keys():
        #A_t_prev_i = A_t_prev[i]
        sum_A_t_prev_i = sum(A_t_prev[seg_number].values())
        if sum_A_t_prev_i != 0.0:
            A_t_prev[seg_number] = {j: A_t_prev[seg_number][j]/sum_A_t_prev_i for j in A_t_prev[seg_number].keys()}
        else:
            print("Zero sum A_t at time = "+str(emission_times[t])) 
    """

"""
PAYOFF: FIGURE OUT WHAT THE MOST LIKELY PATH WAS
"""

# Find the state (segment) associated with the highest probability for the last step.
z_pred[nT-1] = np.argmax([ P[k,nT-1] for k in range(0,nSeg) ])
#z_pred[nT-1] = np.argmax([ P[nT-1][k] for k in P[nT-1].keys() ])

# Find the states (segments) associated with the most likely path leading to this final state.
for t in range(1,nT)[::-1]:
    z_pred[t-1] = back[z_pred[t],t]
    #z_pred[t-1] = back[t][z_pred[t]]

"""
PLOT THE ROUTE PATH, THE GPS PATH, AND THE CORRECTED GPS PATH.
"""
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

plt.plot(emission_latlons_unreduced[:,0],emission_latlons_unreduced[:,1],'b-',lw=1.0)
plt.plot(emission_latlons_reduced[:,0],emission_latlons_reduced[:,1],'g-',lw=2.0)
plt.plot(node_latlons[:,0],node_latlons[:,1],'k-')

max_dist_err = 0.0
for t in range(0,nT):
    seg_number_selected = z_pred[t]
    
    """
    route_latlon_i = candidates[t,seg_number_selected]
    plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'rx')
    """
 

    #for seg_id in range(0,nSeg):
    for seg_number in candidates[t]:
        route_latlon_i = candidates[t][seg_number]['latlon']
        if seg_number==seg_number_selected:
            dist_err = get_dist_lat_lon_m(emission_latlons_reduced_skipped[t],route_latlon_i)
            if dist_err > max_dist_err:
                max_dist_err = dist_err

            plt.plot([emission_latlons_reduced_skipped[t][0],route_latlon_i[0]],[emission_latlons_reduced_skipped[t][1],route_latlon_i[1]],'r-',lw=1.0)
            #plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'rx')
            #plt.plot(emission_latlons[t][0],emission_latlons[t][1],'rx')
            #plt.plot(route_latlon_i[0],route_latlon_i[1],'rx')
        #else:
            #plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'b--',lw=1.0)
            #plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'b--x')
            #plt.plot([emission_latlons[t][0],route_latlon_i[0]],[emission_latlons[t][1],route_latlon_i[1]],'bx')
ratio = (1.0/math.cos(math.radians(node_latlons[0][1])))
ax.set_aspect(ratio)
plt.show()

"""
PLOT HORIZ SPEED vs TIME, VERT SPEED vs TIME, ELEVATION vs DISTANCE
All using the route data. Final data should not use the gps points, except as comparison.
""" 
emission_dists = get_distance_from_start_m(emission_latlons_unreduced)
#plt.plot(emission_times,

# STATISTICS
dts = []
ratios = []
dist_errs = []
#for t in range(0,len(emission_latlons_reduced)-1):
for t in range(0,nT-1):
    dist_between_emissions = dists_gps_cum[emission_times_reduced_skipped[t+1]] \
                           - dists_gps_cum[emission_times_reduced_skipped[t]] 
    #dist_between_emissions = get_dist_lat_lon_m(emission_latlons_reduced[t],emission_latlons_reduced[t+1])
    dist_route = get_route_dist(candidates,node_latlons,dists_cum,z_pred[t],z_pred[t+1],t,t+1)
    dts.append(abs(dist_between_emissions - dist_route))
    ratios.append(abs(dist_between_emissions - dist_route)/dist_between_emissions)
    dts.append(abs(dist_between_emissions - dist_route))
    dist_errs.append(get_dist_lat_lon_m(emission_latlons_reduced_skipped[t],candidates[t][z_pred[t]]['latlon']))
print("Route differences:")
print("Maximum difference between two gps points' distance and the chosen route distance: "+str(max(dts))+" m")
print("for a ratio of: "+str(ratios[np.argmax(dts)]))
print("Maximum ratio of route distance difference to gps distance: "+str(max(ratios)))
print("Beta (approx): "+str(np.median(dts)/math.log(2)))

print("Location errors:")
print("Maximum distance between a gps point and its matching route point: "+str(max(dist_errs))+" m")
print("sigma_z (approx): "+str(1.4826*np.median(dist_errs))+" m")

