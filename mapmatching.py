from utils.hikehelper import get_dist_lat_lon_m, get_distance_from_start
import numpy as np
import math

def find_closest_point(emission_latlon,node_latlons):
    node_lons = [node_latlons[0][0], node_latlons[1][0]]
    node_lats = [node_latlons[0][1], node_latlons[1][1]]
    ye = emission_latlon[1]
    #xe = emission_latlon[0]*math.cos(math.radians(0.5*(node_lons[0]+node_lons[1])))
    xe = emission_latlon[0]

    if node_lons[0] == node_lons[1]: # Special case: vertical line, or no movement.
        index_1 = node_lats.index(min(node_lats)) # the segment is perfectly N-S
        index_2 = 1 - index_1
        x = node_lons[index_1]
        y1 = node_lats[index_1]
        y2 = node_lats[index_2]
        if ye <= y1:
            return [x,y1]
        elif ye >= y2:
            return [x,y2]
        else:
            return [x,ye]
    else: # find point along this horizontalish line
        index_1 = node_lons.index(min(node_lons)) # the segment is not perfectly N-S
        index_2 = 1 - index_1
        x1 = node_lons[index_1]
        x2 = node_lons[index_2]
        y1 = node_lats[index_1]
        y2 = node_lats[index_2]
        #x1 = lon1*math.cos(math.radians(0.5*(y1+y2)))
        #x2 = lon2*math.cos(math.radians(0.5*(y1+y2)))

        m = (y2 - y1) / (x2 - x1)  # slope of segment
        if m != 0.0:
            m_perp = -1/m              # slope of perpendicular line
            b = y1 - m*x1
            b_perp = ye - m_perp*xe

            x_perp = (b - b_perp) / (m_perp - m)
            y_perp = m*x_perp + b
        else: # completely horizontal line
            x_perp = xe
            y_perp = y1

        #lon_perp = x_perp/math.cos(math.radians(0.5*(y1+y2)))

        if x_perp <= x1:
            return [x1, y1]
        elif x_perp >= x2:
            return [x2, y2]
        else:
            return [x_perp, y_perp]

def get_route_dist(candidates,node_latlons,dists_cum,seg_number_start,seg_number_end,time_start,time_end):
    candidate_ti = candidates[time_start][seg_number_start]['latlon']
    candidate_tj = candidates[time_end][seg_number_end]['latlon']
    dist_along_ti = get_dist_lat_lon_m(candidate_ti,node_latlons[seg_number_start]) # distance from segment start node to candidate point
    dist_along_tj = get_dist_lat_lon_m(candidate_tj,node_latlons[seg_number_end])
    return abs(dists_cum[seg_number_end] + dist_along_tj - dists_cum[seg_number_start] - dist_along_ti)

