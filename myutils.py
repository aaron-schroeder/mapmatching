import numpy as np
import math
from geopy.point import Point
from geopy.distance import great_circle as distance

# Maybe this should become a method for segments
def find_closest_point(emission_point,segment):
    lat_emission = emission_point.latitude
    lon_emission = emission_point.longitude
    node_start = segment.node_start
    node_end = segment.node_end
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
        # find point along this horizontal or sloped line
        node_start_is_first = (0==np.argmin([node_start.longitude, node_end.longitude]))
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

def get_route_dist(candidates,segments,dists_cum,seg_number_start,seg_number_end,time_start,time_end):
    candidate_ti = candidates[time_start][seg_number_start]
    candidate_tj = candidates[time_end][seg_number_end]
    # distance from segment start node to candidate point (which lies along segment)
    dist_along_ti = distance(
        candidate_ti.point,
        segments[seg_number_start].node_start).meters
    dist_along_tj = distance(
        candidate_tj.point,
        segments[seg_number_end].node_start).meters
    return abs(dists_cum[seg_number_end] + dist_along_tj - dists_cum[seg_number_start] - dist_along_ti)
