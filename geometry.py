##########################################################################################
##########################################################################################
####
####             Convex hull, with qhull, quick hull algorithm
#### 
##########################################################################################

 
import math
import itertools as it
import operator
import numpy as np
import random
import pylab as pl

import sys

import mathlib

## implementation of quickhull algorithm based on
## The Quickhull algorithm for Convex Hulls

## implementation of Bresenham algorithm 'stay as close as possible to the straight line'

def convexHull(points_set):
    '''
    Returns the convex hull for a set of points points_set
    using quickhull algoritm
    '''
    convex_hull = [] # here you store the hull

    ## min and max (abscissa) points in set
    left_pt = points_set[min(enumerate(points_set),key = operator.itemgetter(1))[0]]
    right_pt = points_set[max(enumerate(points_set),key = operator.itemgetter(1))[0]]

    convex_hull.append(left_pt)
    convex_hull.append(right_pt)
    
    above_set, below_set = pointsPartition(points_set, left_pt, right_pt)
    
    if(above_set != []):

        convex_hull = findHullAbove(above_set, left_pt, right_pt, convex_hull)

    if(below_set != []):
        
        convex_hull = findHullBelow(below_set, left_pt, right_pt, convex_hull)
        
    return clockWiseHull(convex_hull)


def pointsPartition(points_set, left_pt, right_pt):
    '''
    Return set of all points that are above line joining left_pt and right_pt and
    set of all points that are below line joining left_pt and right_pt 
    '''
    above_set, below_set = [],[]
    #print ' left_pt ' + str(left_pt)+ ' right_pt ' + str(right_pt)
    
    if(left_pt[0] == right_pt[0]):
        #print 'partition gives nothing '
        return above_set, below_set

    above_set, below_set = [],[]
    
    for p in points_set:
        
        x_p, y_p = p[0], p[1]
        x_bound = x_p
        y_bound = (right_pt[1] - left_pt[1])/float(right_pt[0] - left_pt[0]) * (x_bound - right_pt[0]) + right_pt[1]

        if(y_p >= y_bound and ((x_p, y_p) not in [left_pt, right_pt])):
            above_set.append(p)
            
        elif(y_p < y_bound and ((x_p, y_p) not in [left_pt, right_pt])):
            below_set.append(p)

    return above_set, below_set

 
def findHullAbove(part_set, left_pt, right_pt, convex_hull):
    '''
    Return points on convex hull from set of points part_set that are on the right side
    of oriented line from left_pt to right_pt
    '''
    if(len(part_set) == 0):
        return convex_hull
 
    left_pt = [(x,y) for x,y in [left_pt]][0]
    right_pt = [(x,y) for x,y in [right_pt]][0]
    
    if(left_pt[0] == right_pt[0]):

        convex_hull.append(left_pt)
        convex_hull.append(right_pt)
        
        convex_hull = [list(X) for X in np.unique([(x,y) for x,y in convex_hull])]

    else:
         
        furthest = [(x,y) for x,y in [ findFurthestPoint(part_set, left_pt, right_pt, 'a')] ][0]
        
        convex_hull.append(furthest)

        convex_hull = [list(X) for X in np.unique([(x,y) for x,y in convex_hull])]
        
        if(left_pt[1] == right_pt[1]):
            
            above_left_set = [p for p in part_set if (p[1] > left_pt[1])]
            below_left_set = []
 
        else:
            
            above_left_set, below_left_set = pointsPartition(part_set, left_pt, furthest)

        if(above_left_set != []):
           
            above_left_set = [list(X) for X in np.unique([(x,y) for x,y in above_left_set])]
            
            convex_hull = findHullAbove(above_left_set, left_pt, furthest, convex_hull)


        # find points above right 
        above_right_set, below_right_set = pointsPartition(part_set, furthest, right_pt)

        if(above_right_set != []):
            
            above_right_set = [list(X) for X in np.unique([(x,y) for x,y in above_right_set])]
            convex_hull = findHullAbove(above_right_set, furthest, right_pt, convex_hull)

    return convex_hull

        
def findHullBelow(part_set, left_pt, right_pt, convex_hull):
    '''
    Return points on convex hull from set of points part_set that are on the right side
    of oriented line from left_pt to right_pt
    '''
    if(len(part_set) == 0):
        return convex_hull

    left_pt = [(x,y) for x,y in [left_pt]][0]
    right_pt = [(x,y) for x,y in [right_pt]][0]

    if(left_pt[0] == right_pt[0]):
    
        convex_hull.append(left_pt)
        convex_hull.append(right_pt)

        convex_hull = [list(X) for X in np.unique([(x,y) for x,y in convex_hull])]

        if(left_pt[1] == right_pt[1]):
            
            above_left_set = []
            below_left_set = [p for p in part_set if (p[1] < left_pt[1])]
 
        else:
            
            above_left_set, below_left_set = pointsPartition(part_set, left_pt, furthest)

    else:
        
        furthest = [(x,y) for x,y in [findFurthestPoint(part_set, left_pt, right_pt,'b')]][0]

        convex_hull.append(furthest)

        convex_hull = [list(X) for X in np.unique([(x,y) for x,y in convex_hull])]

        # find points above left
        above_left_set, below_left_set = pointsPartition(part_set, left_pt, furthest)
        
        if(below_left_set != []):
            
            below_left_set = [list(X) for X in np.unique([(x,y) for x,y in below_left_set])]
            
            convex_hull = findHullBelow(below_left_set, left_pt, furthest, convex_hull)

        # find points above right
        above_right_set, below_right_set = pointsPartition(part_set, furthest, right_pt)

        if(below_right_set != []):
            
            below_right_set = [list(X) for X in np.unique([(x,y) for x,y in below_right_set])]
            convex_hull = findHullBelow(below_right_set, furthest, right_pt, convex_hull)
    
    return convex_hull




def clockWiseHull(hull_points):
    '''
    Return points in hull_points, but ranged in
    clockWise logic
    '''
    leftmost, rightmost = hull_points[min(enumerate(hull_points),key = operator.itemgetter(1))[0]], hull_points[max(enumerate(hull_points),key = operator.itemgetter(1))[0]]

    # partition hull points, above and below
    above, below = pointsPartition(hull_points, leftmost, rightmost)

     
    sorted_above = [above[idx] for idx in [i[0] for i in sorted(enumerate(above),key=operator.itemgetter(1))]]

    sorted_below = [below[idx] for idx in [i[0] for i in sorted(enumerate(below), key=operator.itemgetter(1))][::-1]]
    
    hull = [leftmost] + sorted_above + sorted_below
    
    return hull



def full_hull(hull_points):
    '''
    From the output of quick hull algorithm, complete to a full hull
    '''
    full_hull =[]
    
    last = hull_points[0]
        
    for current in hull_points[1:]:

        full_hull = full_hull + buildLine(last, current)

        last = current

    full_hull = full_hull + buildLine(last, hull_points[0])

    full_hull = [list(X) for X in np.unique([(x,y) for x,y in full_hull])]

    return full_hull



def findFurthestPoint(part_set, left_pt, right_pt, where):
    '''
    Return point in part_set furthest from left_pt and right_pt
    '''
    distances = []

    max_d = 0
    
    if(where == 'a'): # above
        max_p = (sys.maxint, -1 * sys.maxint)
    elif(where == 'b'): # below
        max_p = (- 1 * sys.maxint, sys.maxint)
    else:
        return
    
    for p in part_set:

        d = mathlib.euclidean_distance(left_pt,p) + mathlib.euclidean_distance(right_pt,p)

        if(d > max_d):
            max_d = d
            max_p = p
            
        elif(d == max_d):
            
            if(where == 'a' and (p[0] < max_p[0] or p[1] > max_p[1])):
                max_d = d
                max_p = p
                
            elif(where == 'b' and (p[0] > max_p[0] or p[1] < max_p[1])):
                max_d = d
                max_p = p
            
    return max_p

    
def closest(point, lst):

    lst = [(i,j) for i,j in [X for X in lst]]
    
    return min(enumerate([mathlib.euclidean_distance(X,point) for X in lst]),key = operator.itemgetter(1))[0]

    
def buildLine(left_pt, right_pt):
    '''
    Return set of points in segment from left_pt to right_pt
    so that points in set are as close as possible to segment
    Computed with Bresenham algorithm
    '''
    points = []

    if(left_pt[0] == right_pt[0]):
        
        for y in range(min(right_pt[1],left_pt[1]),max(right_pt[1],left_pt[1])):

            points = points + [(left_pt[0], y)]
            
    else:

        slope = (right_pt[1] - left_pt[1])/float(right_pt[0] - left_pt[0])
        
        for x in range(min(right_pt[0],left_pt[0]),max(right_pt[0],left_pt[0])):

            y = int(math.floor(slope * (x - left_pt[0]) + left_pt[1]))

            points = points + [(x,y)]

    points = points + [right_pt]
    
    return points



 
