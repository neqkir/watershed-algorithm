##########################################################################################
##########################################################################################
####
####             WATERSHED Algorithm, peak detection
#### 
##########################################################################################
####
####             Label points as background or blob 
####            
##########################################################################################

## V6 adds some strategy to deal with coeluted zones

## based on Lindenberg Watershed algorithm

## Connected-component labeling preprocessing with Union Find data structure


 
import math
import itertools as it
import operator
from operator import itemgetter
import numpy as np
import random
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import pylab as pl
from pylab import *
import os, sys
 
import CC_labeling_8, CC_labeling_4, mathlib, peakDetect,filters,watershedVizu, geometry, io, shape_V2



####################################### DISJOINT SET

class disjointSet:

    def __init__(self, n):
        self.parent = [0]*n
        self.rank = [0]*n

    def makeSet(self,x):
        self.parent[x] = x
        self.rank[x] = 0

    def union(self,x,y):
        xRoot = self.find(x)
        yRoot = self.find(y)
        if(xRoot == yRoot):
            return

        if(self.rank[xRoot] < self.rank[yRoot]):
            self.parent[xRoot] = yRoot
            
        elif(self.rank[xRoot] > self.rank[yRoot]):
            self.parent[yRoot] = xRoot
            
        else:
            self.parent[yRoot] = xRoot
            self.rank[xRoot] += 1
            
    def find(self, x):
        if (self.parent[x] != x):
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    

#######################################  PRIORITY QUEUE
        
import heapq

class PriorityQueue(object):
    
    def __init__(self, key=lambda x: x[0]):
        self.l = []
        self.key = key

    def __len__(self):
        return len(self.l)

    def push(self, obj):
        heapq.heappush(self.l, (self.key(obj), obj[1]))

    def push_unique(self, obj):
        if(self.find(obj)==False):
            heapq.heappush(self.l, (self.key(obj), obj [1]))

    def pop(self):
        return heapq.heappop(self.l)[-1]

    def empty(self):
        if (self.l == []):
            return True
        return False

    def reset(self):
        resets = []
        while(self.l != []):
            resets.append(heapq.heappop(self.l)[-1])
        return resets            
        
    def find(self,item):
        #print ' self. l ' + str(self.l)
        for x in self.l:
            if (x[1] == item[1]):
                #print 'true ' 
                return True
        return False

        
###########################################################

class WATERSHED:

    matGC2D = [] # original matrix with abundance (elevation matrix)
    matRetT = [] # original retention time matrix
    GC2D = 0 # GC2D visualization object
    blobs = [] # list of blobs
    # a blob in blobs is a list of (i,j) points coordinates 
    blobs_labels = [] # matrix with labels 0: unlabeled ; i in 1:n_blobs: blob i ; -1: background  
    queueList = [] # a list, with queues: one blob one queue
    cclabeling = 0
    cclabeling4 = 0
    cells_regions = []
    volumes = []
    widths = []
    peaks = []
    
    def __init__(self,pathA=None,pathRetT=None, mod=None, imin = 0, imax=None):
        queue = []
        
        if(pathA==None and pathRetT==None and mod==None):
            return
        self.GC2D = see2D.GC2D(pathA, pathRetT, mod)
        self.modulate(imin, imax)
        self.blobs_labels = [] 
        self.queueList = []
        self.cells_regions = []
        self.volumes = []
        self.widths = []
        self.peaks = []

    def connectedComponentLabeling(self,mat):
        '''
        Turn matrix into an array with regions of connected pixels with same intensity
        '''
        self.cclabeling = CC_labeling_8.CC_lab(mat)
        self.cclabeling.connectedComponentLabel()
        
        self.cclabeling4 = CC_labeling_4.CC_lab(mat)
        self.cclabeling4.connectedComponentLabel()

        self.blobs_labels = []  

 
    def modulate(self,imin,imax):
        '''
        Modulate 1D data to 2D data between imin and imax
        '''
        self.GC2D.modulate(imin, imax)
        self.matGC2D = self.GC2D.matGC2D
        self.matRetT = self.GC2D.matRetT
        self.matGC2D = self.matGC2D[:-1] 
        self.matRetT = self.matRetT[:-1]



    def watershed(self, mat, mat_original = None, maxRegions = None, regions = True, refine = True):
            '''
            Returns list of peaks detected by WATERSHED algorithm for given elevation matrix mat
            . mat is matrix on which watershed should operate (some filtered, enhanced matrix)
            . mat_original is synthetic matrix before any filtering or enhancement
            . maxRegions is list of maxima regions, will be calculated if not given as an input
            . regions is true if input list maxRegions contains labels for regions, if not, pixels are first converted
            to regions labels
            '''
            self.connectedComponentLabeling(mat) # convert to regions
            self.blobs_labels = [0] * self.cclabeling.forest.utils
 
            N, M = len(mat), len(mat[0])

            if(maxRegions == None):
                #maxima = peakDetect.selectMarkers(N,M,mathlib.matrixToVector(mat))   # markers, sorted
                labelng, filtr, maxRegions = filters.maximum_filter_regions(mat, self.cclabeling)
                
            elif(regions == False): # this is pixel-wise watershed
                maxRegions = [self.cclabeling.labels[i][j] for i,j in maxRegions]

            self.queueList = [] # a list, with queues: one blob one queue
            
            for rLabel in maxRegions: # create a list of Priority queues, one for each maximum
   
                pq = PriorityQueue()
                pq.push_unique([(-1*self.cclabeling.labelToElevation(rLabel)),rLabel]) # initially, regions are unlabeled
                self.queueList.append(pq)  # priority in queue is opposite of elevation
                
            completed_blobs = [0]*len(maxRegions) 
            iter_ = 0
            
            while(True): # while some priority queue for a maximum has some elements

                if(self.emptyQueues() == True or min(completed_blobs) == 1):
                    break

                iter_ += 1             
                blob = 0
                
                for pq in self.queueList:

                    insert_queue_current_blob = True
                        
                    blob += 1  # each list is for another blob
 
                    if(completed_blobs[blob-1] != 1 and self.checkIfIdentified(blob, mat) == True): 
                        resets = pq.reset()
                        completed_blobs[blob-1] = 1

                        if(min(completed_blobs) == 1): # if last completed blob have still some unlabel in its pq, they are background

                            for r in resets:
                                self.blobs_labels[r] = -1
                                
                    if(min(completed_blobs) == 1):
                        break
                    
                    if(pq.empty()):
                        #print 'pq is empty ' 
                        continue  # for loop 
                     
                    item = pq.pop() # label for next region to consider, from current queue list

                    if(self.blobs_labels[item] != 0): # if labeled already
                        #print 'WARNING: Region labeled already '  
                        a=0
                        
                    else:
                        
                        higher_neighbors, unlabeled_neighbors = self.neighborhood(item, mat)  # retrieve information for neighbors

                        if (higher_neighbors == []):
                            self.blobs_labels[item] = blob
 
                        else:  # some higher level neighbors
           
                            first_neighbor_label = max([self.blobs_labels[higher_neighbors[nei]] for nei in range(len(higher_neighbors))]) # first higher neighbor labels (for blobs)

                            update_label = True # boolean, update current pixel label with current blob                            
                         
                            for rlabel_n in higher_neighbors:  # region label rlabel_n
 
                                blob_label = self.blobs_labels[rlabel_n] # blob label for this region (0,-1,1:N_blobs)
    
                                if(blob_label == -1):   # if one higher neighbor is background
                                    self.blobs_labels[item] = -1 # then it is background
                                    update_label = False
                                    break # no need to consider more neighbor for this labelling
                                    
                                elif(blob_label != first_neighbor_label and blob_label != 0):  # if two neighbors belong to two different blobs
                                    self.blobs_labels[item] = -1 # then it is background
                                    update_label = False
                                    break # no need to consider more neighbor for this labelling
 
                            if(update_label == True): # all higher neighbors are labeled and with same blob label
                                if(completed_blobs[first_neighbor_label-1] == 0 and self.isNeighborForBlob(first_neighbor_label, item, mat)):
                                    # and self.isNeighborForBlob(first_neighbor_label, item, mat)
                                    self.blobs_labels[item] = first_neighbor_label

                        for unlblob_label in unlabeled_neighbors:  # add all unlabeled point's neighbors into queue
                            pq.push_unique([(-1*self.cclabeling.labelToElevation(unlblob_label)),unlblob_label])

            labels_mat = self.matrixFromLabelsList(mat)

##            ## find coeluted blobs
##            for blob in range(len(maxRegions)+1)[1:]:
##                
##                isCoeluted, family = self.check_if_coeluted(blob, mat)

##                if(isCoeluted == True): # 1 if blob is coeluted (to be examined in second pass)
##                    coeluted_blobs[blob-1] = 1
##
##                coeluted_neighbors[blob-1] = family
##
##            zones = disjointSet(len(maxRegions)+1)
##            # make sets
##            for k in range(len(coeluted_blobs)):
##                if(coeluted_blobs[k] == 1):
##                    zones.makeSet(k+1)
##            
##            ## build coeluted zones
##            i = 0
##            seen = []
##            for cn in coeluted_neighbors:
##                 
##                if(cn != []):
##                    
##                    for blob in cn:
##
##                        zones.union(blob,i+1)
##
##                i += 1
## 
##            ## build self.cells_regions
##            labels_mat = self.matrixFromLabelsList(mat)
##
##            io.matToFile('gaussian2_V6_w.txt', labels_mat)
##                
##            watershedVizu.plot_heatmap(labels_mat)
##            
##            ## deal with coeluted zones
##            self.refineCoelutedZone_as(zones, coeluted_blobs, mat_original, labels_mat)
 
            return labels_mat
     

    def find_coeluted_zones(self, mat, mat_original, maxRegions, labels_mat = None):
        
        coeluted_blobs, coeluted_neighbors = [0]*len(maxRegions), [0]*len(maxRegions)

        if(labels_mat == None):
            labels_mat = self.matrixFromLabelsList(mat)
        
        for blob in range(len(maxRegions)+1)[1:]:
                
            isCoeluted, family = self.check_if_coeluted(blob, mat)
            print 'blob ' + str(blob) + ' is coeluted ' + str(isCoeluted) + ' , family ' + str(family)
                
            if(isCoeluted == True): # 1 if blob is coeluted (to be examined in second pass)
                coeluted_blobs[blob-1] = 1

            coeluted_neighbors[blob-1] = family

        zones = disjointSet(len( maxRegions ) + 1)
        
        # make sets
        for k in range(len( coeluted_blobs )):
            if( coeluted_blobs[k] == 1 ):
                zones.makeSet( k + 1 )
            
        ## build coeluted zones
        i = 0
        seen = []
        for cn in coeluted_neighbors:        
            if(cn != []):   
                for blob in cn:
                    zones.union(blob,i+1)
            i += 1
 
        
        return zones, coeluted_blobs, coeluted_neighbors

    
        
    def refine_clustering(self, zones, coeluted_blobs, coeluted_neighbors, mat, mat_original, maxRegions, labels_mat = None):
 
        ## build self.cells_regions
        if(labels_mat == None):
            labels_mat = self.matrixFromLabelsList(mat)

        io.matToFile('gaussian2_V6_w.txt', labels_mat)
                
        watershedVizu.plot_heatmap(labels_mat)
            
        ## deal with coeluted zones
        self.refineCoelutedZone_as(zones, coeluted_blobs, mat_original, labels_mat)
 
        return self.matrixFromLabelsList(mat)



    def emptyQueues(self, listQ = None):
        '''
        Return True if all queues are empty,
        Removes empty queues from list
        '''
        toRemove = []

        if(listQ == None):
            listQ = self.queueList
            
        noQueues = True
        for pq in listQ:
            if(pq.empty()):
                #print 'pq empty ' 
                toRemove.append(pq)
            else:
                #print 'pq n empty '
                noQueues = False
        
        return noQueues


    def isNeighborForBlob(self, label_color, region, mat):
        '''Return True if indeed region region is in the neighborhood of some region labeled as
        label_color
        '''
        if(self.blobs_labels == []):
            return
        #print 'isNeighborForBlob, is region in the neighborhood of some region colored as label_color '
        #print 'label_color ' + str(label_color)
        #print 'region ' + str(region)
         
        neighbors = self.cclabeling.neighborRegions(region,mat) # 4-connectivity neighbors for region

        for n in neighbors:
            if(self.blobs_labels[n] == label_color): # a neighbor is colored with label_color
                return True

        return False
                
        
    def checkIfIdentified(self, label, mat):
        '''
        Return True if zone connected component region with label label
        is surronded exclusively with background
        '''
        if(self.blobs_labels == []):
            return
        
        #print 'checkIfIdentified'
        iblob = 0 # region label in self.cclabeling
        count_blob, count_not_coeluted = 0, 0 # number of region labeled with this blob label at current stage
        #print 'self.blobs_labels ' + str(self.blobs_labels)
        
        isIdentified = True 
        
        for blob in self.blobs_labels:
            
            if(blob == label):
                count_blob += 1
 
                neighbors = self.cclabeling.neighborRegions(iblob,mat)
                
                for n in neighbors: # examine neighbors looking for unlabeled
                    cclabel = self.cclabeling.forest.find(n)[0] 
                    
                    if(self.blobs_labels[cclabel] == 0):
                        #print 'examining neighbors, return False because of this color ' + str(self.blobs_labels[cclabel])
                        return False # zone is not totally resolved
 
            iblob += 1
    
        if(count_blob == 0):
            print 'this blob did not begin ' 
            return False 
        
        else:
            print 'blob ' + str(label) + ' is complete ' 
            return True

        
    def check_if_coeluted(self, label, mat):
        '''
        Return True if blob with label label in [1:nbBlobs] is in some coelution zone
        Return list of blobs that are coeluted with blobs label and direct neighbors
        '''
        
        if(self.blobs_labels == []):
            return
        
        iblob = 0 # region label in self.cclabeling
        count_blob, count_not_coeluted = 0, 0 # number of region labeled with this blob label at current stage

        isCoeluted, coeluted = True, []
        
        for region in self.blobs_labels:
            
            if(region == label):
                count_blob += 1
                #print 'current region, belonging to current blob ' + str(region)
                #print 'index for region in cclabeling ' + str(iblob)
  
                neighbors = self.cclabeling.neighborRegions(iblob,mat)
                #print 'neighbors ' + str(neighbors)
                neighbors_in_background = True
                
                for n in neighbors: # examine neighbors looking for unlabeled
                    #print 'n ' + str(n)
                    cclabel = self.cclabeling.forest.find(n)[0] 

                    if(self.blobs_labels[cclabel] != label and self.blobs_labels[cclabel] > 0):
                       
                        coeluted.append(self.blobs_labels[cclabel]) # add this blob label to the blobs in this coelution zone

                        neighbors_in_background = False

                if(neighbors_in_background == True):

                    count_not_coeluted += 1
            
            iblob += 1

        #print 'count_blob ' + str(count_blob)

        #print 'count_not_coeluted ' + str(count_not_coeluted)
        
        if(count_blob == count_not_coeluted):
            isCoeluted = False

        #print 'blob ' + str(label) + ' is coeluted '
        
        return isCoeluted, list(np.unique(coeluted))


  
    def refineCoelutedZone_ch(self, zones, coeluted_blobs, mat, mat_labels ):
        '''
        Fill coeluted zone with appropriate blob colors
        Deal with coelution
        Version with points inside convex hull 
        '''
        
        pointsToHull = [[]] * (len(coeluted_blobs)+1) ## sets of points which will be to be convex-hulled
        parents_zones = [] # number of coeluted zones

        for k in range(len(coeluted_blobs)):

            if(coeluted_blobs[k] == 1):
                print 'here '
                parent_zone = zones.find(k + 1)
                parents_zones = parents_zones + [parent_zone]
                pointsToHull[parent_zone] = pointsToHull[parent_zone] + self.cells_regions[k + 1] ## add points in blob
                
        parents_zones = list(np.unique(parents_zones))

        
        print 'parent for zones ' +str(parents_zones) 
        
        # hull the points to hull in each zone
        for parent in parents_zones:

            print 'len for things in pointsToHull ' + str(len(pointsToHull[parent]))
            
            hull = geometry.convexHull(pointsToHull[parent])

            print 'hull terminated' 
             
            
            valley, full = self.findValley(hull, mat, mat_labels)

            pl.scatter([X[0] for X in pointsToHull[parent]], [X[1] for X in pointsToHull[parent]],s=50)
            pl.scatter([X[0] for X in hull], [X[1] for X in hull],c='r',s=50)
            pl.scatter([X[0] for X in full], [X[1] for X in full],c='k',s=50)
            pl.scatter([X[0] for X in valley], [X[1] for X in valley],c='g',s=70)
            pl.show()

    def refineCoelutedZone_as(self, zones, coeluted_blobs, mat, mat_labels ):
        '''
        Fill coeluted zone with appropriate blob colors
        Deal with coelution
        Version with points inside alpha shape
        '''
        
        pointsToHull = [[]] * (len(coeluted_blobs)+1) ## sets of points which will be to be convex-hulled
        parents_zones = [] # number of coeluted zones

        for k in range(len(coeluted_blobs)):

            if(coeluted_blobs[k] == 1): # curren blob belongs to some coelution zone
                parent_zone = zones.find(k + 1)
                parents_zones = parents_zones + [ parent_zone ]
                pointsToHull[ parent_zone ] = pointsToHull[ parent_zone ] + self.cells_regions[k + 1] ## add points in blob
                
        parents_zones = list(np.unique(parents_zones))

        print 'parent for zones ' +str(parents_zones) 
        
        # hull the points to hull in each zone
        for parent in parents_zones:

            print 'len for things in pointsToHull ' + str(len(pointsToHull[parent]))

##            print pointsToHull[parent]

            sh = shape_V2.alphaShape( pointsToHull[parent] )
            
            alpha_neighbors = sh.alpha_neighbors() # compute alpha neighbors

            # alpha_shape_hull = self.alpha_shape_hull(mat, alpha_neighbors)

            # print 'shape terminated' 
             
            # valley, full = self.findValley(hull, mat, mat_labels)

            # pl.scatter([X[0] for X in pointsToHull[parent]], [X[1] for X in pointsToHull[parent]],s=50)
            # pl.scatter([X[0] for X in hull], [X[1] for X in hull],c='r',s=50)
            # pl.scatter([X[0] for X in full], [X[1] for X in full],c='k',s=50)
            # pl.scatter([X[0] for X in valley], [X[1] for X in valley],c='g',s=70)
            # pl.show()



    def alpha_shape_hull(self, mat, alpha_neighbors):
        '''
            Minimum spanning alpha shape hull
            From points identified as alpha_neighbors for the alpha shape of coelution zone,
            return the complete hull for the alpha shape (the crust)
            This done by minimizing cumulative distance between points in alpha shape
        '''
        #print 'MINIMUM SPANNING VALLEY '
        
        d_set = disjointSet(len(alpha_neighbors))
        combinations = [X for X in it.combinations(alpha_neighbors, 2)]
        
        print 'combinations ' +str(combinations)

        id_in_ash_combinations = [idx for idx in it.combinations(range(len(alpha_neighbors)),2)]
        #print 'id_in_an_combinations ' + str(id_in_an_combinations)

        combinations_vertices = [ c.vertices() for c in combinations ]
        print 'combination vertices ' + str(combinations_vertices)

        edges = [ mathlib.euclidean_distance(X,Y) for X,Y in combinations_vertices ]

        #print 'edges ' + str(edges)

        # build forest with single-point-trees
        for k in range(len(alpha_neighbors)):
            
            d_set.makeSet(k)
        
        complete_as_hull = [X for X in alpha_neighbors]
   
        parents = [d_set.find(i) for i in range(len(alpha_neighbors))]
        
        
        # while forest not spanning and not all edges have been examined
        while(min(parents) != max(parents) and (min(edges) < sys.maxint)):
            
            #print 'while loop, spanning ' + str(spanning) + ' edges ' + str(edges) 

            idx,m = min(enumerate(edges), key = operator.itemgetter(1)) # find min weight edge
    
            #print 'idx ' + str(idx) + ' , m ' + str(m)
            
            id_in_ash = id_in_ash_combinations[idx] # indices in valley for vertices in this edge
            
            if(d_set.find(id_in_ash[0]) != d_set.find(id_in_ash[1])):

                line = geometry.buildLine(combinations[idx][0],combinations[idx][1])

                cpt_background = 0
                
                for lpx, lpy in line:

                    if(labels_mat[lpx][lpy] != -1): # line to bridge the gap is unfortunately crossing some not background point
                        
                        cpt_background += 1

                if(cpt_background < 2):    

                    complete_as_hull = complete_as_hull + line
                    
                    #print 'complete valley ' +str(complete_valley)

                    d_set.union(id_in_ash[0], id_in_ash[1])

            edges[idx] = sys.maxint # edge has been examined (either added or discarded)

            parents = [d_set.find(i) for i in range(len(alpha_neighbors))]
                    
        return complete_as_hull

    
            
    def findValley(self, hull_points, elevations_mat, labels_mat):
        '''
        Find valleys inside zone demarcated by convex hull.
        Valley is made of local minima along lines joining 
        '''
        
        if(len(hull_points) < 3):
            print 'WARNING: not enough points in hull'
            return

        valley = []
        
        full_hull = []

        print 'full hull ' + str(full_hull)
        
        # up-sample convex hull
 
        last = hull_points[0]
        
        for current in hull_points[1:]:

            print 'current in hull clock wise ' + str(current)
            
            full_hull = full_hull + geometry.buildLine(last, current)

            last = current

        full_hull = full_hull + geometry.buildLine(last, hull_points[0])

        full_hull = [list(X) for X in np.unique([(x,y) for x,y in full_hull])]

        print 'full hull terminated'
        
        # find valley

        for h_pt1, h_pt2 in it.combinations(hull_points,2): # evaluate min on each segment joining two points in full hull
            
            segment = geometry.buildLine(h_pt1, h_pt2)

            valley = valley + self.find_local_minima(segment, elevations_mat, labels_mat)

        valley = [list(X) for X in np.unique([(x,y) for x,y in valley])]
        print 'terminated valley ' + str(valley)
        
        valley = self.bridge_valley(valley, labels_mat)
            
        return valley, full_hull


    def bridge_valley(self, valley, labels_mat):
        '''
        Minimum spanning valley
        '''
        
        #print 'MINIMUM SPANNING VALLEY '
        
        d_set = disjointSet(len(valley))

        combinations = [X for X in it.combinations(valley,2)]
        
        #print 'combinations ' +str(combinations)

        id_in_valley_combinations = [idx for idx in it.combinations(range(len(valley)),2)]

        #print 'id_in_valley_combinations ' + str(id_in_valley_combinations)
        
        edges = [mathlib.euclidean_distance(X,Y) for X,Y in combinations]

        #print 'edges ' + str(edges)

        # build forest with single-point-trees
        for k in range(len(valley)):
            
            d_set.makeSet(k)
        
        complete_valley = [X for X in valley]
 
            
        parents = [d_set.find(i) for i in range(len(valley))]
        
        
        # while forest not spanning and not all edges have been examined
        while(min(parents) != max(parents) and (min(edges) < sys.maxint)):
            
            #print 'while loop, spanning ' + str(spanning) + ' edges ' + str(edges) 

            idx,m = min(enumerate(edges), key = operator.itemgetter(1)) # find min weight edge
    
            #print 'idx ' + str(idx) + ' , m ' + str(m)
            
            id_in_valley = id_in_valley_combinations[idx] # indices in valley for vertices in this edge
            
            if(d_set.find(id_in_valley[0]) != d_set.find(id_in_valley[1])):

                print 'find parent in tree, current point ' + str(id_in_valley[0]) + ' in valley ' + str(valley[id_in_valley[0]]) + ' parent is ' + str(d_set.find(id_in_valley[0]))
                
                print 'find parent in tree, current point ' + str(id_in_valley[1]) + ' in valley ' + str(valley[id_in_valley[1]]) + ' parent is ' + str(d_set.find(id_in_valley[1]))
                
                line = geometry.buildLine(combinations[idx][0],combinations[idx][1])

                cpt_background = 0
                
                for lpx, lpy in line:

                    if(labels_mat[lpx][lpy] != -1): # line to bridge the gap is unfortunately crossing some not background point
                        
                        cpt_background += 1

                if(cpt_background < 2):    

                    complete_valley = complete_valley + line
                    
                    #print 'complete valley ' +str(complete_valley)

                    d_set.union(id_in_valley[0],id_in_valley[1])

            edges[idx] = sys.maxint # edge has been examined (either added or discarded)

            parents = [d_set.find(i) for i in range(len(valley))]
                    

        return complete_valley




##    def bridge_valley(self, valley, labels_mat):
##
##        combinations = [X for X in it.combinations(valley,2)]
##        print 'combinations ' + str(combinations)
##        complete_valley = [X for X in valley]
##        
##        for pt_val in valley:
##            
##            bounded = False
##
##            if(combinations == []):
##                break
##            
##            candidate_combinations = [p for p in combinations if (p[0] == pt_val or p[1] == pt_val)]
##            distances = [mathlib.euclidean_distance(X,Y) for X,Y in  candidate_combinations]
##            print 'candidate_combinations ' + str(candidate_combinations)
##            if(candidate_combinations == []):
##                continue
##            
##            while(bounded == False and min(distances) < sys.maxint):
##                print 'distances ' + str(distances)
##                idx,m = min(enumerate(distances), key = operator.itemgetter(1))
##
##                line = geometry.buildLine(candidate_combinations[idx][0],candidate_combinations[idx][1])
##
##                cpt, cpt_background = 0,  0
##                
##                for lpx, lpy in line:
##
##                    if(labels_mat[lpx][lpy] != -1): # line to bridge the gap is unfortunately crossing some not background point
##                        cpt_background += 1
##                        if(cpt_background > 1):
##                            distances[idx] = sys.maxint
##
##                    else:
##                        cpt = cpt + 1
##
##                if((cpt + cpt_background) == len(line) and cpt_background <= 1 ):    
##
##                    complete_valley = complete_valley + line
##                    #print 'complete valley ' +str(complete_valley)
##                    combinations = [X for X in combinations if X != candidate_combinations[idx]]
##                    bounded = True
##
##        return complete_valley
            
            
        
    def find_local_minimum(self, segment_points, mat, labels_mat):
        '''
        Return global minimum in unlabeled zone on one segment made of points
        in segment_points. Minimum 'on segment', with 2-connectivity.
        '''
        if(len(segment_points) < 3):
            return []

        unlabeled = [X for X in segment_points if (labels_mat[X[0]][X[1]] == -1)]
 
        if(unlabeled == []):
            return []
        
        #print '[mat[X[0]][X[1]] for X in unlabeled] ' + str([mat[X[0]][X[1]] for X in unlabeled])

        idx, m = min(enumerate([mat[X[0]][X[1]] for X in unlabeled]), key=operator.itemgetter(1)) 
        print 'here is local min ' + str(idx ) + ' m ' + str(m)
        return [unlabeled[idx]]



    def find_local_minima(self, segment_points, mat, labels_mat):
        '''
        Return local minima in unlabeled zone on one segment made of points
        in segment_points, in sense of 8-connectivity
        '''
        if(len(segment_points) < 3):
            return []

        print 'find local minima ' 
        
        unlabeled = [(i,segment_points[i]) for i in range(len(segment_points)) if (labels_mat[segment_points[i][0]][segment_points[i][1]] == -1)]

        print 'unlabeled ' + str(unlabeled)
        
        if(len(unlabeled) < 3):
            return []

        minima = []

        N,M = len(mat), len(mat[0])
        u_is_minimum = False

        for k in range(len(unlabeled)-1):

            if(k == 0):
                continue
            
            last_i, last_j = unlabeled[k-1][1]
            last_idx = unlabeled[k-1][0]
            last_val = mat[last_i][last_j]

            u_i, u_j = unlabeled[k][1]
            next_i, next_j = unlabeled[k+1][1]

            u_idx = unlabeled[k][0]
            next_idx = unlabeled[k+1][0]
            
            u_val = mat[u_i][u_j]

            print 'u_val ' + str(u_val) + ', ui ' + str(u_i) + ', uj ' + str(u_j)
            
            next_val = mat[next_i][next_j]

            
            
            if(u_idx != last_idx + 1): # new zone

                print 'NEW ZONE  :  u_idx ' + str(u_idx) + ' next val ' + str(next_val) + ' u_val ' + str(u_val)
                
                if(next_val > u_val): # just look at next value

                    u_is_minimum = True
                    minima = minima + [(u_i, u_j)]

                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
                    

                else:

                    u_is_minimum = False


            elif((u_idx != next_idx - 1) or (u_is_minimum and last_val == u_val)): # new zone

                print 'END OF ZONE  :  u_idx ' + str(u_idx) + ' last_idx ' + str(last_idx) + ' u_val ' + str(u_val)
                
                if(last_val > u_val): # just look at next value

                    u_is_minimum = True
                    minima = minima + [(u_i, u_j)]

                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
                    

                else:

                    u_is_minimum = False

                    
            else: # any other zone

    
                if((last_val > u_val and next_val > u_val) or (u_is_minimum and last_val == u_val)):  

                    print 'SAME ZONE  :  next val ' + str(next_val) + ' u_val ' + str(u_val) + ' last val ' + str(last_val)
                    
                    u_is_minimum = True

                    print 'same zone added ' + str(u_i) + ', ' + str(u_j)
                    
                    minima = minima + [(u_i, u_j)]
                
                else:
                
                    u_is_minimum = False


        return minima



##    def canny_valley_detector(self, mat, ):
##        '''
##        Examine points in background, and returns points that are belonging to valley in coeluted zones
##        '''
##        if(len(segment_points) < 3):
##            return []
##
##        print 'find local minima ' 
##        
##        unlabeled = [(i,segment_points[i]) for i in range(len(segment_points)) if (labels_mat[segment_points[i][0]][segment_points[i][1]] == -1)]
##
##        print 'unlabeled ' + str(unlabeled)
##        
##        if(len(unlabeled) < 3):
##            return []
##
##        minima = []
##
##        N,M = len(mat), len(mat[0])
##        u_is_minimum = False
##
##        for k in range(len(unlabeled)-1):
##
##            if(k == 0):
##                continue
##            
##            last_i, last_j = unlabeled[k-1][1]
##            last_idx = unlabeled[k-1][0]
##            last_val = mat[last_i][last_j]
##
##            u_i, u_j = unlabeled[k][1]
##            next_i, next_j = unlabeled[k+1][1]
##
##            u_idx = unlabeled[k][0]
##            next_idx = unlabeled[k+1][0]
##            
##            u_val = mat[u_i][u_j]
##
##            print 'u_val ' + str(u_val) + ', ui ' + str(u_i) + ', uj ' + str(u_j)
##            
##            next_val = mat[next_i][next_j]
##
##            
##            
##            if(u_idx != last_idx + 1): # new zone
##
##                print 'NEW ZONE  :  u_idx ' + str(u_idx) + ' next val ' + str(next_val) + ' u_val ' + str(u_val)
##                
##                if(next_val > u_val): # just look at next value
##
##                    u_is_minimum = True
##                    minima = minima + [(u_i, u_j)]
##
##                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
##                    
##
##                else:
##
##                    u_is_minimum = False
##
##
##            elif((u_idx != next_idx - 1) or (u_is_minimum and last_val == u_val)): # new zone
##
##                print 'END OF ZONE  :  u_idx ' + str(u_idx) + ' last_idx ' + str(last_idx) + ' u_val ' + str(u_val)
##                
##                if(last_val > u_val): # just look at next value
##
##                    u_is_minimum = True
##                    minima = minima + [(u_i, u_j)]
##
##                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
##                    
##
##                else:
##
##                    u_is_minimum = False
##
##                    
##            else: # any other zone
##
##    
##                if((last_val > u_val and next_val > u_val) or (u_is_minimum and last_val == u_val)):  
##
##                    print 'SAME ZONE  :  next val ' + str(next_val) + ' u_val ' + str(u_val) + ' last val ' + str(last_val)
##                    
##                    u_is_minimum = True
##
##                    print 'same zone added ' + str(u_i) + ', ' + str(u_j)
##                    
##                    minima = minima + [(u_i, u_j)]
##                
##                else:
##                
##                    u_is_minimum = False
##
##
##        return minima

    
        

    def find_local_minima_hessian(self, segment_points, mat, labels_mat):
        '''
        Return local minima in unlabeled zone
        For each point in unlabeled zone find direction with steepest descent
        And find minimum along this direction 
        '''
        if(len(segment_points) < 3):
            return []

        print 'find local minima ' 
        
        unlabeled = [(i,segment_points[i]) for i in range(len(segment_points)) if (labels_mat[segment_points[i][0]][segment_points[i][1]] == -1)]

        print 'unlabeled ' + str(unlabeled)
        
        if(len(unlabeled) < 3):
            return []

        minima = []

        N,M = len(mat), len(mat[0])
        u_is_minimum = False

        for k in range(len(unlabeled)-1):

            if(k == 0):
                continue
            
            last_i, last_j = unlabeled[k-1][1]
            last_idx = unlabeled[k-1][0]
            last_val = mat[last_i][last_j]

            u_i, u_j = unlabeled[k][1]
            next_i, next_j = unlabeled[k+1][1]

            u_idx = unlabeled[k][0]
            next_idx = unlabeled[k+1][0]
            
            u_val = mat[u_i][u_j]

            print 'u_val ' + str(u_val) + ', ui ' + str(u_i) + ', uj ' + str(u_j)
            
            next_val = mat[next_i][next_j]
            
            if(u_idx != last_idx + 1): # new zone

                print 'NEW ZONE  :  u_idx ' + str(u_idx) + ' next val ' + str(next_val) + ' u_val ' + str(u_val)
                
                if(next_val > u_val): # just look at next value

                    u_is_minimum = True
                    minima = minima + [(u_i, u_j)]

                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
                    

                else:

                    u_is_minimum = False


            elif((u_idx != next_idx - 1) or (u_is_minimum and last_val == u_val)): # new zone

                print 'END OF ZONE  :  u_idx ' + str(u_idx) + ' last_idx ' + str(last_idx) + ' u_val ' + str(u_val)
                
                if(last_val > u_val): # just look at next value

                    u_is_minimum = True
                    minima = minima + [(u_i, u_j)]

                    print 'new zone added ' + str(u_i) + ', ' + str(u_j)
                    

                else:

                    u_is_minimum = False

                    
            else: # any other zone

    
                if((last_val > u_val and next_val > u_val) or (u_is_minimum and last_val == u_val)):  

                    print 'SAME ZONE  :  next val ' + str(next_val) + ' u_val ' + str(u_val) + ' last val ' + str(last_val)
                    
                    u_is_minimum = True

                    print 'same zone added ' + str(u_i) + ', ' + str(u_j)
                    
                    minima = minima + [(u_i, u_j)]
                
                else:
                
                    u_is_minimum = False


        return minima

    
    
            
    def neighborhood(self,x_label,mat=None):
        '''
        Return two lists of coordinates for points in neighborhood
        of region in mat corresponding in CC_labeling_8 to label x_label
            - a list with regions with elevation in mat more than this current region only
            - a list with unlabeled regions only
        Lists contains labels for regions
        '''
        if(self.cclabeling == 0):
            return

        if(self.blobs_labels == []):
            self.blobs_labels = [0]*self.cclabeling.forest.utils
            
        if(mat == None):
            mat = self.matGC2D
        
        labeled_higher_n, unlabeled_n = [], []
        
        i,j = self.cclabeling.forest.find(x_label)[1] # representative cell for x_label

        neighbors = self.cclabeling.neighborRegions(x_label,mat) # list of labels for neighboring regions
 
        for label_n in neighbors:

            if(label_n < 0 or label_n >= self.cclabeling.forest.utils):  # out of limits
                continue 

            nx, ny = self.cclabeling.forest.find(label_n)[1]  # representative for region labeled by label_n

            if(self.blobs_labels[label_n] == 0): # append only unlabeled points to unlabeled_n
                unlabeled_n.append(label_n)
     
            # append only higher labeled elevation points to labeled_higher_n
            if(mat[nx][ny] >= mat[i][j] and self.blobs_labels[label_n] != 0): 
                labeled_higher_n.append(label_n)
   
        return labeled_higher_n, unlabeled_n 


 
    def sortWithElevation(self, positions, mat):
        '''
        positions: array with positions in the vector corresponding to matrix mat
        mat: elevation matrix
        Returns positions sorted with respect to elevation
        '''
        N,M = len(mat), len(mat[0])
        
        elevation = [] # elevation for points in position
        
        for i in position:
            x,y = mathlib.mapToPosition(i,N,M)
            elevation.append(x,y)
            
        idx, sortList = mathlib.sortWithIndex(elevation)
        
        return [position[x] for x in idx]

 
    
    def matrixFromLabelsList(self, mat = None):
        ''' Fill in matrix with labels from list of blobs:
        1. blobs_labels contains labels by connected component regions 
        2. each cell in toPlot is watershed classification label for the
        region (i,j) belongs to 
        '''
        if(self.blobs_labels == []):
            return
        
        if(mat == None):
            mat = mat = self.matGC2D

        N,M = len(mat), len(mat[0])

        self.cells_regions = [[]] * max(np.unique(self.blobs_labels) + 1) # first element with cells in background
        # element i contains cells for ith blob label (watershed label)
        
        toPlot = []

        for i in range(N):
            
            tmp = []

            for j in range(M):

                wtsd_label = self.blobs_labels[self.cclabeling.labels[i][j]]
                
                tmp.append(wtsd_label) # watershed label for the region (i,j) belongs to

                if(wtsd_label > 0):
                    
                    self.cells_regions[wtsd_label] = self.cells_regions[wtsd_label] + [(i,j)] # add i,j to region with label wtsd_label
                    #print 'wtsd_label > 0, self.cells_regions ' + str(self.cells_regions)
                elif(wtsd_label == -1):

                    self.cells_regions[0] = self.cells_regions[0] + [(i,j)] # add i,j to region with label -1 (bakcgoround)
                    
            toPlot.append([x for x in tmp])

        #print 'self.cells_regions ' + str(self.cells_regions)
        
        return toPlot 
        
 

    def elevationMapFromLabels(self, mat, labels_mat):
        '''After processing, transform matrix with labels for blobs into elevation map with orignal elevation in mat
        for maximum in this region
        Return coordinates for maximum TODO
        '''
        elevationsMat = mathlib.sameDimZeroMatrix(labels_mat)

        N,M = len(labels_mat), len(labels_mat[0])

        labels_cclabeling = CC_labeling_8.CC_lab(labels_mat) # component labeling sur matrix with labels
    
        self.peaks = [0] * len(np.unique(self.blobs_labels)+1) # position and value for maximum cell in blob
        print 'np.unique(self.blobs_labels) ' + str(np.unique(self.blobs_labels))
        #print 'self.blobs_labels ' + str(self.blobs_labels)
        
        labels_cclabeling.connectedComponentLabel(N * M)


        for wtsh_label in np.unique(self.blobs_labels): # iter over watershed_labels label for regions 

            if(wtsh_label > 0):
                
                neighbors = self.cells_regions[wtsh_label]
                #print 'neighbors ' + str(neighbors) 
                max_idx, max_flood_fill = max(enumerate([mat[ni][nj] for ni,nj in neighbors]), key=operator.itemgetter(1)) 
                #print 'max_flood_fill '+ str(max_flood_fill)
                

                for i,j in neighbors:
                    
                    elevationsMat[i][j] = max_flood_fill

                print 'wtsh_label ' + str(wtsh_label)
                self.peaks[wtsh_label] = [max_flood_fill, neighbors[max_idx]]
                    

            elif(wtsh_label == -1):

                neighbors = self.cells_regions[0]

                for i,j in neighbors:
                    
                    elevationsMat[i][j] = -1

        return elevationsMat





    def integration(self, mat):
        '''
        Compute volumes for each blob region
        '''
        
        self.volumes = [0] * len(np.unique(self.blobs_labels)) # blobs volumes
        self.widths  = [0] * len(np.unique(self.blobs_labels)+1) # support widths

        surf_cell = 1
        
        for wtsh_label in np.unique(self.blobs_labels):
                
            if(wtsh_label > 0):

                neighbors = self.cells_regions[wtsh_label]

                for i,j in neighbors:
                    #print 'i,j for neighbors in integration '  + str(i) + ', ' + str(j)
                    self.volumes[wtsh_label-1] = self.volumes[wtsh_label-1] + mat[i][j] * surf_cell

                min_x, max_x = min([i for i,j in neighbors]), max([i for i,j in neighbors])
                min_y, max_y = min([j for i,j in neighbors]), max([j for i,j in neighbors])

                self.widths[wtsh_label] = (max_x - min_x + 1, max_y - min_y + 1)



                    
    def getMaxBlobSupport(self, mat):
        ''' Return rectangle support for maximum support blob
        '''
        if(self.widths == []):
            return

        max_support_surf = - 1
        idx = -1
        
        for k in range(len(self.widths)):
            
            if(self.widths[k] == 0):
                continue

            else:
                x,y = self.widths[k]
                if(x*y > max_support_surf):
                    max_support_surf = x*y
                    max_support = self.widths[k]
                    idx = k

        return max_support, idx




