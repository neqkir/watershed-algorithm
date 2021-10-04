##########################################################################################
##########################################################################################
####
####             Connected-Component labeling Two-Pass algorithm  (4-connectivity)
#### 
##########################################################################################
####
####             Cluster elevation matrix - iso-elevation regions 
####            
##########################################################################################

import mathlib
import numpy as np
## based on http://en.wikipedia.org/wiki/Connected-component_labeling

## based on
## Connected Component Labeling Algorithms for Gray-Scale Images 
## and Evaluation of Performance using Digital Mammograms, R. Dharshana Yapa
## and  K. Harada

#######################################################################################################

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

class disjointRegions(disjointSet):

    def __init__(self, n):
        self.parent = [0]*n
        self.rank = [0]*n
        self.cell = [0]*n # top left hand corner cellfor regions
        self.utils = 0 # effective number of new created label 
        self.neighbors = [[]]*n

    def makeSet(self,x,i,j):
        self.parent[x] = x
        self.rank[x] = 0
        self.cell[x] = (i,j)
        self.utils += 1

    def union(self,x,y):
        xRoot,xCell = self.find(x)
        yRoot,yCell = self.find(y)
        if(xRoot == yRoot):
            return

        if(self.rank[xRoot] < self.rank[yRoot]):
            self.parent[xRoot]  = yRoot 
        elif(self.rank[xRoot] > self.rank[yRoot]):
            self.parent[yRoot]  = xRoot 
        else:
            self.parent[yRoot] = xRoot
            self.rank[xRoot] += 1

    def find(self, x):
        if (self.parent[x] != x):
            self.parent[x] = self.find(self.parent[x])[0]

        return self.parent[x], self.cell[self.parent[x]]  # return representative for region: label, cell

    def topDown(self,x):
        '''
        Return all labels whose parent is x's parent
        '''
        tDlabel  = [] 
        for i in range(self.utils):
            if(self.find(i)[0] == self.find(x)[0]):
                tDlabel.append(i)
                 
        return tDlabel

##########################################################################################################

class CC_lab:

    def __init__(self,mat):
        self.labels = []
        self.forest = 0
        self.mat = mat
            
    def connectedComponentLabel(self,n_clusters=0):
            '''
            Label regions belonging to the same color level
            Algorithm uses disjoint-set data structure (equivalence table) and
            it is based on
            A Run-Based Two-Scan Labeling Algorithm, L. H, Y. Chao, K. Suzuki, IEEE
            '''
            N, M = len(self.mat), len(self.mat[0])
                 
            if(n_clusters==0):
                n_clusters = N * M
           
            
            self.labels = sameDimZeroMatrix(self.mat) # build label matrix with appropriate dimensions
            label = 0 # next label to be assigned

            self.forest = disjointRegions(n_clusters) # forest to record equivalences
            
            for i in range(N):
                
                for j in range(M):

                    neighbors = self.connectedNeighbors(i,j) # neighbors with same value
                    
                    #print 'neighbors ' + str(neighbors) +' ' + str(i) + ' '+ str(j)
                    if (neighbors == [[],[]]): # no neighbors at all
                         
                        self.labels[i][j] = label # new label
                        self.forest.makeSet(label,i,j)
                        label += 1
                        
                        if(i > 0):
                                self.forest.neighbors[self.labels[i][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j]]] + [self.labels[i-1][j]])
                                self.forest.neighbors[self.labels[i-1][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i-1][j]]] + [self.labels[i][j]])

                        if(j > 0):
                                self.forest.neighbors[self.labels[i][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j]]] + [self.labels[i][j-1]])
                                self.forest.neighbors[self.labels[i][j-1]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j-1]]] + [self.labels[i][j]])
                                
                    else:
                        ##find minimum neighbor
                        if(neighbors[0] != [] and neighbors[1] == []): # some west neighbor only

                            x,y = neighbors[0][0][0],neighbors[0][0][1]
                            self.labels[i][j] = self.labels[x][y]

                            if(i > 0):
                                self.forest.neighbors[self.labels[i][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j]]] + [self.labels[i-1][j]])
                                self.forest.neighbors[self.labels[i-1][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i-1][j]]] + [self.labels[i][j]])
                                #print ' added to region vicinity ' + str(self.labels[i-1][j]) + ' for label ' + str(self.labels[i][j])
                                #print 'self.forest.neighbors ' + str(self.forest.neighbors)

                        elif(neighbors[1] != [] and neighbors[0] == []): # some north neighbor only
                            #print 'adopted north' 
                            x,y = neighbors[1][0][0],neighbors[1][0][1]
                            self.labels[i][j] = self.labels[x][y]

                            if(j > 0):
                                self.forest.neighbors[self.labels[i][j]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j]]] + [self.labels[i][j-1]])
                                self.forest.neighbors[self.labels[i][j-1]] = np.unique([x for x in self.forest.neighbors[self.labels[i][j-1]]] + [self.labels[i][j]])
                                #print ' added to region vicinity ' + str(self.labels[i][j-1]) + ' for label ' + str(self.labels[i][j])
                                #print 'self.forest.neighbors ' + str(self.forest.neighbors)

                        else: # some west and north neighbors
                            
                            x1,y1,x2,y2 = neighbors[0][0][0],neighbors[0][0][1],neighbors[1][0][0],neighbors[1][0][1]
                            l1,l2 = self.labels[x1][y1],self.labels[x2][y2]
                            self.labels[i][j] = min(l1,l2)
                            self.forest.union(l1, l2)  # union labels, maybe they are different
 
            # second pass
            for i in range(N):

                for j in range(M):

                    self.labels[i][j] = self.forest.find(self.labels[i][j])[0]

                  
    def neighborRegions(self,labelij,mat):
        '''
        Return neighbor regions for the region containing cell (i,j) with label labelij
        '''
        neighbors = []
        R,C = len(mat), len(mat[0])
 
        parentij = self.forest.find(labelij)[0] # representative for the region which the cell i,j belongs to
        areas = self.forest.topDown(parentij)  # list of areas in the region which the cell i,j belongs to
        #print 'areas ' + str(areas)
        M = max(areas)
         
        #print 'representative cells ' + str([self.forest.cell[s] for s in areas])
        ## look at neighbors for top left hand corners of areas in region containing (i,j)
        for k,l in [self.forest.cell[s] for s in areas]: # iter over top left hand corner cells for areas in region
            
            if(l > 0 and self.labels[k][l-1] not in areas):
                neighbors.append(self.forest.find(self.labels[k][l-1])[0]) # look north
                #print 'added north to neighbors ' + str(self.forest.find(self.labels[k][l-1])[0])
                
            if(k > 0 and self.labels[k-1][l] not in areas):
                neighbors.append(self.forest.find(self.labels[k-1][l])[0]) # look west
                #print 'added west to neighbors ' + str(self.forest.find(self.labels[k-1][l])[0])
                 

        ## look in neighbors areas for recorded neighbor regions
        for area in areas:
            neighboor_n_w = self.forest.neighbors[area]
            for k in range(len(neighboor_n_w)):
                neighbors.append(self.forest.find(neighboor_n_w[k])[0])

             
        return np.unique(neighbors)


    def labelToElevation(self, label_i):
        '''
        From label, retrieve elevation in matrix self.mat
        '''
        x,y = self.forest.cell[label_i]
        return self.mat[x][y]

        
        
    def connectedNeighbors(self,i, j): 
        '''
        Return coordinates for neighbors of pixel i, j that have same value
        as pixel i, j
        '''
        neighbors = []
        neighbors.append([])
        neighbors.append([])

        N, M = len(self.mat), len(self.mat[0])
        
        if(i >= N or j >= M or i < 0 or j < 0): # exceed dimensions
            return []

        val = self.mat[i][j]
        
        if(i == 0 and j == 0): ## top left-hand corner, no neighbors
            return neighbors

         
        # i or j is not zero
        if(j > 0):
            
            if(val == self.mat[i][j-1]):
                neighbors[0].append((i, j-1)) ## west
            
                
        if(i > 0):
            
            if(val == self.mat[i-1][j]):
                neighbors[1].append((i-1, j)) ## north

        return neighbors


     
def sameDimZeroMatrix(mat):
    '''
    Create matrix with same dimensions as mat with zeros
    '''
    zeros = []
    
    for i in range(len(mat)):
        tmp = []
        for j in range(len(mat[i])):
            tmp.append(0)
        zeros.append([x for x in tmp])
        
    return zeros
