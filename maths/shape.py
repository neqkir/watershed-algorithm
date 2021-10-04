##########################################################################################
##########################################################################################
####
####             Geometry: alpha shape in the plane
#### 
##########################################################################################
## implementation of Watson Bowyer algorithm for Delaunay triangulation
## implementation of Voronoi diagram algorithm
## implementation of Edelsbrunner algorithm for alpha-shapes computation

import math
import itertools as it
import operator
import numpy as np
import random
import pylab as pl
import scipy
import os, sys

from scipy import stats
import mathlib, CC_labeling_8, filters, geometry

##########################################################################################

## Half-edge data structure
## based on http://www.flipcode.com/archives/The_Half-Edge_Data_Structure.shtml

class HE_edge:

    pair = []
    face = []
    next_edge = []
    end_vert = []

    def __init__(self, vert_, pair_ = None, face_ = None, next_ = None):

        self.end_vert = vert_ # vertex at the end of the half-edge
        self.pair = pair_ # oppositely oriented half edge
        self.face = face_ # face the half-edge borders
        self.next_edge = next_ # next-edge around the face in anticlockwise order
        
        if(self.end_vert.edge == None):
            self.end_vert.edge = self.pair
            

    def vertices(self):
        ''' query vertices for this edge '''
        vert1 = self.end_vert # vertex at the end of the half-edge 
        vert2 = self.pair.end_vert # vertex at the start of the half-edge

        return vert1, vert2


    def faces(self):
        ''' query faces for this edge'''
        face1 = self.face # face the half-edge borders
        face2 = self.pair.face # face the pair half-edge borders

        return face1, face2


    def edges_around_face(self):
        ''' query all edges around the face that this half-edge borders'''
        
        if(self.face == None):
            
            #print 'edges around face is None' + str(self.end_vert.x) + ', ' + str(self.end_vert.y)
            return []
        
        edge = self.face.edge
        
        return [edge, edge.next_edge, edge.next_edge.next_edge]



    def edges_around_vertex(self):
        ''' query all edges around the vertex that this half-edge has as an end'''
        edge = self.end_vert.edge

        edges = [edge]
        contnd = True

        while(contnd):
            edge = edge.next_edge
            edges = edges + [edge]
            contnd = (self.end_vert.edge != edge)

        return edges 


    def faces_around_vertex(self):
        ''' query all faces around the vertex that this half-edge has as an end'''
        edge = self.end_vert.edge
        face = edge.face
        
        faces = [face]
        contnd = True

        while(contnd):
            edge = edge.next_edge
            face = edge.face
            faces = faces + [face]
            contnd = (self.end_vert.edge != edge)

        return faces


    def compare(self, e):
        '''
        Compare half-edge (same side)
        '''
        
        if(self.end_vert == [] or self.pair.end_vert == [] or e.end_vert == [] or e.pair.end_vert == []):        
            return False
            
        if(self.end_vert.x == e.end_vert.x and self.end_vert.y == e.end_vert.y and self.pair.end_vert.x == e.pair.end_vert.x and self.pair.end_vert.y == e.pair.end_vert.y):
            #print 'compare true ' + str(e.end_vert.x) + ', ' + str(e.end_vert.y)
            return True

        return False

    def half_compare(self,e):

        if(self.end_vert == []  or e.end_vert == []):        
            return False

        if(self.end_vert.x == e.end_vert.x and self.end_vert.y == e.end_vert.y and self.pair.end_vert.x == e.pair.end_vert.x and self.pair.end_vert.y == e.pair.end_vert.y):
            return True

        return False
        
        
class HE_vert:

    edge = []
    x, y = 0, 0

    def __init__(self, x, y, edge_ = None):

        self.edge = edge_ # edge which has this vertex as a start
        self.x = x
        self.y = y


class HE_face:

    edge = []
    ccc, ccr = 0, 0 # circumcenter, circumradius
    
    def __init__(self, edge_ ):

        self.edge = edge_

    def get_vertices_around_face(self):
        '''query all edges around the face that this half-edge borders'''
        
        edges = self.edge.edges_around_face() # edges around this face (edges around the face bordered by self.edge)
        #print 'edges around face ' + str(edges)
        vertices = []

        if(edges == []):
            return []
        
        for e in edges:
            vertices = vertices + [e.end_vert]

        return vertices


    def getCircumCircle(self):
        ''' compute center and radius for circumcircle '''
        
        edges_around = self.edge.edges_around_face()
   
        if(edges_around == []):
            #print 'in get circum got no edges' 
            return 0, 0
        
        ax, ay, bx, by, cx, cy = edges_around[0].end_vert.x, edges_around[0].end_vert.y, edges_around[1].end_vert.x, edges_around[1].end_vert.y, edges_around[2].end_vert.x, edges_around[2].end_vert.y
        #print 'ax, ay, bx, by, cx, cy ' + str(ax) + ', ' + str(ay) + ', ' + str(bx) + ', ' + str(by) + ', ' + str(cx) + ', ' + str(cy)  
        Dval = 2.0 * (ax*by - ax*cy - bx*ay + bx*cy + cx*ay - cx*by)
        #print Dval
        cccx = (ax**2 + ay**2) * (by - cy) - (bx**2 + by**2) * (ay - cy) + (cx**2 + cy**2) * (ay - by)
        cccy = -1 * (ax**2 + ay**2) * (bx - cx) + (bx**2 + by**2) * (ax - cx) - (cx**2 + cy**2) * (ax - bx)
            
        cccx, cccy = cccx / float(Dval), cccy / float(Dval)

        self.ccc = HE_vert(cccx,cccy)
        self.ccr = math.sqrt((ax - cccx)**2 + (ay - cccy)**2)
        #print 'in get circum got ' + str(self.ccc.x) + ', ' + str(self.ccc.y) + ', ' + str(self.ccr)
    

    def isInCircumCircle(self, vertex):
        ''' True if vertex is in circumcenter '''
        if(self.ccc == 0 or self.ccr == 0):
            self.getCircumCircle()
            
        ccx = self.ccc.x 
        ccy = self.ccc.y
        
        return (math.sqrt((vertex.x - ccx)**2 + (vertex.y - ccy)**2) <= self.ccr)

    
        
##########################################################################################

##  Delaunay triangulation
        
class delaunayTriangulation:

    points_set = [] # original set of vertices
    edges = [] # edges in the delaunay triangulation
    faces = []
    faces_flags = [] # flags as complete or not for triangles
     

    def __init__(self, points_set):

        self.points_set = [points_set[i] for i in [k[0] for k in sorted(enumerate([X[0] for X in points_set]), key = operator.itemgetter(1))]]
        #print self.points_set
        self.edges = []
        self.faces = []
        self.faces_flags = []
        max_x,min_x = max([X[0] for X in self.points_set]),min([X[0] for X in self.points_set])
        delta_x = max_x - min_x
        max_y,min_y = max([X[1] for X in self.points_set]),min([X[1] for X in self.points_set])
        delta_y = max_y - min_y

        self.x0, self.y0 = min_x - 2*delta_x, max_y + 2*delta_y
        self.x1, self.y1 = min_x + delta_x/2, min_y - 2*delta_y
        self.x2, self.y2 = max_x + 2*delta_x, max_y + 2*delta_y
 
    def setSuperTriangle(self):
        '''
        Take points far away from points in set and
        form a super-triangle out of them
        '''
         

##        a0 = HE_vert(0,-16384)#16384  32768
##        a1 = HE_vert(-32768,16384)
##        a2 = HE_vert(32768,16384)

        a0, a1, a2 = HE_vert(self.x0,self.y0),HE_vert(self.x1,self.y1),HE_vert(self.x2,self.y2)
         

        # build new triangle, out of nowhere

        triangle = self.add_face(a0, a1, a2)

       
        

    def removeSuperTriangle(self):

        self.edges = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]

        toDelete = []
        
        for f in self.faces:

            edge = f.edge
            
            if (self.isInSuperTriangle(edge)):
                toDelete = toDelete + [f]

        #print 'toDelete' + str(len(toDelete))
        for f in toDelete:
            self.remove_face_from_list(f)
 
        
        toDelete = []
        for e in self.edges:
            
            if((e.end_vert.x, e.end_vert.y)  in [(self.x0,self.y0),(self.x1,self.y1),(self.x2,self.y2)]) or ((e.pair.end_vert.x, e.pair.end_vert.y) in [(self.x0,self.y0),(self.x1,self.y1),(self.x2,self.y2)]):
                toDelete = toDelete + [e]

        for e in toDelete:
            self.remove_half_edge(e)
            

    def isInSuperTriangle(self, edge):

        edges = edge.edges_around_face()

        
        for e in edges:

            if(e.end_vert == [] or e.pair.end_vert == []):
                continue
            
            if(((e.end_vert.x, e.end_vert.y)  in [(self.x0,self.y0),(self.x1,self.y1),(self.x2,self.y2)]) or ((e.pair.end_vert.x, e.pair.end_vert.y) in [(self.x0,self.y0),(self.x1,self.y1),(self.x2,self.y2)])):
                return True

        return False
                 

    def triangulate(self):
        '''
        Triangulate and return set of edges
        '''
        if(len(self.points_set) < 3):
            print 'WARNING: not enough points in data set '
            return

        self.setSuperTriangle()
        
        ######################## PLOT
        pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
##        for e in self.edges:
##            pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
##        pl.show()
        ########################
                
        # draw first star
        p = self.points_set[0]
        new_vert = HE_vert(p[0], p[1])

        for ek in self.faces[0].edge.edges_around_face() :
            self.add_face(ek.end_vert, ek.pair.end_vert, new_vert)

        self.faces = self.faces[1:]
        
##        pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
##        for e in self.edges:
##            pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
##        pl.show()
                     
        for p in self.points_set[1:]: # iterate over points
            ##print ' p' 
            new_vert = HE_vert(p[0], p[1])
            edges_in_circumcircling_triangles = [] # list of edges to be updated with edges in triangles that circumcicle the new vertex
            remove_triangles = [] # keep track of those triangles you will have to delete

            for k in range(len(self.faces)): # iterate over triangles
                ##print 'iteration over faces ' + str(k) + ', len(faces)  ' + str(len(self.faces))
                triangle = self.faces[k] # current triangle
                
                if(triangle.isInCircumCircle(new_vert)):
                    ## print 'new vertex is in triangle ' 
                    ## retrieve edges in triangle (and their pairs)
                    edges_in_triangle = triangle.edge.edges_around_face()

##                    print 'edges_in_triangle '
##                    for e in edges_in_triangle:
##                        print str(e.end_vert.x) + ', ' + str(e.end_vert.y)  

                    copy_edges = [HE_edge(HE_vert(X.end_vert.x, X.end_vert.y)) for X in edges_in_triangle] # copy edges
                    copy_pairs = [HE_edge(HE_vert(X.pair.end_vert.x, X.pair.end_vert.y)) for X in edges_in_triangle] # copy pairs
                    
                    for e_i in range(len(copy_edges)): # pair pairs
                        copy_edges[e_i].pair = copy_pairs[e_i] 
                    
                    edges_in_circumcircling_triangles = edges_in_circumcircling_triangles + copy_edges # add edges in this triangle to list of edges to examine
                    
                    #print 'we will remove triangle ' + str(self.faces[k]) + ' with edge ' + str(self.faces[k].edge.end_vert.x) + ', ' + str(self.faces[k].edge.end_vert.y)

                    remove_triangles = remove_triangles + [self.faces[k]] # (plan to) delete this triangle
            
##            print 'all these triangles I have to remove ' + str(remove_triangles) + ', ' + str(len(self.faces))
##
##            print 'edges to be examined (edges_in_circumcircling_triangles)'
##            for e in edges_in_circumcircling_triangles:
##                print e.end_vert.x, e.end_vert.y, e.pair.end_vert.x, e.pair.end_vert.y
##            
            # delete the triangle you have to delete
            for tr in remove_triangles:

##                print 'remove triangle ' + str(tr) + ' with edges '
##                for e in tr.edge.edges_around_face():
##                    print e.end_vert.x, e.end_vert.y
##                    
                self.remove_face_from_list(tr)
            
##            print 'self.edges after removing triangles '
##            for e in self.edges:
##                print 'x ' + str(e.end_vert.x) + ', y ' + str(e.end_vert.y) + ', pair x ' + str(e.pair.end_vert.x) + ', pair y ' + str(e.pair.end_vert.y)
##
##            print 'stored edges are # ' + str(len(edges_in_circumcircling_triangles))
##            for e in edges_in_circumcircling_triangles:
##                print 'x ' + str(e.end_vert.x) + ', y ' + str(e.end_vert.y) + ', pair x ' + str(e.pair.end_vert.x) + ', pair y ' + str(e.pair.end_vert.y)

            remove_cavity_interior = []
##            print 'before remove edges in cavity ' + str(len(edges_in_circumcircling_triangles))
##            for e in edges_in_circumcircling_triangles:
##                print e.end_vert.x, e.end_vert.y, e.pair.end_vert.x, e.pair.end_vert.y
            # examine all edges belonging to triangles that circumcircle the new vertex
            for k in range(len(edges_in_circumcircling_triangles)):
                
                ek = edges_in_circumcircling_triangles[k]

                for l in range(len(edges_in_circumcircling_triangles)):
                    
                    if(l == k):
                        continue

                    el = edges_in_circumcircling_triangles[l]
##                    print 'el ' + str(el.end_vert.x) + ', ' + str(el.end_vert.y) + ', ' + str(el.pair.end_vert.x) + ', ' + str(el.pair.end_vert.y)
##                    print 'ek ' + str(ek.end_vert.x) + ', ' + str(ek.end_vert.y) + ', ' + str(ek.pair.end_vert.x) + ', ' + str(ek.pair.end_vert.y)
                    # same direction
                    if(ek.end_vert.x == el.pair.end_vert.x and ek.end_vert.y == el.pair.end_vert.y and ek.pair.end_vert.x == el.end_vert.x and ek.pair.end_vert.y == el.end_vert.y):
                        remove_cavity_interior = remove_cavity_interior + [k] # (plan to )remove ek from edges_in_circumcircling_triangles
                    # other direction
                    elif(ek.pair.end_vert.x == el.end_vert.x and ek.pair.end_vert.y == el.end_vert.y and ek.end_vert.x == el.pair.end_vert.x and ek.end_vert.y == el.pair.end_vert.y):
                        remove_cavity_interior = remove_cavity_interior + [k] # (plan to )remove ek from edges_in_circumcircling_triangles

##            print 'len cavity interior ' + str(len(remove_cavity_interior))
            for k_e in remove_cavity_interior:
                
                e = edges_in_circumcircling_triangles[k_e]
##                print 'remove edges in cavity ' + str(e.end_vert.x) + ', ' + str(e.end_vert.y) + ' ; ' + str(e.pair.end_vert.x) + ', ' + str(e.pair.end_vert.y) 
                self.remove_half_edge(e)

##              ############   PLOT
##                X_remain = [X.end_vert.x for X in edges_in_circumcircling_triangles]
##                Y_remain = [X.end_vert.y for X in edges_in_circumcircling_triangles]
##                pl.scatter(X_remain, Y_remain, c = 'k', s = 100)
##                pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##                pl.scatter(p[0],p[1],c='r')
##
##                self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
##                for e in self.edges:
##                    pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
##                pl.show()
##              ##########
                          
            edges_in_circumcircling_triangles = [edges_in_circumcircling_triangles[k] for k in range(len(edges_in_circumcircling_triangles)) if(k not in remove_cavity_interior)] # not removed edges
            # around removed triangles
                
            for k in range(len(edges_in_circumcircling_triangles)):

                 ek = edges_in_circumcircling_triangles[k]

##                 if(ek.end_vert.x == ek.pair.end_vert.x and  ek.end_vert.x == new_vert.x):
##                     continue
##                 elif(ek.end_vert.y == ek.pair.end_vert.y  and  ek.end_vert.y == new_vert.y):
##                     continue
##                 elif((ek.end_vert.x !=  ek.pair.end_vert.x and ek.end_vert.x != new_vert.x) and self.areColinear(ek.end_vert, ek.pair.end_vert, new_vert)):
##                     continue
                    
                 if(ek.end_vert != [] and ek.pair.end_vert != []): # check whether it has been nullified in the head and tail
                     
                     self.add_face(ek.end_vert, ek.pair.end_vert, new_vert)
                     
                    # print 'I just added this face in mesh (' + str(ek.end_vert.x) + ', ' + str(ek.end_vert.y) + '), (' + str(ek.pair.end_vert.x) + ', ' + str(ek.pair.end_vert.y) + '), (' + str(new_vert.x) + ', ' + str(new_vert.y) + ')'
            
##            ##########   PLOT
##            X_remain = [X.end_vert.x for X in edges_in_circumcircling_triangles]
##            Y_remain = [X.end_vert.y for X in edges_in_circumcircling_triangles]
##            pl.scatter(X_remain, Y_remain, c = 'k', s = 100)
##            pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##            pl.scatter(p[0],p[1],c='r')
##
##            self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
##            for e in self.edges:
##                pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
##            pl.show()
##            ##########
##            
##            pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##            pl.scatter(p[0],p[1],c='r')

            self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
            
        
        #print self.edges
        
##        X_remain = [X.end_vert.x for X in edges_in_circumcircling_triangles]
##        Y_remain = [X.end_vert.y for X in edges_in_circumcircling_triangles]
##        pl.scatter(X_remain, Y_remain, c = 'k', s = 100)
##        pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
##        pl.scatter(p[0],p[1],c='r')
##
##        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
##        for e in self.edges:
##            pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
##        pl.show()
        
        self.removeSuperTriangle()

        self.drawEdges()

        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
        

    def areColinear(self,a1,a2,a3):
        
        if((a2.y - a1.y)/(float)(a2.x - a1.x) == (a3.y - a1.y)/(float)(a3.x - a1.x)):
            return True
        return False

    def add_face(self, v0, v1, v2):
        '''
        Add triangle to list of triangles and update list of edges
        '''
        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
        
        # counter clock wise everything
        v0, v1, v2 = self.counterClockWise(v0,v1,v2)
        
        # create edges
        edge0,edge1,edge2 = HE_edge(v1),HE_edge(v2),HE_edge(v0)
         
        edge0.next_edge = edge1
        edge1.next_edge = edge2
        edge2.next_edge = edge0
        
        # pairs
        pair0 = HE_edge(v0)
        pair1 = HE_edge(v1)
        pair2 = HE_edge(v2)
        
        # pair pairs
        edge0.pair, edge1.pair, edge2.pair = pair0,pair1,pair2
        pair0.pair, pair1.pair, pair2.pair = edge0,edge1,edge2

        face = HE_face(edge0) # this is the face of any half-edge surronding the face
        
        # face      
        edge0.face = face
        edge1.face = face
        edge1.face = face
        
        self.add_half_edge(edge0)
        self.add_half_edge(edge1)
        self.add_half_edge(edge2)
        self.add_half_edge(pair0)
        self.add_half_edge(pair1)
        self.add_half_edge(pair2)

        face.getCircumCircle()
        
        # append to list of faces
        self.faces = self.faces + [face]
        self.faces_flags = self.faces_flags + [0]
        #print 'in add_face ' + str(face.ccc.x)
        return face

    
 
    def add_half_edge(self, edge_add):
        '''
        Add edge to self list of edges
        '''
        
        # check if already stored
        for i in range(len(self.edges)):    
            edge = self.edges[i]

            # if end and start vertex are same, and face is same
            if(edge.end_vert.x == edge_add.end_vert.x and edge.end_vert.y == edge_add.end_vert.y
               and edge.pair.end_vert.x == edge_add.pair.end_vert.x and edge.pair.end_vert.y ==  edge_add.pair.end_vert.y):
                
                #print 'omit to add ' + str(edge_add.end_vert.x) + ', ' + str(edge_add.end_vert.y) + ', with next '
##                if(edge_add.next_edge != None):
##                    print str(edge_add.next_edge.end_vert.x) + ', ' + str(edge_add.next_edge.end_vert.y)
##                else:
##                    print 'None' 
                return
            
        # else, store it now
 
        self.edges = self.edges + [edge_add]
        
##        print 'added ' + str(edge_add.end_vert.x) + ', ' + str(edge_add.end_vert.y) + ', with next '
##        if(edge_add.next_edge != None):
##            print str(edge_add.next_edge.end_vert.x) + ', ' + str(edge_add.next_edge.end_vert.y)
##        else:
##            print 'None'

        #print 'from add half edge ' + str(edge_add.end_vert.x) + ', ' + str(edge_add.end_vert.y)
        #for e in self.edges:
            #print str(e.end_vert.x) + ', ' + str(e.end_vert.y)
 


##    def remove_face(self, tr_del, triangles_to_remove):
##        ''' Remove this triange, tr_del, from mesh
##        Given that you will have to withdraw all triangles in triangles_to_remove as well
##        '''
##        # remove edges
##        edge, n_edge, n2_edge = tr_del.edge, tr_del.edge.next_edge, tr_del.edge.next_edge.next_edge
##
##        remove  = []
##
##        if(edge.end_vert != []):
##            if(self.belongToFace(edge, triangles_to_remove) == False):
##                print 'not belong to face edge ' + str(edge.end_vert.x) + ', ' + str(edge.end_vert.y)
##                remove = remove + [0]
##            else:
##                print 'belong to face edge ' + str(edge.end_vert.x) + ', ' + str(edge.end_vert.y)
##
##        if(n_edge.end_vert != []):
##            if(self.belongToFace(n_edge, triangles_to_remove) == False):
##                print 'not belong to face n_edge ' + str(n_edge.end_vert.x) + ', ' + str(n_edge.end_vert.y)
##                remove = remove + [1]
##            else:
##                print 'belong to face edge ' + str(n_edge.end_vert.x) + ', ' + str(n_edge.end_vert.y)
##
##        if(n2_edge.end_vert != []):   
##            if(self.belongToFace(n2_edge, triangles_to_remove) == False):
##                print 'not belong to face n2_edge ' + str(n2_edge.end_vert.x) + ', ' + str(n2_edge.end_vert.y)
##                remove = remove + [2]
##            else:
##                print 'belong to face edge' + str(n2_edge.end_vert.x) + ', ' + str(n2_edge.end_vert.y)
##
##
##        # remove triangle
##        remove_face_from_list(tr)
##
##        # remove edges
##        if(0 in remove):
##            self.remove_half_edge(edge)
##            
##        if(1 in remove):
##            self.remove_half_edge(n_edge)
##
##        if(2 in remove):
##            self.remove_half_edge(n2_edge)


    def remove_face_from_list(self, tr_del):

        edges_del = tr_del.edge.edges_around_face()

        #print 'check remove from list '

        #for e in edges_del:
            #print e.end_vert.x, e.end_vert.y, e.pair.end_vert.x, e.pair.end_vert.y
        
        # remove triangle
        for i in range(len(self.faces)):
            #print 'face i ' + str(i)
            tr = self.faces[i]
            edges = tr.edge.edges_around_face()

            #for e in edges:
                #print e.end_vert.x, e.end_vert.y, e.pair.end_vert.x, e.pair.end_vert.y

            # if end and start vertex are same, and face is same
            if(edges_del[0].compare(edges[0]) and edges_del[1].compare(edges[1]) and edges_del[2].compare(edges[2])):
                
                self.faces = self.faces[:i] + self.faces[i+1:] # remove this face for this edge, not for its pair

                break # found it

        #print 'after ' + str(len(self.faces))
##        for tr in self.faces:
##            edges = tr.edge.edges_around_face()
##            for e in edges:
                #print e.end_vert.x, e.end_vert.y, e.pair.end_vert.x, e.pair.end_vert.y
        
 
    def belongToSeveralFace(self, e):
        '''
        Return true if edge belongs to some face
        '''
        in_face = 0
        #print 'find if belongs to some face '
        for face in self.faces:
            
            e1,e2,e3 = face.edge, face.edge.next_edge, face.edge.next_edge.next_edge 
       
            if(e1.compare(e) or e2.compare(e) or e3.compare(e)):
                
                in_face = in_face + 1
                
                print 'in_face ' + str(in_face)

                if(in_face > 1):
                    
                    return True
                       
        return False

 
        
    def remove_half_edge(self,edge_del):
        '''
        Remove this edge from list of edges
        !! When remove is called, all half-edges with no pair are deleted
        !! one half_edge is deleted not its half-edge
        '''
 
        if(edge_del.end_vert == [] or edge_del.pair.end_vert == []):
                edge_del.pair.end_vert = []
                edge_del.end_vert == []
                self.edges = [X for X in self.edges if(X.end_vert != [] and X.pair != None)]
                return

        self.edges = [X for X in self.edges if(X.end_vert != [] and X.pair != None)]
  
        for i in range(len(self.edges)):
            
            edge = self.edges[i]

            if(edge_del.end_vert == []):
                edge_del.pair.end_vert == []
                break # already found

            #print 'remove half edge debug ' + str(edge_del.end_vert) + ', ' + str(edge.end_vert) + ', ' + str(edge.pair)
            # if end and start vertex are same, and same for pair
            if(edge_del.compare(edge) or edge_del.compare(edge.pair)):
                
                self.edges[i].pair.end_vert = []
                self.edges[i].end_vert = []

        self.edges = [X for X in self.edges if(X.end_vert != [] and X.pair != None)]



        

    def voronoi(self):
        '''
        Return the Voronoi diagram that is dual to this Delaunay triangulation
        '''
        voronoi = []

##        print len(self.faces)
##        for f in self.faces:
##            print 'ef ' + str(f.edge) + ', ' + str(f.edge.end_vert.x) + ', ' + str(f.edge.end_vert.y)

        circums_edges = [] # list of circumcenters for faces around vertex for half_edge
        faces_around_vtx_edges = [] # list of faces around vertex for half-edge, same order
        
        for e in self.edges:
            
            #print 'e ' + str(e)

            vtx = e.end_vert
            
            faces_around_vtx, circums = [], []

            #print 'vtx ' + str(vtx.x) + ', ' + str(vtx.y)
            
            for f in self.faces:
                
                edges = f.edge.edges_around_face()
                #print 'edges around ' + str(edges)
                
                for ef in edges:
                    
                    if(ef.end_vert.x == vtx.x and ef.end_vert.y == vtx.y):

                        faces_around_vtx = faces_around_vtx + [f]

                        circums = circums + [[f.ccc.x,f.ccc.y]]

            circums = [list(X) for X in list(np.unique([(x,y) for x,y in circums]))]
            
            idx, circums = self.clockWisePolygon(circums)

            faces_around_vtx_edges = faces_around_vtx_edges + [[faces_around_vtx[i] for i in idx]]

            circums_edges = circums_edges + [circums]
            
##            print 'circums_edges ' + str(circums_edges) + ' faces around vertex ' + str(faces_around_vtx_edges)

        self.drawVoronoi( faces_around_vtx_edges, circums_edges )
        
        return faces_around_vtx_edges, circums_edges




    def drawVoronoi(self, faces = None, circumcenters = None):

        if( faces == None or circumcenters == None ):
            faces, circumcenters = self.voronoi()

        for e_circum in circumcenters:

            pl.plot([X[0] for X in e_circum], [X[1] for X in e_circum], c='k')

        for e in self.edges:
            pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')

        pl.scatter([X[0] for X in self.points_set],[X[1] for X in self.points_set], c='b')

        pl.show()            
                        

    def clockWisePolygon(self, points):
        
        if(len(points) == 1):
            return [0], points
        
        idx_ordered = [k[0] for k in sorted(enumerate(points), key = operator.itemgetter(1))]
        ordered = [points[i] for i in idx_ordered] # x ordered
        
        left_most, right_most = ordered[0], ordered[len(ordered)-1]

        if(right_most[0] == left_most[0]): # a column of more than one points
            idx_sorted = [k[0] for k in sorted(enumerate([X[1] for X in points]), key = operator.itemgetter(1))]
            return idx_sorted, [points[i] for i in idx_sorted]
        
        a_compare = (right_most[1] - left_most[1])/float(right_most[0] - left_most[0])
        above, below, above_indices, below_indices = [], [], [], []
        
        for i in range(len(ordered)):

            p = ordered[i]

            if(p == left_most or p == right_most):
                continue
            
            y_compare = a_compare * (p[0] - left_most[0]) + left_most[1]

            if(p[1] >= y_compare):
                above = above + [p]
                above_indices = above_indices + [i]
                
            else:
                below = below + [p]
                below_indices = below_indices + [i]

        counterCW = [left_most]
        indices = [0]
        
        for k in range(len(above)):
            counterCW = counterCW + [above[k]]
            indices = indices + [above_indices[k]]

        counterCW = counterCW + [right_most]
        indices = indices + [len(ordered)-1]

        
        for k in range(len(below)):
            counterCW = counterCW + [below[::-1][k]]
            indices = indices + [below_indices[k]]
        #print 'counterCW ' + str(counterCW)
        return indices, counterCW 
        

    def counterClockWise(self, v0, v1, v2):
        ''' Order vertices in a CCW manner
        '''
        v = [v0, v1, v2]
        set_v = [(v0.x,v0.y),(v1.x,v1.y),(v2.x,v2.y)]
        set_v_idx = [i for i in [k[0] for k in sorted(enumerate([X[0] for X in set_v]), key = operator.itemgetter(1))]]

        ccw_v0 = v[set_v_idx[0]]
        
        if(set_v[set_v_idx[0]][0] == set_v[set_v_idx[2]][0]):
            if(set_v[set_v_idx[1]][1] > set_v[set_v_idx[0]][1]):
                ccw_v1 = v[set_v_idx[2]]
                ccw_v2 = v[set_v_idx[1]]
            elif(set_v[set_v_idx[1]][1] <= set_v[set_v_idx[0]][1]):
                ccw_v1 = v[set_v_idx[1]]
                ccw_v2 = v[set_v_idx[2]]
        else:

            a_compare = (set_v[set_v_idx[2]][1] - set_v[set_v_idx[0]][1])/float(set_v[set_v_idx[2]][0] - set_v[set_v_idx[0]][0])
            y_compare = (set_v[set_v_idx[1]][0] - set_v[set_v_idx[0]][0]) * a_compare + set_v[set_v_idx[0]][1]

            if(set_v[set_v_idx[1]][1] == y_compare):
                ccw_v1 = v[set_v_idx[1]]
                ccw_v2 = v[set_v_idx[2]]

            elif(set_v[set_v_idx[1]][1] > y_compare):
                ccw_v1 = v[set_v_idx[2]]
                ccw_v2 = v[set_v_idx[1]]

            elif(set_v[set_v_idx[1]][1] < y_compare):
                ccw_v1 = v[set_v_idx[1]]
                ccw_v2 = v[set_v_idx[2]]
 
        return ccw_v0, ccw_v1, ccw_v2


    def drawEdges(self):
        if(self.points_set != []):
            pl.scatter([X[0] for X in self.points_set], [X[1] for X in self.points_set])
        self.edges  = [X for X in self.edges if (X.end_vert != [] and X.pair != None and X.pair.end_vert != [])]
        for e in self.edges:
            pl.plot([e.end_vert.x, e.pair.end_vert.x],[e.end_vert.y, e.pair.end_vert.y], c = 'b')
        pl.show()
            
##########################################################################################

##  Alpha shape

class alphaShape:

    points_set = []
    convex_hull = []
    DT = []
    edges = []
    faces = []
    faces_edges = []
    circums_edges = []
   


    def __init__(self, points_set, DT=None):

        self.points_set = points_set

        if(DT == None):
            self.DT = delaunayTriangulation(points_set)
            self.DT.triangulate()
        else:
            self.DT = DT

        self.edges = self.DT.edges
        dt_faces = self.DT.faces

        self.faces_edges, self.circums_edges = self.DT.voronoi() # voronoi faces and circumcenter for each vertex at the end of an edge in self.edge





    def alpha_neighbors(self):
        '''
        '''
        alpha_neighbors = []
        concave_hull = []

        hull = geometry.convexHull(self.points_set)
        full_hull = geometry.full_hull(hull) # full convex hull

        alpha_mins, alpha_maxs = [], []

        for k in range(len(self.edges)):

            e = self.edges[k]

            bounds = self.alpha_bounds(k, full_hull)

            if(bounds == False):
                alpha_mins = alpha_mins + ['']
                alpha_maxs = alpha_maxs + ['']
                continue
            else:
                alpha_min, alpha_max = bounds

            alpha_mins = alpha_mins + [alpha_min]
            alpha_maxs = alpha_maxs + [alpha_max]
        print "alpha_min " + str(alpha_mins)
        print "alpha_max " + str(alpha_maxs)
        #alpha_max = (np.mean([X for X in alpha_mins if X != '']) + np.mean([X for X in alpha_maxs if X != '']))/float(2)
        #alpha = np.mean([X for X in alpha_mins if X != ''])
        alpha = np.mean([X for X in alpha_maxs if X != ''])
        alpha = list(scipy.stats.mstats.mquantiles([X for X in alpha_maxs if X != ''],[0.80]))[0]
        #print alpha_maxs, alpha_mins

        
        for k in range(len(self.edges)):
            
            e = self.edges[k]

            if((alpha_maxs[k] != '' and  alpha_mins[k] <= alpha  and alpha <= alpha_maxs[k]) or (alpha_maxs[k] == '' and alpha_mins[k] <= alpha)):
                
                alpha_neighbors = alpha_neighbors + [e] 

##        not_in_concave_hull = []
##        for k in range(len(alpha_neighbors)):
##
##            ek = alpha_neighbors[k]
##
##            for l in range(len(alpha_neighbors)):
##
##                if(l != k):
##                            
##                    el = alpha_neighbors[l]
##
##                    if((el.pair.end_vert.x == ek.end_vert.x and el.pair.end_vert.y == ek.end_vert.y and el.end_vert.x == ek.pair.end_vert.x and el.end_vert.y == ek.pair.end_vert.y) # not equal to pair
##                       or (el.end_vert.x == ek.end_vert.x and el.end_vert.y == ek.end_vert.y and el.pair.end_vert.x == ek.pair.end_vert.x and el.pair.end_vert.y == ek.pair.end_vert.y)): # not equal to edge
##
##                        not_in_concave_hull = not_in_concave_hull + [k]
##                        continue


##        concave_hull = [alpha_neighbors[k] for k in range(len(alpha_neighbors)) if k not in not_in_concave_hull]
        
        return alpha_neighbors, alpha_neighbors
            



            
    def alpha_bounds(self, index, convex_hull):
        '''
        Return alpha_min and alpha_max for a given edge at position index in self.edges
        There is some distinction in bounds finding, whether points belongs or not to convex hull
        '''
        edge = self.edges[index]

        circums_pair = self.circums_edges[self.index_edge(edge.pair)]
        #print 'circums pair ' + str(circums_pair)
        circums_edge = self.circums_edges[index]
        #print 'circums edge ' + str(circums_edge)

        intersect = [X for X in circums_pair if(X in circums_edge)] # voronoi edge dual to delaunay edge at index
        
        if(len(intersect) < 2):
            return False

        v1, v2 = intersect[0], intersect[1]

        if(v1 == v2):
            return False
        
        #print 'intersect '
        #print v1, v2

        v = math.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2)

        sigma_1 = math.sqrt((edge.end_vert.x - v1[0])**2 + (edge.end_vert.y - v1[1])**2) # first circumradius
        sigma_2 = math.sqrt((edge.end_vert.x - v2[0])**2 + (edge.end_vert.y - v2[1])**2) # second circumradius

        alpha_min = min(math.sqrt(sigma_1**2 + ((sigma_1**2 - sigma_2**2 + v**2)/float(2 * v))**2),min(sigma_1, sigma_2))

        alpha_max = max(sigma_1, sigma_2)
        
        #alpha_min = min(sigma_1, sigma_2)

        
        #print ( v2[0] - v1[0] )*(edge.pair.end_vert.x - edge.end_vert.x) + ( v2[1] - v1[1] )*(edge.pair.end_vert.y - edge.end_vert.y) # check scalat product
##        pl.plot([edge.end_vert.x, edge.pair.end_vert.x], [edge.end_vert.y, edge.pair.end_vert.y], c='r')
##        pl.plot([v1[0], v2[0]] , [v1[1], v2[1]])  
##        pl.show()

##        if([edge.end_vert.x, edge.end_vert.y] in convex_hull):
##            
##            return alpha_min, ''
        
        return alpha_min, alpha_max

        


    def index_edge(self, edge):
        '''
        Return index in self.edges for some edge you don't know the index
        '''
        for k in range(len(self.edges)):
            
            e = self.edges[k]

            if(e.end_vert.x == edge.end_vert.x and e.end_vert.y == edge.end_vert.y and
               e.pair.end_vert.x == edge.pair.end_vert.x and e.pair.end_vert.y == edge.pair.end_vert.y):

                return k

        return -1 # not found
                

        

        
    
        

 
