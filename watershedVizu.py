##########################################################################################
##########################################################################################
####
####             Vizualization tools -- random gaussian maps and watershed algorithm
#### 
##########################################################################################
 
import math
import itertools as it
import numpy as np
import random
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import pylab as pl
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3

import scipy
import os, sys

import mathlib, filters

def gradientImage(mat):
    '''
    Turn image into gradient image
    Derivatives are turned in discrete
    centered differences along rows and columns
    '''
    N = len(mat) # first dimension
    M = len(mat[0]) # secondary chromatogram
    
    gradX,gradY,gradMagnitude = [],[],[]

    for i in map(lambda x:x+1, range(N-2)):
        tmp = [] # new line
        for j in map(lambda x:x+1, range(M-2)):
            gradx = (mat[i+1][j] - mat[i-1][j])/float(2)
            grady = (mat[i][j+1] - mat[i][j-1])/float(2)
            gradMag = np.sqrt(gradx**2+grady**2)
            tmp.append([gradx, grady, gradMag])
        gradX.append([x[0] for x in tmp])
        gradY.append([x[1] for x in tmp])
        gradMagnitude.append([x[2] for x in tmp])

    return gradX, gradY, gradMagnitude



def plot3D(mat):
    """
    3d plot signal  
    """
    fig=pl.figure()
    ax = p3.Axes3D(fig)

    nrows = len(mat)
    ncols = len(mat[0])

    mat = mathlib.adjustMat(nrows, ncols,mat)
        
    x = arange(0, nrows, 1)
    y = arange(0, ncols, 1)
    X,Y = meshgrid(x,y)
     
    z = [tuple(mat[i]) for i in range(nrows)]
    zzip = zip(*z)
    
    ax.plot_surface(X,Y,zzip,rstride = 1, cstride = 1,cmap=cm.jet)
    pl.colorbar()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    pl.show()
    
  


def plotHeatMap(mat):
    '''
    Plot heatmap with intensity matrix mat (color depends on intensity)
    '''
     
    nrows = len(mat)
    ncols = len(mat[0])

    mat = mathlib.adjustMat(nrows, ncols,mat)
        
    x = arange(0, nrows+1, 1)
    y = arange(0, ncols+1, 1)
    X,Y = meshgrid(x,y)
     
    z = [tuple(mat[i]) for i in range(nrows)]
    zzip = zip(*z)

    pl.pcolor(X,Y,zzip)
    pl.colorbar()
    pl.axis([0,nrows,0,ncols]) ## modified
    pl.show()



def createBlob(x, y, col, mat, decr, spread):
    '''
    Add one squared blob with witdth spread to elevation matrix mat in (x,y)
    and return modified matrix.
    Blob's apex has color related to elevation col.
    '''
    dimx, dimy = len(mat), len(mat[0])
    mat[x][y] = col
    
    for incr in range(1,int(math.floor(spread/2))):
        
        col = col - 1.2 * decr
        visit = list(it.product((-incr,incr),range(-incr,incr+1))) + list(it.product(range(-incr,incr+1),(-incr,incr)))
        visit += [(0,incr)] + [(0,-incr)]+ [(incr, 0)] + [(-incr, 0)]
        
        for i,j in visit:
            
            if(j != 0 or i != 0):
                
                if(x+i >=0 and y+j >= 0 and x+i < dimx and y+j < dimy):
                    
                    mat[x+i][y+j] = col

    return mat

def heatMapGenerator(dimx, dimy, n_blobs, spread=None, randomBackground=False):
    '''
    Compute map with n_blobs squared blobs
    for a map at least 20x20 pixels
    '''
    if(dimx < 20 or dimy < 20):
        print ('WARNING: please draw a larger map') 
        return
    
    mapMat = []
    col = 0
    dx = 1/float(dimx*dimy)

    if (spread == None):
        spread = min([dimx, dimy])/float(5)
    
    # initialize map with background
    for i in range(dimx):
        tmp = []
        for j in range(dimy):
            tmp.append(col)
            col += dx
        mapMat.append(tmp)
        
    # randomize background
    if(randomBackground):
        mapMat = randomizeBackground(mapMat)
    rnd = np.random.random_sample(n_blobs)
    # add blobs
    rdx = [(int) (math.floor(x * dimx)) for x in np.random.random_sample(n_blobs)]
    rdy = [(int) (math.floor(x * dimy)) for x in np.random.random_sample(n_blobs)]
    points = zip(rdx, rdy)
    i = 0
    slack = math.floor(min([dimx, dimy])/float(20))
    dblob = min(0.1,1/float(spread))
    for x,y in points:
        spread = spread + np.random.randint(-slack,slack)
        col = 10 + i * dblob 
        matMap = createBlob(x,y,col,mapMat,dblob,spread)
        i += 1  

    return mapMat, points

def randomizeBackground(mat,noise_level,n_pts = None):
    '''
    Add random background color. n_pts points are colored in matrix
    if n_pts is given. Else, all points are colored.
    '''
    if(mat == []):
        return
    
    dimx, dimy = len(mat), len(mat[0])

    if(n_pts != None): # color n_pts points in matrix

        rdx = [(int) (math.floor(x * dimx)) for x in np.random.random_sample(n_pts)]
        rdy = [(int) (math.floor(x * dimy)) for x in np.random.random_sample(n_pts)]
        points = zip(rdx, rdy)

        for x,y in points: # add noise
            mat[x][y] = mat[x][y] + max(0, noise_level * (1 + np.random.randint(-1,1)))
            
    else: # color all points in matrix
        for i in range(dimx):
            for j in range(dimy):
                mat[i][j] = mat[i][j] + max(0, noise_level * (1 + np.random.randint(-1,1)))
                
    return mat




def GCxGCGenerator(dimx, dimy, n_blobs, elevation, addWhiteNoise = False, noise_level = 0.01, oriented = False):
    '''
    Compute map with n_blobs blobs for a map at least 20x20 pixels
    If addWhiteNoise is true, add white noise with mean 0 and variance noise_level, default 0.01.
    If oriented is True, blobs are randomly oriented, and not necessarily parallel to axis.
    '''
    if(dimx < 20 or dimy < 20):
        print ('WARNING: please draw a larger map') 
        return
    
    Xin, Yin = mgrid[0:dimx, 0:dimy]

    mapMat = []

    # initialize map with 0 background
    for i in range(dimx):
        tmp = []
        for j in range(dimy):
            tmp.append(0.)
        mapMat.append(tmp)
        
    # add blobs
    rdx = [(int) (math.floor(x * (dimx-1))) for x in np.random.random_sample(n_blobs)]
    rdy = [(int) (math.floor(y * (dimy-1))) for y in np.random.random_sample(n_blobs)]
    points = zip(rdx, rdy)

    volumes = []
    for j in range(len(points)):
        volumes.append(0.)
        
    k = 0
    
    for center_x, center_y in points:
        
        width_x = math.floor(min([dimx, dimy])/float(10*np.random.randint(2,5))) * (1 + np.random.random()) 
        width_y = math.floor(min([dimx, dimy])/float(10*np.random.randint(2,5))) * (1 + np.random.random())
        col = elevation * (1 + (np.random.random()))
        
        if(oriented):
            
            gauss_ = mathlib.gaussian_oriented(col, center_x, center_y, width_x, width_y)(Xin,Yin)
  
        else:
            
            gauss_ = mathlib.gaussian(col, center_x, center_y, width_x, width_y)(Xin,Yin)
            
        for i in range(dimx):
            for j in range(dimy):
                mapMat[i][j] += gauss_[i][j]
                volumes[k] += gauss_[i][j]

        k += 1
               
    if(addWhiteNoise == True):
        mapMat = filters.add_gaussian_noise(mapMat, 0, noise_level)
        
    return mapMat, points, volumes




def plot_heatmap(mat):
    pl.matshow(mat, cmap=cm.jet)
    pl.colorbar()
    pl.show()


    
def plotWatershedCarpetPlot(points, mat, cells_regions, N_COLOR):
    '''
    Plot carpet plot, elevation is given in mat
    whereas label for region is in labels_mat
    '''
    fig=pl.figure()
    ax = p3.Axes3D(fig)

    nrows = len(mat)
    ncols = len(mat[0])

    mat = mathlib.adjustMat(nrows, ncols,mat)

    x = arange(0, nrows, 1)
    y = arange(0, ncols, 1)
    X,Y = meshgrid(x,y)

    z = [tuple(mat[i]) for i in range(nrows)]
    zzip = zip(*z)
    
    cm = pl.get_cmap('gist_rainbow')
    color = []
    COLORS = np.empty(X.shape, dtype=str)
    for i in range(N_COLOR):
        color = color + [cm(1.*i/N_COLOR)]  # color will now be an RGBA tuple

    COLORS[:,:] = color[0]
    
    color = color[1:]

    for c in len(cells_regions):

        cells = cells_regions[c]
        
        for xc,yc in cells:
            COLORS[xc,yc] = color[c]

    ax.plot_surface(X,Y,zzip,rstride = 1, cstride = 1,cmap=cm.jet)
    
    Xp,Yp = [x[0] for x in points], [x[1] for x in points]

    Zp = []
    
    for k in range(len(points)):
        Zp.append(mat[Xp[k]][Yp[k]])
    ax.scatter(Xp,Yp,Zp,c='r')
  
    pl.show()
    


 
            
            

    
