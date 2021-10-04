##########################################################################################
##########################################################################################
####
####             WATERSHED Algorithm, peak detection  -- test
#### 
##########################################################################################
import sys
sys.path.append('C:\\Users\\') ##todo

import math
import numpy as np
#import pylab as pl
import matplotlib
import scipy
#from scipy import ndimage
#import mpl_toolkits.mplot3d.axes3d as p3
#from pylab import *

import scipy

from operator import itemgetter


import os

DATA_DIR = os.path.abspath(os.path.join(os.getcwd(), 'IMAGE'))
sys.path.append(DATA_DIR)

import  peakDetect, watershed_V6, watershedVizu, filters, io


### 1.
# Generate matrix, enhancement, filter for maxima and vizualize matrix with maxima
# compare list of synthetic maxima with list of resolved maxima
mat = io.fileToMatrix('C:\\Users\\') ##todo

watershedVizu.plot_heatmap(mat) 

peakdetectV6 = watershed_V6.WATERSHED()

filtered = peakDetect.gaussian_enhancementV1(mat,1) # noise = 0.01 filtering, 1/200 of points to build noise variance

watershedVizu.plot_heatmap(filtered)

cclabeling, filtr, maxima = filters.maximum_filter_regions(filtered)

toCartesian = [ cclabeling.forest.cell[X] for X in maxima]

print ("Maxima "  + str(toCartesian))
 
labels_to_matrix = cclabeling.matrixFromLabelsList(filtr, len(mat), len(mat[0]))

pl.matshow(labels_to_matrix, cmap=cm.gray_r)
pl.colorbar()
pl.show()
 
### 2.
### Watershed, identification of coeluted zones, clusering refinement
mat = io.fileToMatrix('F:\\') ##todo
watershedVizu.plot_heatmap(mat) 
peakdetectV6 = watershed_V6.WATERSHED()
filtered = peakDetect.gaussian_enhancementV1(mat,1) # noise = 0.01 filtering, 1/200 of points to build noise variance
watershedVizu.plot_heatmap(filtered)

# watershed is performed on filtered matrix, but original matrix is needed
labels = peakdetectV6.watershed( filtered, mat )
watershedVizu.plot_heatmap( labels )
cclabeling, filtr, maxima = filters.maximum_filter_regions(filtered)

# coeluted zones 
zones, coeluted_blobs, coeluted_neighbors = peakdetectV6.find_coeluted_zones( filtered, mat, maxima )

 
### 3.
### refinement
peakdetectV6.refine_clustering( zones, coeluted_blobs, coeluted_neighbors, filtered, mat, maxima )
