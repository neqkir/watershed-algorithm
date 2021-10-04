##########################################################################################
##########################################################################################
####
####             Filters for image processing
#### 
##########################################################################################


#####  FILTERS
## - maximum filter
## - Sobel filter
## - enghancement filter (Lagrangian)
## - Gaussian filter at given scale 
## - LEE's adaptative filter 


##### TOOLS
## - padding (value padding, copy padding)
## - maximum filter working for connected regions
## - tools to add noise (addidtive, multiplicative, both)


import math
import itertools as it
import numpy as np
import random
import pylab as pl

import mathlib, CC_labeling_8 , sys


def maximum_filter(mat, str_elem = 0):
    ''' Return
    . output of maximum filter with connectivity defined by structuring
    element str_elem applied to matrix mat
    . list of maxima
    1. str_elem is a 0's and 1's matrix with 1's for neighboring points
    str_elem is allowed to have only odd x and y dimensions
    2. output is 1 if pixel value is more or equal to all neighboring pixels
    '''
 
    if(str_elem == 0): # default is height-connectivity
        
        str_elem = []

        for i in range(3):
            tmp = []
            for j in range(3):
                tmp.append(1)  
            str_elem.append([x for x in tmp])
 
    dimx_se, dimy_se = len(str_elem), len(str_elem[0])
    
    if(dimx_se%2 == 0 or dimy_se%2 == 0):
        print 'WARNING: structuring element should have odd x and y dimensions'
        return
    
    filtered_img, maxima = [], []
    
    N,M = len(mat), len(mat[0])
     
    for i in range(N):
        tmp = []
        for j in range(M):

            if(apply_maximum_filter(mat,str_elem,i,j) == 1):
                maxima.append((i,j))
                tmp.append(1)
            else:
                tmp.append(0)

        filtered_img.append([x for x in tmp])

    return filtered_img, maxima


def apply_maximum_filter(mat,str_elem,i,j):
    '''Apply maximum filter at i, j, with structuring element str_elem
    '''
    dimx_se, dimy_se = len(str_elem), len(str_elem[0])

    if(dimx_se%2 == 0 or dimy_se%2 == 0):
        print 'WARNING: structuring element should have odd x and y dimensions'
        return

    N,M = len(mat), len(mat[0])

    if(i < 0 or i >= N or j < 0 or j >= M):
        return

    var_x = math.floor(dimx_se/2)
    var_y = math.floor(dimy_se/2)

    val = mat[i][j]

    is_max = True

    for k in range(int(i - var_x), int(i + var_x + 1)):
        if(is_max == False):
            break
        for l in range(int(j - var_y), int(j + var_y + 1)):
            
            if(k >= 0 and k < N and l >= 0 and l < M):

                if(str_elem[int(k - i + var_x)][int(l - j + var_y)] == 1):
                    if(mat[k][l] > val):
                        is_max = False
                        break
                
    if(is_max == True):
        return 1
    
    else:
        return 0


def maximum_filter_regions(mat, cclabeling = None):
    '''
    Maximum filter a matrix that has previously been connected-component labeled
    '''
    N,M = len(mat),len(mat[0])

    if(cclabeling == None): ## connected component labeling with 8-connectivity
        cclabeling = CC_labeling_8.CC_lab(mat)  
        cclabeling.connectedComponentLabel(N*M) # at most N*M regions

    regions_mat, regions = cclabeling.labels, np.unique(cclabeling.forest.parent)

    ## find maximum regions 
    filtered_regions = [0] * max( cclabeling.forest.parent )
    
    max_regions = []
    i = 0
    
    for r in regions: # iter over regions
        
        isMax = True

        ri,rj = cclabeling.forest.cell[r]
        valr = mat[ri][rj]
        neighbors = cclabeling.neighborRegions(r, mat)

        if(neighbors == []):
            break

        for label_n in neighbors:

            if(label_n < 0 or label_n >= cclabeling.forest.utils):  # out of limits
                continue 

            nx, ny = cclabeling.forest.cell[label_n]  # representative for region labeled by label_n

            if(mat[nx][ny] > valr): # append only higher elevation points to higher_n
                isMax = False
                break
        
        if(isMax == True):
            print 'r' + str(r)
            filtered_regions[r] = 1
            max_regions.append(r)

    return cclabeling, filtered_regions, max_regions



def gaussian_maximum_filter(mat, scale = 0.01, dimx_kernel = 0, dimy_kernel = 0, padding = True, isLocal = False):
    '''local gaussian maximum filter with connectivity defined by dimx_kernel and dimy_kernel, applied to matrix mat. It is simply gaussian filtering followed
    by local maximum filtering. Filtering at scale t operated with a gaussian kernel with variance t along x and y axis
    OUTPUT
    . matrix with 0's and 1's, 1 for maxima after gaussian blur
    . list of maxima
    INPUT
    . scale is blur variance along both dimensions
    . padding True means all cells in mat are examined. output matrix is same size as mat. padding
    False means kernel is adjusted to fit corners and bordures  
    . isLocal True means that variance (scale) is adaptative. isLocal False means that variance is scale
    '''
    if(dimx_kernel == 0 and dimy_kernel == 0):
        dimx_kernel, dimy_kernel = 5, 5
        
    # first locally blur, with padding
    g_filtered = gaussian_filter(mat,dimx_kernel,dimy_kernel, scale, padding)
    
    # apply maximum filter with 8-connectivity
    return maximum_filter(g_filtered)



def gaussian_maximum_filter_regions(mat, scale = 0.01, dimx_kernel = 0, dimy_kernel = 0, padding = True, isLocal = False):
    '''
    Maximum filter a matrix that has previously been connected-component labeled
    '''
    if(dimx_kernel == 0 and dimy_kernel == 0):
        dimx_kernel, dimy_kernel = 5, 5
        
    # first locally blur, with padding
    g_filtered = gaussian_filter(mat,dimx_kernel,dimy_kernel, scale, True)

    # apply maximum filter with regions-like-connectivity
    return maximum_filter_regions(g_filtered)


    
def gaussian_kernel(sigma, dimx_filter, dimy_filter):
    '''
    Gaussian kernel for filter with input given dimension, and for given variance sigma
    '''
    # dimx and dimy for filter should be odd
    if(dimx_filter%2 == 0 or dimy_filter%2 == 0):
        return
    
    Xin, Yin = pl.mgrid[0:dimx_filter, 0:dimy_filter]
    
    center_x, center_y = math.floor(dimx_filter/2), math.floor(dimy_filter/2)

    kernel = mathlib.gaussian_kernel(center_x, center_y, sigma)(Xin,Yin)
    
    return kernel


def gaussian_filter(mat, dimx_filter, dimy_filter, sigma, padding = False, padding_copy = False):
    '''
    Perform gaussian filtering on matrix mat, for dimx_filter,dimy_filter dimensions filter
    with variance w_x, w_y alongx and y and with amplitude height
    '''
    N, M = len(mat),len(mat[0])
    
    # create gaussian filter
    kernel = gaussian_kernel(sigma, dimx_filter, dimy_filter)
 
    # convolution with mat
    filtered = convolution(mat, kernel, padding, padding_copy)
     
    return filtered

def convolution(img, kernel, padding = False, padding_copy = False, pad_value = 0):
    '''
    Convolute image img with kernel. Options for padding are available: padding uses padding with a constant value pad_value
    and padding_copy uses padding with same value as edge value (in the corner, value is zero)
    '''
    N,M,n,m = len(img),len(img[0]),len(kernel),len(kernel[0])
    #print range(int(m/2),M-int(m/2)+1)

    conv_mat = []
    #print 'M ' + str(M) + ' N ' + str(N) + ' m ' + str(m) + ' n ' + str(n)

    if(padding):
        img = add_padding(img,N,M,n,m,pad_value)
        N = N + 2 * int(math.floor(n/2))
        M = M + 2 * int(math.floor(m/2))

    if(padding_copy):
        img = add_padding_copy(img,N,M,n,m)
        N = N + 2 * int(math.floor(n/2))
        M = M + 2 * int(math.floor(m/2))
        
    for i in range(int(n/2),N-int(n/2)):
        tmp = []
        for j in range(int(m/2),M-int(m/2)):
            tmp.append(apply_filter(img, kernel,i,j,N,M,n,m))
        conv_mat.append([x for x in tmp])

    return conv_mat


def add_padding(image, N, M, n, m, pad_value = 0):
    '''
    Helper function to add padding  
    '''
    padded_img  = []
 
    pad_x, pad_y = int(math.floor(n/2)) , int(math.floor(m/2)) # bordure extension lengths
    
    for i in range(pad_x): # add lines
        tmp = []
        for j in range(2*pad_y + M):
            tmp.append(pad_value)
        padded_img.append([x for x in tmp])

    for i in range(N):
        tmp = []
        for j in range(pad_y):
            tmp.append(pad_value) # bordure extension
        for j in range(M):
            tmp.append(image[i][j]) # copy image
        for j in range(pad_y):
            tmp.append(pad_value) # bordure extension
        padded_img.append([x for x in tmp])
            
    for i in range(pad_x): # add lines
        tmp = []
        for j in range(2*pad_y + M):
            tmp.append(pad_value)
        padded_img.append([x for x in tmp])

    return padded_img



def add_padding_copy(image, N, M, n, m):
    '''
    Add a padding copying range at the edge of matrix
    '''
    padded_img  = []
 
    pad_x, pad_y = int(math.floor(n/2)) , int(math.floor(m/2)) # bordure extension lengths
    
    for i in range(pad_x): # add lines
        tmp = []
        for j in range(pad_y):
            tmp.append(0) # bordure extension, 0 in corner
        for j in range(M):
            tmp.append(image[0][j]) # bordure extension, copy first row
        for j in range(pad_y):
            tmp.append(0) # bordure extension, 0 in corner
        padded_img.append([x for x in tmp])

    for i in range(N):
        tmp = []
        for j in range(pad_y):
            tmp.append(image[i][0]) # bordure extension
        for j in range(M):
            tmp.append(image[i][j]) # copy image
        for j in range(pad_y):
            tmp.append(image[i][M-1]) # bordure extension
        padded_img.append([x for x in tmp])
            
    for i in range(pad_x): # add lines
        tmp = []
        for j in range(pad_y):
            tmp.append(0) # bordure extension, 0 in corner
        for j in range(M):
            tmp.append(image[N-1][j]) # bordure extension, copy last row
        for j in range(pad_y):
            tmp.append(0) # bordure extension, 0 in corner
        padded_img.append([x for x in tmp])

    return padded_img




def apply_filter(img, kernel,i,j,N,M,n,m):
    '''
    Apply filter to image
    '''
    res = 0
 
    for k in range(-int(n/2),int(n/2)+1):
        for l in range(-int(m/2),int(m/2)+1):
             
            res += img[i+k][j+l]*kernel[k+int(n/2)][l+int(m/2)]
            
    return res




def sobel_filter(mat):
    '''
    Convolue image with Sobel filter
    '''
    Gx_kernel =[[-1,0,1],[-2,0,2],[-1,0,1]]
    Gy_kernel =[[-1,-2,-1],[0,0,0],[1,2,1]]

    Gx = convolution(mat, Gx_kernel, False, True)
    Gy = convolution(mat, Gy_kernel, False, True)

    NG, MG = len(Gx), len(Gx[0])
    magnitude, theta = [],[]
    
    for i in range(NG):

        tmp_m, tmp_t = [],[]

        for j in range(MG):
         
            tmp_m.append(math.sqrt(Gx[i][j]**2 + Gy[i][j]**2))
            tmp_t.append(math.atan2(Gx[i][j],Gy[i][j]))

        magnitude.append(tmp_m)
        theta.append(tmp_t)

    return magnitude, theta


    
def enhancement_filter(mat):
    '''
    Return enhanced image, obtained computing convolution with Lagrangien filter and substracting the result to original image
    Image is padded with 0 before convolution with enhancement filter.
    '''
    kernel = [[0,-1,0],[-1,5,-1],[0,-1,0]]
    
    return convolution(mat, kernel, False, True)



def gaussian_enhancement(mat, scale, dimx_filter = 5, dimy_filter = 5):
    '''
    Return image after smoothing at given scale followed by enhancement (convolution to gaussian filter)
    '''
    gaussian_filtered = gaussian_filter(mat, dimx_filter, dimy_filter, scale, False, True)

    return enhancement_filter(gaussian_filtered)

    
######## Adaptative filtering

def adaptative_filter_local_stats(mat):
    '''
    Filter image using local statistics, after
    Digital image Enhancement and Noise Filtering by Use of Local Statistics, J.S. Lee, March 1980
    '''
    # filter
    filtered_img = []
    
    N,M = len(mat), len(mat[0])
    
    img = add_padding_copy(mat,N,M,9,9)
    
    N = N + 2 * int(math.floor(9/2))
    M = M + 2 * int(math.floor(9/2))
     
    for i in [X+int(math.floor(9/2)) for X in range(N-int(math.floor(9/2)))]:
        
        tmp = []

        for j in [X+int(math.floor(9/2)) for X in range(M-int(math.floor(9/2)))]:

            var, mean, min_var, min_m  = local_signal_stats(img,i,j,N,M)

            if(min_m != 0):
                xij_m = (mean - min_m)/float(min_m)
            else:
                xij_m = 0

            if(min_var**2 + min_m**2 != 0):
                Qij = (var + mean**2)/float(min_var**2 + min_m**2) - xij_m**2 - min_var**2
            else:
                Qij = 0

            if(xij_m**2 * min_var**2 + min_m**2 * Qij + min_var**2 != 0):
                kij = (min_m * Qij)/float(xij_m**2 * min_var**2 + min_m**2 * Qij + min_var**2)
            else:
                kij = 0

            sij_hat = xij_m + kij * (img[i][j] - xij_m * min_m - min_var)

            tmp.append(sij_hat)
           
        filtered_img.append([x for x in tmp])

    return filtered_img




def local_signal_stats(mat, i, j,N,M):
    '''
    Estimate local mean and variance 
    '''
    m, n = 9, 9 # 9 x 9 window
    
 
    if(i < 0 or i >= N or j < 0 or j >= M):
        print 'here'
        return

    var_x = math.floor(n/2)
    var_y = math.floor(m/2)

    mean, var = 0,0
    
    # mean
    for k in range(int(i - var_x), int(i + var_x + 1)):

        for l in range(int(j - var_y), int(j + var_y + 1)):
            
            if(k >= 0 and k < N and l >= 0 and l < M):

                mean = mean + mat[k][l]

    mean = mean/float(n*m)

    # variance
    for k in range(int(i - var_x), int(i + var_x + 1)):

        for l in range(int(j - var_y), int(j + var_y + 1)):
            
            if(k >= 0 and k < N and l >= 0 and l < M):

                var = var + (mat[k][l] - mean)**2

    var = var/float(n*m-1)

    # noise estimation, subwindow with lowest elevation
    n, m = 3, 3
    var_x = math.floor(n/2)
    var_y = math.floor(m/2)
    
    min_var = sys.maxint
    min_m = sys.maxint
    
    for u, v in [(i-3,j-3),(i,j),(i+3,j+3),(i-3,j),(i-3,j+3),(i,j-3),(i,j+3),(i+3,j-3),(i+3,j)]:
        
        mean_w, var_w = 0, 0 

        for k in range(int(u - var_x), int(u + var_x + 1)):

            for l in range(int(v - var_y), int(v + var_y + 1)):
            
                if(k >= 0 and k < N and l >= 0 and l < M):

                        # mean
                        mean_w = mean_w + mat[k][l]

        mean_w = mean_w / (float) (n*m)

        for k in range(int(u - var_x), int(u + var_x + 1)):

            for l in range(int(v - var_y), int(v + var_y + 1)):
            
                if(k >= 0 and k < N and l >= 0 and l < M):

                        # mean
                        var_w = var_w + (mat[k][l] - mean_w)**2

        var_w = var_w / float (n*m-1)

        if(mean_w < min_m):
            
            min_m = mean_w

        if(var_w < min_var):
            
            min_var = var_w

    print var, mean, min_var, min_m
    
    return var, mean, min_var, min_m

  
 
######## Noise generation

def add_gaussian_noise(mat, mean, sigma):
    '''
    Add gaussian noise with constant mean and sigma  
    '''
    noised = mathlib.deepCopy(mat)
    
    for i in range(len(mat)):

        for j in range(len(mat[0])):

            noised[i][j] += np.random.randn() * sigma + mean

    return noised


def mult_gaussian_noise(mat, mean, sigma):
    '''
    Multiply by gaussian noise with constant mean and sigma 
    '''
    noised = mathlib.deepCopy(mat)
    
    for i in range(len(mat)):

        for j in range(len(mat[0])):

            noised[i][j] = noised[i][j] * (np.random.randn() * sigma + mean)

    return noised


def mult_add_gaussian_noise(mat, mean1, sigma1, mean2, sigma2):
    '''
    Modify matrix with additive and multiplicative gaussian noise 
    '''
    noised = mathlib.deepCopy(mat)
    
    for i in range(len(mat)):

        for j in range(len(mat[0])):

            noised[i][j] = noised[i][j] * (np.random.randn() * sigma1 + mean1) + (np.random.randn() * sigma2 + mean2)

    return noised
    
