##########################################################################################
##########################################################################################
####
####             Initial peak detection to seed WATERSHED algorithm with markers
#### 
##########################################################################################
 
 

 
import math
 
import numpy as np
import random
import itertools as it
import operator
from operator import itemgetter
import os, sys
DATA_DIR = os.path.abspath(os.path.join(os.getcwd(), 'IMAGE'))
sys.path.append(DATA_DIR)
import CC_labeling_8, mathlib, filters


 
#############################  First method for rough peak detection:: one D peak detection selecting local max without taking noise into account
    
def roughPeakDetect(size1,size2,signal):
    """
    Returns list of roughly identified peaks in signal 
    """

    i = 1
    peaks = []

    if(signal[0] > max(signal[1], signal[2])):
        peaks.append(0)
    
    size = len(signal)
    
    while(i < size-2):

        if(signal[i-1] < signal[i] and signal[i+1] < signal[i] ): # apex
            jmax, j, jminleft = i, i-2, i-1

            # look for min left
            while((j>0) and (signal[j] < signal[jminleft]+ 10**(-5))): # + thshld * localNoiseMin)):
                if(signal[j] < signal[jminleft]):
                    jminleft = j
                j = j-1

            j,jminright = i+2, i+1

            # look for min right
            while((j < size-1) and (signal[j] < signal[jminright]+ 10**(-5))): #+ thshld * localNoiseMin)):
                if(signal[j] < signal[jminright]):
                    jminright = j
                j = j+1
     
            # valley - valley
            valleysDiff = abs(signal[jminleft] - signal[jminright])
            # minimum max - valley difference
            maxDiff =  min(abs(signal[jmax] - signal[jminright]), abs(signal[jmax] - signal[jminleft]))

            if(valleysDiff < 2.0 * maxDiff):
                peaks.append(jmax)
                
            i = jminright + 1

        else:
            i = i+1

            
    return peaks




def selectMarkers(size1,size2,signal):
    """
    Return refined list of peaks in signal from list obtained with roughPeakDetect
    Examination of first derivative
    """
    # print signal
    rPeaks = roughPeakDetect(size1,size2,signal) # first list of potential peaks (indices)
    i = 1 
    peaks = []
    #print 'rpeaks ' + str(rPeaks)
    while(i < len(rPeaks) - 2):
    
        if(signal[rPeaks[i-1]] < signal[rPeaks[i]] and signal[rPeaks[i+1]] < signal[rPeaks[i]] ): # apex
            peaks.append(rPeaks[i])
        i = i + 1         
    return [mathlib.mapToPosition(p,size1,size2) for p in peaks]


  
def maximum_filter_maxima(mat):
    """
    Return maxima as obtained with maximum filter
    """
    # apply maximum_filter to all cells in matrix
    filtered_img, maxima = filters.maximum_filter(mat)
        
    return maxima     


 
def getValley(mat, list_peaks, peak):
    '''
    Find max and mean valley for peak. Five closest neighbors of peak in list_peaks are selected.
    Minimum is looked after in rectangle regions separating peak and a nearest peak.
    Return highest valley and mean valley over five values.
    '''
    nearest_neighbors = []
    distances = []

    for npeak in list_peaks:
        distances.append(mathlib.euclidean_distance(npeak, peak))

    nearest_neighbors = [list_peaks[k] for k in mathlib.minElemsVector(distance,5)]

    valleys = []

    px, py = peak

    # valley
    for nx,ny in nearest_neighbors:
        
        valleys.append(minMatrix(mat[min(px,nx):(min(px,nx)+math.abs(nx-px)+1)][min(py,ny):(min(py,ny)+math.abs(ny-py)+1)])) 
        
    return min(valleys), np.mean(valley)

    
 
 

def removeBackgroundNoise(mat, dimx_support, dimy_support, k, list_peaks):
    '''
    Background removal, V1.
    At each window position, compute local noise and record remaining peaks. 
    Windows have k times the dimensions of support. Support is the minimum size rectangle containing the largest blob.
    Background removal occurs after first watershed determined some blob statistics, as blob support.
    For final decision, np.unique(remaining_peaks) should be kept
    . listPeaks contains list of coordinates for roughly detected peaks
    . mat is elevation matrix
    . window for local noise calculation is k * (dimensions for support)
    '''
    N,M = len(mat), len(mat[0])
    
    window_x, window_y = int(math.floor(k * dimx_support)), int(math.floor(k * dimy_support))

    

    peaks_pos = []
    for p in list_peaks:
        if(p != 0):
            peaks_pos.append(p[1])
         
    print ('peaks_pos ' + str(peaks_pos))
            
    localNoise = [] # store local noise
    remaining_peaks = [] # store remaining peaks
    print ('window_x ' + str(window_x), ', window_y ' + str(window_y))
    
    wi = 0
    while(wi < N - window_x - 1):

        wj = 0
        
        while(wj < M - window_y - 1):

            # at given sub-window, compute local noise
            
            pixels_noise_compute = [] # pixels to be taken into account for noise evaluation
   
            sub_window = [X[wj:int(wj + window_y + 1)] for X in mat[wi:int(wi + window_x + 1)]]

            #print 'sub_window ' + str(sub_window)
            
            five_smallest = mathlib.minElemsMatrix(sub_window, 5)
    
            #print 'five_smallest ' + str(five_smallest)

            for k,l in five_smallest: # add pixels to compute local noise
                
                if(k > 0):
                    pixels_noise_compute = pixels_noise_compute + [(k-1,l)]
                    
                if(k < window_x-1):
                    pixels_noise_compute = pixels_noise_compute + [(k+1,l)]
                    
                if(l > 0):
                    pixels_noise_compute = pixels_noise_compute + [(k,l-1)] 
                    
                if(l < window_y-1):
                    pixels_noise_compute = pixels_noise_compute + [(k,l+1)]
                    
            if(pixels_noise_compute == []): # something went wrong
                return
            
            vals = [sub_window[x][y] for x,y in pixels_noise_compute]
            
            mean_noise, var_noise = np.mean(vals), np.std(vals)

            localNoise.append((mean_noise, var_noise))
            
            # find list of peaks in this sub-window
             
            lpoints_window = list(it.product(range(wi,(wi + window_x + 1)),range(wj,(wj + window_y + 1))))
            
            lpeaks_window = set(lpoints_window).intersection(set(peaks_pos))
            
            #print 'lpeaks window ' + str(lpeaks_window)
            
            for pi,pj in lpeaks_window:
                
                if(mat[pi][pj] >= mean_noise + 3 * var_noise):
                    
                    remaining_peaks.append((pi,pj))
                    
            #print 'local Noise ' + str(localNoise)
            
            #print 'remaining peaks ' + str(remaining_peaks)

            wj = wj + int(math.floor(window_y/2)) # half window overlap

        wi = wi + int(math.floor(window_x/2)) # half window overlap
   
    return localNoise, remaining_peaks


def gaussian_enhancementV1(mat, noise = None, nb_background_points=None):
    '''
    Gaussian enhancement with a variance for gaussian smoothing that is roughly evaluated
    for the whole matrix. this noise corresponds approximately to a background noise.
    . nb_background_points: number of points to take into account to find background noise mean
    and variance
    . if noise is None, filter is adaptative
    '''
    N,M = len(mat), len(mat[0])

    filtered = mathlib.deepCopy(mat)
    
    ## Roughly estimate background noise
    
    pixels_noise_compute = [] # pixels to be taken into account for noise evaluation
    
    if(nb_background_points == None):
        int(math.floor((N * M)/200)) # 1/200 smallest points used for background noise rough evaluation
        
    smallest = mathlib.minElemsMatrix(mat, nb_background_points) # pick up 

    for k,l in smallest: # add pixels to compute local noise
        
        if(k > 0):
            pixels_noise_compute = pixels_noise_compute + [(k-1,l)]
            
        if(k < N-1):
            pixels_noise_compute = pixels_noise_compute + [(k+1,l)]
            
        if(l > 0):
            pixels_noise_compute = pixels_noise_compute + [(k,l-1)] 
            
        if(l < M-1):
            pixels_noise_compute = pixels_noise_compute + [(k,l+1)]
            
    if(pixels_noise_compute == []): # something went wrong
        return
    
    vals = [mat[x][y] for x,y in pixels_noise_compute]
    
    mean_noise, var_noise = np.mean(vals), np.std(vals)

    ## all points less than mean_noise are background
    for i in range(N):
        for j in range(M):
            if(mat[i][j] < mean_noise):
                mat[i][j] = 0
                
    ## Enhance and sharpenize

    if(noise == None):
        noise = var_noise
        
    filtered = filters.gaussian_enhancement(mat, noise)

    ## all points less than mean_noise are background
    for i in range(N):
        for j in range(M):
            if(filtered[i][j] < mean_noise):
                filtered[i][j] = 0
                
    return filtered



##def gaussian_enhancementV2(mat):
##    '''
##    Gaussian enhancement with a variance for gaussian smoothing that is roughly evaluated
##    locally for subwindows in the input matrix.
##    this noise corresponds approximately to a local background noise.
##    '''
##    N,M = len(mat), len(mat[0])
## 
##    pixels_noise_compute = [] # pixels to be taken into account for noise evaluation
##
##    smallest = mathlib.minElemsMatrix(sub_window, 5)
##
##    print 'five_smallest ' + str(five_smallest)
##
##    for k,l in five_smallest: # add pixels to compute local noise
##        
##        if(k > 0):
##            pixels_noise_compute = pixels_noise_compute + [(k-1,l)]
##            
##        if(k < window_x-1):
##            pixels_noise_compute = pixels_noise_compute + [(k+1,l)]
##            
##        if(l > 0):
##            pixels_noise_compute = pixels_noise_compute + [(k,l-1)] 
##            
##        if(l < window_y-1):
##            pixels_noise_compute = pixels_noise_compute + [(k,l+1)]
##            
##    if(pixels_noise_compute == []): # something went wrong
##        return
##    
##    vals = [sub_window[x][y] for x,y in pixels_noise_compute]
##    
##    mean_noise, var_noise = np.mean(vals), np.std(vals)
##
##    localNoise.append((mean_noise, var_noise))
##    
##    # find list of peaks in this sub-window
##     
##    lpoints_window = list(it.product(range(wi,(wi + window_x + 1)),range(wj,(wj + window_y + 1))))
##    
##    lpeaks_window = set(lpoints_window).intersection(set(peaks_pos))
##    
##    print 'lpeaks window ' + str(lpeaks_window)
##    
##    for pi,pj in lpeaks_window:
##        
##        if(mat[pi][pj] >= mean_noise + 3 * var_noise):
##            
##            remaining_peaks.append((pi,pj))
##            
##    print 'local Noise ' + str(localNoise)
##    
##    print 'remaining peaks ' + str(remaining_peaks)
##
##    wj = wj + int(math.floor(window_y/2)) # half window overlap
##
##        wi = wi + int(math.floor(window_x/2)) # half window overlap
##   
##    return localNoise, remaining_peaks
