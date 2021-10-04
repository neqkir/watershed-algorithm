###########  MATH helper methods

import numpy as np
import math
import pylab as pl
import os, sys
DATA_DIR = os.path.abspath(os.path.join(os.getcwd(), 'IMAGE'))
sys.path.append(DATA_DIR)
import filters

def matrixToVector(mat):
    '''
    Build matrix from vector
    Zero pad if signal matrix has second dimension
    vector unequally lenghted
    '''
    matTmp = adjustMat(len(mat),len(mat[0]),mat)
     
    vector = np.asarray(matTmp).ravel()
     
    return vector.tolist()



def deepCopy(mat):
    '''
    Deep copy for matrix
    '''
    cpy = []
    for i in range(len(mat)):
        tmp = []
        for j in range(len(mat[i])):
            tmp.append(mat[i][j])
        cpy.append([x for x in tmp])

    return cpy


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


            
def mapToPosition(k, size1, size2):
    '''
    Return coordinates in a size1 x size2 matrix for a point
    located at index k in vector
    '''
    row = (int)(math.floor(k/size1))
    col = k%size2
    return (row,col)

   
def adjustMat(size1, size2, mat):
    '''
    Shrink matrix to desired size
    Zero padding
    '''
    if(len(mat) >= size1):
        mat = mat[:size1]
    else: # zero padding
        for i in range(len(mat), size1):
            mat.append([0]*size2)
                
    matTmp = []
    for i in range(size1):
        if(len(mat[i])==size2):
            matTmp.append([x for x in mat[i]])
        elif(len(mat[i])>size2):
            tmp = []
            for j in range(size2):
                tmp.append(mat[i][j])
            matTmp.append([x for x in tmp])
        else:
            tmp = []
            for j in range(len(mat[i])):
                tmp.append(mat[i][j])
            matTmp.append([x for x in tmp]+[0]*(size2-len(mat[i])))

    return matTmp


            
def sortWithIndex(list_):
    idx = [i[0] for i in sorted(enumerate(list_), key=lambda x:x[1])]
    sortedList = [i[1] for i in sorted(enumerate(list_), key=lambda x:x[1])]
    return idx, [sortedList[j] for j in idx]
    
def sort(list_):
    idx = [i[0] for i in sorted(enumerate(list_), key=lambda x:x[1])]
    sortedList = [i[1] for i in sorted(enumerate(list_), key=lambda x:x[1])]
    return [sortedList[j] for j in idx]
    

def planeHalton(n, p2=2):
    '''
    Halton implementation of Quasi-Random Number Generator
    in 2 dimension
    see http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoint.pdf
    Sampling with Hammersley and Halton points,
    Tien-Tsin Wong,Wai-Shing Luk,Pheng-Ann Heng
    '''
    points = []
    pos = 0
    for k in range(n):
        u, p, kk = 0, 0.5, k
        while(kk == 0):
            #print 'first while ' + str(kk)
            if(kk & 1):
                u += p
            p *= 0.5
            kk >>= 1
        v = 0
        ip = 1/float(p2)
        p, kk = ip, k
        while(kk == 0):
            print ' second while ' 
            a = kk % p2
            if(a != 0):
                v += a * p
            p *= ip
            kk = (int) (kk/p2)
        
        points[pos] = u
        pos += 1
        points[pos] = v
        pos += 1

    return points



def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def gaussian_kernel(center_x, center_y, sigma):
    """Returns a gaussian function with the given parameters"""
    sigma = float(sigma)
    return lambda x,y: (1/(2*math.pi*sigma**2)) * np.exp(-1*(1/(2*sigma**2))*((center_x-x)**2+(center_y-y)**2))

def gaussian_oriented(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2+np.random.random()*np.random.randint(-1,2)*((center_x-x)/width_x)*((center_y-y)/width_y))/2)

def minMatrix(mat):
    '''Return indices of the first minimum value in a list of lists.'''
    return min((n, i, j)for i, L2 in enumerate(mat)for j, n in enumerate(L2))[1:]


def minElemsMatrix(mat, nb_elems):
    '''Return indices of the nb_elem first minimum values in a list of lists.'''
    vect = matrixToVector(mat)

    N,M = len(mat), len(mat[0])
    
    sortIdx = [i[0] for i in sorted(enumerate(vect), key=lambda x:x[1])]
     
    return [mapToPosition(i,N,M) for i in sortIdx[:nb_elems]]

def minElemsVector(mat, nb_elems):
    '''Return indices of the nb_elem first minimum values in a list.'''
    N,M = len(mat), len(mat[0])
    
    sortIdx = [i[0] for i in sorted(enumerate(vect), key=lambda x:x[1])]
     
    return [i for i in sortIdx[:nb_elems]]



def euclidean_distance(p1, p2):
 
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
