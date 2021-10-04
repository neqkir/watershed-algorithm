##########################################################################################
##########################################################################################
####
####             IO
#### 
##########################################################################################
 
import numpy as np
import mathlib

def imread(path):
    
    try:
        from PIL import Image
    except ImportError:
        raise ImportError("Could not import the Python Imaging Library (PIL)"
                          " required to load image files.  Please refer to"
                          " http://pypi.python.org/pypi/PIL/ for installation"
                          " instructions.")
    
    fp = open(path, "rb")
    im = Image.open(fp)
 
    result = np.array(im)
    fp.close()
    return result  



def matToFile(path, mat):
    '''Write matrix to file
    '''
    N, M = len(mat), len(mat[0])

    vec = mathlib.matrixToVector(mat)
    
    f = open(path, 'w')
    
    f.write(str(N) + "\n")
    f.write(str(M) + "\n")
    
    for v in vec:
        f.write(str(v)+"\n")


def vecToFile(path, vec):
    '''Write matrix to file
    '''
    N = len(vec) 
 
    f = open(path, 'w')
    
    f.write(str(N) + "\n")
    
    for v in vec:
        f.write(str(v)+"\n")

        
def fileToMatrix(path):
    '''Build matrix from data in file
    '''
    f = open(path, 'r')
    
    N = int(f.readline().split("\n")[0])
    M = int(f.readline().split("\n")[0])
    
    print N
    print M

    vec = []

    for line in f.readlines():
        vec.append(float(line.split("\n")[0]))
        
    vec = list(vec)    

    mat = np.reshape(vec,(N,M))

    return mat.tolist()


def fileToVector(path):
    '''Build matrix from data in file
    '''
    f = open(path, 'r')
    
    N = int(f.readline().split("\n")[0])
 
    print N
 
    vec = []

    for line in f.readlines():
        vec.append(float(line.split("\n")[0]))

    return vec  
