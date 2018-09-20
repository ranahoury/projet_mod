import numpy as np


def save_nuc(lis,name="list_of_all.txt"):
    np.savetxt(name, np.ndarray.flatten(np.array(lis)), delimiter = ",")
    return
    

def load_nuc(name="list_of_all.txt"):    
    lisflat=np.genfromtxt(name, delimiter=",")
    lis=lisflat.reshape((int(len(lisflat)/(35*2*2)),35,2,2))
    return (lis)

def unlist(L):
    M=[]
    for i in L:
        for j in i:
            M+=[j]
    return (M)