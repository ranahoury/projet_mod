from random import*
import numpy as np
from pylab import*
import matplotlib.pyplot as plt

class Model(object):
    def __init__(self,x=0.008,y=0.008,z=0.008,a=0.008,g=0.002,coopa=0.008,coopm=0.008,eps=0.6,Pdn=10**-4,C=10**-5,K=150):
        self.x=x
        self.y=y
        self.z=z
        self.a=a
        self.g=g
        self.coopa=coopa
        self.coopm=coopm
        self.eps=eps
        self.Pdn=Pdn
        self.C=C
        self.K=K
        self.compute_matrices()
    
    def compute_matrices(self):
        self.M0=np.array([[1-self.x, self.x, 0], [self.y+0.5*self.g, 1-self.y-self.z-self.g, self.z+0.5*self.g], [0, self.a, 1-self.a]])
        self.MT=np.array([[-self.eps, 0, self.eps], [0 ,-self.eps ,self.eps], [0, 0, 0]])
        self.Mcoopa=np.array([[0,0,0],[1,-1,0],[0,1,-1]])
        self.Mcoopm=np.array([[-1,1,0],[0,-1,1],[0,0,0]])

    def evol(self,nuc,T,Pn):
        for i in range(2):
            if T==0:
                r=random()
                if r<=Pn:
                    nuc[i,1]=1
            if T==1:
                r=random()
                if r<=self.Pdn:
                    nuc[1]=0
            l=[0,1]
            l.remove(i)
            j=l[0]
            M=self.M0+nuc[i,1]*self.MT
            if nuc[j,0]==2:
                M=self.M0+nuc[i,1]*self.MT+self.coopm*self.Mcoopm
            if nuc[j,0]==0:
                M=self.M0+nuc[i,1]*self.MT+self.coopa*self.Mcoopa
            r=random()
            if r<=M[nuc[i,0],0] :
                nuc[i,0]=0
            elif r<=(M[nuc[i,0],0]+M[nuc[i,0],1]) :
                nuc[i,0]=1
            else :
                nuc[i,0]=2
        return(nuc)

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

def activation(etat):
    etat=np.array(unlist(etat))
    etat=np.ndarray.tolist(etat)
    act=0
    for nucleation in [0,1]:
        act=act+etat.count([0, nucleation])
    return act

def simulation(winter, spring,n=35): #winter et srping en jours
    model=Model()
    L=np.array([[[0,0],[0,0]] for i in range(n)])
    liste=[np.copy(L)]
    gene=[]
    T=0
    for i in range(winter*1440):
        Pn=(model.C*(i**2)/(model.K*(1440**2)+i**2))
        for j in range(len(L)):
            L[j]= model.evol(L[j],T,Pn)
        if i%14400==0:                      #7200 minutes =5 jours
            liste.append(np.copy(L))
        if i%60==0:
            gene.append(activation(L))
    T=1
    for i in range(spring*1440):
        for j in range(len(L)):
            L[j]= model.evol(L[j],T,Pn)
        if i%14400==0:
            liste.append(np.copy(L))
        if i%60==0:
            gene.append(activation(L))
    return(liste,gene)