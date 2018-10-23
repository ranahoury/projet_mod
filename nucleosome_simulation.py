from random import*
import numpy as np
from pylab import*
import matplotlib.pyplot as plt

class Model(object):
    def __init__(self,x=0.002,y=0.002,z=0.002,a=0.002,g=0.002,coopa=0.008,coopm=0.008,eps=0.6,Pdn=10**-4,C=10**-5,K=150):
        """
        initialization function that determines the parameters then computes the mmatrices
        Parameters
        __________
         
        x = base probability to switch from A to U
        y = base probability to switch from U to A
        z = base probability to switch from U to M
        a = base probability to switch from M to U
        g = additional probability to the U->A and U->M transitions
        coopa = the probability we add to the transition M->U  and U->A when neigbourh is acetylated
        coopm = the probability we add to the transition A->U  and U->M when neigbourh is acetylated 
        eps = additional probability for A->M and U->M when there is nucleation
        Pdn = probability in summer for the histone to be nucleated, this parameter has a fixed value given in the article
        C = the maximum probability per sweep with which a locus can become competent to nucleate
        K = effective "dissociation constant" for the time Ͳ dependent probability for a locus to become competent to nucleate
        
        Returns 
        _________
        None
        _________________________
    
        Achille
        """
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
        """Creation of the stochastic matrices of the markov chain
        
        Parameters
        __________
         
        M0 = the matrix for histones with no nucleation and non-depending on the other histones state
        MT = the matrix for histones when we take into account their nucleation
        Mcoopa = the matrix of histones when one or two of their neighboors are 𝖺𝖼𝗍𝗂𝗏𝖺𝗍𝖾𝖽
        Mcoopm is the matrix of histones when one or two of their neighbors are  𝗆𝖾𝗍𝗁𝗒𝗅𝖺𝗍𝖾𝖽
        
        Returns 
        _________
        None
        _________________________
    
        Achille
        """
        self.M0=np.array([[1-self.x, self.x, 0], [self.y+0.5*self.g, 1-self.y-self.z-self.g, self.z+0.5*self.g], [0, self.a, 1-self.a]])
        self.MT=np.array([[-self.eps, 0, self.eps], [0 ,-self.eps ,self.eps], [0, 0, 0]])
        self.Mcoopa=np.array([[0,0,0],[1,-1,0],[0,1,-1]])
        self.Mcoopm=np.array([[-1,1,0],[0,-1,1],[0,0,0]])
        
    def evol(self,nuc,T,Pn):
        """Evolution of the nucleation during the time and depending on the temperature
        
        The first part of evol function allows changing of the nucleation probability during the time (0 enucleated, 1 nucleated), depending on the cold (T==0) or the warm (T==1)
        
        Parameters
        __________
         
        T = the temperature (0 or 1 for cold or warm)
        r = random number between 0 and 1
        Pdn = probability in summer for the histone to be nucleated, this parameter has a fixed value given in the article
        Pn = probability in winter for the histone to be nucleated, its value changes depending on the function "simulation" (depending on the time), explained later.
         
        Returns 
        _________
        None
        _________________________
    
        Cécile
        """
        """Impact of the neighbors on the methylation
        
            The second part of this function changes the probability of methylation/acetylation depending on the neighbors methylation/acetylation : if one or two neighbor(s) of one histones are methylated, this histone has more changes to be methylated. 

            Parameters
            __________
            l = a list used to know the neighbouring histone's rank by removing the current histone'srank from l
            M = the Stochastic matrix for state changes in the current situation
            coopa = the probability we add to the transition M->U  and U->A when neigbourh is acetylated
            coopm = the probability we add to the transition A->U  and U->M when neigbourh is acetylated
            
            Returns 
            _________
            nuc = the histone modified
            _________________________
    
            Cécile & Achille
            """
        for i in range(2):
            if T==0:
                r=random()
                if r<=Pn:
                    nuc[i,1]=1
            if T==1:
                r=random()
                if r<=self.Pdn:
                    nuc[i,1]=0
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
    """ saves a list of states in a text file
        Parameters
        __________
         
        lis =the list
         
        Returns 
        _________
        none
        _________________________
    
        Achille
    """
    np.savetxt(name, np.ndarray.flatten(np.array(lis)), delimiter = ",")
    return
    
def load_nuc(name="list_of_all.txt"):
    """loads a text file generated by save_nuc
        Parameters
        __________
         
        name=the file name
         
        Returns 
        _________
        a list of state
        _________________________
    
        Achille
    """
    lisflat=np.genfromtxt(name, delimiter=",")
    lis=lisflat.reshape((int(len(lisflat)/(35*2*2)),35,2,2))
    return (lis)

def activation(etat):
    """ returns the number of activated histones
        Parameters
        __________
         
        etat =the state on wich we will calculate the number of activated histones 
         
        Returns 
        _________
        act= the number of activated histones 
        _________________________
    
        Achille
    """
    def unlist(L):
        """
        pick up the histones into one list
        """
        M=[]
        for i in L:
            for j in i:
                M+=[j]
        return (M)
    etat=np.array(unlist(etat))
    etat=np.ndarray.tolist(etat)
    act=0
    for nucleation in [0,1]:
        act=act+etat.count([0, nucleation])
    return act

def simulation(winter, spring, pas=10, n=35,L=None,model=Model()): #winter et srping en jours
    """Apply the program for each day and each nucleosome
        
        All the functions above will be applied for each nucleosome, for each winter day (T=0) and then each spring day (T=1). It will first change the probability for histones to be nucleated depending on the time (normally each minute, with 1440 minutes in one day). Then it changes the nucleosomes state at each minute. Finally it will make 2 graphs : one graph with the picture of the gene each "pas" days (i%(1440*pas)), and another for the number of nucleosomes activated each hour (i%*60)
        
        Parameters
        __________
         
        winter = number of days in winter
        spring = number of days in spring 
        n = number of nucleosomes
        pas = to represent the nucleosomes each "pas" days only
        L = represents the status of the histones with their nucleosomes inside
        Pn = probability in winter for the histone to be nucleated, its value changes depending on the time i, C : the maximum probability per sweep with which a locus can become competent to nucleate, and K : effective "dissociation constant" for the time Ͳ dependent probability for a locus to become competent to nucleate
        
         
        Returns 
        _________
        None
        _________________________
    
        Cécile
    """
    if L is None:
        L=np.array([[[0,0],[0,0]] for i in range(n)])
    liste=[np.copy(L)]
    gene=[]
    T=0
    for i in range(winter*1440):
        Pn=(model.C*(i**2)/(model.K*(1440**2)+i**2))
        for j in range(len(L)):
            L[j]= model.evol(L[j],T,Pn)
        if i%(1440*pas)==0:                      #1440 minutes =1 jours
            liste.append(np.copy(L))
        if i%60==0:
            gene.append(activation(L))
    T=1
    for i in range(spring*1440):
        for j in range(len(L)):
            L[j]= model.evol(L[j],T,Pn)
        if i%(1440*pas)==0:
            liste.append(np.copy(L))
        if i%60==0:
            gene.append(activation(L))
    return(liste,gene)