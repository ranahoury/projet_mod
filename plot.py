import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

def plot_nuc(L,abscisse,figure):
    """Visualization of H3 histones as 2 half-circles
    
    Function allows to visualize the proportion of acetylated (activated)/methylated (suppressed)/unmodified nucleosomes of the gene of interest             (FLC, in our case)
    
    Parameters
    __________
    
    param1 : L
             A list of histones/nucleosomes containing the information about epigenetic marks (0 - acetylated, 1 - unmodified, 2 - methylated), as well     as whether a histone=>nucleosome is nucleated or not (0 - non-nucleated, 1 - nucleated).
             
    param2 : abscisse
             A value to be specified to determine the heigth of 'the nucleosome line'.
             
    param3 : figure
             The parameter have to be associated with the plotting function => #1.
             
    Returns
    _______
    
    None
    _________________________
    Вивчайте українську мову :)
    Кирило.
    """
    
    plt.rcParams['figure.figsize']=10,10
    #1)figure =plt.figure(0)
    #2)figure.clf()

    for i_n,nucleosome in enumerate(L):
        for i_h,histone in enumerate(nucleosome):
            if i_h==0:
                th1=0
                th2=180
            elif i_h==1:
                th1=180
                th2=360
            if histone[0]==0:
                fc='red'
            elif histone[0]==1:
                fc='grey'
            elif histone[0]==2:
                fc='royalblue'

            w= Wedge(center=[i_n,abscisse], r=.5, theta1=th1, theta2=th2, facecolor=fc)
            figure.gca().add_patch(w)

            if histone[1]==1:
                if i_h==0:
                    c = Circle(xy=[i_n,abscisse+0.5/2],radius=0.2,color='black')
                elif i_h==1:
                    c = Circle(xy=[i_n,abscisse-0.5/2],radius=0.2,color='black')
                figure.gca().add_patch(c)



    figure.gca().axis('equal')
    
    return(None)





        
def histo_nuc(p,fig):   
    """ Plot histogram  
    
    
    
    
    """
    p=np.ndarray.tolist(p)# liste pour count pour histo
    count ={}
    for state in [0,1,2]:
        for nucleation in [0,1]:
            count[(state, nucleation)]= p.count([state, nucleation])
    h = [1,2,3,4,5,6]
    BarName = ['[0, 0]','[0, 1]','[1, 0]','[1, 1]','[2, 0]','[2, 1]']
    fig.clf()
    for j,(state, nucleation) in enumerate(count.keys()):
        height=count[(state, nucleation)]
        absc=h[j]
        if state==0:
            fc='red'
        elif state==1:
            fc='grey'
        elif state==2:
            fc='royalblue'
        ha = "" if nucleation == 0 else "*"
        plt.bar(absc, height, 1.0, color=fc, hatch=ha)
    plt.xticks(h, BarName, rotation=40)
    plt.ylim(0,70)
    fig.canvas.draw()        
    return 



def draw(L):
    p=[]
    fig =plt.figure(0)
    for i,l in enumerate(L):
        p=unlist(l)
        p=np.array(p)
        #np.savetxt( "p"+str(i)+".txt",p, delimiter=",")
        histo_nuc(p,fig)
        fig.savefig('SimpleBar'+str(i)+'.png')
    return
    