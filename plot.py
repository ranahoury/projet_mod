import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

def plot_nuc(L):
    plt.rcParams['figure.figsize']=10,10
    figure =plt.figure()
    figure.clf()

    for i_n,nucleosome in enumerate(L):
        for i_h,histone in enumerate(nucleosome):
            if i_h==0:
                th1=0
                th2=180
            elif i_h==1:
                th1=180
                th2=360
            if histone[0]==0:
                fc='royalblue'
            elif histone[0]==1:
                fc='grey'
            elif histone[0]==2:
                fc='red'

            w= Wedge(center=[i_n,0], r=.5, theta1=th1, theta2=th2, facecolor=fc)
            figure.gca().add_patch(w)

            if histone[1]==1:
                if i_h==0:
                    c = Circle(xy=[i_n,0.5/2],radius=0.2,color='black')
                elif i_h==1:
                    c = Circle(xy=[i_n,-0.5/2],radius=0.2,color='black')
                figure.gca().add_patch(c)



    figure.gca().axis('equal')
    
    return(None)