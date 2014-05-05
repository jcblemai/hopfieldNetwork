from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
plot_dic={'cmap':cm.gray,'interpolation':'nearest'}

tmax = 100

class hopfieldNetwork:
    def __init__(self,N):
        self.N = N
    
    def makePattern(self,P=1,ratio=0.5):
        """
        """
        self.pattern = -np.ones((P,self.N))
        idx = int(ratio*self.N)
        for i in range(P):
            self.pattern[i,:idx] = 1
            self.pattern[i] = np.random.permutation(self.pattern[i])
        
        self.weight = np.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                sum = 0
                for k in range(P):
                    sum += self.pattern[k,i]*self.pattern[k,j]
                self.weight[i,j] = 1./self.N * sum
                

    def overlap(self, mu):
        return ((1./self.N)*sum(self.x*self.pattern[mu]))

    def pixDistance(self, mu) : 
	return 1-overlap(self, mu)

    def grid(self,mu=None):
        if mu is not None:
            x_grid = reshape(self.pattern[mu],(self.N,1))
        else:
            x_grid = reshape(self.x,(self.N,1))
        return x_grid
    
    def dynamic(self,i):
        """
        """

        h = sum(self.weight[i]*self.x)
        self.x[i] = sign(h)
    
    
    def run(self,mu=0,flip_ratio=0):
        """
        """
        try:
            self.pattern[mu]
        except:
            raise IndexError, 'pattern index too high'
        
        # set the initial state of the net
        self.x = copy(self.pattern[mu])
        flip = permutation(arange(self.N))
        idx = int(self.N*flip_ratio)
        self.x[flip[0:idx]] *= -1
        
        t = [0]
        overlap = [self.overlap(mu)]
        print overlap
        # prepare the figure
        figure()
        
        # plot the current network state
        subplot(221)
        g1 = imshow(self.grid(),**plot_dic)# we keep a handle to the image
        axis('off')
        title('x')
        
        # plot the target pattern
        subplot(222)
        imshow(self.grid(mu=mu),**plot_dic)
        axis('off')
        title('pattern %i'%mu)
        
        # plot the time course of the overlap
        subplot(212)
        g2, = plot(t,overlap,'k',lw=2) # we keep a handle to the curve
        axis([0,tmax,-1,1])
        xlabel('time step')
        ylabel('overlap')
        
        # this forces pylab to update and show the fig.
        draw()
        x_old = copy(self.x)
        
        for i in range(tmax):
            print i
            # run a step
            flip = permutation(arange(self.N))
            for k in flip:
                self.dynamic(k)
            t.append(i+1)
            overlap.append(self.overlap(mu))
            
            # update the plotted data
            g1.set_data(self.grid())
            g2.set_data(t,overlap)

            # update the figure so that we see the changes
            draw()

            # check the exit condition
            i_fin = i+1
            if sum(abs(x_old-self.x))==0:
                break
            x_old = copy(self.x)
            time.sleep(1)
            #os.system("pause")
            print i_fin
        print 'pattern recovered in %i time steps with final overlap %.3f'%(i_fin,overlap[-1])
        
def test():
    h = hopfieldNetwork(200)
    h.makePattern(P=5)
    h.run(flip_ratio=0.2)
        
#if __name__ == "__main__":

