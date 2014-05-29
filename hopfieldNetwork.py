from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
from scipy import stats
import cProfile

plot_dic={'cmap':cm.gray,'interpolation':'nearest'}

tmax = 100
frac = 0.5 #fraction of neurons in state +1 in a pattern  
class hopfieldNetwork:
    def __init__(self,N):
        self.N = N
    
    def makePattern(self,P,ratio):
        self.pattern = np.random.choice([1, -1], size=(P,self.N), p=[ratio, 1 - ratio])
        self.makeWeight(P)
    
    def makeWeight(self, P):
        self.weight = 1./self.N * (np.fromfunction(lambda i, j: (sum(self.pattern[k,i]*self.pattern[k,j] for k in range(P))), shape=(self.N,self.N), dtype = int).astype(double))

    def overlap(self, mu):
        return ((1./self.N)*sum(self.x*self.pattern[mu]))

    
    def pixelDistance(self, mu):
        return 100*(1-((1./self.N)*sum(self.x*self.pattern[mu])))
    

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
        pixDist = [self.pixelDistance(mu)]
        
        x_old = copy(self.x)
        
        for i in range(tmax):
            #print i
            # run a step
            flip = permutation(arange(self.N))
            for k in flip:
                self.dynamic(k)
            t.append(i+1)
            overlap.append(self.overlap(mu))
            pixDist.append(self.pixelDistance(mu))
            
            # check the exit condition
            i_fin = i+1
            if sum(abs(x_old-self.x))==0:
                return self.pixelDistance(mu)
                break
            x_old = copy(self.x)

        return -1
        
    def runAndPlot(self,mu=0,flip_ratio=0):
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
        pixDist = [self.pixelDistance(mu)]
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
            #print i
            # run a step
            flip = permutation(arange(self.N))
            for k in flip:
                self.dynamic(k)
            t.append(i+1)
            overlap.append(self.overlap(mu))
            pixDist.append(self.pixelDistance(mu))
            
            # update the plotted data
            g1.set_data(self.grid())
            g2.set_data(t,overlap)

            # update the figure so that we see the changes
            draw()
            
            #print "Overlap : " , self.overlap(mu) , ", PixDist : " , self.pixelDistance(mu)

            # check the exit condition
            i_fin = i+1
            if sum(abs(x_old-self.x))==0:
                return self.pixelDistance(mu)
                break
            x_old = copy(self.x)
            time.sleep(2)
            #os.system("pause")
            #print i_fin
        #print 'pattern recovered in %i time steps with final overlap %.3f'%(i_fin,overlap[-1]) 
        return -1

def test():
    h = hopfieldNetwork(200)
    h.makePattern(P=5, ratio = frac)
    h.runAndPlot(flip_ratio=0.3)
    print "test done!"
        
#ex1
def patternRetrieval():
    h = hopfieldNetwork(200)
    file = open('patternRetrieval.dat','w')

    h.makePattern(P=5, ratio = frac)
    for c in np.arange(0.01,0.51,0.01):
        meanError = []
        for i in range(50):
            meanError.append(h.run(flip_ratio=c))
        avr = np.mean(meanError)
        stdDev = stats.sem(meanError)
        strOut = str(c)+"  "+str(avr)+"    "+str(stdDev)+"\n"
        print  c, "  ",avr,"    ",stdDev
        file.write(strOut)
    file.close()

#ex1    
def capacityEstimation():
    file = open('CapacityEstimation.dat','w')
    h = hopfieldNetwork(200)
    for P in range(1,51):
        meanError = []
        for i in range(10):
            h.makePattern(P, ratio = frac)
            meanError.append(h.run(flip_ratio=0.1))
        avr = np.mean(meanError)
        stdDev = stats.sem(meanError)
        strOut = str(P)+"  "+str(avr)+"    "+str(stdDev)+"\n"
        print  P, "  ",avr,"    ",stdDev
        file.write(strOut)
    file.close()

#ex1    
def maxLoad():
    Nlist = [100,250,500]
    for N in Nlist:
        alphaList = []
        h = hopfieldNetwork(N)
        for P in range(1,51):
            avrLoad = []
            for i in range(10):
                meanError = []
                h.makePattern(P, frac)
                for mu in range(P):
                    run = h.run(mu=mu,flip_ratio=0.1)
                    if(run != -1):
                        meanError.append(run)
                    else:
                        print "ERROR NOT CONVERGING"
                avr = np.mean(meanError)
                avrLoad.append(avr)
            if (np.mean(avrLoad) > 2):
                alphaList.append(P/float(N))    
                print "N : ", N, "Pmax : ", P ,"alphaMax : ", P/float(N), "+/- ", stats.sem(alphaList), "\n"
                break
    print "N : ", N, "Alpha_max : ", np.mean(alphaList), "+/- ", stats.sem(alphaList), "\n"
            
    
if __name__ == "__main__":
   cProfile.run('maxLoad()', 'maxLoad.profile')
