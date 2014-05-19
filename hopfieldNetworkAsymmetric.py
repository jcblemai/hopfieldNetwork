from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
from scipy import stats
import cProfile
from hopfieldNetwork import hopfieldNetwork

tmax = 20

def heaviside(x):
    if x == 0:
        return 0.5

    return 0 if x < 0 else 1

class hopfieldNetworkAsymmetric(hopfieldNetwork):
    def __init__(self, N):
        hopfieldNetwork.__init__(self, N)
        
    def makeAsymmetricWeight(self, P, lamda):
        self.AsymWeight = np.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                sum = 0
                for k in range(P-1):
                    sum += self.pattern[k,i]*self.pattern[k+1,j]
                sum += self.pattern[0,i]*self.pattern[P-1,j]        # We loop the last pattern w/ the first one.
                self.AsymWeight[i,j] = float(lamda)/self.N * sum
                
    def makePattern(self,P,ratio, lamda):
        hopfieldNetwork.makePattern(self, P, ratio)
        self.makeAsymmetricWeight(P, lamda)
        
    def initMemory(self):
        self.SPrev = []
        for i in range(self.N):
            self.SPrev.append([1])
            
    def functionG(self, t):
        return (1./self.tau)*heaviside(-t + self.tau)
    
    def dynamic(self,i,t):
        newSbar =  []
        for j in range(self.N):
            Sbar = 0
            for k in range(len(self.SPrev[j])):
                Sbar += self.SPrev[j][t-k]*self.functionG(k)
            newSbar.append(Sbar)
            
        #print newSbar
        
        h = 0
        for l in range(self.N):
            h += self.weight[i][l]*self.x[l]+ self.AsymWeight[i][l]* newSbar[l]
        self.x[i] = sign(h)
        self.SPrev[i].append( self.x[i])
       
        
    def setTau(self, tau):
        self.tau = tau
    
    def run(self,mu=0,flip_ratio=0):
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
            # run a step
            flip = permutation(arange(self.N))
            for k in flip:
                self.dynamic(k,t[-1]) # Syncronous thing TODO
            
            t.append(i+1)
            overlap.append(self.overlap(mu)) #Matrice
            
            for j in range(10):
                print j, " ", self.overlap(j), "; ",
            print ""
            
            pixDist.append(self.pixelDistance(mu))
            
            # check the exit condition
            i_fin = i+1
            x_old = copy(self.x)
            
        #print 'pattern recovered in %i time steps with final overlap %.3f'%(i_fin,overlap[-1]) 
        return -1
    
    
def testRun():
    h = hopfieldNetworkAsymmetric(N = 500)
    h.makePattern(P=10, ratio = 0.5, lamda = 1.5)
    h.setTau(tau = 8)
    h.initMemory()
    h.run(flip_ratio=0.3)
    
#if __name__ == "__main__":
    