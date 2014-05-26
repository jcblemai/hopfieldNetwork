from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
from scipy import stats
import cProfile
from hopfieldNetwork import hopfieldNetwork

tmax = 200

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
                    sum += self.pattern[k+1,i]*self.pattern[k,j]
                sum += self.pattern[0,i]*self.pattern[P-1,j]        # We loop the last pattern w/ the first one.
                self.AsymWeight[i,j] = float(lamda)/self.N * sum
                
    def makePattern(self,P,ratio, lamda):
	self.P = P
        hopfieldNetwork.makePattern(self, P, ratio)
        self.makeAsymmetricWeight(P, lamda)
        
    def initMemory(self):
        self.SPrev = np.ones((self.N, tmax))
            
    def functionG(self, t):
        return (1./self.tau)*heaviside(-t + self.tau)

       
    def dynamicSyncrone(self,t):
        Sbar = np.zeros(self.N)
        for k in range(t):
            Sbar += self.SPrev[:,t-k]*self.functionG(k)

        self.x = sign(np.dot(self.weight, self.x) + np.dot(self.AsymWeight, Sbar))
        
        self.SPrev[:,t] = self.x
        
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
        
        #calculate mean transition time 
        oldFit = -1
        transitionTime = [0]
        
        t = [0]
        overlap = [self.overlap(mu)]
        pixDist = [self.pixelDistance(mu)]
       # x_old = copy(self.x)
        
        for i in range(tmax):
            # run a step
            self.dynamicSyncrone(t[-1])
            
            t.append(i+1)
           # overlap.append(self.overlap(mu)) #Matrice

            for j in range(self.P):
                #print " ", self.overlap(j), "   ",
                if (self.overlap(j) == 1.0):
                    if (oldFit == j):
                        transitionTime[-1] += 1
                    else:
                        oldFit = j
                        #print j,":",transitionTime, " | "                
                        transitionTime.append(0)
            
            pixDist.append(self.pixelDistance(mu))
        return mean(transitionTime)
    
    
def testRun():
    h = hopfieldNetworkAsymmetric(N = 500)
    h.makePattern(P=10, ratio = 0.5, lamda = 1.5)
    h.setTau(tau = 8)
    h.initMemory()
    h.run(flip_ratio=0.1)
    
def varyLambda(lamda):
    h = hopfieldNetworkAsymmetric(N = 500)
    h.makePattern(P=10, ratio = 0.5, lamda=lamda)
    h.setTau(tau = 8)
    h.initMemory()
    h.run(flip_ratio=0)
def transTime(lamba, tu):
	h = hopfieldNetworkAsymmetric(N=500) 
	h.makePattern(P=10, ratio = 0.5, lamda = lamba)
	h.setTau(tau = tu)
	h.initMemory()
	return h.run(flip_ratio=0)	
def varyTau(lamda) : 
	meanT = []
	for t in range(1,20) : 
		meanT.append(transTime(lamda, t) ) 
	print meanT
