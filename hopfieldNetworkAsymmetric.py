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
        self.AsymWeight = float(lamda)/self.N * (np.fromfunction(lambda i, j: (sum(self.pattern[k+1,i]*self.pattern[k,j] for k in range(P-1)) + self.pattern[0,i]*self.pattern[P-1,j]), shape=(self.N,self.N), dtype = int).astype(double))
        
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
            #print ""
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
        return mean(transitionTime[1:]) #We don't take the fist element of the list : because depend of the flip of start
    
    
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
def transTime():
    file = open('transitionTime.dat','w')
    h = hopfieldNetworkAsymmetric(N=500) 
    lamdaList = [1.3, 1.7, 2.2, 2.5]
    tauList = range(1,20)
    for  lamdar in lamdaList:
        for tu in tauList:
            run = []
            for r in range(5): #We to 5 run to be sure
                h.makePattern(P=10, ratio = 0.5, lamda = lamdar)
                h.setTau(tau = tu)
                h.initMemory()
                run.append(h.run(flip_ratio=0))
                print run
            strOut = str(lamdar)+"  "+str(tu)+"    "+str(mean(run))+"\n"
            print  lamdar, "  ",tu,"    ", mean(run)
            file.write(strOut)
        file.write("\n\n")
    
    file.close()

def varyTau(lamda) : 
	meanTPerRun = []
	meanT = []
	for t in range(1,20) : 
	#for each value of tau, run simulation 5 times and take average transition time 
		for r in range(5) : 
			meanTPerRun.append(transTime(lamda, t) ) 
		meanT.append(mean(meanTPerRun)) 
		meanTPerRun = []	
		print t, " ", meanT[t-1]
