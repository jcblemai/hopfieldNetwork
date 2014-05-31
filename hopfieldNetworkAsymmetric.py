from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
from scipy import stats
from hopfieldNetwork import hopfieldNetwork

#Maximum time step
tmax = 400

def heaviside(x):
    """
    Implentation of the mathematical heaviside function
    """
    if x == 0:
        return 0.5
    return 0 if x < 0 else 1

class hopfieldNetworkAsymmetric(hopfieldNetwork):
    """ 
    Implement an assymetric hopfield network, in the scope of a project for the course Biological Modelling 
    of Neural Network by W. Gerstner, at EPFL. This network is used to retrieve sequence of pattern.
    Inherit from a standart neuron network.
    """
    def __init__(self, N):        
        """ 
        Create a network of size N
        """
        hopfieldNetwork.__init__(self, N)
        
    def makeAsymmetricWeight(self, P, lamda):
        """
        Build the N*N weight matrix according to the theorical formula, in a fast way with numpy
        """
        # The lambda expression return the sum we need for the weight. We loop the last pattern with the first one
        #in this expression. It gets in argument the indices i and j by numpy (so we need the matrix
        # to be int). We convert it as a double so we can divide by 1/N. 
        self.AsymWeight = float(lamda)/self.N * (
            np.fromfunction(lambda i, j: (sum(self.pattern[k+1,i]*self.pattern[k,j] for k in range(P-1)) + self.pattern[0,i]*self.pattern[P-1,j]), shape=(self.N,self.N), dtype = int).astype(double)
                                                )
        
    def makePattern(self,P,ratio, lamda):
        """
        Make the stored patterns and the symmetric and assymetric weight matrix.
        """
        self.P = P
        hopfieldNetwork.makePattern(self, P, ratio) #For symmetric weight and the creation of the patterns
        self.makeAsymmetricWeight(P, lamda)
        
    def initMemory(self):
        """
        Initialisation of the memory for the calculation of Sbar : We will memorize
        all the previous states of the network. We initialize with one so we don't 
        influanciate the creation of Sbar (with zeros, the network will stays at zero).
        """
        #fixed size array because usually simuations are done until the end, ie t = tmax
        self.SPrev = np.ones((self.N, tmax))
            
    def functionG(self, t):
        """
        Our filter function
        """
        return (1./self.tau)*heaviside(-t + self.tau)

       
    def dynamicSyncrone(self,t):
        """
         Run a step of the dynamic : We change the state of all the pixel in one
         operation (syncronyous update).
         """
        #Calculate the new Sbar with a sum
        Sbar = np.zeros(self.N)
        for k in range(t):
            Sbar += self.SPrev[:,t-k]*self.functionG(k)
        
        #Vectorialy change the state of the whole network
        self.x = sign(np.dot(self.weight, self.x) + np.dot(self.AsymWeight, Sbar))
        self.SPrev[:,t] = self.x
        
    def setTau(self, tau):
        """
        Set the parameter tau of the filter function
        """
        self.tau = tau
    
    def run(self):
        """
        Run the full time evolution of the network. We start from the first pattern stored.
        We log the data and return the mean transition time
        """

        # set the initial state of the network as pattern 0
        self.x = copy(self.pattern[0])
        

        # Last pattern that was retrieved, to estimate the transition time
        oldFit = 0        
        
        #Logged data : time, and transition time
        t = [0]
        transitionTime = [0]
        
        for i in range(tmax):
            # run a step
            self.dynamicSyncrone(t[-1])
            
            t.append(i + 1)
            print ""
            #Check the overlap for all patterns
            for j in range(self.P):
                print j," ", self.overlap(j), "   ",
                
                #If we are on a pattern
                if (self.overlap(j) == 1.0):
                    #and it's the same than last time : we increase the transitionTime
                    if (oldFit == j):
                        transitionTime[-1] += 1
                    #else we reset the transition time and are now evaluating another pattern
                    else:
                        oldFit = j  
                        transitionTime.append(1)
            
        return mean(transitionTime[1:]) #We don't take the fist element of the list : because depend of the flip of start
    
    
def test():
    """ 
    Test function, used as a test to see if the network still work after changes in the code 
    """
    h = hopfieldNetworkAsymmetric(N = 500)
    h.makePattern(P=5, ratio = 0.5, lamda = 1.5)
    h.setTau(tau = 8)
    h.initMemory()
    h.run()
    
#Ex 2.2
def varyLambda(lamda):
    """ 
    Since the behaviour we want to observe (sequencial evolution) is not easily mesurable,
    this function enable the testing by hand to find the value of lambda for which the network
    work
    """
    h = hopfieldNetworkAsymmetric(N = 500)
    h.makePattern(P=10, ratio = 0.5, lamda=lamda)
    h.setTau(tau = 8)
    h.initMemory()
    h.run()
    
#Ex 2.3
def transTime():
    """
    Automatically evaluate the mean transition time as a function of
    lambda and tau. Store the result in a gnuplot compliant format
    """
    file = open('transitionTime.dat','w')
    h = hopfieldNetworkAsymmetric(N = 500)
    
    # Our data range
    lamdaList = [1.3, 1.7, 2.2, 2.5]
    tauList = range(1,25)
    for  lamdar in lamdaList:
        for tu in tauList:
            #We do 5 run for each set of parameters to have statistical average of the transition time
            run = []
            for r in range(5): 
                h.makePattern(P = 10, ratio = 0.5, lamda = lamdar)
                h.setTau(tau = tu)
                h.initMemory()
                run.append(h.run())
            
            #logging data
            strOut = str(lamdar) + "  " + str(tu) + "    " + str(mean(run)) + "\n"
            print  lamdar, "  ",tu,"    ", mean(run)
            file.write(strOut)
        
        file.write("\n\n")  #to separate run with different lambda
    
    file.close()

#Optional exercise 2.4
def storageCapacity(nbPattern, N = 500):
    """
    Since the behaviour we want to observe (sequencial evolution) is not easily mesurable,
    this function enable the testing by hand to find the maximal value of P for which the
    network behave correctly
    """
    
    h = hopfieldNetworkAsymmetric(N)
    h.makePattern(P=nbPattern, ratio = 0.5, lamda = 1.5)
    h.setTau(tau = 8)
    h.initMemory()
    h.run()
    
    print "This alpha is ", float(nbPattern)/N
