#Authors : Miryam Chaabouni and Joseph Lemaitre
#Inspired from a work found on moodle, by L.Ziegler

from PIL import Image
from copy import copy
import time
from pylab import *
import numpy as np
import os
from scipy import stats

#For profiling to enhance the speed of the code
import cProfile

plot_dic={'cmap':cm.gray,'interpolation':'nearest'}

#Maximum time step
tmax = 100

#fraction of neurons in state (+1) in the created patterns
frac = 0.5 


class hopfieldNetwork:
    """ 
    Implement a standart hopfield network, in the scope of a project for the course Biological Modelling 
    of Neural Network by W. Gerstner, at EPFL.
    """
    def __init__(self,N):
        """ 
        Create a network of size N
        """
        self.N = N
    
    def makePattern(self,P,ratio):
        """
        Make the stored patterns and the symmetric weight matrix.
        """
        #We create all the patterns in one line : it fill a matrix of size (nbPatern * nbNeurons) with a proportion ratio of 1, the others are -1.
        self.pattern = np.random.choice([1, -1], size=(P,self.N), p=[ratio, 1 - ratio])
        self.makeWeight(P)
    
    def makeWeight(self, P):
        """
        Build the N*N weight matrix according to the theorical formula, in a fast way with numpy
        """
        
        # The lambda expression return the sum we need for the weight. It gets in argument the indices i and j by numpy (so we need the matrix
        # to be int). We convert it as a double so we can divide by 1/N. 
        self.weight = 1./self.N * (
            np.fromfunction(lambda i, j: (sum(self.pattern[k,i]*self.pattern[k,j] for k in range(P))), shape=(self.N,self.N), dtype = int).astype(double)
                                  )

    def overlap(self, mu):
        """
        Calculate the overlap between the pattern mu and the current state of the network.
        """
        return ((1./self.N)*sum(self.x*self.pattern[mu]))

    
    def pixelDistance(self, mu):
        """
         Calculate the pixel distance between the pattern mu and the current state of the network.
        """
        return 100*(1-((1./self.N)*sum(self.x*self.pattern[mu])))
    
    def dynamic(self,i):
        """
        Run a step of the dynamic : change the state of one pixel since we do asyncronyous
        update.
        """
        h = sum(self.weight[i]*self.x)
        self.x[i] = sign(h)
    
    
    def run(self,mu=0,flipRatio=0):
        """
        Run the full time evolution of the network. We start from a deviation of the pattern
        mu with flipRatio changes. We log the data and return the pixel distance at convergence.
        """
        #If the pattern exist.
        try:
            self.pattern[mu]
        except:
            raise IndexError, 'pattern index too high'
        
        # set the initial state of the network
        self.x = copy(self.pattern[mu])
        flip = permutation(arange(self.N))  #permutation of number from one to N
        idx = int(self.N*flipRatio)         #We will flip the N*ratio neurons that have indices in the start of the permutation
        self.x[flip[0:idx]] *= -1           #Perform the flip
        
        # Logged data : time, overlap and pixelDistance with pattern mu
        t = [0]
        overlap = [self.overlap(mu)]
        pixDist = [self.pixelDistance(mu)]
        
        
        x_old = copy(self.x)
        
        #Run the network
        for i in range(tmax):
            #We perform the dynamic for all the pixels individually (so changes in other neurons are accounted), in a random order
            flip = permutation(arange(self.N))
            for k in flip:
                self.dynamic(k)
                
            #Update our data list
            t.append(i+1)
            overlap.append(self.overlap(mu))
            pixDist.append(self.pixelDistance(mu))
            
            #Optional printing
            #print "Overlap : " , self.overlap(mu) , ", PixDist : " , self.pixelDistance(mu)
            
            # check the exit condition : There is no change with us and the previous version. We have reach a fixed point.
            if sum(abs(x_old - self.x)) == 0:
                return self.pixelDistance(mu)
                break
            
            # Copy of our state for the convercence criterion
            x_old = copy(self.x)

        # In case of error
        return -1 

def test():
    """ 
    Test function, used as a test to see if the network still work after changes in the code 
    """
    h = hopfieldNetwork(200)
    h.makePattern(P = 5, ratio = frac)
    h.run(flipRatio = 0.3)
    print "test done!"
        
#Ex 1.2
def patternRetrieval():
    """ 
    Log the error in the retriaval of the network for a different changes in the start state : we change
    c, the ratio of pixels flipped for the pattern we want to retrieve (so the distance from the pattern we want to
    retrieve)
    """
    
    file = open('patternRetrieval.dat','w') #We print the results in a gnuplot ready format.
    
    # Creation of our network
    h = hopfieldNetwork(200)
    h.makePattern(P = 5, ratio = frac)
    
    for c in np.arange(0.01, 0.51, 0.01):
        meanError = []
        for i in range(50):                          # The code is fast, so we do the average over 50 run for each c. This gives us a low statistical error
            meanError.append(h.run(flipRatio = c))   #since the run function return the pixel distance at convergence, we log this value
        
        # Statistical values calculation and logging
        avr = np.mean(meanError)
        stdDev = stats.sem(meanError)
        strOut = str(c) + "  " + str(avr) + "    " + str(stdDev) + "\n"
        print  c, "  ",avr,"    ",stdDev
        file.write(strOut)
    
    file.close()

#Ex 1.3 
def capacityEstimation():
    """
    This function generate the data to plot the mean retrieval error against P. It's not
    a report requirement so this graph
    """
    #gnuplot ready logging
    file = open('CapacityEstimation.dat','w')
    h = hopfieldNetwork(200)
    
    for P in range(1, 51):              #We try for different number of pattern stored
        meanError = []
        for i in range(10):             #10 trials to have good statistical error
            h.makePattern(P, ratio = frac)
            meanError.append(h.run(flipRatio=0.1))
        
        #logging
        avr = np.mean(meanError)
        stdDev = stats.sem(meanError)
        strOut = str(P) + "  " + str(avr) + "    " + str(stdDev) + "\n"
        print  P, "  ", avr, "    ", stdDev
        file.write(strOut)
    file.close()

#Ex 1.3    
def maxLoad():
    """
    This code automaticaly find the max load of an hopfield Network, for N = 100, 250 and 500.
    It does the calculation several time by network to have statistically relevant results.
    """
    Nlist = [100, 250, 500]
    # For every network ... 
    for N in Nlist:
        alphaList = []
        h = hopfieldNetwork(N)
        start = int(0.1 * N)                  #We choose a starting point not too far from the true value.
    
        for trials in range(10):            # We try 10 times to have statistically relevant results.
            for P in range(start, 100):     # ... we try to store different number of pattern, up to 100 but we stop before
                meanError = []
                h.makePattern(P, frac)
                
                # We try to retrieve all the patterns
                for mu in range(P):
                    run = h.run(mu = mu,flipRatio = 0.1)
                    if(run != -1):
                        meanError.append(run)           #We have retrieved a pattern
                    else:
                        print "ERROR NOT CONVERGING"    #We never got this error. If you have it then you should discard the data
                
                # We check if they are not correctly retrieved, if so we have our alpha max
                if (np.mean(meanError) > 2):        #Our criterion for the Maximum number of pattern
                    alphaList.append(P/float(N))    
                    print "N : ", N," Run : ", trials, "Pmax : ", P ,"alphaMax : ", P/float(N), "\n"
                    break   #So we don't do unecessary check
            
        # For every size of network we have our Alpha_max averaged over ten try
        print "------> N : ", N, "Alpha_max : ", np.mean(alphaList), "+/- ", stats.sem(alphaList), "\n"
            
    
if __name__ == "__main__":
   cProfile.run('maxLoad()', 'maxLoad.profile')
