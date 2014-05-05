from PIL import Image
from copy import copy
from time import sleep
from pylab import *
import numpy

class hopfieldNetwork:
    def __init__(self,N):
        self.N = N
    def makePattern(self,P=1,ratio=0.5):
        """
        """
        self.pattern = -numpy.ones((P,self.N**2))
        idx = int(ratio*self.N**2)
        for i in range(P):
            self.pattern[i,:idx] = 1
            print(self.pattern)
            self.pattern[i] = numpy.random.permutation(self.pattern[i])
        print(self.pattern)
        self.weight = numpy.zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                sum = 0
                for k in range(P):
                    sum += self.pattern[k,i]*self.pattern[k,j]
                self.weight[i,j] = 1./self.N * sum
        print self.weight
if __name__ == "__main__":
    h = hopfieldNetwork(200)
    h.makePattern()
