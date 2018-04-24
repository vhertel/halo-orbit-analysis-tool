
"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 24.04.2018

OrbitFamily Class and L1Family/L2Family Subclasses
"""

# Imports
import numpy as np
from Utility import OrbitContinuation
from Utility import Plot




# orbitfamily class
class OrbitFamily:

    def __init__(self, x0, orbitDistance, mu):
        self.x0 = x0
        self.mu = mu
        self.orbitDistance = orbitDistance
        self.familyData = np.array([None, None, None, None, None, None, None])

    def getHaloFamily(self, lagrangian):
        outData = OrbitContinuation.natParaConti(self.x0, self.orbitDistance, lagrangian, 8, self.mu)
        self.familyData = outData

    def getNRHOFamily(self):
        pass

    def plot(self, lagrangian, haloFamily="both", background="off"):
        if self.familyData.any() == None:
            OrbitFamily.getHaloFamily(self, lagrangian)
        Plot.plot(self.familyData, self.mu, haloFamily, background)




class L1Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu):
        super().__init__(x0, orbitDistance, mu)
        self.lagrangian = "L1"

    def getHaloFamily(self):
        outData = OrbitFamily.getHaloFamily(self, self.lagrangian)
        self.familyData = outData

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)




class L2Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu):
        super().__init__(x0, orbitDistance, mu)
        self.lagrangian = "L2"

    def getHaloFamily(self):
        outData = OrbitFamily.getHaloFamily(self, self.lagrangian)
        self.familyData = outData

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)
