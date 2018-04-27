"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 24.04.2018

OrbitFamily Class and L1Family/L2Family Subclasses
"""

# Imports
from Utility import OrbitContinuation
from Utility import Plot


# orbitfamily class
class OrbitFamily:

    def __init__(self, x0, orbitDistance, mu, familyData=None):
        self.x0 = x0
        self.mu = mu
        self.orbitDistance = orbitDistance
        self.familyData = familyData

    def getHaloFamily(self, lagrangian):
        # prints status update
        print("STATUS: Generation of a family of Halo Orbits around " + lagrangian + "...\n")
        outData = OrbitContinuation.natParaConti(self.x0, self.orbitDistance, lagrangian, 30, self.mu)
        self.familyData = outData
        print("\nDONE")

    def getNRHOFamily(self):
        pass

    def plot(self, lagrangian, haloFamily="both", background="off"):
        if self.familyData is None:
            OrbitFamily.getHaloFamily(self, lagrangian)
        Plot.plot(self.familyData, self.mu, haloFamily, background)


class L1Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu, familyData=None):
        super().__init__(x0, orbitDistance, mu, familyData)
        self.lagrangian = "L1"

    def getHaloFamily(self):
        OrbitFamily.getHaloFamily(self, self.lagrangian)

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)


class L2Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu, familyData=None):
        super().__init__(x0, orbitDistance, mu, familyData)
        self.lagrangian = "L2"

    def getHaloFamily(self):
        OrbitFamily.getHaloFamily(self, self.lagrangian)

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)
