from __future__ import division
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

k = 10
m = 0.02
N = 200
n = 5714
R = 0.1
dt = R*np.sqrt(m/k)

boundary = False

def init_conditions(N,n):
    ypos = np.zeros((N,n))
    for i in xrange(N):
        ypos[i, 0] = np.sin(7*np.pi*(i/(N-1)))

    ypos[:,-1] = np.copy(ypos[:,0]) #setting the previous time position equal to ynow
    return ypos

class string:
    def __init__(self, y,k,m,N,n,dt, boundary = False):
        self.y = y
        self.k = k
        self.m = m
        self.N = N
        self.n = n
        self.dt = dt
        self.potential_energy = np.zeros((self.N,self.n))
        self.kinetic_energy = np.zeros((N,n))
        self.boundary = boundary

    def solver(self):
        for t in xrange(0,self.n-1):

            """
            These if test determines whether we have open or reflective edges by
            using the boundary and initial condition we have
            """
            if self.boundary == True: #reflective edges
                self.y[0,t+1] = 2*self.y[0,t] - self.y[0,t-1]
                self.y[self.N-1,t+1] = 2*self.y[self.N-1,t] - self.y[self.N-1,
                t-1]

            if self.boundary == False: #open edges
                fac = (self.dt**2*self.k)/self.m
                #left masspoint
                self.y[0,t+1] = 2*self.y[0,t] - self.y[0,t-1] - fac*(
                self.y[0,t] - self.y[1,t]
                )
                #right masspoint

                self.y[self.N-1,t+1] = 2*self.y[self.N-1,t] - self.y[self.N-1,t-1] - fac*(
                self.y[self.N-1,t] - self.y[self.N-2,t])


            #solving the wave equation
            for j in xrange(0,N-2):
                self.y[j+1,t+1] = (2*self.y[j+1,t] - self.y[j+1,t-1] +
                (self.dt**2/self.m)*(self.k*self.y[j,t] +
                self.k*self.y[j+2,t] - 2*self.k*self.y[j+1,t]))

    def calculate_kinetic_energy(self):
        for c in xrange(0,self.n-1):
            for i in xrange(N-1):
                dy = self.y[i,c+1] - self.y[i,c]
                diff_y = dy/dt
                self.kinetic_energy[i,c] = 0.5*self.m*diff_y**2

        self.kinetic_energy = np.sum(self.kinetic_energy, axis = 0)


    def calculate_potential_energy(self):
        for q in xrange(self.N-1):
            self.potential_energy[q,:] = 0.5*self.k*(y[q+1,:] - y[q,:])**2
        self.potential_energy = np.sum(self.potential_energy,axis = 0)


    def plotting(self):
        K = self.kinetic_energy
        V = self.potential_energy
        """plt.subplot(2,1,1)
        plt.plot(np.linspace(0,self.n*self.dt,self.n),K)
        plt.plot(np.linspace(0,self.n*self.dt,slef.n),V)
        plt.xlim(0,25)
        plt.subplot(2,1,2)
        plt.plot(np.linspace(0,self.n*self.dt,self.n),V + K)
        plt.xlim(0,25)
        plt.show()"""


string = string(init_conditions(N,n),k,m,N,n,dt, True)
string.plotting()
#calculate_potential_energy(y,N,m,n,k)
