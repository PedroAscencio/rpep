'''
Script to solve the spherical diffusion equation.
'''

import numpy as np
import matplotlib.pyplot as plt
import time

class Diffusion():
    def __init__(self):
    # Put in constants here.
        self.R = 1e-06 # radius of particle
        self.C_max = 12000 # Cutoff conditions
        self.C_0 = 9500
        self.J_0 = 0.25*9.5e-6
        self.D_s = 1e-16 # Diffusion coefficient
        self.deltat = 10 # Time interval.
        self.Nr = 20 # points radius_interval
        self.T_max = 3600 # End time
        self.Nt = int(np.floor(self.T_max/self.deltat)) # points time interval
        self.t_array = np.linspace(0, self.T_max, self.Nt) # Time array
        self.r_array = np.linspace(0, self.R, self.Nr)/self.R # radius array
        self.deltar =self.r_array[1]-self.r_array[0]# points space interval
        self.C_array = np.empty([self.Nt,self.Nr])
        self.C_array[0,:]=np.ones((1,self.r_array.size))*(self.C_0/self.C_max) # Concentration array, intialised.
        self.aux_01=np.zeros(self.Nr-2)
        self.aux_02=np.zeros(self.Nr-2)
        
    def j(self):
    # Current as function of time. To be implemented later!
        self.j_array= self.J_0*np.sin(6*3.1415*self.t_array/3600)

            
    def boundary_conditions(self,i):
    # Script that defines the boundary conditions.
        N=self.Nr-1
        self.C_array[i,0] = self.C_array[i,1]
        self.C_array[i,N] = self.C_array[i,N-1] - (self.deltar*self.R/(self.C_max*self.D_s)) * self.j_array[i]
        
        
    def f_diff_E(self,i):
        for k in range(1,self.Nr-1):
            aux1=(self.C_array[i,k+1] - self.C_array[i,k-1]) / self.deltar 
            aux2=(self.C_array[i,k+1] - 2 * self.C_array[i,k] + self.C_array[i,k - 1]) / ((self.deltar) ** 2) 
            aux3=self.deltat * self.D_s *(1/(self.R**2)) * ((1 / self.r_array[k])*aux1+aux2)
            self.C_array[i+1,k] = self.C_array[i,k] + aux3
                
            
    def f_diff_CN(self,i):
        for k in range(1,self.Nr-1):
            aux1=(self.C_array[i,k+1] - self.C_array[i,k-1]) / self.deltar 
            aux2=(self.C_array[i,k+1] - 2 * self.C_array[i,k] + self.C_array[i,k - 1]) / ((self.deltar) ** 2) 
            self.aux_01[k-1]=self.deltat * self.D_s *(1/(self.R**2)) * ((1 / self.r_array[k])*aux1+aux2)
            self.C_array[i+1,k] = self.C_array[i,k] + 0.5*(self.aux_01[k-1]+self.aux_02[k-1])
        self.aux_02=self.aux_01
        
    def solver(self):
    # This will solve the diffusion equation.
        self.j()
        
        for i in range(0,self.Nt-1):
            self.boundary_conditions(i)
            self.f_diff_CN(i)
            #print(self.C_array[i+1,self.Nr-1]/self.C_max)
                
            #ax.plot(self.r_array/self.R,self.C_array[i,:]/self.C_max)
            ax.clear()
            ax.set_xlim([0,1])
            ax.set_ylim([0,1])
            ax.plot(self.r_array[:], self.C_array[i,:])
            fig.canvas.draw()   # draw
            fig.canvas.flush_events()

if __name__ == '__main__':
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.ion()

    fig.show()
    fig.canvas.draw()

    class_init = Diffusion()
    solution  = class_init.solver()
    
