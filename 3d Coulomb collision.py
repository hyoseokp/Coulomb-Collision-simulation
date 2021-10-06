import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

class Coulomb_Collsion:

    def set_init_condition(self,
            m1    , m2    ,
            q1    , q2    ,
            x1_0  , y1_0  , z1_0,  # Initial position of particle 1
            x2_0  , y2_0  , z2_0,  # Initial position of particle 2
            vx1_0 , vy1_0 , vz1_0, # Initial velocity of particle 1
            vx2_0 , vy2_0 , vz2_0, # Initial velocity of particle 2
            k, t, time_resolution):
        super().__init__()
        
        self.time_resolution = time_resolution
        
        self.m1    = m1     ;  self.m2    = m2
        self.q1    = q1     ;  self.q2    = q2
        
        self.x1_0  = x1_0   ;  self.y1_0  = y1_0   ;  self.z1_0  = z1_0 
        self.x2_0  = x2_0   ;  self.y2_0  = y2_0   ;  self.z2_0  = z2_0 
        self.vx1_0 = vx1_0  ;  self.vy1_0 = vy1_0  ;  self.vz1_0  = vz1_0 
        self.vx2_0 = vx2_0  ;  self.vy2_0 = vy2_0  ;  self.vz2_0  = vz2_0 
        
        self.t = t
        self.k = k
        
    def Coulomb_force(self, vec, t):
        
        x, y, z, vx, vy, vz = vec
        r = np.linalg.norm([x,y,z])
        M_r = self.m1*self.m2/(self.m1+self.m2)
        
        return [vx, vy, vz, 
                (1/M_r)*self.k*(self.q1*self.q2)*x/r**3, 
                (1/M_r)*self.k*(self.q1*self.q2)*y/r**3,
                (1/M_r)*self.k*(self.q1*self.q2)*z/r**3]

    def Collision(self):

        x_ini  = self.x2_0  -  self.x1_0
        y_ini  = self.y2_0  -  self.y1_0
        z_ini  = self.z2_0  -  self.z1_0
        
        vx_ini = self.vx2_0 - self.vx1_0
        vy_ini = self.vy2_0 - self.vy1_0
        vz_ini = self.vz2_0 - self.vz1_0
        
        time_step = np.linspace(0,self.t, self.t*self.time_resolution)

        vec_ini = [x_ini, y_ini, z_ini, vx_ini, vy_ini, vz_ini]

        vec = odeint(self.Coulomb_force, vec_ini, time_step)


        x1 = (-self.m2*vec[:,0] + (self.m1*self.vx1_0 + self.m2*self.vx2_0)*time_step 
              + self.m1*self.x1_0 + self.m2*self.x2_0)/(self.m1 + self.m2)
        
        y1 = (-self.m2*vec[:,1] + (self.m1*self.vy1_0 + self.m2*self.vy2_0)*time_step 
              + self.m1*self.y1_0 + self.m2*self.y2_0)/(self.m1 + self.m2)
        
        z1 = (-self.m2*vec[:,2] + (self.m1*self.vz1_0 + self.m2*self.vz2_0)*time_step 
              + self.m1*self.z1_0 + self.m2*self.z2_0)/(self.m1 + self.m2)
        
        x2 = (self.m1*vec[:,0] + (self.m1*self.vx1_0 + self.m2*self.vx2_0)*time_step 
              + self.m1*self.x1_0 + self.m2*self.x2_0)/(self.m1 + self.m2)
        
        y2 = (self.m1*vec[:,1] + (self.m1*self.vy1_0 + self.m2*self.vy2_0)*time_step 
              + self.m1*self.y1_0 + self.m2*self.y2_0)/(self.m1 + self.m2)
        
        z2 = (self.m1*vec[:,1] + (self.m1*self.vz1_0 + self.m2*self.vz2_0)*time_step 
              + self.m1*self.z1_0 + self.m2*self.z2_0)/(self.m1 + self.m2)

        return [x1, y1, z1, x2, y2, z2, vec]
    
    def plot_trajectory(self):
        
        x1, y1, z1, x2, y2, z2, vec = self.Collision()
        
        p1 = np.sqrt(vec[self.t*self.time_resolution-1,3]**2 + 
                     vec[self.t*self.time_resolution-1,4]**2 +
                     vec[self.t*self.time_resolution-1,5]**2
                    )
        
        theta1 = np.arctan(y1[self.t*self.time_resolution-1]/x1[self.t*self.time_resolution-1])*180/np.pi
        theta2 = -np.arctan(y2[self.t*self.time_resolution-1]/x2[self.t*self.time_resolution-1])*180/np.pi
        
        plt.figure()
        plt.plot(x1, y1, '-',label='scattering angle='+'%.2f'%(theta1))
        plt.plot(x2, y2, '-',label='rebound angle='+'%.2f'%(theta2))
#        plt.axis('equal')
        plt.title('Trajectory')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.ylim([-3,3])
        plt.xlim([-3,3])
        plt.show()
        
        print('scattering angle :',theta1)
        print('rebound angle :',theta2)
        print('angle sum:',theta1+theta2)
        print('total momentum :' , p1)
        
        return 0
    def ani2d(self):
    
        x1, y1, z1, x2, y2, z2, vec = self.Collision()
        plt.figure()
        for i in range(self.t*self.time_resolution):
            plt.plot(x1[:i],y1[:i], '-')
            plt.plot(x2[:i],y2[:i], '-')
            plt.show(block=False)
            plt.title('Trajectory')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
#            plt.ylim([-3,3])
#            plt.xlim([-3,3])
            plt.pause(0.00001)
            plt.close()
            
    def ani3d(self):

        x1, y1, z1, x2, y2, z2, vec = self.Collision()
        f=100
        for i in range(int(self.t*self.time_resolution/f)):
            graph1 = plt.figure(1)
            graph2 = graph1.add_subplot(111,projection = '3d')
            graph2.plot(x1[:int(f*i)],y1[:int(f*i)],z1[:int(f*i)], 'r')
            graph2.plot(x2[:int(f*i)],y2[:int(f*i)],z2[:int(f*i)], 'b')
            
            plt.show(block=False)
            plt.pause(0.00001)
            plt.close()
            
c = Coulomb_Collsion()

[m1, q2] = 1,-1
[m2, q1] = 100, 100
[x1_0, y1_0, z1_0]    = -1, 1, 0
[x2_0, y2_0, z2_0]    = 0, 0, 0
[vx1_0, vy1_0, vz1_0] = 1, 0, 0
[vx2_0, vy2_0, vz2_0] = 0, 0, 0
k = 1
t = 10
time_resolution = 1000

c.set_init_condition(m1, m2, q1, q2, 
                     x1_0, y1_0, z1_0, 
                     x2_0, y2_0, z2_0, 
                     vx1_0, vy1_0, vz1_0, 
                     vx2_0, vy2_0, vz2_0, 
                     k, t, time_resolution)
#c.plot_trajectory()
#c.ani()
c.ani3d()
