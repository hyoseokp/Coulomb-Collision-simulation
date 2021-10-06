import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

class Coulomb_Collsion:

    def set_init_condition(self,
            m1    , m2    ,
            q1    , q2    ,
            x1_0  , y1_0  , 
            x2_0  , y2_0  , 
            vx1_0 , vy1_0 , 
            vx2_0 , vy2_0 , 
            k, t, time_resolution):
        super().__init__()
        
        self.time_resolution = time_resolution
        self.m1    = m1     ;  self.m2    = m2
        self.q1    = q1     ;  self.q2    = q2
        self.x1_0  = x1_0   ;  self.y1_0  = y1_0
        self.x2_0  = x2_0   ;  self.y2_0  = y2_0
        self.vx1_0 = vx1_0  ;  self.vy1_0 = vy1_0 
        self.vx2_0 = vx2_0  ;  self.vy2_0 = vy2_0
        self.t = t
        self.k = k
        
    def Coulomb_force(self, vec, t):
        
        x, y, vx, vy = vec
        r = np.linalg.norm([x,y])
        M_r = self.m1*self.m2/(self.m1+self.m2)
        
        return [vx, vy, (1/M_r)*self.k*(self.q1*self.q2)*x/r**3, (1/M_r)*self.k*(self.q1*self.q2)*y/r**3]

    def Collision(self):

        x_ini  = self.x2_0  -  self.x1_0
        y_ini  = self.y2_0  -  self.y1_0
        vx_ini = self.vx2_0 - self.vx1_0
        vy_ini = self.vy2_0 - self.vy1_0

        time_step = np.linspace(0,self.t, self.t*self.time_resolution)

        vec_ini = [x_ini, y_ini, vx_ini, vy_ini]

        vec = odeint(self.Coulomb_force, vec_ini, time_step)


        x1 = (-self.m2*vec[:,0] + (self.m1*self.vx1_0 + self.m2*self.vx2_0)*time_step 
              + self.m1*self.x1_0 + self.m2*self.x2_0)/(self.m1 + self.m2)
        
        y1 = (-self.m2*vec[:,1] + (self.m1*self.vy1_0 + self.m2*self.vy2_0)*time_step 
              + self.m1*self.y1_0 + self.m2*self.y2_0)/(self.m1 + self.m2)
        
        x2 = (self.m1*vec[:,0] + (self.m1*self.vx1_0 + self.m2*self.vx2_0)*time_step 
              + self.m1*self.x1_0 + self.m2*self.x2_0)/(self.m1 + self.m2)
        
        y2 = (self.m1*vec[:,1] + (self.m1*self.vy1_0 + self.m2*self.vy2_0)*time_step 
              + self.m1*self.y1_0 + self.m2*self.y2_0)/(self.m1 + self.m2)

        return [x1, y1, x2, y2 ,vec]
    
    def plot_trajectory(self):
        x1, y1, x2, y2, vec = self.Collision()
        plt.figure()
        theta1 = np.arctan(y1[self.t*self.time_resolution-1]/x1[self.t*self.time_resolution-1])*180/np.pi
        theta2 = -np.arctan(y2[self.t*self.time_resolution-1]/x2[self.t*self.time_resolution-1])*180/np.pi
        plt.plot(x1, y1, '-',label='scattering angle='+'%.2f'%(theta1))
        plt.plot(x2, y2, '-',label='rebound angle='+'%.2f'%(theta2))
        plt.axis('equal')
        plt.title('Trajectory')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.ylim([-0.25,0.25])
        plt.xlim([-0.25,0.25])
        plt.show()

        print('scattering angle :',theta1)
        print('rebound angle :',theta2)
        print('angle sum:',theta1+theta2)
        p1 = np.sqrt(vec[self.t*self.time_resolution-1,2]**2 + vec[self.t*self.time_resolution-1,3]**2)
        print('total momentum :' , p1)
        
        return 0
