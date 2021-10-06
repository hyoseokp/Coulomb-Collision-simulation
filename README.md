# Coulomb-Collision-simulation

you can use this program as

choose your initial condition of particles


c = Coulomb_Collsion()

m1    = 1   ; m2    = 1
q1    = 1   ; q2    = 1
x1_0  = -10 ; y1_0  = 0.01
x2_0  = 0   ; y2_0  = -0.01
vx1_0 = 10  ; vy1_0 = 0 
vx2_0 = 0   ; vy2_0 = 0
k = 1
t = 3
time_resolution = 1000

c.set_init_condition(m1, m2, q1, q2, x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0, k, t, time_resolution)
c.plot_trajectory()
