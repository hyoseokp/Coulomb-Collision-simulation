# Coulomb-Collision-simulation

Simulation of the trajectory that is simulated as two particles interact inverse square force (Coulomb Force, Gravitaional Force). : Coulomb collision
(two body problem)

choose your initial condition of particles

<img width="447" alt="image" src="https://user-images.githubusercontent.com/42889193/136187871-78d4c081-5cea-4764-99da-84ab3ac0b110.png">

example initial condition :

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
