import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
plt.style.use('dark_background')

class Particle:
    def __init__(self, color, id=0, charge=1.602E-19, r=np.zeros(2), v=np.zeros(2), rad=0.01, m=1):
        self.id = id
        self.r = r  # Set of the x and y coordinates of the particle
        self.v = v  # Set of the x and y velocities of the particle
        self.rad = rad # Radius of the particle
        self.m = m  # Mass of the particle
        self.charge = charge * (np.random.randint(0, 2) * 2 - 1)  # Charge of the particle
        if self.charge > 0:
            color = "red"
        else:
            color = "blue"
        self.color = color


class Sim:
    X=10
    Y=10
    def __init__(self, dt=0.0001, Np=45):
        self.dt = dt  # Time jump between events
        self.Np = Np  # Number of particles
        self.particles = [Particle(i) for i in range(Np)] # call the particle function Np(number of particles) times to have Np times particles

    def coll_det(self):
        for i, particle1 in enumerate(self.particles):
            x, y = particle1.r
            
                    # Wall collision handling
            if x - particle1.rad < -self.X / 2:
                particle1.r[0] = -self.X / 2 + particle1.rad
                particle1.v[0] *= -1  # reverse x-velocity component
            elif x + particle1.rad > self.X / 2:
                particle1.r[0] = self.X / 2 - particle1.rad
                particle1.v[0] *= -1  # reverse x-velocity component
            if y - particle1.rad < -self.Y / 2:
                particle1.r[1] = -self.Y / 2 + particle1.rad
                particle1.v[1] *= -1  # reverse y-velocity component
            elif y + particle1.rad > self.Y / 2:
                particle1.r[1] = self.Y / 2 - particle1.rad
                particle1.v[1] *= -1  # reverse y-velocity component

            for j, particle2 in enumerate(self.particles[i + 1:], start=i + 1):
                m1, m2, r1, r2, v1, v2 = particle1.m, particle2.m, particle1.r, particle2.r, particle1.v, particle2.v
                dist = np.linalg.norm(r1 - r2)
                if dist <= (particle1.rad + particle2.rad):
                    # Collision handling
                    v1_new = (v1 * (m1 - m2) + 2 * m2 * v2) / (m1 + m2)
                    v2_new = (v2 * (m2 - m1) + 2 * m1 * v1) / (m1 + m2)
                    particle1.v = v1_new
                    particle2.v = v2_new
                    

    # def increment(self):
    #     self.coll_det()
    #     for particle in self.particles:
    #         particle.r += self.dt * particle.v  # Euler's method
    def increment(self):
        self.coll_det()
        for particle in self.particles:
            # #print(type(particle.v))
            # k1 = self.dt * particle.v
            # k2 = self.dt * (particle.v + 0.5 * k1)
            # k3 = self.dt * (particle.v + 0.5 * k2)
            # k4 = self.dt * (particle.v + k3)
    
            # particle.r += (1/6) * (k1 + 2*k2 + 2*k3 + k4) #runge-kutta
            particle.r += particle.v*self.dt
            #self.dt += 0.000001

    def particle_positions(self):
        return [particle.r for particle in self.particles] #returns the positions of the particles

    def particle_color(self):
        return [particle.color for particle in self.particles] # returns the colors of the particles


sim = Sim()

for particle in sim.particles:
    particle.r = np.random.uniform([-sim.X/2, -sim.Y/2], [sim.X/2, sim.Y/2], size=2)
    particle.v = np.array([np.random.randint(-50,50), np.random.randint(-50,50)])
    # above initializes the first positions and velocities
sim.particles[0].color = "purple" # changing the color of first particle just because I can

fig, ax = plt.subplots()
scatter = ax.scatter([], [], s=0.01, alpha=0.5)

def init():
    ax.set_xlim(-sim.X/2, sim.X/2)
    ax.set_ylim(-sim.Y/2, sim.Y/2)
    return scatter,

def update(frame):
    sim.increment()
    positions = sim.particle_positions()
    ax.clear()
    ax.set_xlim(-sim.X / 2, sim.X / 2)
    ax.set_ylim(-sim.Y / 2, sim.Y / 2)
    colors = sim.particle_color()
    ax.scatter(*zip(*positions), c = colors)

# animation = FuncAnimation(fig, update, frames=12000, init_func=init)
# plt.show()
 

animation = FuncAnimation(fig, update, frames=1000, init_func=init)
plt.show()

# ani = FuncAnimation(fig, update, frames=60000, init_func=init)
# ani.save('ani31.gif',writer='pillow',fps=25)

#above is to save the animation as gif instead of straight up playing it

#for now sometimes particle spawn on each other and that causes minor problems