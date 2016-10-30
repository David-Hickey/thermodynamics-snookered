import scipy as sp
import particle as pa
import simulation as si
import animation as an

def init_simulation(sim, particles):
    for particle in particles:
        sim.add_particle(particle)
        
    sim.first_step()



container = pa.Particle2D((0, 0), (0, 0), radius = 100, mass = 1000, immovable = True)

particles = [
    pa.Particle2D((0, 0), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((0, 50), (0, -10), radius = 5, mass = 1),
    pa.Particle2D((0, -12), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((0, -24), (0, 20), radius = 5, mass = 1),
    pa.Particle2D((0, -36), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((20, -48), (0, -7), radius = 5, mass = 1),
    pa.Particle2D((-35, -40), (0, 0), radius = 30, mass = 10),
    pa.Particle2D((20, 0), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((20, 50), (0, -10), radius = 5, mass = 1),
    pa.Particle2D((20, -12), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((20, -24), (0, 20), radius = 5, mass = 1),
    pa.Particle2D((20, -36), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((50, -48), (0, -7), radius = 5, mass = 1),
    pa.Particle2D((40, 0), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((40, 50), (0, -10), radius = 5, mass = 1),
    pa.Particle2D((40, -12), (0, 0), radius = 5, mass = 1),
    pa.Particle2D((40, -24), (0, 20), radius = 5, mass = 1),
    pa.Particle2D((40, -36), (0, 0), radius = 5, mass = 1)
]

sim = si.Simulation(container)
init_simulation(sim, particles)

def next_iteration(animation):
    sim.next_step()

animation = an.Animation(particles + [container], next_iteration)
animation.start()
