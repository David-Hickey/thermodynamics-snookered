#!/usr/bin/python

import scipy as sp
import particle as pa
import simulation as si
import animation as an
import matplotlib.pyplot as pl
import scipy.linalg as spl
import scipy.optimize as spo

# These are constant settings that determine initial conditions and settings.

# Whether the animation is shown or not. If yes, no data is recorded. If ndims
# is >2 and animate is true, it will animate a projection into the xy plane
# (which may give the false impression that balls are passing through one 
# another).
animate = True

# The number of dimensions (2 or 3 supported).
ndims = 2

# The number of balls.
nballs = 40

# The radius of the container (metres).
container_scale = 1.0

# The radius of the larger particle (metres).
brownian_scale = 0.3

# The standard deviation of the normally generated speeds (m/s).
speed_scale = 0.1

# The radius of all the particles to add to the simulation.
radius = 0.05

# The number of frames to run before producing graphs (if animate is False).
iterations = 250

# Settings end.

def next_iteration(animation):
    """
    The next-frame callback function passed to Animation, in order to 
    completely separate simulation and animation code.
    
    Arguments:
        animation: the Animation object representing the animation.
    """
    
    sim.next_step()
    
def plot_walk(positions, velocities):
    def exponential_distribution(x, mfp):
        return sp.exp(-mfp * x) * mfp

    pos_figure = pl.figure("Position")
    
    pos_figure.gca().plot(positions[:, 0], positions[:, 1])
    
    pos_figure.gca().add_artist(an.make_patch(container))
    
    pos_figure.gca().axis([-container.get_radius(),
                           container.get_radius(),
                           -container.get_radius(),
                           container.get_radius()])
                           
    pos_figure.savefig("random_walk%dd.png" % ndims)
    
    # Find where collisions happen.
    velocity_rolled = sp.roll(velocities, 1, axis=0)
    velocity_change = (velocity_rolled - velocities)[1:]
    
    magnitude_change = spl.norm(velocity_change, axis=1)

    collision_positions = positions[1:][magnitude_change != 0]
    collision_positions_rolled = sp.roll(collision_positions, 1, axis=0)
    
    displacements = (collision_positions_rolled - collision_positions)[1:]
    
    distances = spl.norm(displacements, axis=1)
    
    mfp_figure = pl.figure("Mean Free Path")
    freq_dens, bins, patches = mfp_figure.gca().hist(distances, bins=40, normed=True)
    
    bin_width = bins[1] - bins[0]
    bin_centres = (bin_width / 2.0) + bins[:-1]
    
    params, cov = spo.curve_fit(exponential_distribution, 
                                bin_centres, 
                                freq_dens)
                                
    print params, sp.sqrt(cov)

    mfp_figure.gca().plot(bins, exponential_distribution(bins, *params))
    
    
    mfp_figure.gca().legend(["Exponential Distribution", "Path Length"])
    
    mfp_figure.gca().set_xlabel(r"$\rm{Path\ Length\ /\ m}$")
    mfp_figure.gca().set_ylabel(r"$\rm{Probability\ Density\ /\ m^{-1}}$")
    
    mfp_figure.savefig("mean_free_path_%dd" % ndims) 
    

container = pa.Particle(sp.zeros(ndims), 
                        sp.zeros(ndims),
                        radius=container_scale, 
                        mass=0, 
                        immovable=True)
                        
sim = si.Simulation(container, ndims=ndims, log_state=False)

large_particle = pa.Particle(sp.zeros(ndims), 
                             sp.zeros(ndims), 
                             radius=brownian_scale,
                             mass=10)

particles = [large_particle]
sim.add_particle(large_particle)

# Add many particles at random, ensuring there are no overlaps
for i in xrange(nballs):
    velocity = sp.random.normal(size=ndims, scale=speed_scale)
    position = sim.find_unused_position(radius)

    new_particle = pa.Particle(position, velocity, radius=radius, mass=1)

    sim.add_particle(new_particle)
    particles.append(new_particle)

sim.first_step()

velocities = []
positions = []

if animate:
    animation = an.Animation(particles + [container], next_iteration, scale=container.get_radius())
    animation.start(interval=1)
else:
    for i in xrange(iterations):
        sim.next_step()
        
        velocities.append(large_particle.get_current_velocity().copy())
        positions.append(large_particle.get_current_position().copy())
        
        if i % 500 == 0:
            print "Currently at %dth iteration" % i
        
    plot_walk(sp.array(positions), 
              sp.array(velocities))
    
    pl.show()
