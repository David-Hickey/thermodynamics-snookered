import scipy as sp
import scipy.linalg as spl
import particle as pa
import simulation as si
import animation as an
import matplotlib.pyplot as pl
import numpy.lib as li


# These are constant settings that determine initial conditions and settings.

# Whether the animation is shown or not.
animate = True

# The number of balls.
nballs = 40

# The standard deviation of the normally generated speeds (m/s).
speed_scale = 0.05

# The radius of all the particles to add to the simulation.
radius = 0.05

# The number of frames to run before producing graphs.
iterations = 250

# Settings end.

container = pa.Particle((0, 0), (0, 0), radius=1, mass=1e20, immovable=True)
sim = si.Simulation(container, ndims=2)

def make_s_shape(radius=0.068):
    """
    Generate the list of particles that make up the motor.
    
    Arguments:
        radius: the radius of the particles making up the motor.
    """
    
    particles = [
        pa.Particle((0, -0.6),    (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((-0.1, -0.5), (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((-0.2, -0.4), (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((-0.2, -0.2), (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((-0.1, -0.1), (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0, 0),       (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0.1, 0.1),   (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0.2, 0.2),   (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0.2, 0.4),   (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0.1, 0.5),   (0, 0), radius=radius, mass=1e20, immovable=True),
        pa.Particle((0, 0.6),     (0, 0), radius=radius, mass=1e20, immovable=True)
    ]


    return particles
    
def next_iteration(animation):
    """
    The next-frame callback function passed to Animation, in order to 
    completely separate simulation and animation code.
    
    Arguments:
        animation: the Animation object representing the animation.
    """
    
    sim.next_step()
    record_angular_momentum(sim)
    
    if animation.get_current_frame_number() == iterations:
        plot_angular_momentum(angular_momenta)
        
        pl.show()
        
def record_angular_momentum(sim):
    """
    Record the cumlative angular momentum of the motor at current.
    """
    
    total_angular_momentum = sp.zeros(3)

    for particle in motor_particles:
        r = particle.get_current_position()
        p = particle.get_vector_momentum_received()
        
        # Pad the vectors by adding zeroes at the end.
        p_cross = li.pad(p, (0, 1), "constant")
        r_cross = li.pad(r, (0, 1), "constant")
        
        angular_momentum = sp.cross(r_cross, p_cross)
        
        total_angular_momentum += angular_momentum
        
    angular_momenta.append(total_angular_momentum)
        
def plot_angular_momentum(angular_momenta):
    """
    Plot the angular momentum of the motor against time.
    """
    
    angular_momentum_norm = 0.0
    for particle in sim.get_particles_excluding_container():
        r = particle.get_current_position()
        p = particle.get_vector_momentum_received()
        
        # Pad the vectors by adding zeroes at the end.
        p_cross = li.pad(p, (0, 1), "constant")
        r_cross = li.pad(r, (0, 1), "constant")
        
        angular_momentum_norm += spl.norm(sp.cross(r_cross, p_cross))
        
    # This value is useful because it gives context to the angular momentum
    # of the motor.
    print "Total angular momentum norm:", angular_momentum_norm
    
    angular_momenta = sp.array(angular_momenta)
    
    z_component = angular_momenta[:, -1]
    
    times = sp.arange(0, len(z_component), 1) * si.Simulation.time_step
    
    am_figure = pl.figure("Angular Momentum")
    
    am_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    am_figure.gca().set_ylabel(r"$\rm{Angular\ Momentum\ /\ Nms}^{-1}$")        
    
    am_figure.gca().plot(times, z_component)

    am_figure.savefig("angular_momentum_2d.png")
    
    
    
angular_momenta = []

particles = make_s_shape()
motor_particles = particles[:]

for particle in particles:
    sim.add_particle(particle)

velocities = []

# Add many particles at random, ensuring there are no overlaps
for i in xrange(nballs - 1):
    velocity = sp.random.normal(size=2, scale=speed_scale)
    position = sim.find_unused_position(radius)
    
    new_particle = pa.Particle(position, velocity, radius=radius, mass=1)

    particles.append(new_particle)
    sim.add_particle(new_particle)
    
    velocities.append(velocity)
    
# Ensure that the average velocity is zero
position = sim.find_unused_position(radius)
correcting_particle_velocity = -sp.sum(velocities, axis=0)

sim.add_particle(pa.Particle(position,
    correcting_particle_velocity,
    radius=radius,
    mass=1))

sim.first_step()


    
if animate:
    animation = an.Animation(particles + [container], next_iteration, scale=container.get_radius())
    animation.start(interval=1)
else:
    for i in xrange(iterations):
        sim.next_step()
        
        record_angular_momentum(sim)
        
        if i % 500 == 0:
            print "Currently at %dth iteration" % i
        
    plot_angular_momentum(angular_momenta)
    
    pl.show()
    

    

