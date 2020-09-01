#!/usr/bin/env python

import scipy as sp
import scipy.integrate as spi
import scipy.linalg as spl

import particle as pa
import simulation as si
import animation as an

import matplotlib.pyplot as pl

# These are constant settings that determine initial conditions and settings.

# Whether the animation is shown or not. If yes, no data is recorded. If ndims
# is >2 and animate is true, it will animate a projection into the xy plane
# (which may give the false impression that balls are passing through one 
# another).
animate = True

# The number of dimensions (2 or 3 supported).
ndims = 2

# The number of balls.
nballs = 5

# The radius of the container (metres).
container_scale = 1.0

# The standard deviation of the normally generated speeds (m/s).
speed_scale = 0.1

# The radius of all the particles to add to the simulation.
radius = 0.05

# The number of frames to run before producing graphs (if animate is False).
iterations = 250

# Settings end.

container = pa.Particle(sp.zeros(ndims),
                        sp.zeros(ndims),
                        radius=container_scale,
                        mass=1e20, 
                        immovable=True)
    
sim = si.Simulation(container, ndims=ndims)

velocities = []

# Add many particles at random, ensuring there are no overlaps
for i in xrange(nballs - 1):
    velocity = sp.random.normal(size=ndims, scale=0.1)
    position = sim.find_unused_position(radius)

    new_particle = pa.Particle(position, velocity, radius=radius, mass=1)
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

def plot_kinetic_energy(sim):
    """
    Plot the total kinetic energy of all the particles.
    
    Arguments:
        sim: the simulation object.
    """

    states = sim.get_states()

    kinetic_energies = [d["kinetic_energy"] for d in states]
    times = sp.arange(0, len(states), 1) * si.Simulation.time_step

    ke_figure = pl.figure("Kinetic Energy")

    ke_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    ke_figure.gca().set_ylabel(r"$\rm{Kinetic\ Energy\ /\ J}$")        
    
    ke_figure.gca().plot(times, kinetic_energies)

    ke_figure.savefig("kinetic_energy_%dd.png" % sim.get_ndims())

def plot_total_momentum(sim):
    """
    Plot the total momentum of all particles (not conserved, but would expect
    it to stay quite small).
    
    Arguments:
        sim: the simulation object.
    """
    
    states = sim.get_states()

    particle_momentum = sp.array([d["momentum_total_vec"] for d in states])
    container_momentum = sp.array([d["container_momentum_vec"] for d in states])
    
    total_momentum = container_momentum + particle_momentum
    total_momentum = spl.norm(total_momentum, axis=1)
    
    times = sp.arange(0, len(total_momentum), 1) * si.Simulation.time_step

    p_figure = pl.figure("Total Momentum")

    p_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    p_figure.gca().set_ylabel(r"$\rm{Momentum\ /\ Ns}$")        
    
    p_figure.gca().plot(times, total_momentum)

    p_figure.savefig("total_momentum_%dd.png" % sim.get_ndims())

def plot_maxwell_boltzmann(sim, nbins=20):
    """
    Plot a Maxwell Boltzmann distribution histogram and curve.
    
    Arguments:
        sim: the simulation object.
    """

    def maxwell_boltzmann_function(v, T, m, k=1.38e-23):
        return (v ** (sim.get_ndims() - 1)) * sp.exp(- (0.5 * m * v ** 2) 
                                                     / (k * T))

    states = sim.get_states()
    
    # Take the second half of the states, giving the simulation time to
    # reach equilibrium.
    most_recent_states = states[-len(states) / 2:]
    
    # This is fine because we know kinetic energy is conserved, so temperature
    # is the same at all points.
    temperature = most_recent_states[-1]["temperature"]
    
    mass = sim.get_particle_mass()
    
    speeds = []

    for state in most_recent_states:
        speeds.extend(state["speeds"])

    max_speed = max(speeds)

    frequency, bins = sp.histogram(speeds, bins=nbins)

    speed_probabilities = frequency.astype(float) / len(speeds)

    mb_figure = pl.figure("Maxwell Boltzmann Distribution")

    mb_figure.gca().set_xlabel(r"$\rm{Speed\ /\ ms^{-1}}$")
    mb_figure.gca().set_ylabel(r"$\rm{Probability\ Density\ /\ (ms^{-1})^{-1}}$")

    plot_speeds = sp.linspace(0, max_speed, 1000)
    maxwell_boltzmann = maxwell_boltzmann_function(plot_speeds,
                                                   temperature, 
                                                   mass)

    normalisation_constant = spi.quad(maxwell_boltzmann_function, 
                                      0, 
                                      sp.inf, 
                                      args=(temperature, mass))[0]

    maxwell_boltzmann /= normalisation_constant
    
    mb_figure.gca().hist(speeds, bins=nbins, normed=True)
    mb_figure.gca().plot(plot_speeds, maxwell_boltzmann)

    mb_figure.savefig("maxwell_boltzmann_%dd.png" % sim.get_ndims())
    
def plot_pressure(sim):
    """
    Plot the pressure against time.
    
    Arguments:
        sim: the simulation object.
    """
    
    if ndims == 2:
        pressure_units = "Nm^{-1}"
    elif ndims == 3:
        pressure_units = "Nm^{-2}"
    
    states = sim.get_states()
    
    momenta = sp.array([d["container_momentum_total"] for d in states])
    times = sp.arange(1, len(momenta) + 1, 1) * si.Simulation.time_step

    force = momenta / times
    pressure = force / sim.get_surface_area()
    
    p_figure = pl.figure("Pressure")

    p_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    p_figure.gca().set_ylabel(r"$\rm{Pressure\ /\ %s}$" % pressure_units)        
    
    p_figure.gca().scatter(times, pressure)

    p_figure.savefig("pressure_%dd.png" % sim.get_ndims())
    
def plot_ideal_gas(sim):
    """
    Plot the ideal gas law rearranged to be equal to zero.
    
    Arguments:
        sim: the simulation object.
    """

    states = sim.get_states()
    
    enclosed = sim.get_enclosed_volume()
    surface = sim.get_surface_area()
    
    momenta = sp.array([d["container_momentum_total"] for d in states])
    times = sp.arange(1, len(momenta) + 1, 1) * si.Simulation.time_step

    force = momenta / times
    pressure = force / surface
    
    num_particles = sim.get_particle_number()
    temperature = sp.array([d["temperature"] for d in states])
    
    pv = pressure * enclosed
    nkt = num_particles * temperature * 1.38e-23
    
    pvnkt_figure = pl.figure("Ideal Gas Law")
    
    pvnkt_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    pvnkt_figure.gca().set_ylabel(r"$(PV - Nk_BT) \rm{\ /\ J}$")
    
    pvnkt_figure.gca().scatter(times, pv - nkt)
    
    pvnkt_figure.savefig("ideal_gas_law_%dd.png" % sim.get_ndims())    
    
    
def plot_angular_momentum(sim):
    """
    Plot the angular momentum of all balls against time.
    
    Arguments:
        sim: the simulation object.
    """
    states = sim.get_states()
    
    ang_mom = sp.array([d["angular_momentum_total"] for d in states])
    ang_mom = spl.norm(ang_mom, axis=1)
    
    times = sp.arange(1, len(ang_mom) + 1, 1) * si.Simulation.time_step
    
    am_figure = pl.figure("Angular Momentum")
    
    am_figure.gca().set_xlabel(r"$\rm{Time\ /\ s}$")
    am_figure.gca().set_ylabel(r"$\rm{Angular\ Momentum\ /\ Nms}^{-1}$")
    
    am_figure.gca().plot(times, ang_mom)
    
    am_figure.savefig("angular_momentum_particles_%dd.png" % sim.get_ndims())
    

def next_iteration(animation):
    """
    The next-frame callback function passed to Animation, in order to 
    completely separate simulation and animation code.
    
    Arguments:
        animation: the Animation object representing the animation.
    """
    
    sim.next_step()


if animate:            
    animation = an.Animation(sim.get_particles_including_container(), next_iteration, scale=container.get_radius())
    animation.start(interval=1)
else:
    for i in xrange(iterations):
        sim.next_step()
        
        if i % 500 == 0:
            print "Currently at %dth iteration" % i

    plot_total_momentum(sim)
    plot_kinetic_energy(sim)
    plot_maxwell_boltzmann(sim)
    plot_pressure(sim)
    plot_ideal_gas(sim)
    plot_angular_momentum(sim)
    
    pl.show()
