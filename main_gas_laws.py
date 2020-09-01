import scipy as sp
import scipy.optimize as spo

import particle as pa
import simulation as si
import animation as an

import matplotlib.pyplot as pl

# The number of dimensions (2 or 3 supported).
ndims = 3

# The number of balls.
nballs = 60

# The radius of the container (metres).
container_scale = 1.0

# The radius of all the particles to add to the simulation.
radius = 0.05

# The number of frames to run before producing graphs (if animate is False).
iterations = 70

# The number of simulations to try with different kinetic energies.
simulations = 20

# Settings end.

def simulate_once(speed_scale):
    """
    Run a new simulation with the given speed scale and return its final 
    state.
    
    Arguments:
        speed_scale: standard deviation of the gaussian used to generat the
                     velocities.
                     
    Returns: the final state of the simulation as a dictionary.
    """
    
    container = pa.Particle(sp.zeros(ndims),
                            sp.zeros(ndims),
                            radius=container_scale,
                            mass=1e20, 
                            immovable=True)
        
    sim = si.Simulation(container, ndims=ndims)

    velocities = []

    # Add many particles at random, ensuring there are no overlaps
    for i in xrange(nballs - 1):
        velocity = sp.random.normal(size=ndims, scale=speed_scale)
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
    
    for i in xrange(iterations):
        sim.next_step()
        
    states = sim.get_states()    
    
    return states[-1]
 
def plot_pressure_temperature(states):
    """
    Plot the pressure against temperature, and fit the VDW equation to the
    data to determine the constants a and b.
    
    Arguments:
        states: the list of all states.
    """
    
    def fit_vdw(t, a, bN):
        nkt = nballs * 1.38e-23 * t
        v_corrected = enclosed - bN
        
        num_density_sq = (nballs / enclosed) ** 2
        
        return (nkt / v_corrected) - (a * num_density_sq)

    def ideal(t):
        return nkt / enclosed
        

    if ndims == 2:
        enclosed = container_scale ** 2 * sp.pi
        surface = 2 * sp.pi * container_scale
        enclosed_correction = nballs * radius ** 2 * sp.pi
        
        pressure_units = "Nm^{-1}"
    elif ndims == 3:
        enclosed = container_scale ** 3 * 4.0 / 3.0 * sp.pi
        surface = 4 * sp.pi * container_scale ** 2
        enclosed_correction = nballs * radius ** 3 * sp.pi * 4.0 / 3.0
        
        pressure_units = "Nm^{-2}"
    
    momenta = sp.array([d["container_momentum_total"] for d in states])
    
    times = sp.array([d["time_elapsed"] for d in states])

    force = momenta / times
    pressure = force / surface
    
    num_particles = nballs
    temperature = sp.array([d["temperature"] for d in states])
    
    pv = pressure * enclosed
    nkt = num_particles * temperature * 1.38e-23
    
    figure = pl.figure("Pressure Temperature")
    
    params, cov = spo.curve_fit(fit_vdw, temperature, pressure)
    
    fitted_line = figure.gca().plot(temperature, fit_vdw(temperature, *params))
    
    scatter = figure.gca().scatter(temperature, pressure)
    figure.gca().set_xlabel(r"$\rm{Temperature\ /\ K}$")
    figure.gca().set_ylabel(r"$\rm{Pressure}\ /\ %s$" % pressure_units)

    figure.gca().plot(temperature, ideal(temperature))
    
    pl.legend(labels=["Van der Waals Fit", "Ideal Gas", "Raw Data"], 
              loc="upper left")
    
    figure.savefig("pressure_temperature_%dd.png" % ndims)
    
    with open("pressure_temperature_%dd.txt" % ndims, "w") as f:
        for name, value, error in zip(["a", "bN"], params, sp.sqrt(cov.diagonal(0))):
            f.write(name)
            f.write(" = ")
            f.write(str(value))
            f.write(" +/- ")
            f.write(str(error))
            f.write("\n")

    
states = []
for i in xrange(simulations):
    scale = sp.sqrt(i + 1) * 0.01
    
    print "%dth Simulation (scale %g)" % (i, scale)
    
    latest_states = simulate_once(scale)
    
    states.append(latest_states)

plot_pressure_temperature(states)
pl.show()
