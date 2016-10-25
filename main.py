import scipy as sp
import particle as pa
import simulation as si
import matplotlib.pyplot as pl
import matplotlib.animation as an

def init_simulation(sim, particles):
    for particle in particles:
        sim.add_particle(particle)
        
    sim.first_step()

def init_animation():
    main_figure.gca().axis([-100, 100, -100, 100])

def next_step(frame, figure, simulation):
    # Calculate the next simulation step.
    simulation.next_step()

    # Get the axes so we don't have to keep typing figure.gca().
    axes = figure.gca()

    # Clear the figure so previous frames don't persist.
    axes.cla()

    for particle in simulation.get_particles():
        position = particle.get_current_position()
        shape = pl.Circle(position, particle.get_radius(), fill=False)
        axes.add_artist(shape)

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

main_figure = pl.figure()

# We have to assign a variable to this animation or the garbage collector will
# throw it out and nothing will work. :[
sim_animation = an.FuncAnimation(
    main_figure,
    next_step,
    interval=5,
    fargs=(main_figure, sim),
    blit=False,
    init_func=init_animation
)

pl.show()
