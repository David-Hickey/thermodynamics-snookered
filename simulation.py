import scipy as sp
import particle as pa
import matplotlib.pyplot as pl
import matplotlib.animation as an

class Simulation(object):

    time_step = 0.1

    def __init__(self, container):
        self.__container = container
        self.__particles = [self.__container]

        self.__elapsed_time = 0.0

    def add_particle(self, particle):
        self.__particles.append(particle)

    def get_particles(self):
        return tuple(self.__particles)

    def __handle_collision(self):
        queue = self.__particles[:]

        queue.sort(key = lambda p: p.get_next_collision_time())

        #print "HANDLER CALLED"
        
        for particle in queue:
            if particle.get_next_collision_time() - Simulation.time_step < self.__elapsed_time:
                partner = particle.get_next_collision_particle()

                particle.handle_collision(partner)

                pa.calculate_time_to_collision(particle, self.__particles, self.__elapsed_time)


    def first_step(self):
        pa.calculate_time_to_collision_all(self.__particles)

    def next_step(self):
        self.__handle_collision()

        self.__elapsed_time += Simulation.time_step
        
        #print "NEXT STEP"
        #print time.time()
        for particle in self.__particles:
            # Move all the particles forward.
            particle.update_position(Simulation.time_step)
        #print time.time()

container = pa.Particle2D((0, 0), (0, 0), radius = 100, mass = 1000, immovable = True)
test_pa = pa.Particle2D((0, 0), (0, 0), radius = 5, mass = 1)
test_pa_2 = pa.Particle2D((0, 50), (0, -10), radius = 5, mass = 1)
test_pa_3 = pa.Particle2D((0, -12), (0, 0), radius = 5, mass = 1)
test_pa_4 = pa.Particle2D((0, -24), (0, 20), radius = 5, mass = 1)
test_pa_5 = pa.Particle2D((0, -36), (0, 0), radius = 5, mass = 1)
test_pa_6 = pa.Particle2D((20, -48), (0, -7), radius = 5, mass = 1)
big = pa.Particle2D((-35, -40), (0, 0), radius = 30, mass = 10)

sim = Simulation(container)
sim.add_particle(big)
sim.add_particle(test_pa)
sim.add_particle(test_pa_2)
sim.add_particle(test_pa_3)
sim.add_particle(test_pa_4)
sim.add_particle(test_pa_5)
#sim.add_particle(test_pa_6)

sim.first_step()


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

        axes.axis([-100, 100, -100, 100])
        

main_figure = pl.figure()

# We have to assign a variable to this animation or the garbage collector will
# throw it out and nothing will work. :[
sim_animation = an.FuncAnimation(main_figure, next_step, interval=5, fargs=(main_figure, sim), blit=False)
pl.show()

