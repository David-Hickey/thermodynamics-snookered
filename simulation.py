import scipy as sp
import particle as pa
import matplotlib.pyplot as pl
import matplotlib.animation as an

class Simulation(object):

    approximate_time_step = 0.1
    """
    The approximate time step to use. The actual time step will be determined
    such that an integer number of steps occur before the next collision, but
    the step will be as close to this value as possible.
    """

    def __init__(self, container):
        self.__container = container
        self.__particles = [self.__container]

        # time_step is the simulation time between successive iterations. It
        # should be defined such that the ratio of time until next collision
        # and time_step is an integer.
        self.__time_step = sp.nan
        self.__time_until_next_collision = sp.nan
        self.__steps_until_next_collision = sp.nan
        self.__next_colliding_particles = None

    def add_particle(self, particle):
        self.__particles.append(particle)

    def get_particles(self):
        # Think carefully about this - it exposes the particles
        return tuple(self.__particles)

    def __handle_collision(self):
        #print "HANDLING COLLISION IN NEXT"
        # Clean up all the data pertaining to the collision which just
        # passed.
        self.__next_colliding_particles[0].handle_collision(self.__next_colliding_particles[1])

        self.__time_until_next_collision = sp.nan
        self.__time_step = sp.nan
        self.__steps_until_next_collision = sp.nan
        self.__next_colliding_particles = None

    def __setup_for_next_collision(self):
            #print "SETTING UP FOR FUTURE COLLISION"
            # If there are no collisions forecast, calculate when the next one
            # will be.
            self.__time_until_next_collision, self.__next_colliding_particles = pa.calculate_time_to_collision_all(self.__particles)

            # self.__time_until_next_collision, self.__next_colliding_particles

            # until / step = int, but we want step to be close to
            # approximate_time_step. We can work out what the int would be if
            # it wasn't constrained to the integers, then round it to an integer
            # value to find the closest integer to the value we want.
            not_integer = self.__time_until_next_collision / Simulation.approximate_time_step

            #print "not_integer", not_integer
            
            integer = sp.floor(not_integer)

            self.__time_step = self.__time_until_next_collision / integer
            self.__steps_until_next_collision = integer

            print "integer", integer
            print "time_step", self.__time_step
            print "steps", self.__steps_until_next_collision
            print "until", self.__time_until_next_collision

    def next_step(self):
        if sp.isnan(self.__time_until_next_collision):
            self.__setup_for_next_collision()

        for particle in self.__particles:
            # Move all the particles forward.
            particle.update_position(self.__time_step)

        self.__steps_until_next_collision -= 1

        if self.__steps_until_next_collision == 0:
            self.__handle_collision()

            # TODO handle the case where integer == 0 elegantly - this will mean handling > 1 collision / tick.

container = pa.Particle2D((0, 0), (0, 0), radius = 100, mass = 1000, immovable = True)
test_pa = pa.Particle2D((0, 10), (5, 5), radius = 5, mass = 1)
test_pa_2 = pa.Particle2D((90, 0), (0, 7), radius = 5, mass = 1)
test_pa_3 = pa.Particle2D((-50, 0), (0, 7), radius = 5, mass = 1)
test_pa_4 = pa.Particle2D((-90, 0), (0, 7), radius = 5, mass = 1)
test_pa_5 = pa.Particle2D((-30, 0), (0, 7), radius = 5, mass = 1)
test_pa_6 = pa.Particle2D((20, -30), (0, -7), radius = 5, mass = 1)

sim = Simulation(container)
sim.add_particle(test_pa)
sim.add_particle(test_pa_2)
sim.add_particle(test_pa_3)
sim.add_particle(test_pa_4)
sim.add_particle(test_pa_5)
sim.add_particle(test_pa_6)

def next_step(frame, figure, simulation):
    # Calculate the next simulation step.
    sim.next_step()

    # Clear the figure so previous frames don't persist.
    figure.clf()

    # Get the axes so we don't have to keep typing figure.gca().
    axes = figure.gca()

    for particle in simulation.get_particles():
        position = particle.get_current_position()
        shape = pl.Circle(position, particle.get_radius(), fill=False)
        axes.add_artist(shape)

        axes.axis([-100, 100, -100, 100])
        

main_figure = pl.figure()

# We have to assign a variable to this animation or the garbage collector will
# throw it out and nothing will work.
sim_animation = an.FuncAnimation(main_figure, next_step, interval=15, fargs=(main_figure, sim), blit=False)
pl.show()

##import time
##while True:
##    #print "running step"
##    sim.next_step()
##
##    time.sleep(3)
