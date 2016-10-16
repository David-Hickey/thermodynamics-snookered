import scipy as sp
import particle as pa
import matplotlib.pyplot as pl
import matplotlib.animation as an

class Simulation(object):

    approximate_time_step = 0.05
    """
    The approximate time step to use. The actual time step will be determined
    such that an integer number of steps occur before the next collision, but
    the step will be as close to this value as possible.
    """

    acceptable_error = 1e-3

    def __init__(self, container):
        self.__container = container
        self.__particles = [self.__container]

        self.__future_collisions = None
        self.__time_until_next_collision = sp.nan
        self.__steps_until_next_collision = sp.nan
        self.__time_step = Simulation.approximate_time_step
        self.__steps_since_last_collision_calculation = sp.nan

    def add_particle(self, particle):
        self.__particles.append(particle)

    def get_particles(self):
        return tuple(self.__particles)

    def __clean_up_after_collision(self):
        self.__future_collisions = None
        self.__time_until_next_collision = sp.nan
        #self.__steps_until_next_collision = sp.nan
        #self.__time_step = sp.nan
        self.__steps_since_last_collision_calculation = sp.nan

    def __handle_collision(self):
        time_since_last_collision_calculation = self.__steps_since_last_collision_calculation * self.__time_step

        # This is not valid because a collision may change future collisions,
        # which isn't accounted for here.
        #
        # For this reason I propose something be added to setup_for_next that
        # will handle the problem when the number of steps == 0.
        #
        for collision in self.__future_collisions:
            time_until = collision.time_until - time_since_last_collision_calculation
            not_integer_number_of_steps = time_until / Simulation.approximate_time_step
        
            integer_number_of_steps = sp.floor(not_integer_number_of_steps)

            if integer_number_of_steps == 0:
                collision.particle_1.handle_collision(collision.particle_2)
            else:
                break
            

    def __setup_for_next_collision(self):
        self.__future_collisions = pa.calculate_time_to_collision_all(self.__particles)

        self.__time_until_next_collision = self.__future_collisions[0].time_until

        # Calculate the time step length and number of steps.
        not_integer_number_of_steps = self.__time_until_next_collision / Simulation.approximate_time_step
        
        integer_number_of_steps = sp.floor(not_integer_number_of_steps)

        if integer_number_of_steps == 0:
            print self.__time_until_next_collision
            print not_integer_number_of_steps
            print self.__future_collisions

        #self.__time_step = self.__time_until_next_collision / integer_number_of_steps
        #self.__steps_until_next_collision = integer_number_of_steps
        self.__steps_since_last_collision_calculation = 0

    def next_step(self):
        if sp.isnan(self.__time_until_next_collision):
            self.__setup_for_next_collision()

        for particle in self.__particles:
            # Move all the particles forward.
            particle.update_position(self.__time_step)

        #self.__steps_until_next_collision -= 1
        self.__steps_since_last_collision_calculation += 1

        self.__handle_collision()
        self.__clean_up_after_collision()

container = pa.Particle2D((0, 0), (0, 0), radius = 100, mass = 1000, immovable = True)
test_pa = pa.Particle2D((0, 10), (5, 5), radius = 5, mass = 1)
test_pa_2 = pa.Particle2D((90, 0), (0, 7), radius = 5, mass = 1)
test_pa_3 = pa.Particle2D((-50, 0), (0, 7), radius = 5, mass = 1)
test_pa_4 = pa.Particle2D((-90, 0), (0, 7), radius = 5, mass = 1)
test_pa_5 = pa.Particle2D((-30, 0), (0, -7), radius = 5, mass = 1)
test_pa_6 = pa.Particle2D((20, -30), (0, -7), radius = 5, mass = 1)
big = pa.Particle2D((-35, -40), (0, 0), radius = 30, mass = 10)

sim = Simulation(container)
sim.add_particle(big)
sim.add_particle(test_pa)
sim.add_particle(test_pa_2)
sim.add_particle(test_pa_3)
sim.add_particle(test_pa_4)
sim.add_particle(test_pa_5)
sim.add_particle(test_pa_6)


def next_step(frame, figure, simulation):
    # Calculate the next simulation step.
    simulation.next_step()

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
sim_animation = an.FuncAnimation(main_figure, next_step, interval=5, fargs=(main_figure, sim), blit=False)
pl.show()

