import particle as pa

class Simulation(object):

    time_step = 0.07

    def __init__(self, container):
        self.__container = container
        self.__particles = [self.__container]
        #self.__sorted_particles = self.__particles[:]
        #self.__sorted_particles.sort(key = lambda p: p.get_next_collision_time())

        self.__elapsed_time = 0.0

    def add_particle(self, particle):
        self.__particles.append(particle)
        #self.__sorted_particles.append(particle)
        #self.__sorted_particles.sort(key = lambda p: p.get_next_collision_time())

    def get_particles(self):
        return tuple(self.__particles)

    def __handle_collision(self):
        queue = sorted(self.__particles, key = lambda p: p.get_next_collision_time())

        #print queue
        #print "HANDLER CALLED"
        
        for particle in queue:
            if particle.get_next_collision_time() - Simulation.time_step < self.__elapsed_time:
                partner = particle.get_next_collision_particle()

                particle.handle_collision(partner)

                pa.calculate_time_to_collision(particle, self.__particles, self.__elapsed_time)

    def first_step(self):
        pa.calculate_time_to_collision_all(self.__particles)

        #self.__sorted_particles.sort(key = lambda p: p.get_next_collision_time())

    def next_step(self):
        self.__handle_collision()

        self.__elapsed_time += Simulation.time_step
        
        #print "NEXT STEP"
        #print time.time()
        for particle in self.__particles:
            # Move all the particles forward.
            particle.update_position(Simulation.time_step)
        #print time.time()


