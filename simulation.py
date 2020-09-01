import particle as pa
import scipy as sp
import scipy.linalg as spl
import numpy.lib as li

class CollisionIterator(object):
    """
    An iterator representing a collision queue. Was created so that new 
    collisions generated during a tick could be appended to the list of
    collisions handled in a clean and elegant manner.
    """

    def __init__(self, collisions):
        """
        Arguments:
            collisions: an array-like of objects to iterate over.
        """
        self._collisions = list(collisions)
        self._current_index = 0

    def next(self):
        if self._current_index >= len(self._collisions) - 1:
            raise StopIteration()
        
        ret = self._collisions[self._current_index]

        self._current_index += 1

        return ret

    def add_next_collision(self, collision):
        """
        Append an object to the end of the list.
        """
        self._collisions.append(collision)

    def __iter__(self):
        return self
    

class Simulation(object):
    """
    An object representing a collection of particles and offering methods
    to coordinate their interactions and glean data from their behaviour.
    """

    # The time step for the simulation.
    time_step = 0.1
    tol = time_step / 100.0

    def __init__(self, container, ndims, particles=(), log_state=True):
        """ 
        Arguments:
            container: the Particle object containing all other particles.
            ndims: the number of dimensions in which to run the simulation.
            particles: the particles to add to the simulation as an array-like.
            log_state: if True, then the simulation will record data which may
                       be accessed using the get_simulation_states() method.
        """
        
        self._container = container
        self._particles = list(particles) + [self._container]
        self._particles_not_container = list(particles)
        self._elapsed_time = 0.0
        self._ndims = ndims
        self._log_states = log_state
		
        self._simulation_states = []

    def add_particle(self, particle):
        """
        Add a particle to the simulation.
        
        Arguments:
            particle: the Particle object to add.
        """
        
        self._particles.append(particle)
        self._particles_not_container.append(particle)

    def get_particles_including_container(self):
        """
        Get all of the particles in the simulation, including the container.
        
        Returns: a tuple of all particles.
        """
    
        # We return a tuple because they're immutable so the recipient can't
        # change the particle list.
        return tuple(self._particles)
        
    def get_particles_excluding_container(self):
        """
        Get all of the particles in the simulation, excluding the container.
        
        Returns: a tuple of all particles.
        """
        
        # We return a tuple because they're immutable so the recipient can't
        # change the particle list.
        return tuple(self._particles_not_container)
        
    def get_states(self):
        """
        Get the logged states of the simulation. If the log_state option 
        passed to the initialiser is False, then a ValueError will be raised.
        
        Returns: a list of the simulation states, where each state is itself a
                 dictionary.
        """
        
        if not self._log_states:
            raise ValueError("Simulation is not logging states.")

        return self._simulation_states
        
    def get_ndims(self):
        """
        Get the number of dimensions in which the simulation is running.
        
        Returns: the number of dimensions as an int.
        """
        return self._ndims

    def first_step(self):
        """
        Initialise the simulation. This should only be called once, after all
        the particles have been added and the simulation is ready to begin.
        """
        
        self._calculate_time_to_collision_all()
        
    def _prevent_escape(self):
        """
        I've never actually seen a ball escape, but if the simulation runs for
        long enough it sometimes happens.

        This should correct escaping balls without doing anything unphysical.
        """

        container_radius = self._container.get_radius()

        for particle in self._particles_not_container:
            pos = particle.get_current_position()
            
            radius = particle.get_radius()

            distance_sq = sp.dot(pos, pos)

            if distance_sq > (container_radius - radius) ** 2:
                distance = sp.sqrt(distance_sq)
                new_pos = pos * (container_radius - radius) / (distance + container_radius * 1e-6)

                particle.set_position(new_pos)

                particle.handle_collision(self._container)

                self._calculate_time_to_collision_all()
                
    def _handle_collisions_this_tick(self):
        """
        Handles all the collisions in this tick. New collisions arising as a
        result of collisions this tick are also handled if appropriate.
        """
        queue = CollisionIterator(self._particles)
        
        for particle in queue:
            if particle.get_next_collision_time() - Simulation.time_step < self._elapsed_time + Simulation.tol:
                other_particle = particle.get_next_collision_particle()

                particle.handle_collision(other_particle)
                
                particle.calculate_time_to_collision(self._particles, self._elapsed_time, propagate=True)
                other_particle.calculate_time_to_collision(self._particles, self._elapsed_time, propagate=True)

                queue.add_next_collision(particle)
                queue.add_next_collision(other_particle)
                
    def next_step(self):
        """
        Handles the next iteration in the simulation. Updates the position of
        all particles, handles all collisions, and updates the states (if 
        relevant).
        """
        
        if self._log_states:
            self._simulation_states.append(self._find_state_now())

        self._handle_collisions_this_tick()
        
        self._elapsed_time += Simulation.time_step

        for particle in self._particles:
            particle.update_position(Simulation.time_step)
            
        self._prevent_escape()
            
    def get_surface_area(self):
        """
        Get the surface area of the container in 3D, or the circumference in
        2D. Raises a ValueError is the simulation is not 2D or 3D.
        
        Returns: the surface area of the container as a float.
        """
        
        if self._ndims == 2:
            return self._container.get_radius() * 2 * sp.pi
        elif self._ndims == 3:
            return (self._container.get_radius() ** 2) * 4 * sp.pi
            
        raise ValueError("Not supported for number of dimensions")
            
    def get_enclosed_volume(self):
        """
        Get the enclosed volume of the container in 3D, or the enclosed area in
        2D. Raises a ValueError is the simulation is not 2D or 3D.
        
        Returns: the volume enclosed by the container as a float.
        """
    
        if self._ndims == 2:
            return (self._container.get_radius() ** 2) * sp.pi
        elif self._ndims == 3:
            return (self._container.get_radius() ** 3) * (4.0 / 3.0) * sp.pi
            
        raise ValueError("Not supported for number of dimensions")
            
    def get_particle_number(self):
        """
        Get the total number of particles in the simulation.
        
        Returns: the number of particles in the simulation as an int.
        """
        return len(self._particles_not_container)
        
    def get_particle_volume(self):
        """
        Get the total volume of all the particles in the simulation. In 2D, 
        instead returns the total enclosed area of all the particles. Raises
        a ValueError is the simulation is not in 2D or 3D.
        
        Returns: the total volume of all the particles as a float.
        """
        volume = 0.0
        
        for particle in self._particles_not_container:
            if self._ndims == 2:
                volume += (particle.get_radius() ** 2) * sp.pi
            elif self._ndims == 3:
                volume += (particle.get_radius() ** 3) * sp.pi * (4.0 / 3.0)
            else:
                raise ValueError("Not supported for number of dimensions")
                
        return volume
                
    def get_particle_mass(self):
        """
        Get the mass of all the particles in the simulation, assuming that
        they are the same. If the masses are not all the same, then a
        ValueError is rasied.
        
        Returns: the mass of all of the particles as a float.
        """
        mass = self._particles_not_container[0].get_mass()
        
        for particle in self._particles_not_container:
            if abs(particle.get_mass() - mass) > 1e-10:
                raise ValueError("Not all masses are the same")
                
        return mass
            

    def _find_state_now(self):
        """
        Record the current state of the particles in the simulation.
        """
        
        state = {}

        state["time_elapsed"] = self._elapsed_time
        state["container_momentum_total"] = spl.norm(self._container.get_momentum_received())

        momentum_vector = sp.zeros(self._ndims)
        kinetic_energy = 0.0
        velocities = []
        total_angular_momentum = sp.zeros(3)
        
        for particle in self._particles_not_container:
            vel = particle.get_current_velocity()
            pos = particle.get_current_position()
            mass = particle.get_mass()
        
            kinetic_energy += 0.5 * mass * sp.dot(vel, vel)
            momentum_vector += mass * vel
            velocities.append(vel)
            
            # Pad vectors by adding zeroes at the end if necessary.
            if self._ndims == 2:
                mom_3d = mass * li.pad(vel, (0, 1), "constant")
                pos_3d = li.pad(pos, (0, 1), "constant")
            elif self._ndims == 3:
                mom_3d = mass * vel
                pos_3d = pos
            
            total_angular_momentum += sp.cross(pos_3d, mom_3d)

        state["angular_momentum_total"] = total_angular_momentum.copy()
        state["speeds"] = spl.norm(velocities, axis=1)
        state["momentum_total_vec"] = momentum_vector.copy()
        state["container_momentum_vec"] = self._container.get_vector_momentum_received().copy()
        
        state["kinetic_energy_total"] = kinetic_energy
        state["kinetic_energy"] = kinetic_energy / self.get_particle_number()
        state["temperature"] = 2.0 * state["kinetic_energy"] / (1.38e-23 * self._ndims)
        
        return state
           
            
    def find_unused_position(self, radius, attempts=1000):
        """
        Find an unused position in the simulation - that is, a position where a ball
        can be placed without overlapping any other balls.
        
        Arguments:
            radius: the radius of the ball to add as a float.
            attempts: the integer number of attempts to make to find an empty position,
                      after which a ValueError is rasied if no position is 
                      found.
                     
        Returns: the unused position as a numpy ndarray.
        """
        
        container_size = self._container.get_radius()
        
        for i in xrange(attempts):
            position = sp.random.uniform(-container_size, container_size, size=self._ndims)
            
            if radius + sp.sqrt(sp.sum(position ** 2)) >= container_size:
                continue
            
            is_unused_position = True

            for particle in self._particles_not_container:
                displacement = position - particle.get_current_position()
                distance = sp.sqrt(sp.sum(displacement ** 2))

                if distance <= radius + particle.get_radius():
                    is_unused_position = False
                    
                    break

            if is_unused_position:
                return position

        raise ValueError("Can't find an unused position in the simulation")

    def _calculate_time_to_collision_all(self):
        """
        Compare all particles to find the next collisions and store that data in the
        particle objects.
        """
        
        time_now = self._elapsed_time
        
        for particle in self._particles:
            particle.set_next_collision(sp.nan, None)

        for i_1, particle_1 in enumerate(self._particles):
            # We truncate the second list because otherwise we would compare each
            # particle to each other twice. This would be inefficient and would
            # approximately double execution time, so we're not doing that.
            #
            # This has the added benefit of ensuring that no particle is ever
            # compared to itself, so we do not have to check identity explicitly.

            for particle_2 in self._particles[:i_1]:
                time_to_collision = particle_1.calculate_time_to_collision_single(
                    particle_2)
                time_at_collision = time_to_collision + time_now

                if not sp.isnan(time_to_collision):
                    if not particle_1.has_defined_collision() or particle_1.get_next_collision_time() > time_at_collision:
                        particle_1.set_next_collision(time_at_collision, particle_2)
                        
                    if not particle_2.has_defined_collision() or particle_2.get_next_collision_time() > time_at_collision:
                        particle_2.set_next_collision(time_at_collision, particle_1)

