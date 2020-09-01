import scipy as sp
import scipy.linalg as spl

class Particle(object):

    # Global constant, tolerance for floating point comparisons.
    tol = 1e-13

    def __init__(self, initial_position, initial_velocity, mass, radius,
                 immovable=False, colour=None):
        """
        Arguments:
            initial_position: the starting centre position of the particle, as
                              a 1D array-like of numbers.
            initial_velocity: the starting velocity of the particle, as a 1D
                              array-like of numbers. Should have the same
                              dimensions as initial_position.
            mass: the mass of the particle.
            radius: the radius of the particle.
            immovable: whether the particle is immovable or not.
            colour: colour of the particle in matplotlib colour format; 
                    if unspecified colour is set randomly
        """

        self._position = sp.array(initial_position, dtype=float)
        self._velocity = sp.array(initial_velocity, dtype=float)

        if len(self._position) != len(self._velocity):
            raise ValueError("Position and velocity must have same dimensions")
        
        self._mass = float(abs(mass))
        self._radius = float(abs(radius))
        self._immovable = immovable
        self._colour = colour

        self._next_collision_time = sp.nan
        self._next_collision_particle = None

        self._momentum_received = 0.0
        self._vector_momentum_received = sp.zeros(len(self._position))

        if self._immovable:
            self._velocity *= 0


    def get_mass(self):
        """
        Get the mass of the particle.

        Returns: a float representing the mass.
        """
        return self._mass

    def get_radius(self):
        """
        Get the radius of the particle.

        Returns: a float representing the radius.
        """
        return self._radius

    def get_current_position(self):
        """
        Get the current position of the particle.

        Returns: a 1D scipy array of floats representing the position vector.
        """
        return self._position

    def get_current_velocity(self):
        """
        Get the current velocity of the particle.

        Returns: a 1D scipy array of floats representing the velocity vector.
        """        
        return self._velocity

    def set_velocity(self, velocity):
        """
        Set the velocity of the particle. Raises ValueError if the vector
        has the wrong length or if this particle is immovable.
        
        Arguments:
            velocity: an ndarray representing a velocity.
        """
        
        if len(velocity) != len(self._velocity):
            raise ValueError("Velocity must have same size as current")
            
        if self._immovable:
            raise ValueError("Can't change velocity of an immovable particle")
        
        self._velocity = sp.array(velocity, dtype=float)

    def set_position(self, position):
        """
        Set the position of the particle. Raises ValueError if the vector
        has the wrong length or if this particle is immovable.
        
        Arguments:
            velocity: an ndarray representing a velocity.
        """
        
        if len(position) != len(self._position):
            raise ValueError("Position must have same size as current")
            
        if self._immovable:
            raise ValueError("Can't change position of an immovable particle")
        
        self._position = sp.array(position, dtype=float)


    def is_immovable(self):
        """
        Get whether this particle is immovable.
        
        Returns: True if this particle cannot be moved; False otherwise.
        """
        return self._immovable
        
    def get_colour(self):
        """
        Get the colour of this particle, which may be None.
        
        Returns: the colour of this particle as an allowed matplotlib colour.
        """
        return self._colour

    def get_next_collision_time(self):
        """
        Get the time of the next collision for this particle, which may be 
        NaN if this particle has no upcoming collision.
        
        Returns: the time of the next collision as a float.
        """
        return self._next_collision_time

    def get_next_collision_particle(self):
        """
        Get the particle this one will collide with next, which may be None
        if this particle has no upcoming collision.
        
        Returns: the Particle this will next collide with.
        """
        return self._next_collision_particle

    def set_next_collision(self, time, partner):
        """
        Set the time of the next collision and the particle with which this 
        one will collide.
        
        Arguments:
            time: the time of the collision as a float.
            partner: the Particle this will next collide with.
        """
        self._next_collision_time = time
        self._next_collision_particle = partner

    def has_defined_collision(self):
        """
        Get whether this particle has a collision in its forseeable future.
        
        Returns: True if it does; False otherwise.
        """
        return not (sp.isnan(self._next_collision_time)
                    or (self._next_collision_particle is None))

    def intersects(self, other_particle):
        """
        Check whether this particle is intersecting another.
    
        Returns: True if there is an intersection; False otherwise.
        """
        displacement = self._position - other_particle.get_current_position()
        radius_sum = self._radius + other_particle.get_radius()
    
        return sp.dot(displacement, displacement) <= radius_sum ** 2

    def update_position(self, time_step):
        """
        Move this particle forwards in time by time_step. It will move the
        particle forwards by a displacement equal to velocity * time_step.

        Arguments:
            time_step: the time in the future (relative to now) to which the
            particle should be moved.
        """
        self._position += self._velocity * time_step

    def calculate_time_to_collision_single(self, other_particle):
        """
        Find the time until this particle will intersect with other_particle
        using linear extrapolation.

        This is done by analytically solving the equation
        (s_1 - s_2) ** 2 = (r_1 +- r_2) ** 2, which is the point where both
        circles make contact.

        The - in the right hand side of the equation is there for the case where
        one of the circles is currently inside the other.

        Arguments:
            other_particle: another Particle object.

        Returns: a float representing the time until these objects
                 intersect, or NaN if these particles are not on a
                 (future) collision trajectory.
        """

        # Calculate all of the coefficients.
        velocity_diff = self._velocity - other_particle.get_current_velocity()
        position_diff = self._position - other_particle.get_current_position()

        dt_sq_coefficient = sp.dot(velocity_diff, velocity_diff)
        dt_coefficient = 2 * sp.dot(position_diff, velocity_diff)
        constant_coefficient_lhs = sp.dot(position_diff, position_diff)

        # (R1 +- R2) ** 2
        if self.intersects(other_particle):
            rhs = (self._radius - other_particle._radius) ** 2
        else:
            rhs = (self._radius + other_particle._radius) ** 2

        constant_coefficient = constant_coefficient_lhs - rhs

        discriminant = dt_coefficient ** 2 \
                       - (4 * dt_sq_coefficient * constant_coefficient)

        if discriminant < 0:
            return sp.nan


        #sqrt_discriminant = sp.sqrt(discriminant)

        roots = sp.roots([dt_sq_coefficient,
                          dt_coefficient,
                          constant_coefficient])

        if sp.any(roots > Particle.tol):

            # We're only interested in the first collision, that has a dt which
            # is greater than zero, so return the smallest positive root:
            return min(root for root in roots if root > Particle.tol)

        return sp.nan
            
    def calculate_time_to_collision(self, particles, time_now, propagate=True):
        """
        Calculate which of the given particles this will collide with first, 
        and when that will happen.
        
        Arguments:
            particles: the particles to compare this one to.
            time_now: the elapsed time now.
            propagate: whether to also recalculate for the particles this one
                       was originally going to collide with.
        """
        self.set_next_collision(sp.nan, None)
        
        for current_particle in particles:
            if self is not current_particle:
                time_to_collision = self.calculate_time_to_collision_single(
                    current_particle
                )
                
                time_at_collision = time_to_collision + time_now

                if not sp.isnan(time_to_collision):
                    if (not self.has_defined_collision()
                        or self.get_next_collision_time() > time_at_collision):
                        
                        self.set_next_collision(time_at_collision,
                                                    current_particle)
                        
                    if (not current_particle.has_defined_collision()
                        or current_particle.get_next_collision_time() > time_at_collision):

                        current_particle.set_next_collision(time_at_collision,
                                                            self)

        
        if propagate:
            for current_particle in particles:
                current_partner = current_particle.get_next_collision_particle()
                
                current_hits_self = current_particle is self
                self_hits_current = self._next_collision_particle is current_particle
                
                if current_hits_self or self_hits_current:
                    current_particle.calculate_time_to_collision(particles,
                                                time_now,
                                                propagate=False)

    def handle_collision(self, other_particle):
        """
        Update the velocities of these particles as they undergo a collision.
        
        Arguments:
            other_particle: the particle that is colliding with this one.
        """
        
        other_mass = other_particle.get_mass()
        other_velocity = other_particle.get_current_velocity()

        self_momentum_before = self._mass * self._velocity
        other_momentum_before = other_mass * other_velocity

        # Transform into a frame where other is initially at rest. rf_ prefix
        # indicates this frame.
        rf_velocity_self = self._velocity - other_velocity
        
        radial_vector = self._position - other_particle.get_current_position()
        radial_unit_vector = radial_vector / sp.sqrt(sp.dot(radial_vector, radial_vector))

        rf_speed_self_parallel = rf_velocity_self.dot(radial_unit_vector)
        rf_velocity_self_perpendicular = rf_velocity_self - rf_speed_self_parallel * radial_unit_vector

        # This is now a 1D conservation of momentum and energy problem
        if other_particle.is_immovable():
            
            rf_speed_self_parallel_after = -rf_speed_self_parallel
            rf_velocity_self_parallel_after = radial_unit_vector * rf_speed_self_parallel_after
            self_velocity_after = other_velocity + rf_velocity_self_parallel_after + rf_velocity_self_perpendicular
            self._velocity = self_velocity_after

            other_velocity_after = other_particle.get_current_velocity()

            #if 2 * self._position.dot(self._velocity) + self._velocity.dot(self._velocity) > 0:
            #    # Very rarely a small inaccuracy can cause a ball to bounce out
            #    # of the container, which can be fixed by this.
            #    self_velocity_after = other_velocity - rf_velocity_self_parallel_after + rf_velocity_self_perpendicular
            
        elif self._immovable:
        
            rf_speed_parallel_other_after = 2 * rf_speed_self_parallel
            rf_velocity_other_parallel_after = radial_unit_vector * rf_speed_parallel_other_after
            other_velocity_after = other_velocity + rf_velocity_other_parallel_after
            other_particle.set_velocity(other_velocity_after)
        
            self_velocity_after = self._velocity
            
        else:
            rf_speed_other_parallel_after = 2 * self._mass * rf_speed_self_parallel / (self._mass + other_mass)
            rf_speed_self_parallel_after = rf_speed_self_parallel * (self._mass - other_mass) / (self._mass + other_mass)

            rf_velocity_self_parallel_after = radial_unit_vector * rf_speed_self_parallel_after
            rf_velocity_other_parallel_after = radial_unit_vector * rf_speed_other_parallel_after

            other_velocity_after = other_velocity + rf_velocity_other_parallel_after
            self_velocity_after = other_velocity + rf_velocity_self_parallel_after + rf_velocity_self_perpendicular

            self._velocity = self_velocity_after
            other_particle.set_velocity(other_velocity_after)


        # Update the received momentum (for pressure purposes)        
        self_momentum_after = self._mass * self_velocity_after
        self_momentum_change_vec = self_momentum_after - self_momentum_before  
        
        other_momentum_after = other_mass * other_velocity_after
        other_momentum_change_vec = other_momentum_after - other_momentum_before
        
        if self._immovable:
            self.increase_momentum_received(-other_momentum_change_vec)
            other_particle.increase_momentum_received(other_momentum_change_vec)
        else:
            self.increase_momentum_received(self_momentum_change_vec)
            other_particle.increase_momentum_received(-self_momentum_change_vec)

    def increase_momentum_received(self, increase):
        """
        Increment the magnitude of impulse that has been exerted on this particle.
        
        Arguments:
            increase: the increase in impulse, as a numpy ndarray.
        """
        
        if isinstance(increase, sp.ndarray):
            self._vector_momentum_received += increase
            self._momentum_received += sp.linalg.norm(increase)
        else:
            raise ValueError("Increase must be an ndarray")
        
    def get_momentum_received(self):
        """
        Get the total impulse received, as a scalar.
        
        Returns: the impulse received as a float.
        """
        
        return self._momentum_received
        
    def get_vector_momentum_received(self):
        """
        Get the total impulse received, as a vector.
        
        Returns: the impulse received as a numpy ndarray.
        """
        return self._vector_momentum_received

    def __str__(self):
        return ("Particle(radius=%f, mass=%f, pos=%s, vel=%s)"
            % (self._radius,
               self._mass,
               str(self._position),
               str(self._velocity)))
               
