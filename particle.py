import scipy as sp

class Particle2D(object):
    """
    Class representing a generalised 2-dimensional 'spherical' particle.
    """

    def __init__(self, initial_position, initial_velocity=(0, 0), mass=1.0, radius=1.0, immovable=False):
        """
        Create a new particle with the given values.

        Arguments:
            initial_position: the starting centre position of the particle, as
                              a 1D array-like of numbers.
            initial_velocity: the starting velocity of the particle, as a 1D
                              array-like of numbers. Should have the same
                              dimensions as initial_position.
            mass: the mass of the particle
            radius: the radius of the particle
        """

        self.__position = sp.array(initial_position, dtype=float);
        self.__velocity = sp.array(initial_velocity, dtype=float)
        self.__mass = float(abs(mass))
        self.__radius = float(abs(radius))
        self.__immovable = immovable

        self.__next_collision_time = sp.nan
        self.__next_collision_particle = None

        if self.__immovable:
            self.__velocity *= 0

    def get_mass(self):
        """
        Get the mass of the particle.

        Returns: a float representing the mass.
        """
        return self.__mass

    def get_radius(self):
        """
        Get the radius of the particle.

        Returns: a float representing the radius.
        """
        return self.__radius

    def get_current_position(self):
        """
        Get the current position of the particle.

        Returns: a 1D scipy array of floats representing the position vector.
        """
        return self.__position

    def get_current_velocity(self):
        """
        Get the current velocity of the particle.

        Returns: a 1D scipy array of floats representing the velocity vector.
        """        
        return self.__velocity

    def get_next_collision_time(self):
        return self.__next_collision_time

    def get_next_collision_particle(self):
        return self.__next_collision_particle

    def set_next_collision(self, time, partner):
        self.__next_collision_time = time
        self.__next_collision_particle = partner

    def has_defined_collision(self):
        return not (sp.isnan(self.__next_collision_time) or (self.__next_collision_particle is None))

    def calculate_position_later(self, time_step):
        """
        Calculate the position of the particle after a time time_step has
        passed, using linear extrapolation.

        Arguments:
            time_step: the time in the future (relative to now) at which the
            position should be calculated.

        Returns: a one-dimensional scipy array representing the expected
                 position.
        """
        return self.__position + (self.__velocity * time_step)

    def intersects(self, other_particle):
        """
        Check whether this particle is intersecting another.

        Returns: True if there is an intersection; False otherwise.
        """
        return (
            sp.sum((self.__position - other_particle.__position) ** 2)
            < (self.__radius + other_particle.__radius) ** 2
        )

    def update_position(self, time_step):
        """
        Move this particle forwards in time by time_step. It will move the
        particle forwards by a displacement equal to velocity * time_step.

        Arguments:
            time_step: the time in the future (relative to now) to which the
            particle should be moved.
        """
        self.__position = self.calculate_position_later(time_step)

    def calculate_time_to_collision(self, other_particle):
        """
        Find the time until this particle will intersect with other_particle
        using linear extrapolation.

        This is done by analytically solving the equation
        (s_1 - s_2) ** 2 = (r_1 +- r_2) ** 2, which is the point where both
        circles make contact.

        The - in the right hand side of the equation is there for the case where
        one of the circles is currently inside the other. This may be the case
        if the vessel is another instance of circle.

        TODO: if there is time, it should not be difficult to vectorise this
        operation using meshgrid. This will ensure that all of the calculations
        are done very efficiently. If this is done, it will have to be a
        @staticmethod which accepts all the particles in its arguments. Meshgrid
        operations would have to include meshgrid(r1, r1), meshgrid(r1, r2) etc.

        Arguments:
            other_particle: another Particle object.

        Returns: a float representing the time until these objects
                 intersect, or NaN if these particles are not on a
                 (future) collision trajectory.
        """

        # Calculate all of the parameters required and assign to variables.
        #
        # This is significantly more readable than trying to do all of this
        # inline when solving.

        # v1.v1 + v2.v2 - 2v1.v2
        dt_sq_coefficient = sp.sum(
            self.__velocity ** 2
            + other_particle.__velocity ** 2
            - 2 * (self.__velocity * other_particle.__velocity)
        )

        # 2(r1.v1 + r2.v2 - r1.v2 - r2.v1)
        dt_coefficient = 2 * sp.sum(
            (self.__position * self.__velocity)
            + (other_particle.__position * other_particle.__velocity)
            - (self.__position * other_particle.__velocity)
            - (other_particle.__position * self.__velocity)
        )

        # r1.r1 + r2.r2 - 2r1.r2
        constant_coefficient_initial = sp.sum(
            self.__position ** 2
            + other_particle.__position ** 2
            - 2 * (self.__position * other_particle.__position)
        )

        # (R1 +- R2) ** 2
        if self.intersects(other_particle):
            rhs = (self.__radius - other_particle.__radius) ** 2
        else:
            rhs = (self.__radius + other_particle.__radius) ** 2

        constant_coefficient = constant_coefficient_initial - rhs

        # Now we have the coefficients when the LHS is equal to zero, it is
        # possible to solve the equation.
        discriminant = dt_coefficient ** 2 - (4 * dt_sq_coefficient * constant_coefficient)

        #print "#printING SOLUTIONS TO QUADRATIC"
        #print "disc", discriminant
        #print dt_sq_coefficient, dt_coefficient, constant_coefficient
        #print self.__position, other_particle.__position
        #print self.__velocity, other_particle.__velocity
        #print self.__radius, other_particle.__radius

        if discriminant < 0 or dt_sq_coefficient == 0:
            # Particles are not on a collision trajectory, so we can return NaN
            # as described in the docstring.
            return sp.nan
        else:
            sqrt_discriminant = sp.sqrt(discriminant)

            root_1 = (-dt_coefficient + sqrt_discriminant) / (2 * dt_sq_coefficient)
            root_2 = (-dt_coefficient - sqrt_discriminant) / (2 * dt_sq_coefficient)

            # If the roots are both negative there are no collisions in the future.
            if root_1 < 0 and root_2 < 0:
                return sp.nan

            # We're only interested in the first collision, that has a dt which
            # is greater than (or equal to) zero, so return the smallest
            # positive root:

            #print root_1, root_2
            try:
                return min(root for root in (root_1, root_2) if root >= 0)
            except ValueError:
                print sqrt_discriminant, dt_sq_coefficient, dt_coefficient, constant_coefficient

    def handle_collision(self, other_particle):
        """
        Handle the collision between the two particles. The velocity of the
        two particles will be adjusted appropriately.

        Arguments:
            other_particle: another Particle object that this one is colliding
                            with.
        """

        #print "POSITIONS", self.__position, other_particle.__position

        # Find a new set of basis vectors to make the collision maths easier.
        radial_axis = self.__position - other_particle.__position
        radial_magnitude = sp.sqrt(sp.sum(radial_axis ** 2))
        radial_axis /= radial_magnitude
        
        perpendicular_axis = sp.array([-radial_axis[1], radial_axis[0]])


        #print "RADIAL AXIS", radial_axis
        #print "PERPENDICULAR AXIS", perpendicular_axis

        # Use the dot product to work out the velocities in the new basis vectors.
        self_velocity_radial = self.__velocity.dot(radial_axis)
        other_velocity_radial = other_particle.__velocity.dot(radial_axis)

        self_velocity_perpendicular = self.__velocity.dot(perpendicular_axis)
        other_velocity_perpendicular = other_particle.__velocity.dot(perpendicular_axis)

        #print "Other radial", other_velocity_radial
        #print "Other perp  ", other_velocity_perpendicular

        # Solve the equations for elastic collisions and conservation of
        # momentum in the radial direction.
        self_velocity_radial_after = (
            2 * (other_particle.__mass * other_velocity_radial)
            + (self.__mass * self_velocity_radial)
            - (other_particle.__mass * self_velocity_radial)
        ) / (self.__mass + other_particle.__mass)

        #other_velocity_radial_after = (
        #    2 * (self.__mass * self_velocity_radial)
        #    + (other_particle.__mass * other_velocity_radial)
        #    + (self.__mass * other_velocity_radial)
        #) / (self.__mass + other_particle.__mass)

        other_velocity_radial_after = (
            2 * self.__mass * self_velocity_radial
            + other_particle.__mass * other_velocity_radial
            - self.__mass * other_velocity_radial
        ) / (self.__mass + other_particle.__mass)

        #print "Other radial after", other_velocity_radial_after
        #print "Self radial after", self_velocity_radial_after

        # Form vectors for the new velocities in the new basis.
        self_velocity_after = sp.array([self_velocity_radial_after, self_velocity_perpendicular])
        other_velocity_after = sp.array([other_velocity_radial_after, other_velocity_perpendicular])

        # Find the i and j basis vectors in this basis.
        i_vector = sp.array([1, 0])
        j_vector = sp.array([0, 1])

        i_vector_new_basis = sp.array([i_vector.dot(radial_axis), i_vector.dot(perpendicular_axis)])
        j_vector_new_basis = sp.array([j_vector.dot(radial_axis), j_vector.dot(perpendicular_axis)])

        # Convert the velocities in the new basis back to the more useful
        # cartesian form.
        self_velocity_after = sp.array([
            self_velocity_after.dot(i_vector_new_basis),
            self_velocity_after.dot(j_vector_new_basis)
        ])
        other_velocity_after = sp.array([
            other_velocity_after.dot(i_vector_new_basis),
            other_velocity_after.dot(j_vector_new_basis)
        ])

        # Update the velocities of the two particles.
        if not self.__immovable:
            #print "Changing self velocity from", self.__velocity, "to", self_velocity_after
            
            self.__velocity = self_velocity_after
            #self.update_position(0.001)

        if not other_particle.__immovable:
            #print "Changing other velocity from", other_particle.__velocity, "to", other_velocity_after
            
            other_particle.__velocity = other_velocity_after
            #self.update_position(0.001)

        #print "Handled it"


    def __str__(self):
        return ("Particle2D[radius=%f, mass=%f, pos=%s, vel=%s"
            % (self.__radius, self.__mass, str(self.__position), str(self.__velocity)))

    def __repr__(self):
        return self.__str__()

def calculate_time_to_collision_all(particles, time_now=0):
    """
    Compare all particles to find the next collisions and store that data in the
    particle objects.

    Arguments:
        particles: an array-like containing Particle objecs.
        time_now: the time passed when this function is called, defaulting to 0.
    """
    for i_1, particle_1 in enumerate(particles):
        # We truncate the second list because otherwise we would compare each
        # particle to each other twice. This would be inefficient and would
        # approximately double execution time, so we're not doing that.
        #
        # This has the added benefit of ensuring that no particle is ever
        # compared to itself, so we do not have to check identity explicitly.

        for particle_2 in particles[:i_1]:
            time_to_collision = particle_1.calculate_time_to_collision(particle_2)
            time_passed_at_collision = time_to_collision + time_now

            if not sp.isnan(time_to_collision):
                if not particle_1.has_defined_collision() or particle_1.get_next_collision_time() > time_passed_at_collision:
                    particle_1.set_next_collision(time_passed_at_collision, particle_2)
                    
                if not particle_2.has_defined_collision() or particle_2.get_next_collision_time() > time_passed_at_collision:
                    particle_2.set_next_collision(time_passed_at_collision, particle_1)

def calculate_time_to_collision(particle, particles, time_now, propagate=2):
    particle.set_next_collision(sp.nan, None)
    
    for current_particle in particles:
        if particle is not current_particle:
            time_to_collision = particle.calculate_time_to_collision(current_particle)
            time_passed_at_collision = time_to_collision + time_now

            if not sp.isnan(time_to_collision):
                # TODO: These lines hide a bug! Ideally we would have no non-mutual collisions.
                if not particle.has_defined_collision() or particle.get_next_collision_time() > time_passed_at_collision:
                    particle.set_next_collision(time_passed_at_collision, current_particle)
                    
                if not current_particle.has_defined_collision() or current_particle.get_next_collision_time() > time_passed_at_collision:
                    current_particle.set_next_collision(time_passed_at_collision, particle)

    if propagate:
        for current_particle in particles:
            if current_particle.get_next_collision_particle() is particle:
                #current_particle.set_next_collision(sp.nan, None)

                calculate_time_to_collision(current_particle, particles, time_now, propagate - 1)

