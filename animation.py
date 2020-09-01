import matplotlib.pyplot as pl
import matplotlib.animation as an

import random
import string

class Animation(object):
    """
    An object representing the visual animation and offering methods to handle
    most of the animation code.
    """

    def __init__(self, particles, frame_callback, scale):
        """
        Arguments:
            particles: the list of Particle objects to display on the 
            animation.
            frame_callback: the function to call after drawing each new frame.
            scale: half the side length of the animation window (which is 
                   square).
        """
        
        self._particles = particles
        self._frame_callback = frame_callback
        self._scale = scale

        self._frame = 0

        self._patches = map(make_patch, self._particles)
        
    def _init_animation(self):
        """
        Initialise the animation by creating the figures and scaling the axes.
        
        Returns: the list of patches.
        """
    
        self._main_figure.gca().axis(
            [-self._scale, self._scale, -self._scale, self._scale])
        
        self._main_figure.gca().autoscale(False)

        for patch in self._patches:
            self._main_figure.gca().add_artist(patch)

        return self._patches

    def _next_step(self, frame):
        """
        Update the patches to produce the next animation frame.
        
        Arguments:
            frame: an int representing the current frame.
            
        Returns: the list of patches.
        """
        
        # Calculate the next simulation step.
        self._frame_callback(self)

        for particle, patch in zip(self._particles, self._patches):
            patch.center = particle.get_current_position()

        self._frame += 1

        return self._patches

    def get_current_frame_number(self):
        """
        Get the integer number of frames elapsed so far.
        
        Returns: the current frame number as an int.
        """
        
        return self._frame

    def get_pyplot_animation(self):
        """
        Get the matplotlib FuncAnimation object that is dealing with the
        animation.
        
        Returns: the matplotlib FuncAnimation object responsible for the
        animation.
        """
        
        return self._sim_animation
        
    def start(self, block=True, interval=10, size=(10, 10)):
        """
        Begin the animation.
        
        Arguments:
            block: whether to stop anything else executing (other than the 
                   callback function) while the animation runs.
            interval: the time (in milliseconds) between frames.
            size: the size of the plot (in inches, because matplotlib takes
                  inches by default), as a tuple.
        """
        
        self._main_figure = pl.figure("Animation", figsize=size)

        # We have to assign a variable to this animation or the garbage
        # collector will throw it out and nothing will work.
        self._sim_animation = an.FuncAnimation(
            self._main_figure,
            self._next_step,
            interval=interval,
            blit=True,
            init_func=self._init_animation
        )

        pl.show(block)
        
def make_patch(particle):
    """
    Create a patch for a particle.
    
    Arguments:
        particle: a Particle object.
        
    Returns: a Circle patch representing the particle, or a projection of
             the particle onto the xy plane.
    """

    if particle.is_immovable():
        return pl.Circle(
            particle.get_current_position(),
            particle.get_radius(),
            fill=False
        )
    else:
        if particle.get_colour() is not None:
            fill_colour = particle.get_colour()
        else:
            # Generate random colour.
            fill_colour = "#" + string.zfill(
                hex(random.randrange(0, 16777215))[2:], 6)

        return pl.Circle(
            particle.get_current_position(),
            particle.get_radius(),
            color=fill_colour)