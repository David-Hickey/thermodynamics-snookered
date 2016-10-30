import matplotlib

matplotlib.use("TkAgg")

import matplotlib.pyplot as pl
import matplotlib.animation as an

class Animation(object):

    def __init__(self, particles, frame_callback):
        self.__particles = particles
        self.__frame_callback = frame_callback

        self.__patches = map(self.__make_patch, self.__particles)

    @staticmethod
    def __make_patch(particle):
        return pl.Circle(particle.get_current_position(), particle.get_radius(), fill=False)
    
    def __init_animation(self):
        self.__main_figure.gca().axis([-100, 100, -100, 100])
        self.__main_figure.gca().autoscale(False)

        for patch in self.__patches:
            self.__main_figure.gca().add_artist(patch)

        return self.__patches

    def __next_step(self, frame):
        # Calculate the next simulation step.
        self.__frame_callback(self)

        for particle, patch in zip(self.__particles, self.__patches):
            patch.center = particle.get_current_position()

        return self.__patches
        
    def start(self, block=True):
        self.__main_figure = pl.figure()

        # We have to assign a variable to this animation or the garbage collector will
        # throw it out and nothing will work. :[
        self.__sim_animation = an.FuncAnimation(
            self.__main_figure,
            self.__next_step,
            interval=5,
            blit=True,
            init_func=self.__init_animation
        )

        pl.show(block)

    
