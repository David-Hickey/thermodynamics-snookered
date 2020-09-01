I have set up different scripts for all of the different cases. They generally
set up the initial conditions and handle the output data. Each script can be
run without any modifications on the part of the user. The runnable scripts 
can be easily identified because they all begin with "main_".

In each main file, there are settings to customise, defined just below the 
import statements. The comments above each setting explain what it does.

Note that the settings are defined so that they give results in reasonable 
time, rather than taking much longer to give more accurate results.

NB: in most files, setting animate = True will stop graphs being produced. 
This is not a mistake.


=== List of main files ===

main_pressure.py:       produces all of the graphs and extracts useful data 
                        from the simulation. The graphs are produced and then
                        saved in the current working directory.
                                    
main_brownianmotion.py: shows an animation demonstrating Brownian motion. 
                        Plots a graph showing mean free path distribution and 
                        the random walk by the brownian particle.

main_brownianmotor.py:  shows an animation and gets data demonstrating the
                        futility of the Brownian motor. Produces a graph
                        showing angular momentum of the motor against time.
                        It works by considering a motor of extremely high mass
                        such that it does not move over the time scale of the
                        simulation and tries to answer the question of whether
                        it would ever accelerate a significant amount.

main_gas_laws.py:       plots the pressure-temperature graph and derives the
                        constants in the VDW gas law.


=== List of other files ===

simulation.py:          co-ordinates all of the particles and keeps a record
                        of state variables.

particle.py:            contains the Particle class and relevant functions and
                        variables.

animation.py:           handles the graphical frontend.


