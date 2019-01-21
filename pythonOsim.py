# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:50:55 2019

@author: mhmdk
"""

new_states_file = open("corrected_tug_of_war_states.sto", "w")

new_states_file.write("Tug_of_War_Competition\n")
new_states_file.write("nRows = 1\n")
new_states_file.write("nColumns = 7\n")
new_states_file.write("inDegree = no\n")
new_states_file.write("endheader\n")
new_states_file.write("time\tblock_tz\tblock_tz_u\tRightMuscle.activation\tRightMuscle.fiber_length\tLeftMuscle.activation\tLeftMuscle.fiber_length\n")
new_states_file.write("0\t0\t0\t0.01\t1\t0.01\t1")

new_states_file.close()