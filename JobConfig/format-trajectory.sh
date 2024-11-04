#!/bin/bash
#
# This script processes the log file from running m4.fcl on one event to extract and print the trajectory.
#
# Andrei Gaponenko, 2024

awk 'BEGIN{ print "# X Y Z nx ny nz"}; /Pre: /{printf "%6.1f    %6.1f    %6.1f    %6.6f    %6.6f    %6.6f\n", $3, $4, $5, ($6/$9),($7/$9),($8/$9)}' ${1:?Specify an input file}
