# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 18:24:04 2019

@author: MathLab
"""

from math import atan2, degrees
h = 29 # Monitor height in cm
d = 57 # Distance between monitor and participant in cm
r = 1080 # Vertical resolution of the monitor
size_in_px = 240 # The stimulus size in pixels
# Calculate the number of degrees that correspond to a single pixel. This will
# generally be a very small value, something like 0.03.
deg_per_px = degrees(atan2(.5*h, d)) / (.5*r)
deg_per_px
# Calculate the size of the stimulus in degrees
size_in_deg = size_in_px * deg_per_px
(size_in_px, size_in_deg)



size_in_deg = 5.7 # The stimulus size in degrees
# Calculate the number of degrees that correspond to a single pixel. This will
# generally be a very small value, something like 0.03.
deg_per_px = degrees(atan2(.5*h, d)) / (.5*r)
deg_per_px
# Calculate the size of the stimulus in pixels
size_in_px = size_in_deg / deg_per_px
(size_in_px, size_in_deg)