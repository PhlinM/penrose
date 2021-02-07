# example1.py
import math
from penrose import PenroseP3, BtileL, psi

# A simple example starting with a BL tile

scale = 100
# Configuration of the tiling
config = {'draw-arcs': True,
          'normal-arcs': False,
          'arc-colour': 'black',
          'draw-tiles': False,
          'proportion': 0.5}
tiling = PenroseP3(scale, ngen=5, config=config)

# Create the initial tiles, a triangle
theta = 2*math.pi / 5
rot = math.cos(theta) + 1j*math.sin(theta)
A = -scale/2 + 0j
B = scale/2 * rot
C = scale/2 / psi + 0j
tiling.set_initial_tiles([BtileL(A - B.real, B - B.real, C - B.real)])

# Make the tiles up to ngen generations
tiling.make_tiling()

tiling.make_plot()
