# example1.py
import math
from penrose import PenroseP3, BtileL, psi

# A simple example starting with a BL tile

scale = 100
# Configuration of the tiling
config = {'draw-arcs': True,
          'arc-colour': 'black',
          'draw-tiles': False,
          'proportion': 0.5}
tiling = PenroseP3(scale, ngen=7, config=config)

# Create the initial tiles, a triangle
theta = 2*math.pi / 5
rot = math.cos(theta) + 1j*math.sin(theta)
A = -scale/2 + 0j
B = scale/2 * rot
C = scale/2 / psi + 0j
tiling.set_initial_tiles([BtileL(A - B.real, B - B.real, C - B.real)])
tiling.make_tiling()

# Makes the matplotlib figure
tiling.make_plot()


stroke_width = str(psi ** tiling.ngen * tiling.scale *
                   tiling.config['base-stroke-width'])

print('Proportion: {}'.format(tiling.config['proportion']))
print('Generation: {}'.format(tiling.ngen))
print('BSW: {}'.format(stroke_width))
