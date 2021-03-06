# example2.py
import math
import webbrowser
from penrose import PenroseP3, BtileS
from cropSVG import crop

# A "sun", with randomly-coloured tiles and including arc paths|

scale = 100
# Configuration of the tiling
config = {'draw-arcs': True,
          'base-stroke-width': 0.1,
          'arc-colour': '#000',
          'draw-tiles': False,
          'random-tile-colours': True,
          'tile-opacity': 0.6,
          'proportion': 0.2}
tiling = PenroseP3(scale, ngen=7, config=config)

# Create the initial tiles, a triangle
theta = math.pi / 5
#  alpha = math.cos(theta)
rot = math.cos(theta) + 1j*math.sin(theta)
A1 = scale + 0.j
B = 0 + 0j
C1 = C2 = A1 * rot
A2 = A3 = C1 * rot
C3 = C4 = A3 * rot
A4 = A5 = C4 * rot
C5 = -A1
tiling.set_initial_tiles([BtileS(A1, B, C1), BtileS(A2, B, C2),
                          BtileS(A3, B, C3), BtileS(A4, B, C4),
                          BtileS(A5, B, C5)])
tiling.make_tiling()

# Makes the matplotlib figure
tiling.make_plot()

tiling.write_svg('pictures/example2.svg')

crop('example2.svg', 0, 0)

webbrowser.open('C:/Users/flynn/PycharmProjects/penrose/pictures/crop_example2.svg')
