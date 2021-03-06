import math
import random
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.path import Path
import matplotlib.patches as patches

'''
Some maths constants used throughout the code below
'''
# A small tolerance for comparing floats for equality
TOL = 1.e-5
# psi = 1/phi where phi is the Golden ratio, sqrt(5)+1)/2
phi = (math.sqrt(5) + 1) / 2
psi = (math.sqrt(5) - 1) / 2
# phi**2 = phi + 1, phi**4 = 3*phi + 2
phi2 = phi + 1
phi4 = 3*phi + 2
# psi**2 = 1 - psi
psi2 = 1 - psi


# Used when finding the control points for the Bezier curves for RobinsonTriangles
def intersection(start, end, centre):
    """
    Finds the intersection point from between the two tangents
    at start and end on the circle centred at centre.
    """

    centreToStart = start - centre
    centreToEnd = end - centre
    # direction vectors of the tangents
    m1 = centreToStart.imag - centreToStart.real * 1j
    m2 = centreToEnd.imag - centreToEnd.real * 1j

    # Solving for the parameter in equation intersection = start + t * m1
    denominator = m1.imag * m2.real - m1.real * m2.imag
    numerator = (end.imag - start.imag) * m2.real + (start.real - end.real) * m2.imag

    if start == end:
        # Tangents are coincident
        return start
    elif abs(denominator) == 0:
        # tangents are parallel, so no intersection
        raise ValueError

    return start + (numerator / denominator) * m1


class RobinsonTriangle:
    """
    A class representing a Robinson triangle and the rhombus formed from it.
    """

    def __init__(self, A, B, C):
        """
        Initialize the triangle with the ordered vertices. A and C are the
        vertices at the equal base angles; B is at the vertex angle.

        """

        self.A, self.B, self.C = A, B, C

    def centre(self):
        """
        Return the position of the centre of the rhombus formed from two
        triangles joined by their bases.

        """

        return (self.A + self.C) / 2

    # A method for an SVG format for drawing a tile used to draw all of them in make_svg
    def path(self):
        """
        Return the SVG "d" path element specifier for the rhombus formed
        by this triangle and its mirror image joined along their bases. If
        rhombus=False, the path for the triangle itself is returned instead.

        """

        AB, BC = self.B - self.A, self.C - self.B
        xy = lambda v: (v.real, v.imag)
        return 'm{},{} l{},{} l{},{} l{},{}z'.format(*xy(self.A) + xy(AB) + xy(BC) + xy(-AB))

    # The below two methods are used to get SVG commands for the
    # pattern that overlays the penrose tiling, using sections of circular arcs
    def get_arc(self, U, V, W, proportion1=0.5):
        """
        Returns 4 values: radius, start, control point, end.
        These define the circular arc between start and end, with
        a radius related to proportion1. The control point the
        middle point for the Bezier curve that defines the same arc.
        """

        # Centre of the circle
        centre = U + (U - V) * phi
        # Start of the curve
        start = U + (V - U) * proportion1
        # arc radius
        radius = abs(centre - start)
        # How far end is from a vertex U or W
        proportion2 = -phi2 / 2
        if proportion1 < -2.57:
            proportion2 += math.sqrt(phi4 / 4 + proportion1 * (proportion1 + 2 * phi))
        elif proportion1 > -0.666:
            proportion2 += math.sqrt(phi4 / 4 + proportion1 * (proportion1 + 2 * phi))

        if isinstance(self, BtileL):
            # end is on the opposite edge to UV
            end = W + (V - U) * proportion2
        else:
            # end is on UW
            end = U + (W - U) * proportion2

        # ensure we draw the arc for the angular component < 180 deg
        cross = lambda u, v: u.real * v.imag - u.imag * v.real
        US, UE = start - centre, end - centre
        if cross(US, UE) > 0:
            start, end = end, start

        control = intersection(start, end, centre)

        vertices = [
            radius,                        # Radius
            (start.real, start.imag),      # Start
            (control.real, control.imag),  # Control point
            (end.real, end.imag)           # End
        ]

        return vertices

    def arcs(self, proportion=0.5):
        """
        Return the 4 values that specifies the two circular arcs that
        trace over the rhombus defined over a RobinsonTriangle
        """

        # 4th point of the rhombus defined from a RobinsonTriangle
        D = self.A - self.B + self.C

        arc1 = self.get_arc(self.B, self.C, self.A, proportion)
        arc2 = self.get_arc(D, self.C, self.A, proportion)

        return arc1, arc2

    def conjugate(self):
        """
        Return the vertices of the reflection of this triangle about the
        x-axis. Since the vertices are stored as complex numbers, we simply
        need the complex conjugate values of their values.
        """

        return self.__class__(self.A.conjugate(), self.B.conjugate(),
                              self.C.conjugate())


# The below 2 subclasses are the two different fundamental tiles in the penrose tiling.
# They only contain the two different inflation methods each have
class BtileL(RobinsonTriangle):
    """
    A class representing a "B_L" Penrose tile in the P3 tiling scheme as
    a "large" Robinson triangle (sides in ratio 1:1:phi).

    """

    def inflate(self):
        """
        "Inflate" this tile, returning the three resulting Robinson triangles
        in a list.

        """

        # D and E divide sides AC and AB respectively
        D = psi2 * self.A + psi * self.B
        E = psi2 * self.A + psi * self.C
        # Take care to order the vertices here so as to get the right
        return [BtileL(E, D, self.A),
                BtileS(D, E, self.B),
                BtileL(self.C, E, self.B)]


class BtileS(RobinsonTriangle):
    """
    A class representing a "B_S" Penrose tile in the P3 tiling scheme as
    a "small" Robinson triangle (sides in ratio 1:1:psi).

    """

    def inflate(self):
        """
        "Inflate" this tile, returning the two resulting Robinson triangles
        in a list.

        """
        D = psi * self.A + psi2 * self.B
        return [BtileS(D, self.C, self.A),
                BtileL(self.C, D, self.B)]


class PenroseP3:
    """ A class representing the P3 Penrose tiling. """

    def __init__(self, scale=200, ngen=4, config=None):
        """
        Initialise the PenroseP3 instance with a scale determining the size
        of the final image and the number of generations, ngen, to inflate
        the initial triangles. Further configuration is provided through the
        key, value pairs of the optional config dictionary.

        """

        if config is None:
            config = {}
        self.scale = scale
        self.ngen = ngen

        # Default configuration
        self.config = {'width': '100%', 'height': '100%',
                       'stroke-colour': '#fff',
                       'base-stroke-width': 0.05,
                       'margin': 1.05,
                       'tile-opacity': 0.6,
                       'random-tile-colours': False,
                       'Stile-colour': '#08f',
                       'Ltile-colour': '#0035f3',
                       'arc-colour': '#f00',
                       'draw-tiles': True,
                       'draw-arcs': False,
                       'proportion': 0.7,
                       'reflect-x': True,
                       'rotate': 0,
                       'flip-y': False, 'flip-x': False}
        self.config.update(config)
        # And ensure width, height values are strings for the SVG
        self.config['width'] = str(self.config['width'])
        self.config['height'] = str(self.config['height'])

        self.initial_tiles = []
        self.elements = []

        self.crop = False

    def set_initial_tiles(self, tiles):
        self.initial_tiles = tiles
        self.elements = tiles

    def inflate(self):
        """ "Inflate" each triangle in the tiling ensemble."""
        new_elements = []
        for element in self.elements:
            new_elements.extend(element.inflate())
        self.elements = new_elements

    def remove_dupes(self):
        """
        Remove triangles giving rise to identical rhombuses from the
        ensemble.

        """

        # Triangles give rise to identical rhombuses if these rhombuses have
        # the same centre.
        sort_elements = sorted(self.elements,
                               key=lambda e: (e.centre().real, e.centre().imag))
        self.elements = [sort_elements[0]]
        for i, element in enumerate(sort_elements[1:], start=1):
            if abs(element.centre() - sort_elements[i - 1].centre()) > TOL:
                self.elements.append(element)

    def add_conjugate_elements(self):
        """ Extend the tiling by reflection about the x-axis. """

        self.elements.extend([e.conjugate() for e in self.elements])

    def rotate(self, theta):
        """ Rotate the figure anti-clockwise by theta radians."""

        rot = math.cos(theta) + 1j * math.sin(theta)
        for e in self.elements:
            e.A *= rot
            e.B *= rot
            e.C *= rot

    def flip_y(self):
        """ Flip the figure about the y-axis. """

        for e in self.elements:
            e.A = complex(-e.A.real, e.A.imag)
            e.B = complex(-e.B.real, e.B.imag)
            e.C = complex(-e.C.real, e.C.imag)

    def flip_x(self):
        """ Flip the figure about the x-axis. """

        for e in self.elements:
            e.A = e.A.conjugate()
            e.B = e.B.conjugate()
            e.C = e.C.conjugate()

    def make_tiling(self):
        """ Make the Penrose tiling by inflating ngen times. """

        for gen in range(self.ngen):
            self.inflate()

        self.remove_dupes()
        if self.config['reflect-x']:
            self.add_conjugate_elements()
            self.remove_dupes()

        # Rotate the figure anti-clockwise by theta radians.
        theta = self.config['rotate']
        if theta:
            self.rotate(theta)

        # Flip the image about the y-axis (note this occurs _after_ any
        # rotation.
        if self.config['flip-y']:
            self.flip_y()

        # Flip the image about the x-axis (note this occurs _after_ any
        # rotation and after any flip about the y-axis.
        if self.config['flip-x']:
            self.flip_x()

    def get_tile_colour(self, e):
        """ Return a HTML-style colour string for the tile. """

        if self.config['random-tile-colours']:
            # Return a random colour as '#xxx'
            return '#' + hex(random.randint(0, 0xfff))[2:]

        if isinstance(e, BtileL):
            if hasattr(self.config['Ltile-colour'], '__call__'):
                return self.config['Ltile-colour'](e)
            return self.config['Ltile-colour']

        if hasattr(self.config['Stile-colour'], '__call__'):
            return self.config['Stile-colour'](e)
        return self.config['Stile-colour']

    # make_svg and writ_svg are methods that creates and save an
    # SVG style picture from a generated tiling,
    def make_svg(self):
        """ Make and return the SVG for the tiling as a str. """

        xmin = ymin = -self.scale * self.config['margin']
        width = height = 2 * self.scale * self.config['margin']
        viewbox = '{} {} {} {}'.format(xmin, ymin, width, height)

        svg = ['<?xml version="1.0" encoding="utf-8"?>',
               '<svg width="{}" height="{}" viewBox="{}"'
               ' preserveAspectRatio="xMidYMid meet" version="1.1"'
               ' baseProfile="full" xmlns="http://www.w3.org/2000/svg">'.format(self.config['width'],
                                                                                self.config['height'], viewbox),
               '    <defs>',
               '        <clipPath id="window">',
               '            <rect x="0" y="0" width="10" height="10" transform="rotate(0)"/>',
               '        </clipPath>',
               '    </defs>']

        # The tiles' stroke widths scale with ngen
        stroke_width = str(psi ** self.ngen * self.scale *
                           self.config['base-stroke-width'])

        svg.append('    <g stroke-width="{}" stroke-linejoin="round"'
                   ' clip-path="none">'.format(stroke_width))

        proportion = self.config['proportion']

        # Loop over the rhombuses to draw them
        if self.config['draw-tiles']:
            svg.append('        <g stroke="{}" fill-opacity="{}">'
                       .format(self.config['stroke-colour'], self.config['tile-opacity']))
            for e in self.elements:
                svg.append('            <path fill="{}" d="{}"/>'
                           .format(self.get_tile_colour(e),
                                   e.path()))
            svg.append('        </g>')

        # Loop over the rhombuses to draw the arcs
        if self.config['draw-arcs']:
            svg.append('        <g fill="none" stroke="{}" '
                       'stroke-linecap="round">'.format(self.config['arc-colour']))
            for e in self.elements:
                for curve in e.arcs(proportion):
                    s = 'M {} {} A {} {} 0 0 0 {} {}'.format(curve[1][0], curve[1][1],
                                                             curve[0], curve[0],
                                                             curve[3][0], curve[3][1])
                    svg.append('            <path d="{}"/>'.format(s))

            svg.append('        </g>')

        svg.append('    </g>\n</svg>')
        return '\n'.join(svg)

    def write_svg(self, filename):
        """ Make and write the SVG for the tiling to filename. """
        svg = self.make_svg()
        with open(filename, 'w') as fo:
            fo.write(svg)

    # The below two methods are used to create a matplotlib plot,
    # great for playing with the variables, and finding nice patterns
    def make_patch(self, proportion, width, colour):
        """
        Makes a matplotlib patch
        """

        vertices, codes = [], []

        #  Drawing commands
        codes1 = [
            Path.MOVETO,
            Path.CURVE3,
            Path.CURVE3
        ]

        # Iterate over rhombuses of the tiling, calculate
        # vertices for the paths and extend the collections
        for e in self.elements:
            for curve in e.arcs(proportion):
                vertices += curve[1:]
                codes += codes1

        # Create the path & patch object with the parameters
        path = Path(vertices, codes)
        return patches.PathPatch(path, fc='none', ec=colour,
                                 lw=width, capstyle='round')

    def make_plot(self, x_pos=0, y_pos=0, angle=0, ratio=0.3, scale=90):
        """
        Creates a matplotlib figure of the pattern with sliders
        to adjust parameters. Updates the class properties when
        the plot is closed
        """

        prop = self.config['proportion']
        colour = self.config['arc-colour']

        # Configure the figure: setting axis limits and positioning
        x_min = y_min = -self.scale
        x_max = y_max = self.scale

        fig, ax = plt.subplots(figsize=(7, 6))
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0.20)
        plt.axis('equal')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.axis('off')

        # Toggle crop buttons
        rax = plt.axes([0.03, 0.01, 0.1, 0.04])
        Crop = Button(rax, 'Crop', color='beige', hovercolor='0.975')
        # The path that crops the picture
        rectangle = patches.Rectangle((x_pos, y_pos), width=scale,
                                      height=scale * ratio, angle=angle,
                                      transform=ax.transData)

        def crop_toggle(event):
            '''Toggles crop and changes axis limits'''
            self.crop = not self.crop
            if self.crop:
                ax.set_xlim(x_pos, x_pos + scale)
                ax.set_ylim(y_pos, y_pos + scale * ratio)
            else:
                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_min, y_max)

        Crop.on_clicked(crop_toggle)

        # Deflate button
        ax_prev = plt.axes([0.16, 0.01, 0.1, 0.04])
        b_prev = Button(ax_prev, 'Deflate', color='beige',
                        hovercolor='0.975')

        def ngen_prev(event):
            if self.ngen == 0:
                return
            self.ngen -= 1

        b_prev.on_clicked(ngen_prev)

        # Inflate button
        ax_next = plt.axes([0.29, 0.01, 0.1, 0.04])
        b_next = Button(ax_next, 'Inflate', color='beige',
                        hovercolor='0.975')

        def ngen_next(event):
            if self.ngen == 10:
                return
            self.ngen += 1

        b_next.on_clicked(ngen_next)

        # Updates figure when properties are changed
        def change_tiling(event):
            """ Remakes the tiling when ngen changes """
            self.elements = self.initial_tiles
            self.make_tiling()
            ax.axes.texts = []
            ax.annotate('Generation: {}'.format(self.ngen),
                        xy=(0.4, 0.025), xycoords='figure fraction',
                        fontsize=15)

        b_next.on_clicked(change_tiling)
        b_prev.on_clicked(change_tiling)

        # Adds sliders to the figure for proportion and line width
        ax_prop = plt.axes([0.13, 0.06, 0.725, 0.04],
                           facecolor='beige')
        ax_width = plt.axes([0.13, 0.11, 0.725, 0.04],
                            facecolor='beige')

        s_prop = Slider(ax_prop, 'Proportion', -0.666,
                        2, valinit=prop)
        s_width = Slider(ax_width, 'Width', 0.1, 10,
                         valinit=3)

        # Updates figure when properties are changed
        def update(val):
            """ Recalculates the paths making use of make_patch """
            proportion = s_prop.val
            width = s_width.val
            # Clears old patches
            ax.patches = []
            ax.add_patch(self.make_patch(proportion=proportion,
                                         width=width, colour=colour))
            if self.crop:
                ax.patches[0].set_clip_path(rectangle)

        s_prop.on_changed(update)
        s_width.on_changed(update)
        b_next.on_clicked(update)
        b_prev.on_clicked(update)
        Crop.on_clicked(update)

        # Reset button
        reset_ax = plt.axes([0.8, 0.01, 0.1, 0.04])
        button = Button(reset_ax, 'Reset', color='beige',
                        hovercolor='0.975')

        def reset(event):
            self.ngen = 5
            s_prop.reset()
            s_width.reset()

        button.on_clicked(reset)
        button.on_clicked(change_tiling)
        button.on_clicked(update)

        ax.add_patch(self.make_patch(proportion=prop,
                                     width=3, colour=colour))

        ax.annotate('Generation: {}'.format(self.ngen),
                    xy=(0.42, 0.01), xycoords='figure fraction',
                    horizontalalignment='left', verticalalignment='bottom',
                    fontsize=15)

        plt.show()

        # Updates values from slider
        self.config['proportion'] = s_prop.val
        self.config['base-stroke-width'] = s_width.val / (psi ** self.ngen * self.scale)
