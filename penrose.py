import math
import random
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.path import Path
import matplotlib.patches as patches

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
    A class representing a Robinson triangle and the rhombus formed from it.f

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

    def path(self, rhombus=True):
        """
        Return the SVG "d" path element specifier for the rhombus formed
        by this triangle and its mirror image joined along their bases. If
        rhombus=False, the path for the triangle itself is returned instead.

        """

        AB, BC = self.B - self.A, self.C - self.B
        xy = lambda v: (v.real, v.imag)
        if rhombus:
            return 'm{},{} l{},{} l{},{} l{},{}z'.format(*xy(self.A) + xy(AB) + xy(BC) + xy(-AB))
        return 'm{},{} l{},{} l{},{}z'.format(*xy(self.A) + xy(AB) + xy(BC))

    def get_arc_d(self, U, V, W, proportion1=0.5, half_arc=False, normal_arcs=True):
        """
        Return the SVG "d" path element specifier for the circular arc between
        start and end, with a radius related to proportion.
        If normal_arcs is True then proportion is 0.5 also if half_arc is True,
        the arc is at the vertex of a rhombus; if half_arc is False, the arc is
        drawn for the corresponding vertices of a Robinson triangle.

        If normal_arcs is False then
        """

        if normal_arcs:
            # centre is at vertex U
            centre = U
            proportion2 = proportion1
        else:
            # centre is along the line UV with a ratio phi
            centre = U + (U - V) * phi
            proportion2 = math.sqrt(phi4 / 4 + proportion1 * (proportion1 + 2 * phi)) - phi2 / 2

        # start in on UV
        start = U + (V - U) * proportion1
        if isinstance(self, BtileL) and not normal_arcs:
            # end is on the opposite edge to UV
            end = W + (V - U) * proportion2
        else:
            # end is on UW
            end = U + (W - U) * proportion2
        # arc radius
        r = abs(centre - start)

        if half_arc:
            # Find the endpoint of the "half-arc" terminating on the triangle
            # base
            UN = V + W - 2 * U
            end = U + r * UN / abs(UN)

        # ensure we draw the arc for the angular component < 180 deg
        cross = lambda u, v: u.real * v.imag - u.imag * v.real
        US, UE = start - centre, end - centre
        if cross(US, UE) > 0:
            start, end = end, start
        return 'M {} {} A {} {} 0 0 0 {} {}'.format(start.real, start.imag,
                                                    r, r, end.real, end.imag)

    def arcs(self, proportion=0.5, half_arc=False, normal_arcs=True):
        """
        Return the SVG "d" path element specifiers for the two circular arcs
        about vertices A and C. If half_arc is True, the arc is at the vertex
        of a rhombus; if half_arc is False, the arc is drawn for the
        corresponding vertices of a Robinson triangle.

        """

        D = self.A - self.B + self.C
        if normal_arcs:
            arc1_d = self.get_arc_d(self.A, self.B, D, 0.5, half_arc, normal_arcs)
            arc2_d = self.get_arc_d(self.C, self.B, D, 0.5, half_arc, normal_arcs)
        elif isinstance(self, BtileS):
            arc1_d = self.get_arc_d(self.B, self.C, self.A, proportion, half_arc, normal_arcs)
            arc2_d = self.get_arc_d(D, self.C, self.A, proportion, half_arc, normal_arcs)
        elif isinstance(self, BtileL):
            arc1_d = self.get_arc_d(self.B, self.C, self.A, proportion, half_arc, normal_arcs)
            arc2_d = self.get_arc_d(D, self.C, self.A, proportion, half_arc, normal_arcs)
        else:
            raise ValueError
        return arc1_d, arc2_d

    def get_curve(self, U, V, W, proportion1):
        """
        Takes 3 vertices of a rhombus and calculates the start, control
        and end points for a quadratic bezier curve
        """

        # Centre of the circle
        centre = U + (U - V) * phi
        # Start of the curve
        start = U + (V - U) * proportion1
        # How far end is from a vertices
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

        control = intersection(start, end, centre)

        vertices = [
            (start.real, start.imag),      # Start
            (control.real, control.imag),  # Control point
            (end.real, end.imag)           # End
        ]

        return vertices

    def curves(self, proportion):
        """
        Finds 2 sets of triple points to define 2 Bezier
        curves defined over the rhombus A, B, C, D
        """

        # 4th point of the rhombus
        D = self.A - self.B + self.C
        #
        curve1 = self.get_curve(self.B, self.C, self.A, proportion)
        curve2 = self.get_curve(D, self.C, self.A, proportion)
        return curve1, curve2

    def conjugate(self):
        """
        Return the vertices of the reflection of this triangle about the
        x-axis. Since the vertices are stored as complex numbers, we simply
        need the complex conjugate values of their values.
        """

        return self.__class__(self.A.conjugate(), self.B.conjugate(),
                              self.C.conjugate())


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
                       'normal-arcs': True,
                       'reflect-x': True,
                       'draw-rhombuses': True,
                       'rotate': 0,
                       'flip-y': False, 'flip-x': False,
                       }
        self.config.update(config)
        # And ensure width, height values are strings for the SVG
        self.config['width'] = str(self.config['width'])
        self.config['height'] = str(self.config['height'])

        self.elements = []

    def set_initial_tiles(self, tiles):
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
        sort_elements = sorted(self.elements, key=lambda e: (e.centre().real, e.centre().imag))
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
        if self.config['draw-rhombuses']:
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
        draw_rhombuses = self.config['draw-rhombuses']
        normal_arcs = self.config['normal-arcs']

        # Loop over the rhombuses to draw them
        if self.config['draw-tiles']:
            svg.append('        <g stroke="{}" fill-opacity="{}">'
                       .format(self.config['stroke-colour'], self.config['tile-opacity']))
            for e in self.elements:
                svg.append('            <path fill="{}" d="{}"/>'
                           .format(self.get_tile_colour(e),
                                   e.path(rhombus=draw_rhombuses)))
            svg.append('        </g>')

        # Loop over the rhombuses to draw the arcs
        if self.config['draw-arcs']:
            svg.append('        <g fill="none" stroke="{}" '
                       'stroke-linecap="round">'.format(self.config['arc-colour']))
            for e in self.elements:
                arc1_d, arc2_d = e.arcs(proportion, half_arc=not draw_rhombuses, normal_arcs=normal_arcs)
                svg.append('            <path d="{}"/>'.format(arc1_d))
                svg.append('            <path d="{}"/>'.format(arc2_d))
            svg.append('        </g>')

        svg.append('    </g>\n</svg>')
        return '\n'.join(svg)

    def write_svg(self, filename):
        """ Make and write the SVG for the tiling to filename. """
        svg = self.make_svg()
        with open(filename, 'w') as fo:
            fo.write(svg)

    def make_patch(self, ax, proportion, line_width, colour):
        """
        Makes a matplotlib patch and adds it to ax, part of a figure
        """

        # Clears old paths
        ax.patches = []
        vertices, codes = [], []

        #  Drawing commands
        codes1 = [
            Path.MOVETO,
            Path.CURVE3,
            Path.CURVE3,
        ]

        # Iterate over rhombuses of the tiling, calculate
        # vertices for the paths and extend the collections
        for e in self.elements:
            vertices1, vertices2 = e.curves(proportion)
            vertices += vertices1 + vertices2
            codes += codes1 + codes1

        # Create the path & patch object with the parameters
        path = Path(vertices, codes)
        patch = patches.PathPatch(path, fc='none', ec=colour,
                                  lw=line_width, capstyle='round')
        ax.add_patch(patch)

    def make_plot(self):
        """
        Creates a matplotlib figure of the pattern with sliders to adjust parameters
        """

        initial_tiles = self.elements
        self.make_tiling()

        prop = self.config['proportion']
        colour = self.config['arc-colour']
        gen = self.ngen

        # Configure the figure
        x_min = y_min = -self.scale * self.config['margin']
        x_max = y_max = self.scale * self.config['margin']

        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        plt.axis('equal')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        self.make_patch(ax, proportion=prop, line_width=3, colour=colour)

        # Adds sliders to the figure for proportion and line width
        ax_gen = plt.axes([0.13, 0.05, 0.725, 0.03],
                          facecolor='beige')
        ax_prop = plt.axes([0.13, 0.1, 0.725, 0.03],
                           facecolor='beige')
        ax_width = plt.axes([0.13, 0.15, 0.725, 0.03],
                            facecolor='beige')

        s_gen = Slider(ax_gen, 'Generation', 0, 10, valinit=gen, valstep=1)
        s_prop = Slider(ax_prop, 'Proportion', 0, 2, valinit=prop)
        s_width = Slider(ax_width, 'Width', 0.1, 10, valinit=3)

        # Updates figure when properties are changed
        def update1(val):
            self.ngen = int(s_gen.val)
            self.set_initial_tiles(initial_tiles)
            self.make_tiling()

        s_gen.on_changed(update1)

        def update2(val):
            proportion = s_prop.val
            width = s_width.val
            self.make_patch(ax, proportion=proportion, line_width=width, colour=colour)
            fig.canvas.draw_idle()

        s_gen.on_changed(update2)
        s_prop.on_changed(update2)
        s_width.on_changed(update2)

        # Reset button
        reset_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(reset_ax, 'Reset', color='beige', hovercolor='0.975')

        def reset(event):
            s_gen.reset()
            s_prop.reset()
            s_width.reset()

        button.on_clicked(reset)

        plt.show()

        # Updates values from slider
        self.config['proportion'] = s_prop.val
        self.config['base-stroke-width'] = s_width.val / (psi ** self.ngen * self.scale)
