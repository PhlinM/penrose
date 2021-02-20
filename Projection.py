import math
from random import random, seed

"""
Generates a tiling by choosing a plane in
n-dimensional space and projecting selected points
form a hypercube onto the plane.
"""


class Geometry:

    def __init__(self):
        self.fillStyle = 'none'
        self.strokeStyle = 'none'
        self.Path = None

    def begin_path(self):
        # Start a list of path commands
        self.Path = []

    def move_to(self, x, y):
        self.Path.append('M {} {} '.format(x, y))

    def line_to(self, x, y):
        self.Path.append('L {} {} '.format(x, y))

    def close_path(self):
        self.Path.append('Z')

    def arc(self, d):
        self.Path = d


class DeBruijn:

    def __init__(self, n, width=1024, height=768, init=None):

        if init is None:
            # random choice for init
            init = random()
        # Sets init as the seed
        print('Seed: {}'.format(init))
        seed(init)

        # Set limits and make choice for the dimension
        self.min_dimension = 4
        self.max_dimension = 12
        self.doPenrose = False

        diff = self.max_dimension - self.min_dimension + 1
        if n > 1:
            self.dim_n = n
            if self.dim_n == 5:
                # Draw arcs for Penrose tiling
                self.doPenrose = True
        else:
            self.dim_n = self.min_dimension + math.floor(diff * random())

        self.centre = (0, 0)
        self.height = height
        self.width = width

        self.decorationPalette = ['#ff0000', '#00ff00', '#ffffff', '#cccccc']
        self.tile_elements = []
        self.arc_elements = []

        phase = 6 * random()
        sA = math.pi / self.dim_n
        theta = (self.dim_n % 2 + 1) * sA
        norm = math.sqrt(2 / self.dim_n)

        # Compute the necessary margin around the image, by seeing how far the projection
        # of a hypercube vertex can be from the projection of its centre.
        self.margin = 0.5 * norm / math.sin(sA / 2)
        # Colour and increment size
        self.colour_fac = random()
        dpr = 1  # device pixel ratio
        coordWidth = 20
        # coordinate increment for 1 pixel
        self.increment = coordWidth / (math.sqrt(self.dim_n) * 300 * dpr)

        # Generate grid parameters
        # (Grid['x'][i], Grid['y'][i]) is the normal to the lines of grid i, which have equations:
        #
        # (x, y) . (Grid['x'][i], Grid['y'][i]) = C - 0.5 - Grid['d'][i], for C an integer

        # These grid lines are the intersections of the hyperplanes that form the
        # boundaries of the lattice's dual fundamental domains (i.e. hypercubes
        # centred on lattice points), x_i = C - 0.5, with the plane
        # Grid['d'][ ] + x Grid['x'][ ] + y Grid['y'][ ]

        self.Grid = {'x': [norm * math.cos(phase + theta * k) for k in range(self.dim_n)],
                     'y': [norm * math.sin(phase + theta * k) for k in range(self.dim_n)],
                     'd': [random() for _ in range(self.dim_n - 1)]}
        self.Grid['d'].append(self.dim_n / 2 - sum(self.Grid['d']))

        self.Pairs = []
        for fD in range(self.dim_n):
            self.Pairs.append([])
            for eD in range(self.dim_n):
                xf = self.Grid['x'][fD]
                yf = self.Grid['y'][fD]
                xe = self.Grid['x'][eD]
                ye = self.Grid['y'][eD]
                t = math.acos((xf * xe + yf * ye) /
                              math.sqrt((xf * xf + yf * yf) * (xe * xe + ye * ye))) / math.pi

                self.Pairs[fD].append({
                    'shape': 2 * min(t, 1 - t),
                    'orientation': abs(math.atan2(xf, yf) + math.atan2(xe, ye)) / (2 * math.pi),
                    'determinant': ye * xf - xe * yf
                })

    def calculate(self):
        # Clear any elements previous calculated
        self.tile_elements = []
        self.arc_elements = []

        wc = self.centre[0] * self.increment
        half_width = self.width * self.increment * 0.5 + self.margin
        x_margin = [wc - half_width,
                    wc + half_width]
        hc = self.centre[1] * self.increment
        half_height = self.height * self.increment * 0.5 + self.margin
        y_margin = [hc - half_height,
                    hc + half_height]

        x_ep = [None, None]
        y_ep = [None, None]

        # Find all intersections of grid lines, and hence identify
        # all lattice point for which a dual fundamental domain intersects the plane

        # For each family of lines ...
        for eD in range(self.dim_n):
            xe = self.Grid['x'][eD]
            ye = self.Grid['y'][eD]
            de = self.Grid['d'][eD] + 0.5

            # Identify those grid lines that lie partly inside the margins
            lower1 = math.inf  # max number?
            upper1 = -math.inf  # min number?
            for j in range(2):
                for i in range(2):
                    k = math.ceil(xe * x_margin[i] + ye * y_margin[j] + de)
                    if k < lower1:
                        lower1 = k
                    if k - 1 > upper1:
                        upper1 = k - 1

            # For each relevant grid line in the eD family ...
            for eV in range(lower1, upper1 + 1):
                ae = eV - de
                nn = 0

                # Find (x, y) values of endpoints of this grid line, where it hits the margins
                for ii in range(2):
                    x = x_margin[ii]
                    y = (ae - self.Grid['x'][eD] * x) / self.Grid['y'][eD]
                    if y_margin[0] <= y <= y_margin[1]:
                        x_ep[nn] = x
                        y_ep[nn] = y
                        nn += 1
                        if nn == 2:
                            break
                    y = y_margin[ii]
                    x = (ae - self.Grid['y'][eD] * y) / self.Grid['x'][eD]
                    if x_margin[0] <= x <= x_margin[1]:
                        x_ep[nn] = x
                        y_ep[nn] = y
                        nn += 1
                        if nn == 2:
                            break

                # For each other family of grid lines ...
                for fD in range(eD + 1, self.dim_n):
                    pair = self.Pairs[eD][fD]
                    xf = self.Grid['x'][fD]
                    x_fea = xf * ae
                    yf = self.Grid['y'][fD]
                    y_fea = yf * ae
                    df = self.Grid['d'][fD] + 0.5

                    # Identify lowest and highest fD coordinates for current eD line
                    lower2 = math.ceil(xf * x_ep[0] + yf * y_ep[0] + df)
                    upper2 = math.floor(xf * x_ep[1] + yf * y_ep[1] + df)

                    if upper2 < lower2:
                        lower2, upper2 = upper2 + 1, lower2 - 1

                    # For each relevant grid line in the fD family
                    for fV in range(lower2, upper2 + 1):
                        # Find intersection with the eD grid line
                        af = fV - df

                        determinant = pair['determinant']
                        if determinant == 0:
                            print('Determinant was 0')
                            raise ValueError

                        x = (y_fea - ye * af) / determinant
                        y = (xe * af - x_fea) / determinant

                        # Identify fixed coordinates in other dimensions
                        # and perform partial calculations for projection onto plane.
                        xs = ys = 0
                        for m in range(self.dim_n):
                            if m != eD and m != fD:
                                gxm = self.Grid['x'][m]
                                gym = self.Grid['y'][m]
                                gdm = self.Grid['d'][m]
                                z = math.floor(gxm * x + gym * y + gdm + 0.5) - gdm
                                xs += z * gxm
                                ys += z * gym

                        # For odd n, get the sum of the lattice coordinates for the
                        # lowest-sum corner of the selected square
                        c_sum = i_mean = 0
                        if self.doPenrose:
                            for m in range(self.dim_n):
                                gxm = self.Grid['x'][m]
                                gym = self.Grid['y'][m]
                                gdm = self.Grid['d'][m]
                                rev = gxm * x + gym * y + gdm
                                if m == eD or m == fD:
                                    c_sum += round(rev - 0.5)
                                else:
                                    c_sum += round(rev)

                        # Produce points with +/- coordinates in intersection dimensions
                        # to complete projection, and convert to screen coordinates
                        xp, yp, index = [], [], []
                        for m in range(4):
                            e0 = math.floor(m / 2)
                            f0 = math.floor(((m + 1) % 4) / 2)
                            ee = ae + 0.5 - e0
                            ff = af + 0.5 - f0
                            xa = xs + ee * xe + ff * xf
                            ya = ys + ee * ye + ff * yf
                            xp.append(wc + xa / self.increment)
                            yp.append(hc + ya / self.increment)

                            if self.doPenrose:
                                index.append(c_sum + 2 - e0 - f0)
                                i_mean += index[m]

                        i_mean //= 4

                        # Determine tile colour
                        col = self.decorationPalette[i_mean % 4]

                        # Draw tile
                        self.draw_tile(xp, yp, col)

                        if self.doPenrose:
                            self.draw_arc(xp, yp, pair, i_mean, index)

    def draw_tile(self, xp, yp, col):
        tile = Geometry()
        tile.begin_path()
        for m in range(4):
            if m == 0:
                tile.move_to(xp[m], yp[m])
            else:
                tile.line_to(xp[m], yp[m])
        tile.close_path()
        tile.fillStyle = col
        tile.strokeStyle = '#000000'

        self.tile_elements.append(tile)

    def draw_arc(self, xp, yp, pair, i_mean, index):

        paths = []

        thin = pair['shape'] < 0.5
        dx = xp[1] - xp[0]
        dy = yp[1] - yp[0]
        edge2 = dx * dx + dy * dy
        edge = math.sqrt(edge2)
        rgi = [5 - i_mean, 3 * i_mean - 5]

        for k in range(2):
            path = Geometry()
            path.strokeStyle = self.decorationPalette[k]
            radius = 0.3 * edge if thin or k == 1 else 0.7 * edge
            for q in range(4):
                if index[q] == rgi[k]:
                    qp = (q + 1) % 4
                    qm = (q + 3) % 4
                    dxm = xp[qm] - xp[q]
                    dym = yp[qm] - yp[q]
                    dxp = xp[qp] - xp[q]
                    dyp = yp[qp] - yp[q]
                    cp = dxp * dym - dxm * dyp
                    a = math.acos((dxp * dxm + dyp * dym) / edge2)
                    if cp > 0:
                        a0 = math.atan2(dyp, dxp)
                    else:
                        a0 = math.atan2(dym, dxm)

                    start_x = str(xp[q] + radius * math.cos(a0 + a))
                    start_y = str(yp[q] + radius * math.sin(a0 + a))

                    end_x = str(xp[q] + radius * math.cos(a0))
                    end_y = str(yp[q] + radius * math.sin(a0))

                    largeArcFlag = '0' if a <= math.pi else '1'

                    d = ' '.join(["M", start_x, start_y, "A", str(radius),
                                  str(radius), '0', largeArcFlag, '0', end_x, end_y])

                    path.begin_path()
                    path.arc(d)
                    paths.append(path)

                    break

        self.arc_elements.append(paths)

    def make_svg(self, window_width=1000):
        """ Make and return the SVG for the tiling as a str. """

        width = self.width
        height = self.height
        x, y = self.centre

        viewbox = '{} {} {} {}'.format(x - 0.5 * width, y - 0.5 * height, width, height)

        svg = ['<?xml version="1.0" encoding="utf-8"?>',
               '<svg width="{}" height="{}" viewBox="{}"'
               ' xmlns="http://www.w3.org/2000/svg">'
               ''.format(window_width, window_width * height / width, viewbox)]

        # The tiles' stroke widths scale with ngen
        stroke_width = 2

        svg.append('    <g stroke-width="{}" stroke-'
                   'linejoin="round">'.format(stroke_width))

        # Loop over the rhombuses to draw them
        svg.append('        <g>')
        for tile in self.tile_elements:
            string = '          <path '
            string += 'stroke="{}" '.format(tile.strokeStyle)
            string += 'fill="{}" '.format(tile.fillStyle)
            string += 'd="'
            for command in tile.Path:
                string += command
            string += '"/>'
            svg.append(string)
        svg.append('        </g>')

        svg.append('        <g>')
        for arcs in self.arc_elements:
            for arc in arcs:
                string = '        <path '
                string += 'fill="{}" '.format(arc.fillStyle)
                string += 'stroke="{}" d="'.format(arc.strokeStyle)
                string += arc.Path
                string += '"/>'
                svg.append(string)
        svg.append('        </g>')

        svg.append('    </g>\n</svg>')
        return '\n'.join(svg)

    def write_svg(self, filename):
        """ Make and write the SVG for the tiling to filename. """
        svg = self.make_svg()
        with open(filename, 'w') as fo:
            fo.write(svg)
