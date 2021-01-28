"""
Changes the properties of a SVG, such as different
stroke thickness or cropping
"""
# Name of the picture to be edited within the directory pictures/
file = 'example2.svg'

# Cropping rectangle configuration
x_pos, y_pos, angle, ratio, scale = -20, -13.5, 7, 0.3, 90
thickness = 1

margin = scale*0.02
height = ratio * scale
viewBox = '{} {} {} {}'.format(x_pos - margin, y_pos - margin,
                               scale + 2*margin, height + 2*margin)

# Opens file and copies all the line to a list
picture = open("pictures/" + file)
svg_element_list = picture.readlines()
picture.close()

# Replace the relevant lines with new configuration
svg_element_list[1] = '<svg width="{}" height="{}" viewBox="{}" ' \
                      'preserveAspectRatio="xMidYMid meet" version="1.1" baseProfile="full"' \
                      ' xmlns="http://www.w3.org/2000/svg">'.format(scale + 2*margin, height + 2*margin, viewBox)
svg_element_list[4] = '            <rect x="{}" y="{}"' \
                      ' width="{}" height="{}" transform=' \
                      '"rotate({})"/>\n'.format(x_pos, y_pos, scale, height, angle)
svg_element_list[7] = '    <g style="stroke: white; stroke-width: {}; ' \
                      'stroke-linejoin: round;" clip-path="url(#window)" ' \
                      'transform="rotate(-{})">\n'.format(thickness, angle)

# Open the file and write its new file contents
picture = open("pictures/crop_" + file, "w")
new_file_contents = "".join(svg_element_list)
picture.write(new_file_contents)
picture.close()
