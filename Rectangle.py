"""
Adds a rectangle to a SVG picture
"""

# Rectangle configuration
file = 'example1.svg'
x_pos, y_pos, angle, ratio, scale = -45, -13.5, 0, 0.3, 90

# Opens file and copies all the line to a list
picture = open("pictures/" + file)
svg_element_list = picture.readlines()

picture.close()

# Save last two lines, adds the rectangle then appends the last lines again
holder = svg_element_list[-2:]
svg_element_list[-2] = '<g transform="translate({},{})' \
                       ' rotate({})">\n'.format(x_pos, y_pos, angle)
svg_element_list[-1] = '<rect x="0" y="0" width="{}" height="{}"' \
                       ' fill="none" stroke="red"/>\n'.format(scale, ratio*scale)
svg_element_list.append(holder[0])
svg_element_list.append(holder[0])
svg_element_list.append(holder[1])

# Open the file and write its new file contents
picture = open("pictures/rect_" + file, "w")
new_file_contents = "".join(svg_element_list)
picture.write(new_file_contents)
picture.close()
