def rectangle(file, x_pos, y_pos, angle=0, ratio=0.3, scale=90):
    """
    Adds a rectangle to a SVG picture
    """

    # Opens file and copies all the line to a list
    picture = open("pictures/" + file)
    svg_element_list = picture.readlines()

    picture.close()

    # Save last two lines, adds the rectangle then appends the last lines again
    holder = svg_element_list[-2:]
    svg_element_list[-2] = '<g transform="translate({},{})' \
                           ' rotate({})">\n'.format(x_pos, y_pos, angle)
    svg_element_list[-1] = '<rect x="0" y="0" width="{}" height="{}"' \
                           ' fill="none" stroke="red"/>\n'.format(scale, ratio * scale)
    svg_element_list.append(holder[0])
    svg_element_list.append(holder[0])
    svg_element_list.append(holder[1])

    # Open the file and write its new file contents
    picture = open("pictures/rect_" + file, "w")
    new_file_contents = "".join(svg_element_list)
    picture.write(new_file_contents)
    picture.close()


def crop(file, x_pos, y_pos, angle=0, ratio=0.3, scale=90, thickness=0.5):
    """
    Changes the properties of a SVG, such as different
    stroke thickness or cropping
    """

    margin = scale * 0
    height = ratio * scale
    viewBox = '{} {} {} {}'.format(x_pos - margin, y_pos - margin,
                                   scale + 2*margin, height + 2*margin)

    # Opens file and copies all the line to a list
    picture = open("pictures/" + file)
    svg_element_list = picture.readlines()
    picture.close()

    # Replace the relevant lines with new configuration
    svg_element_list[1] = '<svg width="{}" height="{}" viewBox="{}"' \
                          ' preserveAspectRatio="xMidYMid meet" ' \
                          'version="1.1" baseProfile="full"' \
                          ' xmlns="http://www.w3.org/2000/svg">' \
                          '\n'.format(1024, 1024 * ratio, viewBox)
    svg_element_list[4] = '            <rect x="{}" y="{}"' \
                          ' width="{}" height="{}" transform="rotate({})"/>' \
                          '\n'.format(x_pos, y_pos, scale, height, angle)
    svg_element_list[7] = '    <g stroke-width="{}" stroke-linejoin="round" ' \
                          'clip-path="url(#window)" transform="rotate(-{})">' \
                          '\n'.format(thickness, angle)

    # Open the file and write its new file contents
    picture = open("pictures/crop_" + file, "w")
    new_file_contents = "".join(svg_element_list)
    picture.write(new_file_contents)
    picture.close()
