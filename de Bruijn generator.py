from Projection import DeBruijn
import webbrowser

n = 5

something = DeBruijn(n)

something.calculate()

something.write_svg('Pattern.svg')
webbrowser.open('C:/Users/flynn/PycharmProjects/penrose/Pictures/Pattern.svg')
