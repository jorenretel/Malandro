'''Everything that has to do with color in the GUI.
   The qualitative 12 colors are directly taken from
   colorbrewer: http://colorbrewer2.org/
'''


def pick_color_by_percentage(percentage, asHex=True):
    '''Get a color from a predefined scale for a percentage.
           args: percentage: the percentage on the scale
                 asHex:      boolean: if true return hex (default)
                                      else: return rgb tuple
    '''

    if not 0 <= percentage <= 100:
        raise ValueError('percentage should be between 0 and 100')

    RGB = pick_color_from_scale(percentage,
                                _standardBackGroundColorRGB,
                                _niceGreenRGB)

    if asHex:
        return rgb_to_hex(RGB)
    else:
        return RGB


def pick_color_from_scale(percentage, minRGB, maxRGB):
    '''Pick color from scale.
           args: percentage: percentage on scale
                 minRGB:     color corresponding to 0%
                 maxRGB:     color corresponding to 100%
           returns: rgb tuple
    '''

    percentage = float(percentage)
    newRGB = []

    for minimal, maximal in zip(minRGB, maxRGB):

        colorPart = int(minimal + (percentage / 100.0) * (maximal - minimal))
        newRGB.append(colorPart)

    return tuple(newRGB)


def pickColorByPercentage(percentage):
    '''deprecated'''
    percentage = float(percentage)
    if percentage < 1:

        return 'grey83'

    if percentage > 80:

        red = 0
        green = int(percentage / 100.0 * 255.0)
        blue = int(255.0 - (percentage / 100.0 * 255.0))

    elif percentage > 50:

        green = 0
        blue = int(percentage / 80.0 * 255.0)
        red = int(255.0 - (percentage / 80.0 * 255.0))

    else:

        red = 255
        green = 0
        blue = 0

    red = rgb_to_hex(red)
    green = rgb_to_hex(green)
    blue = rgb_to_hex(blue)

    color = '#' + red + green + blue
    return color


def rgb_to_hex(rgb):
    '''Returns hex for rgb iterable'''

    return '#' + ''.join([format(x, '02x') for x in rgb])


def grey_scale(rgb):
    '''Get grey scale for rgb tuple or list'''
    return tuple([int(((rgb[0] * 299) + (rgb[1] * 587) + (rgb[2] * 114)) / 1000)] * 3)


highLightRed = rgb_to_hex((227, 26, 28))
highLightYellow = rgb_to_hex((253, 191, 111))
_standardBackGroundColorRGB = (212, 212, 212)  # This is grey83
_niceGreenRGB = (0, 130, 60)

#_niceGreenRGB = (51,160,44)

# http://colorbrewer2.org/ qualitative 12 colors
colorSeries = ['#a6cee3',
               '#1f78b4',
               '#b2df8a',
               '#33a02c',
               '#fb9a99',
               '#e31a1c',
               '#fdbf6f',
               '#ff7f00',
               '#cab2d6',
               '#6a3d9a',
               '#ffff99',
               '#b15928']
