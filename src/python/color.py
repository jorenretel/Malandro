

def pick_color_by_percentage(percentage, asHex=True):
  
  standardBackGroundColorRGB = (212,212,212) #This is grey83
  #niceGreenRGB = (49,163,84)
  #niceGreenRGB = (0,69,41)
  niceGreenRGB = (0,130,60)
  RGB = pick_color_from_scale(percentage, standardBackGroundColorRGB, niceGreenRGB)
  
  #RGB = grey_scale(RGB)
  
  if asHex:
    return rgb_to_hex(RGB)
  else:
    return RGB
  
  
def pick_color_from_scale(percentage, minRGB, maxRGB):
  
  percentage = float(percentage)
  
  newRGB = []
  
  for minimal, maximal in zip(minRGB, maxRGB):
  
    colorPart = int(minimal + (percentage / 100.0) * (maximal - minimal))
    
    newRGB.append(colorPart)

  return tuple(newRGB)
  

def pickColorByPercentage(percentage):

  percentage = float(percentage)
  
  if percentage < 1 :
    
    return 'grey83'
  
  if percentage > 80 :

    red = 0
    green = int(percentage/100.0 * 255.0)
    blue = int(255.0 - (percentage/100.0 * 255.0))
    

  elif percentage > 50 :

    green = 0 
    blue = int(percentage/80.0 * 255.0)
    red = int(255.0 - (percentage/80.0 * 255.0))

  else :

    red = 255
    green = 0
    blue = 0

  red = self.rgbToHex(red)
  green = self.rgbToHex(green)
  blue = self.rgbToHex(blue)

  color = '#' + red + green + blue
  return color
  
def rgb_to_hex(rgb):
  
  return '#' + ''.join([format(x,'02x') for x in rgb])
      
def grey_scale(rgb):
  
  return tuple([int(((rgb[0] * 299) + (rgb[1] * 587) + (rgb[2] * 114)) / 1000)]*3)