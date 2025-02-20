#!/usr/bin/env python

import cairo

# Set up an image surface (PNG format, 300x200 pixels)
width, height = 300, 200
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
context = cairo.Context(surface)

# Set background color (white)
context.set_source_rgb(1, 1, 1)  
context.paint()

# Draw a red rectangle
context.set_source_rgb(1, 0, 0)  # RGB (red)
context.rectangle(140, 25, 100, 150)  # x, y, width, height
context.move_to(70, 25)
context.line_to(70, 175)
context.stroke()

# Save to file
surface.write_to_png("number_10_pycairo.png")

print("Image saved as number_10_pycairo.png")

