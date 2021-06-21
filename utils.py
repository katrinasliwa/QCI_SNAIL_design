import numpy as np
import os
import json
import re
import math
from collections import OrderedDict

import gdspy
from pint import UnitRegistry
ureg = UnitRegistry()

DEFAULT_UNIT = 'um'
LAYERS_FILE = os.path.join(os.path.dirname(__file__), 'configs', 'layers.json')

class Constants:
    """ Base unit is um """   
    # nm = 0.001
    # mm = 1000.0
    # um = 1.0
    # mil = 25.4
    # inches = mil*1000
    nm = ureg('nm')
    mm = ureg('mm')
    um = ureg('um')
    mil = ureg('mil')
    inches = ureg('inches')

def load_json(fn):
    """ Load json from file. Return as dictionary. """
    with open(fn) as f:
        json_file = json.load(f, object_pairs_hook=OrderedDict)
    return json_file

def load_layers(fn=None):
    # try:
    #     print ("attempting to load from local config directory")
    #     if fn is None:
    #         fn = os.path.join('.', 'configs', 'layers.json')
    #     return load_json(fn)
    # except:
    #     if fn is None:
    #         fn = LAYERS_FILE
    #     print("Using default layers mapping from %s" % (LAYERS_FILE))
    #     return load_json(fn)
    if fn is None:
        raise ValueError("specify layers config file")
        # fn = LAYERS_FILE
    return load_json(fn)

def parse_value_unit(val_str, default_unit=DEFAULT_UNIT):
    """ Convert dimensioned value string to float """
    return ureg(val_str).to(DEFAULT_UNIT).magnitude

# def translate_dict(in_dict):
#     out_dict = OrderedDict()
#     for k, v in in_dict.items():
#         # if re.search(r'\d', v): # check if there are any numbers in the string
#         if v[0].isdigit():
#             uv = ureg(v)
#             if type(uv) in [int, float]: 
#                 out_dict[k] = uv
#             else:
#                 out_dict[k] = parse_value_unit(v)
#         else:
#             out_dict[k] = v
#     return out_dict

def translate_dict(in_dict):
    out_dict = OrderedDict()
    for k, v in in_dict.items():
        # attempt to convert numbers using pint
        try:
            out_dict[k] = parse_value_unit(v)
        except:
            out_dict[k] = v
    return out_dict

def process_json(fn):
    json_file = load_json(fn)
    return translate_dict(json_file)

################################################

def generate_dose_file(layers_dose, fn, relative_dose=100.0):
    layers_dose_list = [(k, float(v)/relative_dose) for k,v in layers_dose.items()]
    np.savetxt(fn, layers_dose_list, fmt='%u(0), %.3f')

################################################

def get_bounding_box_list(polygon_list):
    """ Return the bounding box for a list of Polygon objects """
    [[global_x_min, global_y_min], [global_x_max, global_y_max]] = [[0, 0], [0, 0]]
    for p in polygon_list:
        [[x_min, y_min], [x_max, y_max]] = p.get_bounding_box()
        
        global_x_min = x_min if x_min < global_x_min else global_x_min
        global_y_min = y_min if y_min < global_y_min else global_y_min
        global_x_max = x_max if x_max > global_x_max else global_x_max
        global_y_max = y_max if y_max > global_y_max else global_y_max
    return [[global_x_min, global_y_min], [global_x_max, global_y_max]]

def get_cell_center(cell):
    [[x_min, y_min], [x_max, y_max]] = cell.get_bounding_box()
    cu_x, cu_y = (0.5*(x_max + x_min), 0.5*(y_max + y_min))
    return cu_x, cu_y

def get_cell_size(cell):
    [[x_min, y_min], [x_max, y_max]] = cell.get_bounding_box()
    sizex, sizey = ((x_max - x_min), (y_max - y_min))
    return sizex, sizey

################################################

def unite(p1, p2, layer=0): 
    # TODO make layer selection more intuitive
    return gdspy.fast_boolean(p1, p2, 'or', layer=layer, precision=1e-8)

def unite_list(ps, layer=0):
    assert (len(ps) > 1)
    out = ps[0]
    for p in ps[1:]:
        out = unite(out, p, layer=layer)
    return out

def intersect(p1, p2, layer=0): 
    return gdspy.fast_boolean(p1, p2, 'and', layer=layer, precision=1e-8)

def subtract(p1, p2, layer=0): 
    return gdspy.fast_boolean(p1, p2, 'not', layer=layer, precision=1e-8)

################################################

def draw_rectangle(center, xsize, ysize, layer=0):
    x0, y0 = center
    lower_left = (x0-0.5*xsize, y0-0.5*ysize)
    upper_right = (x0+0.5*xsize, y0+0.5*ysize)
    return gdspy.Rectangle(lower_left, upper_right, layer=layer)

def draw_circle(center, radius, layer=0):
    return gdspy.Round(center, radius, layer=layer, number_of_points=50)

#"""Funnel subroutines, taken from https://gist.github.com/Alquimista/1274149"""
def binomial(i, n):
#    """Binomial coefficient"""
    return math.factorial(n) / float(
        math.factorial(i) * math.factorial(n - i))

def bernstein(t, i, n):
#    """Bernstein polynom"""
    return binomial(i, n) * (t ** i) * ((1 - t) ** (n - i))

def bezier(t, points):
#    """Calculate coordinate of a point in the bezier curve"""
    n = len(points) - 1
    x = y = 0
    for i, pos in enumerate(points):
        bern = bernstein(t, i, n)
        x += pos[0] * bern
        y += pos[1] * bern  
    return x, y

def bezier_curve_range(n, points):
#    """Range of points in a curve bezier"""
    for i in np.arange(n):
        t = i / float(n - 1)
        yield bezier(t, points)

def funnel(start_width, end_width,layer=0, length = 500,steepness = 0.5,steps =50):
    sw = start_width
    ew = end_width

    steps = 50
    x = np.zeros(steps*2+1)
    y = np.zeros(steps*2+1)

    control_points = [(0,sw/2),
                      (steepness*length,sw/2),
                      ((1-steepness)*length,ew/2),
                      (length,ew/2)]

    #Draw the first curve
    for idx,pt in enumerate(bezier_curve_range(steps, control_points)):
        x[idx] = pt[0]
        y[idx] = pt[1]

    #Draw the second curve (backwards, to keep the points in the right order)
    for idx,pt in enumerate(bezier_curve_range(steps, control_points)):
        x[2*steps-idx-1] = pt[0]
        y[2*steps-idx-1] = -pt[1]

    #Finish the plg.
    x[-1] = x[0]
    y[-1] = y[0]

    pts = list(np.zeros(2*steps+1))
    for idx in range(steps*2+1):
        pts[idx] = (x[idx],y[idx])

    plg = gdspy.Polygon(pts,layer=layer)
    return plg

