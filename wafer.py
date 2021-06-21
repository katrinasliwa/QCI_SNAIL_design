import numpy as np
import json
import gdspy

from .import utils
from .design import Design

class Wafer(Design):

    def create(self, **cell_kwargs):
        self.make_cell(**cell_kwargs)
        cell = self.cell

        # extract attributes
        wafer_size = self.get_property('wafer_size')
        flat_length_primary = self.get_property('flat_length_primary')
        flat_length_secondary = self.get_property('flat_length_secondary')
        frame_width = self.get_property('frame_width')

        layer = self.get_layer_number('layer')

        # calculate flat distance
        if flat_length_primary != 0:
            flat_distance_primary = np.sqrt(wafer_size**2-flat_length_primary**2)/2 
        else:
            flat_distance_primary = wafer_size

        if flat_length_secondary != 0:
            flat_distance_secondary = np.sqrt(wafer_size**2-flat_length_secondary**2)/2
        else:
            flat_distance_secondary = wafer_size

        # create frame outline with flat
        frame_outer = gdspy.Round([0,0], wafer_size*0.5, number_of_points = 144)
        frame_outer = utils.intersect(frame_outer, 
                                gdspy.Rectangle((-wafer_size, -flat_distance_primary),
                                                (wafer_size, flat_distance_secondary)))

        frame_inner = gdspy.Round([0,0], (wafer_size*0.5-frame_width), number_of_points = 144)
        frame_inner = utils.intersect(frame_inner, 
                            gdspy.Rectangle((-wafer_size, -(flat_distance_primary-frame_width)),
                                            (wafer_size, (flat_distance_secondary-frame_width))))

        frame = utils.subtract(frame_outer, frame_inner, layer=layer)

        cell.add(frame)

