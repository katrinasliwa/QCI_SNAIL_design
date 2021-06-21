import os
import copy
import json
from collections import OrderedDict
import numpy as np
import gdspy

from . import utils
from .layout import LayoutGroup
from .design import DesignProperties
from .sweeper import Sweeper, SweeperProperties

gds_name = 'module_transmons'

working_dir = os.path.dirname(__file__)
config_dir = os.path.join(working_dir, 'configs')

def extract_design_properties(fn):
    params = utils.process_json(fn)
    return DesignProperties(params)

def extract_sweep_properties(fn, sweep="sweep"):
    params = utils.process_json(fn)

    sweep_params_properties = params[sweep]["properties"]
    sweep_params_layers = params[sweep]["layers"]
    layout_properties = params["layout"]
    sweep_properties = SweeperProperties(parameters=sweep_params_properties, layers=sweep_params_layers)

    return sweep_properties, layout_properties

def extract_design_sweep_properties(fn, sweep="sweep", layers_fn=None):
    params = utils.process_json(fn)

    design_properties = params["design"]
    sweep_params_properties = params[sweep]["properties"]
    sweep_params_layers = params[sweep]["layers"]
    layout_properties = params["layout"]

    design_properties = DesignProperties(design_properties, layers_fn=layers_fn)
    sweep_properties = SweeperProperties(parameters=sweep_params_properties, layers=sweep_params_layers, layers_fn=layers_fn)

    return design_properties, sweep_properties, layout_properties

def generate_sweep_layout(name, design, fn, **kwargs):
    """ Helper method to generate and layout.

    Dictionary format: 
    {
        "design": {
            ...
        },
        "sweep": {
            "properties": {
                ...
            },
            "layers": {
                ...
            }
        },
        "layout": {
            ...
        }
    }
    
    Parameters
    ----------
    name : [str]
        [description]
    design : [Design]
        [description]
    fn : [str]
        filename to configuration
    
    Returns
    -------
    [Sweeper]
        Returns Sweeper
    """

    design_properties, sweep_properties, layout_properties = extract_design_sweep_properties(fn, **kwargs)

    layers_fn = kwargs.pop('layers_fn', None)

    sweep = Sweeper(name, design, design_properties, layers_fn=layers_fn)
    sweep.generate_variations(sweep_properties)
    sweep.layout_variations(layout_properties)

    return sweep

def generate_process_monitor_squares(name, fn, **kwargs):
    from mask_maker_v3.process_monitor_squares import ProcessMonitor_Squares
    return generate_sweep_layout(name, ProcessMonitor_Squares, fn, **kwargs)

def generate_transmons(name, fn, **kwargs):
    from mask_maker_v3.transmon import Transmon
    return generate_sweep_layout(name, Transmon, fn, **kwargs)

def generate_transmon_chips(name, fn, **kwargs):
    from mask_maker_v3.transmon_chip import TransmonChip
    return generate_sweep_layout(name, TransmonChip, fn, **kwargs)

def generate_transmon_chips_filter(name, fn, **kwargs):
    from mask_maker_v3.transmon_chip import TransmonChip_Filter
    return generate_sweep_layout(name, TransmonChip_Filter, fn, **kwargs)

def generate_alignment_marks(name, fn, **kwargs):
    from mask_maker_v3.alignment_marks import AlignmentMarks
    return generate_sweep_layout(name, AlignmentMarks, fn, **kwargs)

def generate_transmon_witness(name, fn_chip, fn_witness, **kwargs):
    from mask_maker_v3.transmon import Transmon

    layers_fn = kwargs.pop('layers_fn', None)

    # for this one we need update the transmon configuration with the chip layout
    design_properties, sweep_properties, _ = extract_design_sweep_properties(fn_chip, layers_fn=layers_fn, **kwargs)

    # extract only junction sweep, ignore everything else
    sweep_properties = sweep_properties.get_component('transmon')
    sweep_properties.properties = OrderedDict()

    # get witness junction DesignProperties
    witness_params = utils.process_json(fn_witness)
    design_properties_witness = DesignProperties(witness_params['design'], layers_fn=layers_fn)

    layout_properties = witness_params['layout']

    sweep = Sweeper(name, Transmon, design_properties_witness, layers_fn=layers_fn)
    sweep.generate_variations(sweep_properties)
    sweep.layout_variations(layout_properties)

    return sweep

def get_transmon_chip_spacing(fn_chip, layers_fn=None):
    """ Get distance between two chips (sum of chip_width and the street_width) """
    design_properties, _, layout_properties = extract_design_sweep_properties(fn_chip, layers_fn=layers_fn)
    chip_width = design_properties.get_property("chip_width")
    street_width = layout_properties['layout']['pad_y']
    offset = chip_width + street_width
    return offset

def update_layout_spacing(fn_chip, fn_witness, layers_fn=None):
    """ Update the layout spacing for witness junctions according to chip and street width.
        Used for witness junctions located between chips """

    offset = get_transmon_chip_spacing(fn_chip, layers_fn=layers_fn)
    witness_json = utils.load_json(fn_witness)
    witness_json['layout']['layout']['dy'] = offset
    with open(fn_witness, 'w') as fp:
        json.dump(witness_json, fp, indent=2)

def generate_dose_test_1d(names, fns, dx=0, dy=800, layout_name='dose_test', layers_fn=None):

    dose_test_list = []
    for name, fn in zip(names, fns):
        dose_test_list.append(generate_transmons(name, fn, layers_fn=layers_fn))
    
    dose_test_layout = LayoutGroup(layout_name)
    for idx, d in enumerate(dose_test_list):
        dose_test_layout.add(d.layout, (dx*idx, dy*idx))
    dose_test_layout.create()

    return dose_test_layout

def generate_wafer_outline(name, fn, **kwargs):
    from mask_maker_v3.wafer import Wafer
    wafer_params = utils.load_json(fn)
    wafer = Wafer('wafer', wafer_params, **kwargs)
    wafer.create()
    return wafer

def generate_rotated_wafer(layout, suffix="rotated"):
    # rotate by 90 deg for the EBPG
    name = layout.name + "_%s" % (suffix)
    top_cell_rotated = gdspy.Cell(name)
    top_cell_rotated.add(gdspy.CellReference(layout.cell, rotation=-90))
    return top_cell_rotated

def generate_layer_doses(layout, required_layers=['fine_bias_layer'], fn=None, layers=None):

    layer_doses = layout.get_layer_doses()
    if layers is None:
        layers = utils.load_layers(fn=fn)

    for required_layer in required_layers:
        layer_number, dose_number = layers[required_layer]
        layer_doses[layer_number] = dose_number

    # fix layer order
    layer_doses_ordered = OrderedDict()
    for layer_number in sorted(layer_doses.keys()):
        layer_doses_ordered[layer_number] = layer_doses[layer_number]

    return layer_doses_ordered
