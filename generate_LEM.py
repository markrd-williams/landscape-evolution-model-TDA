"""
This file runs the fastscape models prepared, and saves the run to disk.
"""

import numpy as np
import xsimlab as xs

import fastscape
from fastscape.models import (sediment_model, basic_model)
from fastscape.processes import (LinearDiffusion, StreamPowerChannel,
                                 TotalErosion, UniformRectilinearGrid2D)

import warnings

warnings.filterwarnings("ignore", category=PendingDeprecationWarning)

mountain_in_ds = xs.create_setup(
        model=sediment_model,
        clocks={
            'time': np.arange(0, 2e7 + 1, 2e4),
            'out': np.arange(0, 2e7 + 1, 2e4),
            },
        master_clock='time',
        input_vars={
            'grid__shape': [201, 201],
            'grid__length': [1e5, 1e5],
            'boundary__status': ['fixed_value', 'fixed_value', 'fixed_value', 'fixed_value'],
            'flow__slope_exp': 1.,
            'spl': {
                'k_coef_bedrock': 1e-4,
                'k_coef_soil': 1.5e-4,
                'g_coef_bedrock': 1.,
                'g_coef_soil': 1.,
                'area_exp': 0.4,
                'slope_exp': 1.
                },
            'uplift__rate': 1e-3,
            'diffusion': {
                'diffusivity_bedrock': 1e-2,
                'diffusivity_soil': 1.5e-2
                }
            },
        output_vars={
            'topography__elevation': 'out',
            'drainage__area': 'out',
            'erosion__rate': 'out'
            }
        )


@xs.process
class DippingDyke:
    """Mimics the effect on erosion rates of a dyke dipping at
    a given angle, that is buried beneath the landscape and that is
    progressively exhumed by erosion.
    
    {{attribute_section}}
    
    """
    x_position = xs.variable(description='initial x-position of exposed dyke')
    width = xs.variable(description='dyke fixed width')
    angle = xs.variable(description='dyke dipping angle in degrees')
    
    grid_shape = xs.foreign(UniformRectilinearGrid2D, 'shape')
    x = xs.foreign(UniformRectilinearGrid2D, 'x')
    
    etot = xs.foreign(TotalErosion, 'cumulative_height')
    
    k_coef = xs.foreign(StreamPowerChannel, 'k_coef', intent='out')
    diffusivity = xs.foreign(LinearDiffusion, 'diffusivity', intent='out')
    
    def run_step(self):
        cotg = 1. / np.tan(np.radians(self.angle))
        
        dyke_pos = self.x - self.x_position - self.etot * cotg
        
        in_dyke = (dyke_pos - self.width) * (dyke_pos + self.width) <= 0
        
        self.k_coef = np.where(in_dyke, 1e-5, 2e-5)
        self.diffusivity = np.where(in_dyke, 0.05, 0.1)
        
        
dyke_model = basic_model.update_processes({'dyke': DippingDyke})
dyke_in_ds = xs.create_setup(
        model=dyke_model,
        clocks={
            'time': np.arange(0, 2e7 + 1, 2e4),
            'out': np.arange(0, 2e7 + 1, 2e4),
            },
        master_clock='time',
        input_vars={
            'grid__shape': [201, 201],
            'grid__length': [1e5, 1e5],
            'dyke': {
                'x_position': 1e4,
                'width': 2e3,
                'angle': 30.
                },
            'uplift__rate': 1e-3,
            },
        output_vars={
            'topography__elevation': 'out',
            'terrain__slope': 'out',
            }
        )

mountain_out_ds = mountain_in_ds.xsimlab.run(model=sediment_model)
mountain_out_ds['border'] = mountain_out_ds.border.astype("S6")
mountain_out_ds.to_netcdf(f"./LEM/mountain.nc")
mountain_in_ds.to_netcdf(f"./LEM/info/mountainparams.nc")

dyke_out_ds = dyke_in_ds.xsimlab.run(model=dyke_model)
dyke_out_ds['border'] = dyke_out_ds.border.astype("S6")
dyke_out_ds.to_netcdf(f"./LEM/dyke.nc")
dyke_in_ds.to_netcdf(f"./LEM/info/dykeparams.nc")
