import sys
from pathlib import Path

try:
    from fractaldna.dna_models import dnachain as dna
except (ImportError, ModuleNotFoundError):
    sys.path.append(str(Path.cwd().parent.parent.parent))
    from fractaldna.dna_models import dnachain as dna

import matplotlib.pyplot as plt
from mayavi import mlab

# Disable this option for interactive rendering
mlab.options.offscreen = True

# Enable this option for an interactive notebook
# mlab.init_notebook()
# Make a straight solenoid in a 750 Å box
chain = dna.Solenoid(voxelheight=750)
# MayaVI plots are best for visualisation here
plot = chain.to_strand_plot()

# Save the figure
plot.scene.save_jpg("single_solenoid_strand_plot.jpg")

# In an interactive notebook, you can refer to the figure to
# interact with it.
# plot
plot = chain.to_line_plot()
plot.scene.save_jpg("single_solenoid_line_plot.jpg")