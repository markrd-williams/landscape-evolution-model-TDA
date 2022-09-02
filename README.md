# Analysing Landscape Evolution Models with TDA and vineyards

This is the code for "Using topological data analysis and persistence vineyards as a descriptor for evolving landscapes".

## Dependencies: 
To install the required conda packages to use the environment.yml file included.

```
conda env create --file environment.yml
conda activate LEMvineyard
```

There are some dependencies which are not available on conda.
A custom fork of [Ripser.py](https://github.com/ggrindstaff/ripser.py) is required, maintained by Dr. Gillian Grindstaff. To install, ensure first the conda environment from the previous command is active. Then run the following:

```
git clone https://github.com/ggrindstaff/ripser.py.git
cd ripser.py
pip install -e .
```

To run the Dionysus vineyards code, [Dionysus 1](https://www.mrzv.org/software/dionysus/) (Not version two!), is required.
Specifically, once Dionysus is built, the file `/path/to/dionysus/build/examples/pl-functions/pl-vineyard` must be added to `PATH`.

## Running Guide:

First LEM must be generated. This process is random but results will be similar each run given initial conditions. 
To do this run `generate_LEM.py`.

The code to generate the vineyards and draw the plots can then be run inside `landscape_statistics.ipynb`
