
The visco-frictional earthquake cycle model described in 'Loading stress controls the statistics of modelled earthquake cycles on visco-frictional fault zones',
by Adam Beall, Martijn van den Ende, Jean-Paul Ampuero and Ake Fagereng. These scripts are input files for the earthquake cycle simulator QDYN (https://github.com/ydluo/qdyn) and post-processing. The model incorporates a 1D fault which can deform by Newtonian viscous creep or rate-and-state friction, embedded in elastic half-spaces.


The script 'visc.py' runs the earthquake cycle model, where the inputs can be varied to reproduce any of the heterogenous fault models in the study.
Use 'plot.py' to colalte the synthetic earthquake catalogue and visualise the model output.

The script 'visc.py' can also be used to run the single asperity test model, by setting 'visc_dist' to 'single'.


Requirements:
- QDYN (and associated requirements), with the viscosity branch:
git clone https://github.com/ydluo/qdyn.git -b feature/viscosity
- Python
- Matplotlib
- Numpy
- Scipy
