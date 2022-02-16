"""

Visco-frictional earthquake cycle simulation, using QDYN

The key input paramaters are:
    shear-zone width - W
    max. viscosity contrat - visc_contrast
    viscosity distribution type - visc_dist

The output is then collated and visualised using plot.py


Written by Adam Beall and Martijn van Den Ende

"""

import matplotlib.pyplot as plt
import numpy as np

from pyqdyn import qdyn
import plot_functions as qdyn_plot

def load():

    params = {}
    # ---- Key paramaters
    # Shear-zone thickness
    params['W'] = 100.0
    params['visc_contrast'] = 100.0
    # seed for patch size and viscosity random generation 
    np.random.seed(3)  # model-set A
    #np.random.seed(5)  # model-set B
    #np.random.seed(10)  # model-set C


    # Model duration
    params['maxtime'] = 1000 * 3600 * 24 * 365.0    # model time in seconds

    # options are log, powerlaw or bimodal
    params['visc_dist'] = 'log' 




    # Instantiate the QDYN class object
    p = qdyn()
    # Get the settings dict
    set_dict = p.set_dict

    # Global simulation parameters
    set_dict["MESHDIM"] = 1        # Simulation dimensionality (1D fault in 2D medium)
    set_dict["FINITE"] = 1         # Finite fault
    set_dict["TMAX"] = params['maxtime']     # Maximum simulation time [s]
    set_dict["NTOUT"] = 50       # Save output every N steps
    set_dict["NXOUT"] = 1          # Snapshot resolution (every N elements)
    set_dict["V_PL"] = 1e-9        # Plate velocity [m/s]
    set_dict["MU"] = 3e10          # Shear modulus [Pa]
    set_dict["SIGMA"] = 100e6        # Effective normal stress [Pa]
    set_dict["ACC"] = 1e-7         # Solver accuracy
    set_dict["SOLVER"] = 2         # Solver type (Runge-Kutta
    set_dict["FAULT_TYPE"] = 2



    params['res'] = 1024     # Mesh resolution
    params['L'] = 30e3      # Fault length
    set_dict["N"] = params['res']       
    set_dict["L"] = params['L']            

    # Location on the fault (middle) for time series output
    set_dict["IC"] = params['res'] // 2

    # RSF parameter values
    set_dict["SET_DICT_RSF"]["RNS_LAW"] = 2
    set_dict["SET_DICT_RSF"]["A"] = 0.009    # Direct effect (will be overwritten later)
    set_dict["SET_DICT_RSF"]["B"] = 0.02      # Evolution effect
    set_dict["SET_DICT_RSF"]["DC"] = 1e-2     # Characteristic slip distance
    set_dict["SET_DICT_RSF"]["V_SS"] = set_dict["V_PL"]    # Reference velocity [m/s]
    set_dict["SET_DICT_RSF"]["V_0"] =  1e-7         #1e-5
    set_dict["SET_DICT_RSF"]["TH_0"] =  set_dict["SET_DICT_RSF"]["DC"] / set_dict["V_PL"]    # Initial state [s]
    set_dict["SET_DICT_RSF"]["MU_SS"] = 0.6   


    # Set parameter values and generate mesh
    p.settings(set_dict)
    p.render_mesh()

    return params,p



if __name__ == '__main__':

    params,p = load()
    # Override mesh values

    # Generate random patch sizes

    # power-law exponent
    D = -1.0

    wl_min = 100   # min. patch width
    wl_max = 1e3   # max. patch width

    # redundant random number generator, but needed to get same seed as in paper
    np.random.random(int(1e5))


    # array of random patch widths
    xrand = np.random.random(200)
    arrW = (wl_min**D + (wl_max**D - wl_min**D)*xrand)**(1./D)


    # Generate random viscosities

    eta_min = 18   # min. log(eta)
    eta_max = eta_min + np.log10(params['visc_contrast'])   # max. log(eta)


    xrand = np.random.random(200)

    # Choose viscosity distribution
    if params['visc_dist'] == 'log':
        arrStressRand = 10. ** ( xrand * (eta_max-eta_min) + eta_min)
    elif params['visc_dist'] == 'bimodal':
        arrStressRand = np.tile([0.,1.],len(xrand)//2) 
        arrStressRand *= 10.**eta_max - 10.**eta_min
        arrStressRand += 10.**eta_min
    elif params['visc_dist'] == 'powerlaw':
        arrStressRand = (eta_min**D + (eta_max**D - eta_min**D)*xrand)**(1./D)
    else:
        raise NameError('No viscosity distribution chosen')


    # map patches to fault x co-ordinates
    arrWCumul = np.cumsum(arrW)

    arrX = np.linspace(0,params['L'],params['res'])
    arrPert = np.zeros(len(arrX))
    arrMaterial = np.zeros(arrX.size)

    pCount = 0
    for i,xi in enumerate(arrX):
        if xi > arrWCumul[pCount+1]:
            pCount += 1
        material = pCount % 2
        arrMaterial[i] =  arrStressRand[pCount+1]
            
    arrW = arrW[:pCount]


    # Set viscosity
    p.mesh_dict["INV_VISC"] = 1. / arrMaterial * params['W']



    # Write input to qdyn.in
    p.write_input()
    # On Windows 10, set this flag to True 
    # (see http://ydluo.github.io/qdyn/getting_started.html#additional-notes-for-windows-10-users)
    p.W10_bash = False


    # Plot the viscosity distribution

    plt.plot(arrX, params['W']/(p.mesh_dict["INV_VISC"]))
    plt.axhline(60e6/1e-9*params['W'], ls=":", c="k")
    plt.xlabel("Position (m)")
    plt.ylabel("Viscosity (Pa s)")
    plt.tight_layout()
    plt.yscale('log')
    plt.ylim(10.**eta_min,10.**eta_max)
    plt.savefig('eta_dist_%.2f.pdf' %params['W'])
    np.savetxt('eta_dist_%.2f.txt' %params['W'],np.vstack([arrX/1e3,params['W']/p.mesh_dict["INV_VISC"]]).T)


    # Run model
    p.run()

