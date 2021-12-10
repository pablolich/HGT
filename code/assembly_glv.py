#!/usr/bin/env python3

__appname__ = '[assembly.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from functions import *

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Number of metabolites 
    m = 200
    #Number of species
    s = 5
    #Preference matrix
    P = np.zeros(shape=(s, m))
    #Sample a matrix of metabolic preferences with norm 1
    for i in range(s):
        P[i,:] = sample_preferences(m)
    #Create interaction matrix (for now it is just one species)
    A = interaction_matrix(P)
    #Assign growth rates
    r = np.ones(s)
    #Find equilibrium
    x_star = find_equilibria(r, A)
    ind_present = np.where(x_star > 0)[0]
    #Initialize number of species vector
    n_sp_vec = list()
    #Start iteration
    while len(ind_present) <= round(0.2*m):
        invasive = False
        while not invasive:
            #Select randomly a species from the community
            sp_ind = np.random.choice(ind_present)
            #Create a new mutated species by inducing a perturbation of 3% in 
            #one of the resource consumption rates of the selected species
            v_new = mutate(P[sp_ind,:], 0.03)
            #Add v_new, the new consumer vector to matrix P
            P_new = np.vstack([P, v_new])
            #Create new matrix of interactions
            A_new = interaction_matrix(P_new)
            #Check weth or not A_new is singular 
            invasive = invasion_criteria(1, A_new, x_star)
        #The mutant can invade and matrix a is non-singular; add it to the 
        #community
        P = P_new
        A = A_new
        #Update number of species
        s += 1
        #Find new equilibrium
        x_star_new = find_equilibria(np.ones(s), A)
        #Check feasibility
        feasible = np.all(x_star_new >= 0) 
        if feasible:
            #Assuming that feasibility imply local stability
            x_star = x_star_new
        else:
            #Keep integrating and pruning until a feasible (and thus,  
            #stable) equilibrium is reached
            #Introduce invader (that we know can invade) at low initial 
            #abundance
            x0 = np.hstack([x_star, 1e-5])
            t_span = [0, 2000]
            #Integrate until solutions are constant 
            sol = integrate_n(lotka_volterra, 10, t_span, x0, A, s, 
                              tol=1e-6)
            #Check if the invader is still went extinct
            if sol.y[-1, -1] < 0:
                #If invader goes extinct, terminate this loop
                break
            x_star = sol.y[:, -1]
            s = len(x_star)
            #Get indices of species that are present
            ind_present = np.where(x_star > 0)[0]
            #Feasibility and (local) stability are ensured by 
            #integration. 
            #Record species number
            n_sp_vec.append(len(ind_present))
            print('Number of coexisting species: ', len(ind_present),
                  end = '\r')
        #We pruned the community succesfully to a feasible and 
        #stable state
    import ipdb; ipdb.set_trace(context = 20)
    plt.plot(np.arange(len(n_sp_vec)), n_sp_vec)
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
