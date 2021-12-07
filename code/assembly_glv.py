#!/usr/bin/env python3

__appname__ = '[assembly.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pylab as plt
from functions import *

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Number of metabolites 
    m = 50
    #Number of species
    s = 5
    #Preference matrix
    P = np.zeros(shape=(s, m))
    #Sample a matrix of metabolic preferences with norm 1
    for i in range(s):
        P[i,:] = sample_preferences(m)
    #Create interaction matrix (for now it is just one species)
    import ipdb; ipdb.set_trace(context = 20)
    A = interaction_matrix(P)
    #Assign growth rates
    r = np.ones(s)
    #Find equilibrium
    x_star = find_equilibria(r, A)
    #Initialize number of species vector
    n_sp_vec = list()
    #Start iteration
    while s <= m:
        print('Number of species: ', s)
        #Select randomly a species from the community
        sp_ind = np.random.randint(s)
        singular = True
        non_invasive = True
        while singular or non_invasive:
            #Create a new mutated species by inducing a perturbation of 3% in 
            #one of the resource consumption rates of the selected species
            v_new = mutate(P[sp_ind,:], 0.03)
            #Add v_new, the new consumer vector to matrix P
            P_new = np.vstack([P, v_new])
            #Create new matrix of interactions
            A_new = interaction_matrix(P_new)
            #Check weth or not A_new is singular 
            singular = check_singularity(A_new)
            non_invasive = bool(invasion_criteria(r, A_new, x_star))
        #The mutant can invade and matrix a is non-singular; add it to the 
        #community
        P = P_new
        A = A_new
        #Update number of species
        s += 1
        #Also on the vector
        n_sp_vec.append(s)
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
            x0 = np.hstack([x_star, 1e-3])
            t_span = [0, 2000]
            #Integrate until solutions are constant 
            sol = interate_n(lotka_volterra, 10, t_span, x0, tol=1e-9)
            #Check if the invader is still went extinct
            if sol[-1, -1] < 0
                #If invader goes extinct, terminate this loop
                break
            x_star = sol[:, -1]
            s = len(x_star)
            #Feasibility and (local) stability are ensured by 
            #integration. 
            #Record species number
            n_sp_vec.append(len(x_star))
            print(n_sp_vec)
        #We pruned the community succesfully to a feasible and 
        #stable state
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
