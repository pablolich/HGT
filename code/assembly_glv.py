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
    s = 2
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
    while len(ind_present) <= round(m):
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
        #feasible = np.all(x_star_new >= 0) 
        #if feasible:
        #    #Assuming that feasibility imply local stability
        #    x_star = x_star_new
        #    #Create matrix of HGT, H
        #    H = compute_h(x_star)
        #    #Create the the matrix of mutated preferences based on HGT
        #    P_hgt = hgt(H, x_star, P)

        #else:
        #Keep integrating and pruning until a feasible (and thus,  
        #stable) equilibrium is reached
        #Introduce invader (that we know can invade) at low initial 
        #abundance
        x0 = np.hstack([x_star, 1e-5])
        t_span = [0, 5000]
        #Integrate until solutions are constant 
        sol = integrate_n(lotka_volterra, 1, t_span, x0, A, s, 
                          tol=1e-6)
        #Check if the invader is still went extinct
        if sol.y[-1, -1] < 0:
            #If invader goes extinct, terminate this loop
            break
        x_star = sol.y[:, -1]
        s = len(x_star)
        #Get indices of species that are present
        ind_present = np.where(x_star > 0)[0]
        #Record species number
        n_sp_vec.append(len(ind_present))
        #Now start HGT
        #Create matrix of HGT, H
        H = compute_h(x_star)
        #Create the the matrix of mutated preferences based on HGT
        P_hgt = hgt(H, x_star, P)
        A_hgt = interaction_matrix(P_hgt)
        #Create megamatrix with original and mutants
        P_tot = np.vstack([P, P_hgt])
        #Compute new interaction matrix
        A_tot = interaction_matrix(P_tot)
        #Only present species give offspirng
        offspring = np.zeros(s)
        x_present_ab = x_star[ind_present]
        offspring[ind_present] = x_present_ab
        #Set initial conditions by introducing mutants at low abundance, but 
        #proportional to the abundance of their parents
        x0_hgt = np.hstack([x_star, offspring])
        #Integrate the system
        sol_hgt = integrate_n(lotka_volterra, 1, t_span, x_star, 
                              A_hgt, s, tol = 1e-6) 
        print("Plotting form assembly")
        for i in range(s):
            if i >= s:
                plt.plot(sol_hgt.t, sol_hgt.y[i,:], linestyle = 'dashed')
            else:
                plt.plot(sol_hgt.t, sol_hgt.y[i,:])

        plt.show()

        print('Number of coexisting species: ', len(ind_present),
              end = '\r')
        #We pruned the community succesfully to a feasible and 
        #stable state
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
