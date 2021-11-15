#!/usr/bin/env python3

__appname__ = '[assembly.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Number of metabolites 
    m = 10
    #Sample a vector of metabolic preferences with norm 1
    v_first = sample_preferences(m)
    #Create interaction matrix (for now it is just one species)
    A = -1
    #Assign growth rates
    r = 1
    #Find equilibrium
    x_star = find_equilibria(r, A)
    #Set number of species
    n_sp = number_sp(x_star)
    #Initialize interaction matrix
    P = v_first
    #Set abundance threshold
    ab_thresh = 1e-6
    #Initialize number of species vector
    n_sp_vec = list()
    #Start iteration
    while n_sp <= m:
        print(n_sp)
        #Select randomly a species from the community
        sp_ind = np.random.randint(n_sp)
        #Create a new mutated species by inducing a perturbation of 3% in one 
        #of the resource consumption rates of the selected species
        v_new = mutate(v_first, 0.2)
        #Add v_new, the new consumer vector to matrix P
        P_new = np.vstack([P, v_new])
        #Create new matrix of interactions
        A_new = interaction_matrix(P_new)
        #Check weth or not A_new is singular 
        singular = check_singularity(A_new)
        while singular:
            v_new = mutate(v_first, 0.03)
            P_new = np.vstack([P, v_new])
            A_new = interaction_matrix(P_new)
            singular = check_singularity(A_new)
        #Check if the new species can invade
        inv = invasion_criteria(r, A_new, x_star)
        print('inv criteria: ', inv)
        if inv > 0:
            #In the case that it can invade, add it to the community 
            P = P_new
            A = interaction_matrix(P)
            #Update number of species
            n_sp += 1
            n_sp_vec.append(n_sp)
            #Find new equilibrium
            x_star_new = find_equilibria(np.ones(n_sp), A)
            #Calculate eigenvalues of A
            eigen_vals_new = find_eigenvals(A)
            #Check feasibility and stability
            try:
                feasible = np.all(x_star_new > 0) 
            except:
                plt.plot(n_sp_vec)
                plt.show()
                raise Exception('A is singular, end loop')
            stable = np.all(eigen_vals_new.real < 0)
            if feasible and stable:
                #Set x_star to the equilibrium reached in the presecne of
                #the invader
                x_star = x_star_new
            else:
                #Keep integrating and pruning until a feasible and stable
                #equilibrium is reached
                #Introduce invader at low initial abundance
                x0 = np.hstack([x_star, 1e-3])
                t_span = [0, 20000]
                while not (feasible and stable):
                    #Integrate until solutions are constant 
                    constant_sol = False
                    #Put a maximum integration times
                    runs = 0
                    while not constant_sol and runs < 10:
                        #Run dynamics until putative equilibrium is reached
                        sol = solve_ivp(lotka_volterra, t_span, x0, 
                                        method = 'BDF', 
                                        args = (A,  np.ones(n_sp)))
                        #Check if solution is constant
                        constant_sol = check_constant(sol.y, tol = 1e-3)
                        #Get rid of extinct species (those below abundance 
                        #threshold)
                        print('Endpoint abundances: ', sol.y[:, -1])
                        x_star, P, A = remove_extinct(sol.y[:,-1], P, A, 
                                                      ab_thresh)
                        print('Dimensions of A: ', A.shape)
                        #Set initial conditions for next integration with the 
                        #endpoints of the previous one
                        x0 = x_star
                        n_sp = len(x_star)
                        runs += 1
                        print('Integration number: ', runs)
                    #Check feasibility and stability
                    feasible = np.all(x_star > 0) 
                    eigen_vals_new = np.linalg.eigvals(A)
                    stable = np.all(eigen_vals_new.real < 0)
                    n_sp_new = len(x_star)
                #We pruned the community succesfully to a feasible and stable
                #state
                #Record species number
                n_sp_vec.append(len(x_star))
                print(n_sp_vec)
        else: 
            #Invasibility criterion is not satisfied, so come up with another 
            #mutant that can invade
            continue
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
