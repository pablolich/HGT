#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys

## CONSTANTS ##


## FUNCTIONS ##

def first_species(m):
    '''
    Generate metabolic preferences of first species

    Parameters:
        m (int): Number of resources available in the environment
    Output:
        pref_vec (1xm): Array with preference vector of species 
    '''
    return(np.random.normal(m))

def mutate(v, p):
    '''
    Generate a single mutation in one of the elements of the preference vector
    by perturbing an element by a percentage
    '''
    #Select randomly an index to perturb
    ind = np.random.randint(len(v))
    #Sample perturbation 
    delta = np.random.uniform(1 - p, 1 + p)
    #Apply perturbation
    v[ind] = delta*v[ind]
    return v

def interaction_matrix(P):
    '''
    Calculate the covariance matrix of P, rescaled such that coefficients
    are in the ingerval [-1, 0]
    '''
    return -(P@P.T + 1)/2

def lotka_volterra(t, x, A, r):
    return np.diag(x)@(r + A@x)


def main(argv):
    '''Main function'''
    #Number of metabolites 
    m = 2
    #Sample first species
    v_first = first_species(m)
    #Create interaction matrix
    A = -1
    #Assign growthrates
    r = 1
    #Find equilibrium
    x_star = -A*r
    #Set number of species
    n_sp = 1
    P = v_first
    while n_sp <= m:
        #Select randomly a species from the community
        sp_ind = np.random.randint(n_sp)
        #Induce a mutation of 3% in selected species
        v_new = mutate(v_first, 0.03)
        #Add v_new to matrix P
        P = np.hstack([P, v_new])
        #Find equilibrium
        x_star = -np.linalg.inv(A)@r
        #Create mutation
        #Run dynamics until equilibrium is reached
        sol = solve_ivp(lotka_voltera, t_span, x0, args = arguments)


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
