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

def invasion_criteria(r_i, A, x):
    #Select interaction coefficients that affect the new species
    a_i_vec = A[-1, 0:-1]
    return r_i + np.dot(a_i_vec, x)

def lotka_volterra(t, x, A, r):
    return np.diag(x)@(r + A@x)

def enzyme_budget():
    return None

def remove_extinct(abundances, preferences, A, tol):
    '''
    Remove extinct species from abundance vector, preferences and interaction
    matrices
    '''
    #Identify indices of extinct species
    ind_ext = np.where(abundances < tol)
    #Remove from abundance vector
    abundances = np.delete(abundances, ind_ext)
    #Remove rows corresponding to extinct species from preference matrix
    preferences = np.delete(preferences, ind_ext, axis = 0)
    #Remove rows and columns from interaction matrix
    A = np.delete(A, ind_ext, axis = 0)
    B = np.delet(A, ind_ext

    return list(abundances, preferences, A)

def main(argv):
    '''Main function'''
    #Number of metabolites 
    m = 2
    #Sample first species
    v_first = first_species(m)
    #Create interaction matrix
    A = -1
    #Assign growth rates
    r = 1
    #Find equilibrium
    x_star = -A*r
    #Set number of species
    n_sp = 1
    #Initialize interaction matrix
    P = v_first
    #Set abundance threshold
    ab_thresh = 1e-6
    #Start iteration
    while n_sp <= m:
        #Select randomly a species from the community
        sp_ind = np.random.randint(n_sp)
        #Create a new mutated species by inducing a perturbation of 3% in one 
        #of the resource consumption rates of the selected species
        v_new = mutate(v_first, 0.03)
        #Add v_new, the new consumer vector to matrix P
        P_new = np.hstack([P, v_new])
        #Create new matrix of interactions
        A_new = interaction_matrix(P_new)
        #Check if the new species can invade
        inv = invasion_criteria(r, A_new, x_star)
        if inv:
            #In the case that it can invade, add it to the community 
            P = P_new
            A = interaction_matrix(P)
            #Find new equilibrium
            x_star_new = -np.linalg.inv(A)@r
            #Check for stability
            eigen_vals_new = np.linalg.eig_vals(A)
            #When the community is feasible and stable, keep adding invaders
            if ( np.all(x_star_new > 0) and np.all(eigen_vals_new.real < 0) ):
                #Set x_star to the equilibrium reached in the presecne of
                #the invader
                x_star = x_star_new
            else:
                #When either of those conditions is not met, numerically 
                #integrate the dynamics to check what is the final community 
                #state after invasion
                #Introduce invader at low initial abundance
                x0 = np.vstack([x_star, 1e-5])
                #Run dynamics until equilibrium is reached
                sol = solve_ivp(lotka_voltera, t_span, x0, args = arguments)
                #Get rid of extinct species (those below abundance threshold)

        else: 
            #Invasibility criterion is not satisfied, so come up with another 
            #mutant that can invade
            continue


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
