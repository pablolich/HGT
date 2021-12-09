import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
%matplitlib inline


### Functions for generalized lotka volterra model ### 

def sample_preferences(m):
    '''
    Sample preferences from a normal distribution and normalize to 1
    '''
    vec = np.random.normal(0, 1, m)
    norm = np.linalg.norm(vec)
    return  vec/norm

def mutate(v, p):
    '''
    Generate a single mutation in one of the elements of the preference vector
    by perturbing an element by a percentage
    '''
    #Sample perturbation 
    delta = np.random.uniform(1 - p, 1 + p, size = len(v))
    #Apply perturbation
    v_mut = v*delta
    #Normalize again
    v_mut_norm = v_mut/np.linalg.norm(v_mut)
    return v_mut_norm

def interaction_matrix(P):
    '''
    Calculate the covariance matrix of P, rescaled such that coefficients
    are in the ingerval [-1, 0]
    '''
    return -(P@P.T + 1)/2

def invasion_criteria(r_i, A, x):
    #Select last row of A: the effect every species (but invader) have on the
    #invader
    effects_on_invader = A[-1, 0:-1] #Note we are leaving out A_ii
    inv = r_i + np.dot(effects_on_invader, x)
    return inv > 0

def find_equilibria(r, A):
    return (-np.linalg.inv(A)@r.reshape(len(r), 1)).T[0]

def number_sp(x):
    try:
        n = len(x)
    except:
        n = 1
    return n

def lotka_volterra(t, x, A, r):
    return (np.diag(x)@(r.reshape(len(r), 1) + A@x.reshape(len(r), 1))).T[0]

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
    A = np.delete(A, ind_ext, axis = 1)
    return abundances, preferences, A

def check_constant(sol_mat, tol):
    '''
    Check if all the solutions have reached steady state (constant)
    '''
    #Get differences between solutions
    diff_sol = sol_mat[:, 1:] - sol_mat[:, 0:-1]
    #Get last 3 timepoints
    last_3 = diff_sol[:, -1:-3:-1]
    #Note that we only impose three because there are no oscillations here. 
    const = np.all(abs(last_3) < tol)
    return const

def check_singularity(A):
    '''
    Check if A is a singular matrix
    '''
    #Calculate determinant
    det_A = abs(np.linalg.det(A))
    if det_A < 1e-16:
        return True
    else:
        return False


### Functions for consumer resource model ### 

def consumer_resource(t, x, C, r, z, K, b):
    N, R = x
    return [np.diag(r)@(C@R - z), K - b*R - np.diag(R)@(C.T@N)]

def cost_model(C):
    '''
    Calculate cost of each species in the community
    '''
    return np.sum(C, axis = 0)

def assembly_cr(C, r, z, K, b):
    '''
    Assembly a community by subsequently adding successful invaders 
    '''
    assembly_cr(C, np.ones(m), z, np.ones(m), np.ones(m))
    
### General functions ### 

def integrate_n(fun, tot_runs, t_span, x0, A, n_sp, tol):
    '''
    Integrate until solutions are constant
    '''
    #Integrate until solutions are constant 
    constant_sol = False
    #Put a maximum integration times
    run_i = 0
    while not constant_sol and run_i < tot_runs:
        #Run dynamics until putative equilibrium is reached
        sol = solve_ivp(fun, t_span, x0, 
                        method = 'BDF', 
                        args = (A,  np.ones(n_sp)))
        #Prepare initial conditions for next integration
        x0 = sol.y[:, -1]
        #Get indices of species with abundance below tolerance
        ext_ind = np.where(x0 < tol)[0]
        #Set to 0 species that are below threshold
        x0[ext_ind] = 0
        #Check if solution is constant
        constant_sol = check_constant(sol.y, 1e-6)
        #Add 1 to the numer of iterations
        run_i += 1
        print('Integration number: ', run_i)
        import ipdb; ipdb.set_trace(context = 20)
        for i in range(n_sp):
            plt.plot(sol.t, sol.y[i, :])
        plt.show()
    return sol
