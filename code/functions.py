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
    #Select randomly an index to perturb
    ind = np.random.randint(len(v))
    #Sample perturbation 
    delta = np.random.uniform(1 - p, 1 + p)
    #Apply perturbation
    v_mut = np.copy(v) 
    v_mut[ind] = delta*v[ind]
    norm = np.linalg.norm(v_mut)
    return v_mut/norm

def interaction_matrix(P):
    '''
    Calculate the covariance matrix of P, rescaled such that coefficients
    are in the ingerval [-1, 0]
    '''
    return -(P@P.T + 1)/2

def invasion_criteria(r_i, A, x):
    #Select interaction coefficients that affect the new species
    try:
        a_i_vec = A[-1, 0:-1]
    except:
        #If A is a scalar because we are in the first iteration
        a_i_vec = A
    inv = r_i + np.dot(a_i_vec, x)
    try:
        if len(inv) > 1:
            import ipdb; ipdb.set_trace(context = 20)
    except:
        inv = r_i + np.dot(a_i_vec, x)
    return inv 

def find_equilibria(r, A):
    if isinstance(A, float) or isinstance(A, int):
        x_star = -1/A*r
    else:
        try:
            x_star = (-np.linalg.inv(A)@r.reshape(len(r), 1)).T[0]
        except: 
            x_star = None
    return x_star

def find_eigenvals(A):
    try:
        eigs = np.linalg.eigvals(A)
    except:
        eigs = A
    return eigs

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
    #Get last 10 timepoints
    last_10 = diff_sol[:, -1:-10:-1]
    const = np.all(last_10 < tol)
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
