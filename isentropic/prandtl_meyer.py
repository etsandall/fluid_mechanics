import numpy as np

def scalar_or_array(func):
  '''Decorator so that func will return a scalar if second argument is a scalar else a numpy array'''
  def wrapper(a, b, *args, **kwargs):
    if np.isscalar(b):
      return func(a, b, *args, **kwargs)
    else:
      return func(a, np.asarray(b), *args, **kwargs)
  return wrapper

############################
# Prandtl-Meyer Relations  #
############################

@scalar_or_array
def M2nu(G, M):
    '''Prandtl-Meyer relation: Mach number -> PM function'''

    # Outputs: nu [Prandtl-Meyer function] in degrees
    # inputs:  G  [specific heat ratio]
    #          M  [Mach number]
    
    return (((np.sqrt((G+1.0)/(G-1.0)))*np.arctan(np.sqrt((G-1.0)*(M**2.0 - 1.0)/(G+1.0))) - 
            np.arctan(np.sqrt(M**2.0 -1.0)))*180.0/np.pi)

@scalar_or_array
def M2dnudM(G, M):
    '''Prandtl-Meyer relation: Mach number -> derivative of PM function'''

    # Outputs: dnu/dM [derivative of Prandtl-Meyer function wrt Mach number]
    # inputs:  G [specific heat ratio]
    #          M [Mach number]
    return (((G-1.0)*np.sqrt((G+1.0)/(G-1.0))*M)/((G+1.0)*np.sqrt((G-1.0)*(M**2.0-1.0)/(G+1.0))*
        ((G-1.0)*(M**2.0-1.0)/(G+1.0)+1.0)) - 1.0/(M*np.sqrt(M**2.0-1.0)))*180.0/np.pi

@scalar_or_array
def nu2M(G, nu):
    '''Prandtl-Meyer relation: PM function -> Mach number'''

    # Outputs: M  [Mach number]
    # inputs:  G  [specific heat ratio]
    #          nu [Prandtl-Meyer function]

    # Newton-Rhapson Method
    if np.size(nu) == 1:
      M = [0.0, 1.5]
      while abs(M[0] - M[1])/M[1] > 1.0e-15:
          M[0] = M[1]
          M[1] = M[0] - (M2nu(G,M[0])-nu)/M2dnudM(G, M[0])
      M2 = M[1]
    else:
      M2 = []
      for n in nu:
        M = [0.0, 1.5]
        while abs(M[0] - M[1])/M[1] > 1.0e-15:
            M[0] = M[1]
            M[1] = M[0] - (M2nu(G,M[0])-n)/M2dnudM(G, M[0])
        M2 += [M[1]]
    return M2

@scalar_or_array
def nu2mu(G, nu):
    '''Prandtl-Meyer relation: PM function -> Mach angle'''

    # Outputs: mu [Mach angle] in degrees
    # inputs:  G  [specific heat ratio]
    #          nu [Prandtl-Meyer function]
    M = nu2M(G,nu)
    return np.arcsin(1.0/M)*180.0/np.pi
