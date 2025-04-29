import numpy as np
from scipy.optimize import fsolve

def M2Aratio(k, M):
  return ((k+1.0)/2.0)**(-(k+1.0)/(2.0*(k-1.0)))*(1.0 + (k-1.0)/2.0*M**2.0)**((k+1.0)/(2.0*(k-1.0)))/M

def Aratio2M(k, Aratio, supersonic=True):
  def f(M):
    return Aratio - M2Aratio(k, M)
  if supersonic:
    M = fsolve(f, 5.0, xtol=1e-12)[0]
    assert M > 1.0, "Mach number is subsonic."
  else:
    M = fsolve(f, 0.1, xtol=1e-12)[0]
    assert M < 1.0, "Mach number is supersonic."
  return M

def M2Pratio(k, M):
  return (1.0 + (k-1.0)/2.0*M**2.0)**(-k/(k-1.0))

def Pratio2M(k, Pratio):
  return np.sqrt((Pratio**((k+1.0)/(-k)) - 1.0)*2.0/(k-1.0))

if __name__ == '__main__':
  gamma = 1.4
  M = 2.0
  Aratio = M2Aratio(gamma, M)
  print(Aratio)
  M2 = Aratio2M(gamma, Aratio)
  print(M2)
  Msub = Aratio2M(gamma, Aratio, False)
  print(Msub)

  Pratio = M2Pratio(gamma, M)
  print(Pratio)
  Mp = Pratio2M(gamma, Pratio)
  print(Mp)
  
