import numpy as np
from scipy.optimize import fsolve

class Shock:
  def __init__(self, gamma, M1, theta=None, beta=None, units='deg'):
    self.gamma = gamma
    self.mach1 = M1

    if beta:
      if units == 'deg':
        beta *= np.pi/180.0
      self.beta = beta
      self.theta = self.beta2theta(gamma, M1, beta)
    elif theta:
      if units == 'deg':
        theta *= np.pi/180.0
      def tbm(b,t,M):
        return np.tan(t) - 2.0/np.tan(b)*((M**2.0*np.sin(b)**2.0 - 1.0)/(M**2.0*(gamma + np.cos(2.0*b)) + 2.0))
      self.beta = fsolve(tbm, theta, (theta, M1), xtol=1e-12)[0]
      self.theta = self.beta2theta(gamma, M1, self.beta)
    else:
      raise SyntaxError('Must specify either flow angle (theta) or shock angle (beta).')

    self.mach2 = 1.0/np.sin(self.beta - self.theta) * np.sqrt((1.0 + (gamma-1.0)/2.0*M1**2.0*np.sin(self.beta)**2.0)/
                                                              (gamma*M1**2.0*np.sin(self.beta)**2.0 - (gamma-1.0)/2.0))
    self.pressure_ratio = 1.0 + 2.0*gamma/(gamma+1.0)*(M1**2.0*np.sin(self.beta)**2.0 - 1.0)
    self.density_ratio = (gamma+1.0)*M1**2.0*np.sin(self.beta)**2.0/((gamma-1.0)*M1**2.0*np.sin(self.beta)**2.0 + 2.0)
    self.temperature_ratio = self.pressure_ratio/self.density_ratio
    self.mach_normal1 = self.mach1*np.sin(self.beta)
    self.mach_normal2 = self.mach2*np.sin(self.beta - self.theta)

    self.stagnation_pressure_ratio = self.pressure_ratio*((1.0 + (gamma-1.0)/2.0*self.mach2**2.0)/(1.0 + (gamma-1.0)/2.0*self.mach1**2.0))**(gamma/(gamma-1.0))

  def beta2theta(self, gamma, M, b):
    return np.arctan(2.0/np.tan(b)*((M**2.0*np.sin(b)**2.0 - 1.0)/(M**2.0*(gamma + np.cos(2.0*b)) + 2.0)))
 
  def __repr__(self):
    return f'''
    Oblique Shock Relations
        gamma = {self.gamma}
        M1 = {self.mach1:.3f}
        theta = {self.theta:.3f} rad = {self.theta*180.0/np.pi:.3f} deg
        beta = {self.beta:.3f} rad = {self.beta*180.0/np.pi:.3f} deg

        M2 = {self.mach2:.3f}
        Mn1 = {self.mach_normal1:.3f}
        Mn2 = {self.mach_normal2:.3f}
        P2/P1 = {self.pressure_ratio:.3f}
        Po2/Po1 = {self.stagnation_pressure_ratio:.3f}
        rho2/rho1 = {self.density_ratio:.3f}
        T2/T1 = {self.temperature_ratio:.3f}
    '''

if __name__ == '__main__':
  print(Shock(1.4, 3.0, 17.2))
