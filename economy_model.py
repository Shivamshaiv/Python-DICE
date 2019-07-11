import numpy as np
import math

class Economy:
    def __init__(self,social_time_prefrence_rate_rho=0.015,elasticity_marginal_utility_alpha=1.5,
    coef_on_damage_exponent_pi2 = 0.0028,damage_exponent_epsilon = 2,exponent_emission_reduction_theta2 = 2.8,
    production_gamma=0.300,saving_rate_s=0.22,preindustrail_carbon_Mpi=592):
        productivity = lambda A,K,L,gamma : A*np.power(K,gamma)*np.power(L,1-gamma)                     # Productivity Function
        abatement_phi = lambda time : 1                                                                # Temp abatement.
        self.time = 1
        self.L = [6411]                                                                                   # World Population in millions
        self.R = [1/(1+social_time_prefrence_rate_rho)]                                                     # Social Time Discount factor
        self.A = [0.0303220]
        self.A_g = [0.16]                   #  Initial growth rate of TFP per decade
        self.K = [137]                      #  2005 capital value
        self.sigma_g = [0.158]              #  Initial rate of decline of carbon intensity per period
        self.sigma = [0.14452]              #  2005 Effective Carbon Intensity
        self.T_at = 0.83                    #  Temperature change from 1900 until 2000                                                           # Temperature change from 1900 until 2000
        self.omega = [1 - (1/(1+(coef_on_damage_exponent_pi2*np.power(self.T_at,damage_exponent_epsilon))))]
        self.BC = [1.26]                    # Cost of backstop technology
        self.mu = [0]
        self.Eind = [(self.sigma[0])*(1-self.mu[0])*productivity(self.A[0],self.K[0],self.L[0],production_gamma)]
        self.lmbda = [(np.power(abatement_phi(self.time),1-exponent_emission_reduction_theta2)*self.BC[0]*
        np.power(self.mu[0],exponent_emission_reduction_theta2)*self.sigma[0])/exponent_emission_reduction_theta2]
        self.Y = productivity(self.A[0],self.K[0],self.L[0],production_gamma)
        #self.Q = [((1 - self.omega)*(1 - self.lmbda))*self.Y]
        self.Q = [55.34]
        self.I = [saving_rate_s * self.Q[0]]
        self.C = [self.Q[0] - self.I[0]]                                                # Total consumption, trillions of 2005 US dollars
        self.c = [self.C[0]/self.L[0]]                                                              # per capita consumption, thosands of 2005 US dollars
        self.U = [(np.power(self.c[0],1-elasticity_marginal_utility_alpha))/(1-elasticity_marginal_utility_alpha) + 1]

        self.E_ind = [self.sigma[0]*(1-self.mu[0])*self.Y]
        self.E_land = [1.1]                                                          # Carbon emissions from land use (ie, deforestation), GtC per period
        self.E = [84.1910]
        self.Tax_tau = self.BC[0]*pow(self.mu[0],exponent_emission_reduction_theta2-1)  # Carbon Tax

        self.M_at = [787]                                                            # Mass of carbon in the atmosphere in 2005
        self.M_up = [1600]                                                           # Mass of carbon in the upper ocean in 2005
        self.M_lo = [10100]                                                          # Mass of carbon in the lower ocean in 2005

        self.T_lo = [0.0068]                                                        #  Temperature change in the lower ocean from 1900 until 2000
        self.F_ex0 = [0.83]                                                         #  Non-CO2 forcings in 2005
        self.F_ex10 = [0.30]                                                        #  Estimate of non-CO2 forcings in 2100

    def abatement_phi(time,phi0=1,phi5=1,phi10=1,phi15=1,phimax=1):
        if time < 5:
            return phi5 + (phi0-phi5)*np.exp(-0.25*time)
        elif time >= 5 and time < 10 :
            return phi10 + (phi5-phi10)*np.exp(-0.25*time)
        elif time >= 10 and time < 15 :
            return phi15 + (phi10-phi15)*np.exp(-0.25*time)
        else:
            return phi15 + (phimax-phi15)*np.exp(-0.25*time)

    def productivity(A,K,L,gamma):
        return A*np.pow(K,gamma)*np.pow(L,1-gamma)

model = Economy()
print(model.__dict__)
