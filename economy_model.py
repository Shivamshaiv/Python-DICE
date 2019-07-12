import numpy as np
from util import compute_external_forcing
from util import update_carbon_masses
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
        self.T_at = [0.83]                    #  Temperature change from 1900 until 2000                                                           # Temperature change from 1900 until 2000
        self.omega = [1 - (1/(1+(coef_on_damage_exponent_pi2*np.power(self.T_at[0],damage_exponent_epsilon))))]
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
        self.F_ex = [0.83]
        self.F_ex0 = self.F_ex[0]                                                         #  Non-CO2 forcings in 2005
        self.F_ex10 = 0.30                                                        #  Estimate of non-CO2 forcings in 2100

        ##### Connstants ############
        self.social_time_prefrence_rate_rho=0.015
        self.elasticity_marginal_utility_alpha=1.5
        self.coef_on_damage_exponent_pi2 = 0.0028
        self.damage_exponent_epsilon = 2
        self.exponent_emission_reduction_theta2 = 2.8
        self.production_gamma = 0.300
        self.saving_rate_s = 0.22
        self.preindustrail_carbon_Mpi = 592
        self.L_Tmax = 8700                                 #  Asymptotic population in the last period
        self.tech_change_decline_deltaa = 0.9
        self.depreciation_technological_change_δK = 0.1
        self.decline_rate_decarbonisation_σd1 = 0.006
        self.temp_increase_doubling_co2 = 3.2              # Temperature increase (◦C) from a doubling of preindustrial CO2.
        self.costdecinline_backstop_tech_percent_BCg = 0.05
        self.fossil_remaining_Ccum = 6000
        self.pop_growth_pop_perdecade_Lg = 0.5
        self.declinerate_growth_productivity_deltab = 0.2

        # Climate Connstants
        self.F2CO2 = 3.8          #  Forcing from a doubling of preindustrial CO2
        self.ξ1 = 0.220           #  Inverse of thermal capacity of the atmosphere and the upper ocean.
        self.ξ2 = 0.310           #  Ratio of the thermal capacity of the deep oceans to the transfer rate from the shallow ocean to the deep ocean
        self.ξ3 = 0.050           #  Transfer rate of heat from the upper ocean to the deep ocean.

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

    def compute_external_forcing(self,time):
        if time > 10:
            return self.F_ex10
        else:
            return self.F_ex0 + 0.1*(self.F_ex0-self.F_ex10)*time

    def update_carbon_masses(self):
        _b11, _b12, _b13 = 0.810712, 0.189288, 0
        _b21, _b22, _b23 = 0.097213, 0.852787, .05
        _b31, _b32, _b33 = 0, 0.003119, 0.996881
        self.carbon_matrix = np.array([
            _b11, _b12, _b13,
            _b21, _b22, _b23,
            _b31, _b32, _b33,
        ]).reshape(3, 3)

        self.mass_matrix = np.array([[self.M_at[-1]],[self.M_up[-1]],[self.M_lo[-1]]])
        self.emmision_matrix = np.array([[self.E[-1]],[0],[0]])

        new_masses = np.matmul(self.carbon_matrix,self.mass_matrix) + self.emmision_matrix
        self.M_at.append(new_masses[0][0])
        self.M_up.append(new_masses[1][0])
        self.M_lo.append(new_masses[2][0])

    def loop(self,t=1):
        for time in range(t):
            self.time = time
            self.L.append(np.sqrt(self.L[-1]*self.L_Tmax))
            self.T_lo.append(self.T_lo[-1]+(self.ξ3*(self.T_at[-1]-self.T_lo[-1])))
            const_lamb = self.F2CO2/self.temp_increase_doubling_co2
            self.F_ex.append(compute_external_forcing(self,time))

            # Update A(t)
            self.A_g.append(self.A_g[0]*np.exp(-self.tech_change_decline_deltaa*time*np.exp(-self.declinerate_growth_productivity_deltab*time)))
            self.A.append((self.A[-1])/(1-self.A_g[-2]))

            # Update K(t)


            # Emisssions
            #self.E_land.append(self.E_land[0]*np.power(0.8,time))

            update_carbon_masses(self)
            F_eta = 1
            self.F = F_eta*(np.log2(self.M_at[-1]/self.preindustrail_carbon_Mpi)) + self.F_ex[-1]
            self.T_at.append(self.T_at[-1]+self.ξ1(self.F[-1]-const_lamb*self.T_at[-1]-self.ξ2(self.T_at[-1]-self.T_lo[-2])))
            self.omega.append(1 - (1/(1+(coef_on_damage_exponent_pi2*np.power(self.T_at[-1],damage_exponent_epsilon)))))


model = Economy()
print(model.__dict__)
model.loop(2)
print(model.__dict__)
