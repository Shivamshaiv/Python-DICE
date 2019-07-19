import numpy as np
from util import compute_external_forcing
from util import update_carbon_masses
from util import productivity
from util import abatement_phi
from util import emmision_control_from_carbontax
import math
import matplotlib.pyplot as plt


class Economy:
    def __init__(self,start_year=2005,social_time_prefrence_rate_rho=0.015,elasticity_marginal_utility_alpha=1.5,
    coef_on_damage_exponent_pi2 = 0.0028,damage_exponent_epsilon = 2,exponent_emission_reduction_theta2 = 2.8,
    production_gamma=0.300,saving_rate_s=0.22,preindustrail_carbon_Mpi=592):
        productivity = lambda A,K,L,gamma : A*np.power(K,gamma)*np.power(L,1-gamma)                     # Productivity Function
        abatement_phi = lambda time : 1                                                                # Temp abatement.
        self.start_year = start_year                                                 # The starting year of the simulation where all values are caliberated.
        self.time = 0
        self.L = [6411]                                                                                   # World Population in millions
        self.L_g = [0.5]
        self.R = [1/(1+social_time_prefrence_rate_rho)]                                                     # Social Time Discount factor
        self.A = [0.0303220]
        self.A_g = [0.16]                   #  Initial growth rate of TFP per decade
        self.K = [137]                      #  2005 capital value
        self.sigma_g = [0.158]              #  Initial rate of decline of carbon intensity per period
        self.sigma = [0.14452]              #  2005 Effective Carbon Intensity
        self.T_at = [0.83]                    #  Temperature change from 1900 until 2000
        self.omega = [1 - (1/(1+(coef_on_damage_exponent_pi2*np.power(self.T_at[0],damage_exponent_epsilon))))]
        self.tipping_omega = [1-np.power((1+np.power(self.T_at[-1]/20.46,2)+np.power(self.T_at[-1]/6.081,6.754)),-1)]
        self.BC = [1.26]                    # Cost of backstop technology
        self.mu = [1]
        self.Eind = [(self.sigma[0])*(1-self.mu[0])*productivity(self.A[0],self.K[0],self.L[0],production_gamma)]
        self.lmbda = [(np.power(abatement_phi(self.time),1-exponent_emission_reduction_theta2)*self.BC[0]*
        np.power(self.mu[0],exponent_emission_reduction_theta2)*self.sigma[0])/exponent_emission_reduction_theta2]
        self.Y = [productivity(self.A[0],self.K[0],self.L[0],production_gamma)]
        #self.Q = [((1 - self.omega)*(1 - self.lmbda))*self.Y]
        self.Q = [55.34]
        self.I = [saving_rate_s * self.Q[0]]
        self.C = [self.Q[0] - self.I[0]]                                                # Total consumption, trillions of 2005 US dollars
        self.c = [self.C[0]/self.L[0]]                                                              # per capita consumption, thosands of 2005 US dollars
        self.U = [(np.power(self.c[0],1-elasticity_marginal_utility_alpha))/(1-elasticity_marginal_utility_alpha) + 1]
        self.test = []

        self.E_ind = [self.sigma[0]*self.Y[0]]
        #self.E_ind = [84.1910-1.1]
        self.E_land = [1.1]                                                          # Carbon emissions from land use (ie, deforestation), GtC per period
        self.E = [self.E_ind[-1]+self.E_land[-1]]
        self.E_cum = [self.E[-1]]
        self.Tax_tau = self.BC[0]*pow(self.mu[0],exponent_emission_reduction_theta2-1)  # Carbon Tax

        self.M_at = [787]                                                            # Mass of carbon in the atmosphere in 2005
        self.M_up = [1600]                                                           # Mass of carbon in the upper ocean in 2005
        self.M_lo = [10100]                                                          # Mass of carbon in the lower ocean in 2005

        self.T_lo = [0.0068]                                                        #  Temperature change in the lower ocean from 1900 until 2000
        self.F_ex = [0.83]
        self.F_ex0 = self.F_ex[0]                                                         #  Non-CO2 forcings in 2005
        self.F_ex10 = 0.30                                                        #  Estimate of non-CO2 forcings in 2100

        self.W = [self.L[-1]*self.U[-1]*self.R[-1]]

        ##### Connstants ############
        self.social_time_prefrence_rate_rho = social_time_prefrence_rate_rho
        self.elasticity_marginal_utility_alpha= elasticity_marginal_utility_alpha
        self.coef_on_damage_exponent_pi2 = coef_on_damage_exponent_pi2
        self.damage_exponent_epsilon = damage_exponent_epsilon
        self.exponent_emission_reduction_theta2 = exponent_emission_reduction_theta2
        self.production_gamma = production_gamma
        self.saving_rate_s = saving_rate_s
        self.preindustrail_carbon_Mpi = preindustrail_carbon_Mpi
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

        self.F = [self.F2CO2*(np.log2(self.M_at[0]/self.preindustrail_carbon_Mpi)) + self.F_ex0]       # Forcing due to CO2


    def loop(self,t,tipping_damage = False,temp_model = "default",carbon_tax = (0,0,0)):
        for time in range(t):
            self.time = time
            self.L_g.append((np.exp(self.L_g[0]*(time+1))-1)/(np.exp(self.L_g[0]*(time+1))-0))
            self.L.append(self.L[-1]*(np.power(self.L_Tmax/self.L[-1],self.L_g[-1])))
            self.R.append(1/np.power(1+self.social_time_prefrence_rate_rho,self.time))
            self.T_lo.append(self.T_lo[-1]+(self.ξ3*(self.T_at[-1]-self.T_lo[-1])))
            const_lamb = self.F2CO2/self.temp_increase_doubling_co2
            self.F_ex.append(compute_external_forcing(self,time))

            # Update A(t)
            self.A_g.append(self.A_g[0]*np.exp(-self.tech_change_decline_deltaa*1*time*np.exp(-self.declinerate_growth_productivity_deltab*time)))
            self.A.append((self.A[-1])/(1-self.A_g[-2]))

            # Update K(t)
            self.K.append(self.I[-1]+(np.power((1-self.depreciation_technological_change_δK),1)*(self.K[-1]*1)))

            #Update productivity
            self.Y.append(productivity(self.A[-1],self.K[-1],self.L[-1],self.production_gamma))

            self.sigma_g.append(self.sigma_g[-1]*(1-self.decline_rate_decarbonisation_σd1)**1)
            self.sigma.append(self.sigma[-1]*(1-self.sigma_g[-2]))

            # Emisssions
            self.E_land.append(self.E_land[0]*np.power(0.8,time))
            self.E_ind.append(self.sigma[-1]*(1-self.mu[-1])*self.Y[-1])
            self.E.append(self.E_land[-1]+self.E_ind[-1])
            self.E_cum.append(np.sum(self.E))

            update_carbon_masses(self)
            F_eta = 3.2    #Forcing of CO2
            self.F.append(F_eta*(np.log(self.M_at[-1]/self.preindustrail_carbon_Mpi)) + self.F_ex[-1])

            if temp_model == "default":
                self.T_at.append(self.T_at[-1]+self.ξ1*(self.F[-1]-const_lamb*self.T_at[-1]-self.ξ2*(self.T_at[-1]-self.T_lo[-2])))
            elif temp_model == "linear":
                #denom = (2*2*self.preindustrail_carbon_Mpi) - (0.94796*self.M_at[0]) - (0.00075*self.M_up[0])
                denom = (2*2*self.preindustrail_carbon_Mpi) - (-0.2*self.M_at[0]) - (0.001*self.M_up[0])
                self.T_at.append(self.T_at[0]+(self.time)*(self.E_cum[-1])*(self.temp_increase_doubling_co2/(denom)))
            self.omega.append(1 - (1/(1+(self.coef_on_damage_exponent_pi2*np.power(self.T_at[-1],self.damage_exponent_epsilon)))))
            self.tipping_omega.append(1-np.power((1+np.power(self.T_at[-1]/20.46,2)+np.power(self.T_at[-1]/6.081,6.754)),-1))

            ### Abetement Function #To do
            self.BC.append(self.BC[0]*np.power(1-self.costdecinline_backstop_tech_percent_BCg,self.time))
            self.mu.append(emmision_control_from_carbontax(self,carbon_tax))
            self.lmbda.append((np.power(abatement_phi(self.time),1-self.exponent_emission_reduction_theta2)*self.BC[-1]*np.power(self.mu[-1],self.exponent_emission_reduction_theta2)*self.sigma[-1])/self.exponent_emission_reduction_theta2)

            # Update Q production Function
            if tipping_damage:
                self.Q.append((1-self.tipping_omega[-1])*(1-self.lmbda[-1])*self.Y[-1])
            else:
                self.Q.append((1-self.omega[-1])*(1-self.lmbda[-1])*self.Y[-1])

            # Update the Investment.
            self.I.append(self.saving_rate_s*self.Q[-1])

            #Update the capital
            self.C.append(self.Q[-1]-self.I[-1])
            self.c.append(self.C[-1]/self.L[-1])

            #Update the utility
            self.U.append((np.power(self.c[-1],1-self.elasticity_marginal_utility_alpha)/(1-self.elasticity_marginal_utility_alpha))+1)

            self.test.append(emmision_control_from_carbontax(self,carbon_tax))


model = Economy()
#print(model.K[-1])
#print(model.__dict__)
model.loop(100)
#print(model.U)
#print(model.__dict__)
#plt.plot(model.A_g)

#plt.show()
