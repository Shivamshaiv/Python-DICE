import numpy as np

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
    return A*np.power(K,gamma)*np.power(L*1,1-gamma)

def compute_external_forcing(self,time):
    if time > 10:
        #return self.F_ex0 + 0.1*(self.F_ex0-self.F_ex10)*time
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

def emmision_control_from_carbontax(self,carbon_tax):
    t5,t10,t15 = carbon_tax
    if self.time <= 10 and self.time > 5:
        output = np.power(t5/self.BC[-1],1/(self.exponent_emission_reduction_theta2-1))
    elif self.time <= 15 and self.time > 10:
        output = np.power(t10/self.BC[-1],1/(self.exponent_emission_reduction_theta2-1))
    elif self.time > 15:
        output = np.power(t15/self.BC[-1],1/(self.exponent_emission_reduction_theta2-1))
    else:
        output = 1

    if output > 1:
        output = 1           # To check that there are no negative emissions.
    return output
