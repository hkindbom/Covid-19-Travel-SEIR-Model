# André Gerbaulet & Robin Sandström
#
# Based on models from "Mathematical Models in Population Biology and Epidemiology" by F. Brauer and C. Castillo-Chavez
#
# The SIR model with disease death will be modelled as N = S + I + R
#
# Assumed that people go from (susceptible state) to (infective state) to (removed state)
#
# N - Population size
# S - The number of people susceptible to the disease
# I - The number of infective people
# R - The number of previously infected people people removed from the disease, either by immunization or death
#
# S, I, R is assumed to follow the differential equation: N = S + I + R
#
# S' = - beta(N) * S * I              People going to infective state
# I' = beta(N) * S * I - alpha * I    People coming from susceptible state & people going to removed state
# N' = - (1 - f) * alpha * I          People leaving the population (death)
#
# beta(N) = C(N) / N,    C(N) = lambda * N^a,    a = 0.05

import numpy as np
import xlrd
import matplotlib.pyplot as plt
from Travelling import Travelling

class SIR_model:
    def __init__(self, beta=0, gamma=0, S=0, I=0, R=0, N=0, dt=1):
        self.beta = beta
        self.gamma = gamma
        self.S_vector = [S]
        self.I_vector = [I]
        self.R_vector = [R]
        self.N_vector = [N]
        self.t_vector = [0]
        self.dt = dt

    @property
    def S(self):
        return self.S_vector[-1]

    @property
    def I(self):
        return self.I_vector[-1]

    @property
    def R(self):
        return self.R_vector[-1]

    @property
    def N(self):
        return self.N_vector[-1]

    @property
    def set_N(self, val):
        self.N_vector[-1] = val

    @property
    def t(self):
        return self.t_vector[-1]

    def add_travellers(self, S_diff, I_diff):
        self.I_vector[-1] = self.I + I_diff
        self.N_vector[-1] = self.N + I_diff + S_diff

    # Adds new state
    def update(self, change_vector):

        self.S_vector.append(self.S + change_vector[0])
        self.I_vector.append(self.I + change_vector[1])
        self.R_vector.append(self.R + change_vector[2])
        self.t_vector.append(self.t + self.dt)
        return

    # Makes one step using Euler method
    def next_step(self):
        dV = [
                (-self.beta * self.I * self.S)/self.N,
                (self.beta * self.I * self.S)/self.N - self.gamma*self.I,
                self.gamma*self.I
             ]

        return dV*self.dt



if __name__ == "__main__":

    travel_rates = [[0, 0], [0, 0]]
    travel_model = Travelling(travel_rates)

    dt = 5

    SIR_sweden = SIR_model(
                            beta=0.85,
                            gamma=0.45,
                            S=10000000,
                            I=10,
                            R=0,
                            N=10000000,
                            dt=dt
                           )

    SIR_nz = SIR_model(
                            beta=0.85,
                            gamma=0.45,
                            S=10000000,
                            I=0,
                            R=0,
                            N=10000000,
                            dt=dt
                           )

    for t in range(0, 100):
        new_state_sv = SIR_sweden.next_step()
        new_state_nz = SIR_nz.next_step()

        S_from_0, I_from_0 = travel_model.generate_travellers_from(SIR_sweden.N,
                                                                   SIR_sweden.I, 0, 1, dt)

        S_from_1, I_from_1 = travel_model.generate_travellers_from(SIR_nz.N,
                                                                   SIR_nz.I, 1, 0, dt)


        SIR_sweden.add_travellers(S_from_1 - S_from_0, I_from_1 - I_from_0)
        SIR_nz.add_travellers(S_from_0 - S_from_1, I_from_0 - I_from_1)

        SIR_sweden.update(new_state_sv)
        SIR_nz.update(new_state_nz)

    fig, axs = plt.subplots(2)
    fig.suptitle('Sweden / NZ')
    axs[0].plot(SIR_sweden.t_vector, SIR_sweden.I_vector, 'g', fillstyle='none', label='Infective')
    axs[0].plot(SIR_sweden.t_vector, SIR_sweden.S_vector, 'r', fillstyle='none', label='Susceptible')
    axs[0].plot(SIR_sweden.t_vector, SIR_sweden.R_vector, 'b', fillstyle='none', label='Recovered')
    axs[0].legend(loc="upper right")

    axs[1].plot(SIR_nz.t_vector, SIR_nz.I_vector, 'g', fillstyle='none', label='Infective')
    axs[1].plot(SIR_nz.t_vector, SIR_nz.S_vector, 'r', fillstyle='none', label='Susceptible')
    axs[1].plot(SIR_nz.t_vector, SIR_nz.R_vector, 'b', fillstyle='none', label='Recovered')
    axs[1].legend(loc="upper left")
    plt.show()
