"""
DTMC SIS model

sources:
- https://www.mdpi.com/2227-7390/6/8/128/htm#B26-mathematics-06-00128
- "An Introduction to Stochastic Epidemic Models" section 3.2
"""
import numpy as np
import matplotlib.pyplot as plt
from Travelling import Travelling

class DTMC:
    def __init__(self, beta, gamma, N, I_start):
        """
        :param beta: contact rate (unit days)
        :param gamma: recovery rate (unit days)
        :param N: population size
        """
        self.beta = beta
        self.gamma = gamma
        self.N = N
        self.I = I_start

    def get_transition_probs(self, delta_t):
        incr_prob = (self.beta/self.N)*self.I*(self.N-self.I)*delta_t
        decr_prob = self.gamma*self.I*delta_t
        stay_prob = 1 - incr_prob - decr_prob
        stay_prob = 0 if stay_prob < 0 else stay_prob

        return [incr_prob, decr_prob, stay_prob]

    def add_travellers(self, S_diff, I_diff):
        self.I += I_diff
        self.N += I_diff + S_diff

    def next_step(self, delta_t):
        transitions = np.array([1, -1, 0])
        transition_probs = self.get_transition_probs(delta_t)
        state_diff = np.random.choice(transitions, 1, p=transition_probs)[0]
        self.I += state_diff

    def get_delta_t(self, max_N_0):
        return 1/(self.gamma*max_N_0) if self.beta <= self.gamma else \
            4*self.beta/((self.beta + self.gamma)**2 * max_N_0)



def plot_result(end_time, I_over_time_0, I_over_time_1):
    plt.plot(I_over_time_0, label='Infectives country 0, T='+str(end_time)+' days')
    plt.plot(I_over_time_1, label='Infectives country 1, T='+str(end_time)+' days')

    plt.xlabel('Time steps')
    plt.ylabel('Infectives')

    plt.legend()
    plt.show()


def main():
    travel_rates = [[0, 100], [100, 0]]
    travel_model = Travelling(travel_rates)

    # Country 0
    beta_0 = 1
    gamma_0 = 0.5
    N_0 = 1000000
    I_start_0 = 1000
    epidemic_model_0 = DTMC(beta_0, gamma_0, N_0, I_start_0)

    # Country 1
    beta_1 = 1
    gamma_1 = 0.5
    N_1 = 1000000
    I_start_1 = 500
    epidemic_model_1 = DTMC(beta_1, gamma_1, N_1, I_start_1)

    end_time = 10
    max_N_0 = max(N_0 + (travel_rates[1][0] - travel_rates[0][1]) * end_time, N_0)
    max_N_1 = max(N_1 + (travel_rates[0][1] - travel_rates[1][0]) * end_time, N_1)

    delta_t = min(epidemic_model_0.get_delta_t(max_N_0), epidemic_model_1.get_delta_t(max_N_1))
    steps = int(end_time // delta_t)
    I_over_time_0 = np.zeros(steps)
    I_over_time_1 = np.zeros(steps)

    for step in range(steps):
        I_over_time_0[step] = epidemic_model_0.I
        I_over_time_1[step] = epidemic_model_1.I

        S_from_0, I_from_0 = travel_model.generate_travellers_from(epidemic_model_0.N,
                                                                   epidemic_model_0.I, 0, 1, delta_t)
        S_from_1, I_from_1 = travel_model.generate_travellers_from(epidemic_model_1.N,
                                                                   epidemic_model_1.I, 1, 0, delta_t)

        epidemic_model_0.add_travellers(S_from_1 - S_from_0, I_from_1 - I_from_0)
        epidemic_model_0.next_step(delta_t)

        epidemic_model_1.add_travellers(S_from_0 - S_from_1, I_from_0 - I_from_1)
        epidemic_model_1.next_step(delta_t)

    plot_result(end_time, I_over_time_0, I_over_time_1)


main()
