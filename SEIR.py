import numpy as np
import matplotlib.pyplot as plt

class SEIR:

    def __init__(self, S0, E0, I0, R0, L, dt, steps, n, beta=0.75, gamma=0.54, alpha=0.2):
        self.S = np.zeros((n, steps))
        self.E = np.zeros((n, steps))
        self.I = np.zeros((n, steps))
        self.R = np.zeros((n, steps))
        self.t = np.zeros((steps))

        self.S[:,0] = S0
        self.E[:,0] = E0
        self.I[:,0] = I0
        self.R[:,0] = R0
        self.t[0] = 0

        self.L = np.array(L)
        self.dt = dt

        self.beta  = beta
        self.gamma = gamma
        self.alpha = alpha

    def next_step(self, t):
        dS = - self.beta * np.multiply(self.S[:,t-1], self.I[:,t-1])
        dE = self.beta * np.multiply(self.S[:,t-1], self.I[:,t-1]) - self.alpha*self.E[:,t-1]
        dI = self.alpha*self.E[:,t-1] - self.gamma*self.I[:,t-1]
        dR = self.gamma*self.I[:,t-1]

        self.S[:,t] = self.S[:,t-1] + dS*self.dt
        self.E[:,t] = self.E[:,t-1] + dE*self.dt
        self.I[:,t] = self.I[:,t-1] + dI*self.dt
        self.R[:,t] = self.R[:,t-1] + dR*self.dt

        self.t[t] = self.dt*t

def plot_SEIR(seir, index):
    fig, axs = plt.subplots(2)
    axs[0].plot(seir.t_vector, seir.I[index], 'g', fillstyle='none', label='Infective')
    return

if __name__ == "__main__":
    steps = 10000
    seir = SEIR(
                [0.99, 0.95], #S0
                [0, 0], #E0
                [0.01, 0.05], #I0
                [0, 0], #R0
                [[0, 0], [0, 0]], #L
                0.01, #dt
                steps, # no of steps
                2, # no of countries
                beta=0.71,
                gamma=0.54,
                alpha=0.54
               )
    for t in range(1, steps):
        seir.next_step(t)
        print(seir.S[0][t] + seir.E[0][t] + seir.I[0][t] + seir.R[0][t])
    #print(seir.I[0], seir.E[0])
    fig, axs = plt.subplots(2)
    axs[0].plot(seir.t, seir.S[0], 'g', fillstyle='none', label='Susceptible')
    axs[0].plot(seir.t, seir.E[0], 'r', fillstyle='none', label='Exposed')
    axs[0].plot(seir.t, seir.I[0], 'b', fillstyle='none', label='Infective')
    axs[0].plot(seir.t, seir.R[0], 'y', fillstyle='none', label='Recovered')
    axs[0].legend(loc="upper right")


    axs[1].plot(seir.t, seir.S[1], 'g', fillstyle='none', label='Susceptible')
    axs[1].plot(seir.t, seir.E[1], 'r', fillstyle='none', label='Exposed')
    axs[1].plot(seir.t, seir.I[1], 'b', fillstyle='none', label='Infective')
    axs[1].plot(seir.t, seir.R[1], 'y', fillstyle='none', label='Recovered')
    axs[1].legend(loc="upper left")
    plt.show()
