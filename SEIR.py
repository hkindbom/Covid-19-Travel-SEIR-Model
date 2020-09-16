import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from Metrics import norm_L

class SEIR:

    def __init__(self, S0, E0, I0, R0, L, dt, steps, n, beta=0.75, gamma=0.54, alpha=0.2):
        self.S = np.zeros((n, steps))
        self.E = np.zeros((n, steps))
        self.I = np.zeros((n, steps))
        self.AI = np.zeros((n, steps))
        self.R = np.zeros((n, steps))
        self.t = np.zeros((steps))

        self.S[:,0] = S0
        self.E[:,0] = E0
        self.I[:,0] = I0
        self.AI[:,0] = I0
        self.R[:,0] = R0
        self.t[0] = 0

        self.L = np.array(L)
        self.dt = dt

        self.beta  = beta
        self.gamma = gamma
        self.alpha = alpha
        self.epsilon = 0.0

    def next_step(self, t):
        dS = -self.epsilon * self.L.dot(self.S[:,t-1]) - self.beta * np.multiply(self.S[:,t-1], self.I[:,t-1])
        dE = -self.epsilon * self.L.dot(self.E[:,t-1]) + self.beta * np.multiply(self.S[:,t-1], self.I[:,t-1]) - self.alpha*self.E[:,t-1]
        dI = -self.epsilon * self.L.dot(self.I[:,t-1]) + self.alpha*self.E[:,t-1] - self.gamma*self.I[:,t-1]
        dR = -self.epsilon * self.L.dot(self.R[:,t-1]) + self.gamma*self.I[:,t-1]

        self.S[:,t] = self.S[:,t-1] + dS*self.dt
        self.E[:,t] = self.E[:,t-1] + dE*self.dt
        self.I[:,t] = self.I[:,t-1] + dI*self.dt
        self.AI[:,t] = self.AI[:,t-1] + (-self.epsilon * self.L.dot(self.I[:,t-1]) + self.alpha*self.E[:,t-1])*self.dt
        self.R[:,t] = self.R[:,t-1] + dR*self.dt
        self.t[t] = self.dt*t



def millions(x, pos):
    x = int(x)
    return '{:,}'.format(x).replace(',', ' ')


def plot_SEIR(confirmed_cases, seir, n):
    '''fig, axs = plt.subplots(n)
    for i in range(n):
        t = [i for i in range(len(confirmed_cases[i]))]
        #axs[i].plot(seir.t, seir.S[i], 'g', fillstyle='none', label='Susceptible')
        #axs[i].plot(seir.t, seir.E[i], 'r', fillstyle='none', label='Exposed')
        axs[i].plot(t, confirmed_cases[i], 'ro', fillstyle='none', label='Confirmed cases')
        axs[i].plot(seir.t, 10000000*seir.AI[i], 'r', fillstyle='none', label='Accumulated Infective')
        #axs[i].plot(seir.t, seir.R[i], 'y', fillstyle='none', label='Recovered')
        #axs[i].legend(loc="upper right")

    plt.show()'''
    fig, axs = plt.subplots(n, 1)
    for i in range(n):
        time_confirmed = [j for j in range(len(confirmed_cases[i]))]
        axs[i].plot(time_confirmed, confirmed_cases[i], 'o', label='Uppmätt')
        #axs[i].plot(x, infected, label='Prognos')
        #axs[1].plot(time_confirmed, deaths_confirmed, 'o', label='Uppmätt')
        #axs[1].plot(x, D, label='Prognos')

        axs[0].get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(millions))
        axs[0].legend(loc="upper left")
        #axs[1].get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(millions))
        #axs[1].legend(loc="upper left")
        #plt.xlabel('Dagar sedan 100 bekräftade fall', fontsize=14)
        #axs[0].set_ylabel('Bekräftat smittade', fontsize=14)
        #axs[1].set_ylabel('Bekräftat döda', fontsize=14)
        # plt.title('Prognos: intensivvårdssökande i Sverige', fontsize=18)
    plt.show()

def import_confirmed(country):
    total_cases = []
    with open('owid-covid-data.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
        for row in spamreader:
            #row[2] = country
            #row[3] = date
            #row[4] = total_cases
            #row[5] = new_cases
            if row[2] == country:
                total_cases.append(row[4])

    return total_cases

if __name__ == "__main__":
    steps = 300
    D = np.array([[10, 0, 0],[0, 10, 0],[0, 0, 10]])
    W = np.array([[0, 5, 5],[5, 0, 5],[5, 5, 0]])
    L_un = D - W
    L = norm_L(L_un, D)
    n = 3
    confirmed_cases = [import_confirmed('Sweden'), import_confirmed('Norway'), import_confirmed('United Kingdom')]
    seir = SEIR(
                [0.95, 1.0, 0.96], #S0
                [0.025, 0, 0.02], #E0
                [0.025, 0, 0.02], #I0
                [0, 0, 0], #R0
                L, #L
                1, #dt
                steps, # no of steps
                n, # no of countries
                beta=0.8,
                gamma=0.54,
                alpha=0.54
               )

    for t in range(1, steps):
        seir.next_step(t)

    plot_SEIR(confirmed_cases, seir, n)
