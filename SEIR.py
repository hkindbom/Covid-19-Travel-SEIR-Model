import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

from Map import plot_graph, plot_map
from Metrics import calc_Lambda, norm_L_norm

class SEIR:

    def __init__(self, S0, E0, I0, R0, L, dt, steps, n, restrictions, beta=0.75, gamma=0.54, alpha=0.2, mobility=0.43, I_trade_off = 0.001):
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

        self.beta = beta
        self.gamma = gamma
        self.alpha = alpha
        self.epsilon = mobility/(365*dt)
        self.beta_per_country = np.array([beta, beta, beta, beta])
        self.restrictions = restrictions
        self.I_trade_off = I_trade_off

    def next_step(self, t):
        self.update_betas(t)
        dS = -self.epsilon * self.L.dot(self.S[:,t-1]) - (self.beta_per_country * np.multiply(self.S[:,t-1], self.I[:,t-1]).T).T
        dE = -self.epsilon * self.L.dot(self.E[:,t-1]) + (self.beta_per_country * np.multiply(self.S[:,t-1], self.I[:,t-1]).T).T - self.alpha*self.E[:,t-1]
        dI = -self.epsilon * self.L.dot(self.I[:,t-1]) + self.alpha*self.E[:,t-1] - self.gamma*self.I[:,t-1]
        dR = -self.epsilon * self.L.dot(self.R[:,t-1]) + self.gamma*self.I[:,t-1]

        self.S[:,t] = self.S[:,t-1] + dS*self.dt
        self.E[:,t] = self.E[:,t-1] + dE*self.dt
        self.I[:,t] = self.I[:,t-1] + dI*self.dt
        self.AI[:,t] = self.AI[:,t-1] + (-self.epsilon * self.L.dot(self.I[:,t-1]) + self.alpha*self.E[:,t-1])*self.dt
        self.R[:,t] = self.R[:,t-1] + dR*self.dt
        self.t[t] = self.dt*t

    def update_betas(self, t):
        if (t-1) % 15 == 0:
            for country_idx in range(len(self.restrictions)):
                if self.I[country_idx, t-1] > self.I_trade_off:
                    self.beta_per_country[country_idx] *= self.restrictions[country_idx]
                #if self.I[country_idx, t-1] < self.I_trade_off:
                #    self.beta_per_country[country_idx] = self.beta


def millions(x, pos):
    x = int(x)
    return '{:,}'.format(x).replace(',', ' ')


def plot_start_values(confirmed_cases, confirmed_recovered_cases, seir, n, countries, populations):
    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=2)
    fig.suptitle('Beta=' + str(seir.beta) + ', gamma=' + str(seir.gamma) + ', alpha=' + str(seir.alpha)) # or plt.suptitle('Main title')
    x = 0
    y = 0
    for i in range(n):
        t = [i for i in range(len(confirmed_cases[i]))]

        axs[x][y].plot(t[:150], confirmed_cases[i][:150], 'rx', fillstyle='none', label='Confirmed cases')
        axs[x][y].plot(seir.t[:150], populations[i]*seir.AI[i][:150], 'r', fillstyle='none', label='Accumulated Infective')

        axs[x][y].plot(t[:150], confirmed_recovered_cases[i][:150], 'bx', fillstyle='none', label='Confirmed recoverd cases')
        axs[x][y].plot(seir.t[:150], populations[i]*seir.R[i][:150], 'b', fillstyle='none', label='Recovered')

        axs[x][y].legend(loc="upper right")
        axs[x][y].title.set_text(countries[i])

        x+=1
        if x > 1:
            y+=1
            x=0

    plt.show()

def plot_traveling(confirmed_cases, confirmed_recovered_cases, seir, n, countries, populations):
    fig, axs = plt.subplots(2, 2)
    fig.tight_layout(pad=2)
    fig.suptitle('Beta=' + str(seir.beta) + ', gamma=' + str(seir.gamma) + ', alpha=' + str(seir.alpha)) # or plt.suptitle('Main title')
    x = 0
    y = 0
    for i in range(n):
        t = [i for i in range(len(confirmed_cases[i]))]
        axs[x][y].plot(t, confirmed_cases[i], 'rx', fillstyle='none', label='Confirmed cases')
        axs[x][y].plot(seir.t, populations[i]*seir.AI[i], 'r', fillstyle='none', label='Accumulated Infective')

        axs[x][y].plot(t, confirmed_recovered_cases[i], 'bx', fillstyle='none', label='Confirmed recoverd cases')
        axs[x][y].plot(seir.t, populations[i]*seir.R[i], 'b', fillstyle='none', label='Recovered')

        axs[x][y].legend(loc="upper right")
        axs[x][y].title.set_text(countries[i])

        x+=1
        if x > 1:
            y+=1
            x=0

    plt.show()

def import_confirmed(country, gamma):
    total_cases = []
    recovered_cases = []
    first_day_over_100 = None
    with open('owid-covid-data.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
        i=0
        for row in spamreader:
            daily_acc = row[4]
            try:
                daily_acc = float(row[4])
            except:
                daily_acc = len(total_cases) -1
            if row[2] == country and daily_acc > 100:
                if first_day_over_100 == None:
                    first_day_over_100 = i
                total_cases.append(daily_acc)
                recovered_cases.append(0)
        i+=1

    for i in range(0, len(total_cases)-int(1/gamma)):
        recovered_cases[i+int(1/gamma)] = total_cases[i]

    return total_cases, recovered_cases


def get_pop(country):
    try:
        xls = pd.ExcelFile('PopulationByCountry.xlsx')
        first_sheet = pd.read_excel(xls)
        row = first_sheet.loc[first_sheet['Country'] == country]
        return int(row['Population'])
    except:
        print('No population found for:', country, '!')
        return None

def get_country_to_sheet(country_sheet_df):
    country_idxs = country_sheet_df.index.tolist()
    country_names = country_sheet_df['GEO/TIME'].tolist()
    country_to_sheet = {}
    for country_idx in country_idxs:
        # to get passengers carried (second half of sheet numbers)
        country_to_sheet[country_names[country_idx]] = 'Data' + str(len(country_names) + 2 + country_idx)
    return country_to_sheet

def import_travel_as_W(countries='all'):
    xls = pd.ExcelFile('Traffic_between_EU_countries.xls')
    # Turkey is skipped since no good data
    first_sheet = pd.read_excel(xls, skiprows = 10, nrows = 34)
    country_to_sheet = get_country_to_sheet(first_sheet)

    if countries == 'all':
        countries = first_sheet['GEO/TIME'].tolist()

    values_missing = False

    W = np.zeros((len(countries), len(countries)))
    for country_idx, country in enumerate(countries):
        country_travel_df = pd.read_excel(xls, country_to_sheet[country], skiprows = 10, nrows = 34)
        other_countries = countries.copy()
        other_countries.remove(country)
        for other_country in other_countries:
            travel_rate = pd.to_numeric(country_travel_df.loc[country_travel_df['GEO/TIME'] == other_country]['2019'], errors ='coerce').fillna(0)

            W[country_idx, countries.index(other_country)] = travel_rate
            if W[country_idx, countries.index(other_country)] < 1:
                values_missing = True

    print('some travel data is 0: ', values_missing)
    return make_symmetric(W)

def make_symmetric(W):
    # Travel from country X to Y in excel sheet for country X may not be same as in sheet for country Y
    dim = W.shape[0]
    for row in range(dim):
        for col in range(dim):
            W[col, row] = W[row, col]
    return W

def get_D_from_W(W):
    row_sums = np.sum(W, axis=1)
    return np.diag(row_sums)

def plot_maps(countries, W, comp_mat, sep_eval, compartment, mobility):
    country_data = {}
    comp_eval = list(comp_mat[:, sep_eval])
    comp_eval = [round(num, 2) for num in comp_eval]
    country_data.update(zip(countries, comp_eval))

    plot_map(country_data, compartment, 'Regional levels of ' + str(compartment) + ' after ' + str(sep_eval) + ' time steps and eps = ' +str(mobility) )
    plot_graph(W, country_data, 'Undirected Travel Graph', 0.001)


def run_model(beta, gamma, alpha, steps, countries, restrictions, S0, E0, I0, R0, n, mobility, I_trade_off, step_eval):

    populations = []
    confirmed_cases = []
    confirmed_recovered_cases = []

    for country in countries:
        conf, rec = import_confirmed(country, gamma)
        confirmed_cases.append(conf)
        confirmed_recovered_cases.append(rec)
        populations.append(get_pop(country))

    total_population = sum(populations)

    W = import_travel_as_W(countries)
    D = get_D_from_W(W)
    L_un = (D - W)
    L = norm_L_norm(L_un, D)

    seir = SEIR(
        S0,  # S0
        E0,  # E0
        I0,  # I0
        R0,  # R0
        L,  # L
        1,  # dt
        steps,  # no of steps
        n,  # no of countries
        restrictions,
        beta=beta,
        gamma=gamma,
        alpha=alpha,
        mobility=mobility,
        I_trade_off=I_trade_off
    )

    for t in range(1, steps):
        seir.next_step(t)

    plot_maps(countries, W, seir.AI, step_eval, "AI", mobility)
    #plot_start_values(confirmed_cases, confirmed_recovered_cases, seir, n, countries, populations)
    return seir

def plot_vectors(x, y, title, countries):
    plt.figure()
    for country_idx, country_name in enumerate(countries):
            plt.plot(x, y[country_idx], label=country_name)
    plt.title(title, fontsize=16)
    plt.xlabel('I0 difference', fontsize=14)
    plt.ylabel('Lambda', fontsize=14)
    plt.legend()
    plt.show()

def main():
    beta = 0.247
    gamma = 0.1056
    alpha = 0.44

    steps = 300
    countries = ['Sweden', 'Denmark', 'Norway', 'Finland']
    n = len(countries)

    restrictions = np.array([0.71, 0.675, 0.625, 0.65])
    I_trade_off = 0.0001
    mobility = 0.43

    """
    S0 = [1] * 4
    E0 = [0] * 4
    I0 = [0.0001, 0.00001, 0.00001, 0.00001]  # [100/populations[i] for i in range(n)]
    R0 = [0] * 4
    for j in range(n):
        E0[j] = I0[j] * 2.5
        S0[j] -= (I0[j] + E0[j])
    """
    other_country_I0 = 0.000001
    SWE_I0s = [0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01]
    step_eval = 100
    Lambda_vecs = np.zeros((len(countries), len(SWE_I0s)))
    diff_I0s_vec = np.array(SWE_I0s) - other_country_I0
    for I0_diff_idx, SWE_I0 in enumerate(SWE_I0s):
        S0 = [1] * 4
        E0 = [0] * 4
        I0 = [SWE_I0, other_country_I0, other_country_I0, other_country_I0]  # [100/populations[i] for i in range(n)]
        R0 = [0] * 4
        for j in range(n):
            E0[j] = I0[j] * beta / gamma
            S0[j] -= (I0[j] + E0[j])

        seir_travel = run_model(beta, gamma, alpha, steps, countries, restrictions, S0, E0, I0, R0, n,
                         mobility, I_trade_off, step_eval)
        seir_no_travel = run_model(beta, gamma, alpha, steps, countries, restrictions, S0, E0, I0, R0, n,
                                0, I_trade_off, step_eval)
        RT = seir_travel.R[:,step_eval]
        RNT = seir_no_travel.R[:,step_eval]

        for country_idx, country_name in enumerate(countries):
            Lambda_vecs[country_idx, I0_diff_idx] = calc_Lambda(RT[country_idx], RNT[country_idx])

    plot_vectors(diff_I0s_vec, Lambda_vecs, 'Lambda showing impact of travelling with SWE as epicentrum & r!=1', countries)

if __name__ == "__main__":
    main()