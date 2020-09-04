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

class SIR_model:
    def __init__(self, covid19_filename, population_filename, country, p=None, steps=100, min_cases=100):
        self.covid19 = self.read_covid19_file(covid19_filename)
        self.populations = self.read_population_file(population_filename)
        self.country = country
        self.steps = steps

        self.cases_confirmed = self.covid19[country][0][self.covid19[country][0] >= min_cases]
        self.deaths_confirmed = self.covid19[country][1][self.covid19[country][0] >= min_cases]
        self.n = np.shape(self.cases_confirmed)[0]    # Number of measurements
        self.x_confirmed = np.arange(self.n)

        if p is None:
            alpha = 0.5
            theta = 0.5
            gamma = 1  # 5 * 10**(-2) # (South Korea)
            a = 0.05
            f = 0.99
            self.par = np.array([alpha, theta, gamma, a, f], dtype=float)
        else:
            self.par = p

        N_start = self.par[2] * self.populations[self.country]
        I_start = min(self.cases_confirmed)
        R_start = min(self.deaths_confirmed)
        self.SIN_first = np.array([[N_start - I_start - R_start],
                                   [I_start],
                                   [N_start - (1 - self.par[4]) * R_start],
                                   [0]])

        I_end = max(self.cases_confirmed)
        R_end = max(self.deaths_confirmed)
        self.SIN_end = np.array([[N_start - I_end - R_end],
                                 [I_end],
                                 [N_start - (1 - self.par[4]) * R_end],
                                 [self.n - 1]])

        self.X_matrix = np.array([np.power(self.x_confirmed, 0),
                                  np.power(self.x_confirmed, 1),
                                  np.power(self.x_confirmed, 2)], dtype=float).transpose()

    def update_initial_condition_for_gamma(self, gamma):
        N_start = gamma * self.populations[self.country]
        I_start = min(self.cases_confirmed)
        R_start = min(self.deaths_confirmed)
        self.SIN_first = np.array([[N_start - I_start - R_start],
                                   [I_start],
                                   [N_start - (1 - self.par[4]) * R_start],
                                   [0]])


    # Returns dictionary on the form {'Sweden' : [np.array(# smittade), np.array(# döda)]}
    @staticmethod
    def read_covid19_file(filename):
        workbook = xlrd.open_workbook(filename)
        worksheet = workbook.sheet_by_index(0)
        current_country = ''
        dictionary = {}
        for row in range(1, worksheet.nrows):
            temp_country = worksheet.cell_value(row, 1)
            if temp_country == current_country:
                cases = np.append(cases, np.array([worksheet.cell_value(row, 2)], dtype=int), 0)
                deaths = np.append(deaths, np.array([worksheet.cell_value(row, 3)], dtype=int), 0)
            else:
                if current_country != '':
                    dictionary[current_country] = [np.cumsum(np.flip(cases)), np.cumsum(np.flip(deaths))]
                cases = np.array([worksheet.cell_value(row, 2)], dtype=int)
                deaths = np.array([worksheet.cell_value(row, 3)], dtype=int)
                current_country = temp_country
        return dictionary

    # Returns dictionary on the form {'Sweden' : int(population size)}
    @staticmethod
    def read_population_file(filename):
        workbook = xlrd.open_workbook(filename)
        worksheet = workbook.sheet_by_index(0)
        dictionary = {}
        for row in range(1, worksheet.nrows):
            dictionary[worksheet.cell_value(row, 0)] = int(worksheet.cell_value(row, 1))
        return dictionary
        return dictionary

    # Used in SIR model: infection rate, higher beta => more infections
    def beta(self, N, p):
        return p[1] * (N**(p[3] - 1))

    # Makes one step using Euler method
    def next_step(self, SIN, dt, p):
        SINp = np.array([[- self.beta(SIN[2], p) * SIN[0] * SIN[1]],
                         [self.beta(SIN[2], p) * SIN[0] * SIN[1] - p[0] * SIN[1]],
                         [- (1 - p[4]) * p[0] * SIN[1]]]).transpose()
        return SIN + dt * SINp

    # Takes M step with Eulers method
    def prediction(self, p_input=None, P=0, n=None, SIN_0=None):
        if p_input is None:
            p = self.par
        else:
            p = p_input
        if SIN_0 is None:
            SIN_0 = self.SIN_first
        if n is None:
            n = self.steps
        if P == 0:
            T = SIN_0[3] + self.n
            x = np.arange(n + 1) / n * T
        else:
            T = P
            x = np.arange(n + 1) / n * T + SIN_0[3]
        dt = T / n
        SIN = np.zeros((3, n + 1))
        SIN[:, 0] = SIN_0[0:3, 0].transpose()
        for i in range(n):
            SIN[:, i+1] = self.next_step(SIN[:, i], dt, p)
        return [x, SIN]

    # Calculates the residuals of fitted curve vs. measured points
    def get_residuals(self, x, y, x_hat, y_hat):
        n = np.shape(x)[0]
        m = np.shape(x_hat)[0]
        r = np.zeros(n, dtype=float)
        for i in range(n):
            for j in range(m-1):
                if x_hat[j] == x[i]:
                    r[i] = y[i] - y_hat[j]
                    break
                elif (x_hat[j] < x[i]) and (x_hat[j+1] > x[i]):
                    w = (x_hat[j+1] - x[i]) / (x_hat[j+1] - x_hat[j])
                    r[i] = y[i] - (w * y_hat[j] + (1 - w) * y_hat[j+1])
                    break
        return r

    def fit_3rd_degree_polynomial(self, y):
        XtX_inv = np.linalg.inv(np.dot(self.X_matrix.transpose(), self.X_matrix))
        parameters = np.dot(np.dot(XtX_inv, self.X_matrix.transpose()), y.transpose())
        return parameters

    def estimate_mean_second_derivative(self, y, dt=1):
        s = np.shape(y)[0]
        sd = np.zeros(s-2)
        for i in range(s-2):
            sd[i] = (y[i+2] - 2 * y[i+1] + y[i]) / (dt**2)
        return np.mean(sd)

    def estimate_theta(self, alpha_input, gamma_input):
        theta_min = 0.5
        d_theta = 0.25
        for j in range(10):
            parameters = np.array([alpha_input, theta_min, gamma_input, self.par[3], self.par[4]])
            prediction = self.prediction(parameters)
            x_prediction = prediction[0]
            SIN = prediction[1]
            S = SIN[0, :]
            N = SIN[2, :]
            residuals = self.get_residuals(self.x_confirmed, self.cases_confirmed, x_prediction, (N - S))
            if np.mean(residuals) == 0:
                break
            elif np.mean(residuals) < 0:
                theta_min -= d_theta
                d_theta *= 0.5
            else:
                theta_min += d_theta
                d_theta *= 0.5
        return theta_min, residuals, parameters


    def estimate_alpha_theta(self, gamma_input):
        alpha_min = .5
        d_alpha = 0.25
        for i in range(10):
            theta_min, residuals, parameters = self.estimate_theta(alpha_min, gamma_input)
            fit = self.fit_3rd_degree_polynomial(residuals)
            if fit[1] == 0:
                break
            elif fit[1] < 0:
                alpha_min -= d_alpha
                d_alpha *= 0.5
            else:
                alpha_min += d_alpha
                d_alpha *= 0.5
        return alpha_min, theta_min, residuals, parameters


    def minimize_parameters(self):
        gamma_values = [1]  # np.power(10, -np.arange(100+1)/20)
        alpha_best = 0
        theta_best = 0
        gamma_best = 0
        rss_best = np.inf
        c = 1
        for gamma_loop in gamma_values:
            if c % 10 == 0:
                print(c)
            c += 1
            self.update_initial_condition_for_gamma(gamma_loop)
            alpha_min, theta_min, residuals, parameters = self.estimate_alpha_theta(gamma_loop)
            current_rss = np.dot(residuals, residuals)  # self.get_rss(p)
            if current_rss < rss_best:
                alpha_best = alpha_min
                theta_best = theta_min
                gamma_best = gamma_loop
                rss_best = current_rss
        return [alpha_best, theta_best, gamma_best]

    def get_rss(self, p):
        prediction = self.prediction(p)
        x_prediction = prediction[0]
        SIN = prediction[1]
        S = SIN[0, :]
        I = SIN[1, :]
        N = SIN[2, :]
        R = N - S - I
        res = self.get_residuals(self.x_confirmed, self.cases_confirmed, x_prediction, I + R)
        rss = np.sum(np.power(res, 2))

if __name__ == "__main__":
    covid19_filename = 'COVID-19-geographic-disbtribution-worldwide-2020-03-13.xls'
    population_filename = 'PopulationByCountry.xlsx'
    country = 'Sweden'
    sir = SIR_model(covid19_filename, population_filename, country)
    par = sir.minimize_parameters()
    # par = [0.14599609375, 0.15673828125, 0.14125375446227545]
    # par = [0.00048828125, 0.11181640625, 0.35481338923357547]  # Sweden
    # par = [0.00048828125, 0.11181640625, 0.35481338923357547] # Sweden 5 % death rate
    # par = [0.44482421875, 0.37353515625, 0.0017782794100389228]  # Italy
    f = 0.95
    a = 0.05
    p = np.array([par[0], par[1], par[2], a, f])

    sir.update_initial_condition_for_gamma(par[2])

    prediction = sir.prediction(p, 0, 500)
    x_prediction = prediction[0]
    SIN = prediction[1]

    '''for t in range(0,100):
        SWEDEN_SIR.next_step()
        UK_SIR.next_step()
        inc_sweden = TRAVEL(UK.state)
        inc_uk = TRAVEL(SWEDEN_SIR.state)
        sweden_sir.add_travels(inc_sweden)
        uk_sir.add_travels(inc_uk)'''

    S = SIN[0, :]
    I = SIN[1, :]
    N = SIN[2, :]
    R = N - S - I
    plt.plot(sir.x_confirmed, sir.cases_confirmed, 'bo', fillstyle='none')
    plt.plot(x_prediction, I)
    R_0 = sir.par[1] * (sir.par[2] * (sir.populations[country])**(sir.par[3] - 1)) / sir.par[0]
    print(R_0)
    # plt.plot(sir.x_confirmed, res, 'bo', fillstyle='none')
    # plt.plot([min(x_prediction), max(x_prediction)], [0, 0])
    plt.show()
