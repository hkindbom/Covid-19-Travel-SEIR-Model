# File containing function for calculating metrics
import numpy as np
from scipy.linalg import fractional_matrix_power

def calc_QoS(Lambdas, rs, epsilon, c):
    """
    :param Lambdas: list of lambdas for each country
    :param rs: list of r for each country
    :param epsilon: travel mobility
    :param c: list of three cost floats
    :return: QoS value (float)
    """
    avg = 0
    N = len(Lambdas)
    for Lambda, r in zip(Lambdas, rs):
        avg += (1/N) * (c[0]*r + c[1]*(1 - Lambda))
    return avg + c[2]*epsilon

def calc_Lambda(RT, RNT):
    """

    :param RT: int nr recovered with travel
    :param RNT: int nr recovered without travel
    :return:
    """
    return (RT - RNT)/RNT

def norm_L(L_un, D):
    D_mod = fractional_matrix_power(D, -0.5)
    return D_mod @ L_un @ D_mod

def norm_L_norm(L_un, D):
    D_mod = fractional_matrix_power(D, -1)
    return D_mod @ L_un

Lambda = [0.6, 0.2, 0.5]
rs = [0.8, 0.2, 0.7]
epsilon = 0.3
c = [1, 4, 2]

RT = 40
RNT = 30



D = np.array([[3000, 0, 0],[0, 5000, 0],[0, 0, 2000]])
W = np.array([[0, 3000, 0],[3000, 0, 2000],[0, 2000, 0]])
L_un = D - W
