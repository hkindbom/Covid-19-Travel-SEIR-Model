# File containing function for calculating metrics

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


Lambda = [0.6, 0.2, 0.5]
rs = [0.8, 0.2, 0.7]
epsilon = 0.3
c = [1, 4, 2]

print(calc_QoS(Lambda, rs, epsilon, c))


RT = 40
RNT = 30
print(calc_Lambda(RT, RNT))
