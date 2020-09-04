# Class for travelling between countries

class Travelling:

    def __init__(self, travel_rates):
        '''
        E.g.
            SWE   ENG
        SWE  0    100
        ENG  200   0
        = [[0, 100], [200, 0]]
        '''
        self.travel_rates = travel_rates # per time unit

    def generate_travellers_from(self, N, I, from_idx, to_idx, delta_t):
        S = N - I
        # nr suseptibles travelling
        S_travel = delta_t * self.travel_rates[from_idx][to_idx] * S / N

        # nr infective travelling
        I_travel = delta_t * self.travel_rates[from_idx][to_idx] * I / N

        return S_travel, I_travel

