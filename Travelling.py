# Class for travelling between countries

class Travelling:

    def __init__(self):
        '''
            SWE   ENG
        SWE  0    100
        ENG  200   0
        '''
        self.travel_rates = [[0, 100],[200, 0]]

    def generate_travellers(self, N, S, I, R, nr_travellers):
        # nr suseptibles travelling
        S_travel = nr_travellers*S/N

        # nr infective travelling
        I_travel = nr_travellers*I/N

        # nr recovered travelling
        R_travel = nr_travellers*R/N

        # Add random generation

        return S_travel, I_travel, R_travel

