import numpy as np

class single_particle(object):
    def sigma_x(self, state, site):
        return state ^ (1 << site)

    def sigma_z(self, state, site):
        return (((state >> site) & 1) << 1) - 1

    def site_shift(self, n_spins, state):
        bulk_bits = (state << 1)
        last_bit = ((1 << (n_spins - 1)) & state) >> (n_spins - 1)
        sub = last_bit << n_spins
        new_state = bulk_bits + last_bit - sub
        return new_state


class hamiltonian(object):
    def __init__(self, n_spins=1, Jx=1., Jy=1., Jz=1., Hx=0., Hz=0., openBC=False):
        mu = (0.927400968 / 1.38064852)  # magnetic moment of an electron in Kelvin/Tesla (assuming s=1/2 and g=2)

        dim = 1 << n_spins
        self.matrix = np.zeros(shape=(dim, dim))
        sp = single_particle()
        #	print 'size f the Hilbert space: %d' %(dim)
        #	sys.stdout.flush()

        n_max = n_spins
        if openBC:
            n_max = n_spins - 1

        # print("check interaction sites")
        for site in range(0, n_max):

            site_next = site + 1
            if (site_next == n_spins): site_next = 0

            # print(site, site_next)

            for state in range(0, dim):
                #		Exchange terms
                sigma_1 = sp.sigma_z(state, site)
                sigma_2 = sp.sigma_z(state, site_next)
                self.matrix[state, state] -= Jz * sigma_2 * sigma_1  # sigma_z_1 sigma_z_2

                middle_state = sp.sigma_x(state, site)
                new_state = sp.sigma_x(middle_state, site_next)
                self.matrix[
                    new_state, state] -= Jx - Jy * sigma_2 * sigma_1  # sigma_x_1 sigma_x_2 + sigma_y_1 sigma_y_2

                #		Zeeman terms along z and along x
                new_state = sp.sigma_x(state, site)
                self.matrix[new_state, state] += mu * Hx
                self.matrix[state, state] += mu * Hz * sigma_1

        if openBC:  # missing Zeeman interaction
            site = n_max
            for state in range(0, dim):
                new_state = sp.sigma_x(state, site)
                sigma_1 = sp.sigma_z(state, site)
                self.matrix[new_state, state] += mu * Hx
                self.matrix[state, state] += mu * Hz * sigma_1

class S_total_x(object):
    def __init__(self, n_spins=1):
        dim = 1 << n_spins
        self.matrix = np.zeros(shape=(dim, dim))
        sp = single_particle()
        for state in range(0, dim):
            for site in range(0, n_spins):
                new_state = sp.sigma_x(state, site)
                self.matrix[new_state, state] += 0.5

class S_total_y(object):
    def __init__(self, n_spins=1):
        dim = 1 << n_spins
        self.matrix = np.zeros(shape=(dim, dim)) + np.zeros(shape=(dim, dim))*1j
        sp = single_particle()
        for state in range(0, dim):
            for site in range(0, n_spins):
                sigma_1 = sp.sigma_z(state, site)
                new_state = sp.sigma_x(state, site)
                self.matrix[new_state, state] += (0. + 0.5*sigma_1*1j)

class S_total_z(object):
    def __init__(self, n_spins=1):
        dim = 1 << n_spins
        self.matrix = np.zeros(shape=(dim, dim))
        sp = single_particle()
        for state in range(0, dim):
            for site in range(0, n_spins):
                self.matrix[state, state] += 0.5 * sp.sigma_z(state, site)



class translation(object):
    def __init__(self, n_spins=1):
        sp = single_particle()
        dim = 1 << n_spins
        self.matrix = np.zeros(shape=(dim, dim))
        for state in range(0, dim):
            new_state = sp.site_shift(n_spins, state)
            self.matrix[new_state, state] = 1

