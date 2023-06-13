from math import sqrt
from random import randrange
import os
from math import exp, pi, acos, atan, sin, cos

import numpy as np
# from ase.cluster.cubic import SimpleCubicFactory
# from ase.lattice.tetragonal import CenteredTetragonalFactory
from mrmc_package import k_range, deltaE_shift, back_k_space, norm_fft

from time import perf_counter as timer


def get_distance(coor):
    return np.array([sqrt((_ ** 2).sum()) for _ in coor])

def sort(xyz, compound):
    distance = get_distance(xyz)
    index = distance.argsort()
    distance.sort()
    coordinate = xyz.take([index], 0)[0]
    element = compound.take([index], 0)[0]
    return coordinate, distance, element

def metropolis(r_0, r, tau=10e-5):
    if r_0 > r:
        return True
    else:
        if tau == 0:
            return False
        met = exp(-(r - r_0) / tau)
        judge = randrange(0, 100)
        if (judge / 100) < met:
            return True
        else:
            return False


class EXP:
    def __init__(self, file_name, k_start, k_end, r_start, r_end):
        self.name = file_name
        self.k_start = k_start
        self.k_end = k_end
        self.r_start = r_start
        self.r_end = r_end
        self.k0 = np.array([])
        self.k = np.array([])
        self.chi = np.array([])
        self.ft = np.array([])
        self.cross = np.array([])
        self.r = np.arange(0, 6 + pi / 102.4, pi / 102.4)
        self.r_range = np.where((r_start < self.r) & (self.r < r_end))[0]
        self.r_cut = np.array([])
        self.ft_cut = np.array([])

        self.init()
        self.max = [np.argmax(np.abs(self.chi)), np.max(np.abs(self.chi))]

    def init(self):
        self.read_exp(self.name)
        self.chi_bottom = np.sum(self.chi ** 2)
        self.ft_bottom = np.sum(self.ft_cut ** 2)
        self.cross_bottom = np.sum(self.cross ** 2)

    def read_exp(self, filename='703K_diff_H2-dry.rex'):
        # reading oscillation data
        with open(filename, 'r') as f_exp:
            k = np.array([])
            xi = np.array([])
            while True:
                lines = f_exp.readline(11)
                if not lines or not lines.find('[ED_BEGIN]') == -1:
                    break
            while True:
                data = f_exp.readline()
                if data.isspace():
                    break
                temp = [i for i in data.split()]
                k = np.append(k, float(temp[0]))
                xi = np.append(xi, float(temp[1]))
            self.k0 = k
            self.process(xi)

    '''def read_bk(self, filename='Tcu.ex3', k_max=20.0):
        # reading oscillation data
        from larch.xafs import autobk
        from larch.io import read_ascii
        f_exp = open(filename, 'r')
        f_bk = open('temp_bk.dat', 'w')
        f_bk.write('#   energy        xmu\n')
        while True:
            lines = f_exp.readline(11)
            if not lines or not lines.find('[EX_BEGIN]') == -1:
                break
        while True:
            data = f_exp.readline()
            if data.isspace():
                break
            f_bk.write(data)
        f_exp.close()
        f_bk.close()
        x_data = read_ascii('temp_bk.dat', labels='energy xmu')
        autobk(x_data.energy, x_data.xmu, group=x_data, rbkg=0.9, kmax=k_max)
        self.k0 = np.asarray(x_data.k)
        chi0 = np.multiply(np.asarray(x_data.chi), self.k0 ** 3)
        self.process(chi0)'''

    def process(self, source):
        self.k, chi = k_range(self.k0, source, self.k_start, self.k_end, False)
        self.chi, self.ft = back_k_space(chi, self.r, self.k.size, self.r_start, self.r_end)
        self.r_cut = self.r[self.r_range]
        self.ft_cut = self.ft[self.r_range]
        self.cross = self.chi * np.transpose([self.ft_cut])

    def r_factor_chi(self, target):
        return np.sum(np.subtract(self.chi, target) ** 2) / self.chi_bottom

    def r_factor_ft(self, target):
        return np.sum(np.subtract(self.ft_cut, target[self.r_range]) ** 2) / self.ft_bottom

    def r_factor_cross(self, target):
        return np.sum(np.subtract(self.cross, target) ** 2) / self.cross_bottom

    def amp_ratio(self, target):
        amp = (self.max[1] / np.abs(target)[self.max[0]]) - 1
        return amp if amp > 0 else 0


class CHI:
    def __init__(self, file_name='cu/chi.dat', k_start=3.0, k_end=12.0, r_start=0.0, r_end=6.0,
                 dE=0.0, s02=1.0, k_size=0):
        self.name = file_name
        self.k_start = k_start
        self.k_end = k_end
        self.dE = dE
        self.s02 = s02
        self.r_start = r_start
        self.r_end = r_end
        self.size = k_size
        self.k0 = np.array([])
        self.k = np.array([])
        self.chi = np.array([])
        self.chi0 = np.array([])
        self.r = np.arange(0, 6 + pi / 102.4, pi / 102.4)

        self.read_file()
        self.ft = np.abs(norm_fft(self.chi, self.r.size))

    def read_file(self):
        with open(self.name, 'r') as target:
            if self.k0.size == 0:
                k = np.array([])
            chi = np.array([])
            if self.name.split('.')[1] == 'dat':
                while True:
                    lines = target.readline()
                    if not lines or not lines.find('phase @#') == -1:
                        break
            if self.k0.size == 0:
                while True:
                    data = target.readline()
                    if not data:
                        break
                    temp = [i for i in data.split()]
                    k = np.append(k, float(temp[0]))
                    chi = np.append(chi, float(temp[1]))
                self.size = k.size if self.size == 0 else self.size
                self.k0 = k.copy()[:self.size]
            else:
                while True:
                    data = target.readline()
                    if not data:
                        break
                    chi = np.append(chi, float(data.split()[1]))
            self.chi0 = np.multiply(chi.copy()[:self.size], self.k0 ** 3)
            self.process(self.s02 * self.chi0)

    def process(self, source):
        chi_shift = deltaE_shift(self.k0, source, dE=self.dE)
        chi_ift = back_k_space(chi_shift, self.size, self.r_start, self.r_end)
        self.k, self.chi = k_range(self.k0, chi_ift, self.k_start, self.k_end, False)

    def r_factor(self, target):
        return np.sum(np.power(np.subtract(self.chi, target), 2)) / np.sum(self.chi ** 2)


'''class TetraFactory(CenteredTetragonalFactory):
    bravais_basis = [[0, 0, 0], [0.5, 0.5, -0.5], [0.5, -0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, -0.5]]
    element_basis = (0, 1, 1, 1, 1)


class OctaFactory(SimpleCubicFactory):
    atomic_basis = np.array([[0., 0., 0.], [.5, 0., 0.], [0., .5, 0.], [0., 0., .5]])
    element_basis = [0, 1, 1, 1]


class CrossFactory(SimpleCubicFactory):
    atomic_basis = np.array([[0., 0., 0.], [.5, 0., 0.], [0., .5, 0.]])
    element_basis = [0, 1, 1]'''


class ATOMS:
    def __init__(self, database='', file=' ', pos=np.array([]), element=np.array([]), spherical=True, random=True,
                 local_range=3.0, surface='', crate_flag=True):
        self.file = file
        self.surface = surface
        self.coordinate_whole = np.array([])
        self.element_whole = np.array([])
        self.cw_temp = np.array([])
        self.ew_temp = np.array([])
        self.distance = np.array([])
        self.coordinate = np.array([])
        self.element = np.array([])
        self.c_temp = np.array([])
        self.e_temp = np.array([])
        self.symbol = np.array([])
        with open(database + r'\table.ini', 'r') as f:
            for i in range(4):
                temp = f.readline().split()[-1]
                if temp == 'Null':
                    break
                self.symbol = np.append(self.symbol, temp)
            self.min_distance = float(f.readline().split()[1])

        self.local_size = element.size
        self.local_range = local_range
        self.center_c = pos[0]
        self.center_e = element[0]
        self.satellite_c = pos[1:]
        self.satellite_e = element[1:]
        self.surface_c = np.array([])
        self.surface_e = np.array([])

        if crate_flag:
            if self.surface == 'TiO2':
                root = self.create_TiO2(random)
                self.deposition(spherical, root)
                self.c_best = self.coordinate_whole.copy()
                self.e_best = self.element_whole.copy()
            elif self.surface == 'Al2O3':
                root = self.create_Al2O3(random)
                self.deposition(spherical, root)
                self.c_best = self.coordinate_whole.copy()
                self.e_best = self.element_whole.copy()
            else:
                self.coordinate = np.vstack((self.center_c, self.satellite_c))
                self.coordinate = np.append(self.center_e, self.satellite_e)
                self.c_best = self.coordinate.copy()
                self.e_best = self.element.copy()
            self.cw_temp = self.coordinate_whole.copy()
            self.ew_temp = self.element_whole.copy()
            self.c_temp = self.coordinate.copy()
            self.e_temp = self.element.copy()


        '''if materials == 'Cu':
            self.min_distance = [2.00, 2.00]
            self.symbol = np.array(['Cu'])
            if crate_flag:
                self.create_cu(layer)
        elif materials == 'Fe':
            self.min_distance = [1.24, 1.24]
            self.symbol = np.array(['Fe'])
            if crate_flag:
                self.create_fe()
        elif materials == 'PtCl4':
            self.min_distance = [2.00, 2.00]
            self.symbol = np.array(['Pt', 'Cl'])
            if crate_flag:
                self.create_ptcl4()
        elif materials == 'PtCl6':
            self.min_distance = [2.00, 2.00]
            self.symbol = np.array(['Pt', 'Cl'])
            if crate_flag:
                self.create_ptcl6()
        elif materials == 'CuOS':
            self.min_distance = [1.5, 1.5, 1.5]
            self.symbol = np.array(['Cu', 'O', 'S'])
            if crate_flag:
                self.create_CuOS(pos)
        elif materials == 'CuON':
            self.min_distance = [1.5, 1.5, 1.5]
            self.symbol = np.array(['Cu', 'O', 'N'])
            if crate_flag:
                self.create_CuON(pos)
        elif materials == 'CuTiO2':
            self.min_distance = [1.5, 1.5, 1.8]
            self.symbol = np.array(['Cu', 'O', 'S', 'Ti'])
            if crate_flag:
                self.create_Cu_TiO2(pos)
            self.adsorption = False
            self.coordinate_best = self.coordinate0.copy()
            self.y1 = -1.87
            self.y2l = -5.28
            self.rx = 4.95
            self.rxl = 5.65
            self.elevation = 43.5 / 180 * pi
            self.eccenter = np.array([-0.33, self.y1])
        elif materials == 'CuAlO':
            self.min_distance = [1.5, 2]
            self.symbol = np.array(['Cu', 'O', 'Al'])
            if crate_flag:
                self.create_CuAlO(pos)
            self.coordinate_best = self.coordinate0.copy()
        if not materials == 'CuTiO2' and not materials == 'CuAlO':
            self.coordinate = self.coordinate0.copy()
            self.coordinate_new = self.coordinate.copy()
            self.coordinate_best = self.coordinate_new.copy()
            self.element0 = self.element.copy()'''

    def create_TiO2(self, ran):
        ele = np.array([])
        coor = np.array([])
        with open(os.getcwd() + r'\TiO2.xyz', 'r') as f:
            f.readline()
            f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                temp = line.split()
                ele = np.append(ele, temp[0])
                coor = np.append(coor, np.array([float(temp[1]), float(temp[2]), float(temp[3])]))
        self.coordinate_whole = np.round(np.reshape(coor, (int(coor.size / 3), 3)), 3)
        self.element_whole = ele.copy()
        self.surface_c = self.coordinate_whole.copy()
        self.surface_e = self.element_whole.copy()
        return randrange(ele.size) if ran else 45

    def create_Al2O3(self, ran):
        ele = np.array([])
        coor = np.array([])
        with open(os.getcwd() + r'\CuAlO_large.xyz', 'r') as f:
            f.readline()
            f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                temp = line.split()
                ele = np.append(ele, temp[0])
                coor = np.append(coor, np.array([float(temp[1]), float(temp[2]), float(temp[3])]))
        self.coordinate_whole = np.round(np.reshape(coor, (int(coor.size / 3), 3)), 3)
        self.element_whole = ele.copy()
        self.surface_c = self.coordinate_whole.copy()
        self.surface_e = self.element_whole.copy()
        adsorb = np.array([], dtype=int)
        for i in range(ele.size):
            if ele[i] == 'O':
                if abs(self.surface_c[i][0]) < 3.5 and abs(self.surface_c[i][1]) < 5 \
                        and abs(self.surface_c[i][2]) > 0.7:
                    adsorb = np.append(adsorb, i)
        return adsorb[randrange(adsorb.size)] if ran else 42

    def deposition(self, sph, root):
        if sph:
            self.center_c = np.array([self.coordinate_whole[root][0] + self.center_c[0]
                                      * sin((180 - self.center_c[2]) / 180 * pi)
                                      * cos((360 - self.center_c[1]) / 180 * pi),
                                      self.coordinate_whole[root][1] + self.center_c[0]
                                      * sin((180 - self.center_c[2]) / 180 * pi)
                                      * sin((360 - self.center_c[1]) / 180 * pi),
                                      self.coordinate_whole[root][2] + self.center_c[0]
                                      * cos((180 - self.center_c[2]) / 180 * pi)])
        else:
            self.center_c += self.coordinate_whole[root]
        self.coordinate_whole = np.vstack((self.coordinate_whole, self.center_c))
        self.element_whole = np.append(self.element_whole, self.center_e)
        for i in range(self.satellite_e.size):
            if sph:
                self.satellite_c[i] = np.array([self.center_c[0] + self.satellite_c[i][0]
                                                * sin((180 - self.satellite_c[i][2]) / 180 * pi)
                                                * cos((360 - self.satellite_c[i][1]) / 180 * pi),
                                                self.center_c[1] + self.satellite_c[i][0]
                                                * sin((180 - self.satellite_c[i][2]) / 180 * pi)
                                                * sin((360 - self.satellite_c[i][1]) / 180 * pi),
                                                self.center_c[2] + self.satellite_c[i][0]
                                                * cos((180 - self.satellite_c[i][2]) / 180 * pi)])
            else:
                self.satellite_c[i] += self.center_c
            self.coordinate_whole = np.vstack((self.coordinate_whole, self.satellite_c[i]))
            self.element_whole = np.append(self.element_whole, self.satellite_e[i])
        distance = get_distance(self.coordinate_whole - self.coordinate_whole[-self.local_size])
        self.coordinate = self.coordinate_whole[-self.local_size:].copy()
        self.element = self.element_whole[-self.local_size:].copy()
        surface_symbol = np.unique(self.surface_e)
        for j in range(surface_symbol.size):
            for i in range(self.surface_e.size):
                if self.surface_e[i] == surface_symbol[j] and distance[i] < self.local_range:
                    self.coordinate = np.vstack((self.coordinate, self.surface_c[i]))
                    self.element = np.append(self.element, self.surface_e[i])
        self.coordinate -= self.coordinate[0]
        self.distance = get_distance(self.coordinate)
        self.t

    '''def create_cu(self, layer):
        from ase.cluster.cubic import FaceCenteredCubic
        surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
        layers = layer
        lc = 3.61
        self.atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)
        coordinate = self.atoms.get_positions()
        self.coordinate0, self.distance, self.element = self.sort(
            coordinate, np.array(self.atoms.get_chemical_symbols()))

    def create_fe(self):
        from ase.cluster.cubic import BodyCenteredCubic
        surfaces = [(1, 0, 0), (0, 1, 0), (1, 1, 1)]
        layers = [2, 2, 3]
        lc = 2.863
        self.atoms = BodyCenteredCubic('Fe', surfaces, layers, latticeconstant=lc)
        coordinate = self.atoms.get_positions()
        self.coordinate0, self.distance, self.element = self.sort(
            coordinate, np.array(self.atoms.get_chemical_symbols()))

    def create_ptcl4(self):
        tetra = TetraFactory()
        self.atoms = tetra(symbol=('Pt', 'Cl'), latticeconstant={'a': 2.65, 'c': 2.65},
                           size=(1, 1, 1), directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]])
        coordinate = self.atoms.get_positions()
        self.coordinate0, self.distance, self.element = self.sort(
            coordinate, np.array(self.atoms.get_chemical_symbols()))

    def create_ptcl4_flat(self):
        cross = CrossFactory()
        self.atoms = cross(('Pt', 'Cl'), surfaces=[(1, 0, 0), (1, 1, 0)],
                           layers=[1, 1], latticeconstant=4.62, center=np.array([0.0, 0.0, 0.0]))
        coordinate = self.atoms.get_positions()
        self.coordinate, self.distance, self.element = self.sort(
            coordinate, np.array(self.atoms.get_chemical_symbols()))

    def create_ptcl6(self):
        octa = OctaFactory()
        self.atoms = octa(('Pt', 'Cl'), surfaces=[(1, 1, 1), (1, 0, 0)],
                          layers=[1, 1], latticeconstant=4.62, center=np.array([0.0, 0.0, 0.0]))
        coordinate = self.atoms.get_positions()
        self.coordinate, self.distance, self.element = self.sort(
            coordinate, np.array(self.atoms.get_chemical_symbols()))

    def create_CuOS(self, init):
        if init.size == 0:
            r1, r2 = 1.84, 2.15
            azi1, azi2 = 45, 225
            ele1, ele2 = 180, 0
        else:
            r1, r2 = init[0][0], init[1][0]
            azi1, azi2 = init[0][1], init[1][1]
            ele1, ele2 = init[0][2], init[1][2]
        azi1, azi2 = randrange(0, 360), randrange(0, 360)
        self.coordinate0 = np.array([[0, 0, 0], [0, 1.3, 1.3], [-1.52, 1.52, 0]])
        self.coordinate0[1][0] = round(r1 * sin(ele1 / 180 * pi) * cos(azi1 / 180 * pi), 3)
        self.coordinate0[1][1] = round(r1 * sin(ele1 / 180 * pi) * sin(azi1 / 180 * pi), 3)
        self.coordinate0[1][2] = round(r1 * cos(ele1 / 180 * pi), 3)
        self.coordinate0[2][0] = round(r2 * sin(ele2 / 180 * pi) * cos(azi2 / 180 * pi), 3)
        self.coordinate0[2][1] = round(r2 * sin(ele2 / 180 * pi) * sin(azi2 / 180 * pi), 3)
        self.coordinate0[2][2] = round(r2 * cos(ele2 / 180 * pi), 3)
        self.distance = np.array([sqrt((_ ** 2).sum()) for _ in self.coordinate0])
        self.element = np.array(['Cu', 'O', 'S'])
        with open(self.file + r'\ini.txt', 'w') as ini:
            ini.write('O(r, azimuth, polar): %.2f %d %d\n' % (r1, azi1, ele1))
            ini.write('S(r, azimuth, polar): %.2f %d %d\n' % (r2, azi2, ele2))
            ini.write(
                'O(x, y, z): %.2f %.2f %.2f\n' % (self.coordinate0[1][0], self.coordinate0[1][1], self.coordinate0[1][2]))
            ini.write(
                'S(x, y, z): %.2f %.2f %.2f\n' % (self.coordinate0[2][0], self.coordinate0[2][1], self.coordinate0[2][2]))

    def create_CuON(self, init):
        if init.size == 0:
            r1, r2 = 1.84, 2.15
            azi1, azi2 = 45, 225
            ele1, ele2 = 180, 0
        else:
            r1, r2 = init[0][0], init[1][0]
            azi1, azi2 = init[0][1], init[1][1]
            ele1, ele2 = init[0][2], init[1][2]
        self.coordinate0 = np.array([[0, 0, 0], [0, 1.3, 1.3], [-1.52, 1.52, 0]])
        self.coordinate0[1][0] = round(r1 * sin(ele1 / 180 * pi) * cos(azi1 / 180 * pi), 3)
        self.coordinate0[1][1] = round(r1 * sin(ele1 / 180 * pi) * sin(azi1 / 180 * pi), 3)
        self.coordinate0[1][2] = round(r1 * cos(ele1 / 180 * pi), 3)
        self.coordinate0[2][0] = round(r2 * sin(ele2 / 180 * pi) * cos(azi2 / 180 * pi), 3)
        self.coordinate0[2][1] = round(r2 * sin(ele2 / 180 * pi) * sin(azi2 / 180 * pi), 3)
        self.coordinate0[2][2] = round(r2 * cos(ele2 / 180 * pi), 3)
        self.distance = np.array([sqrt((_ ** 2).sum()) for _ in self.coordinate0])
        self.element = np.array(['Cu', 'O', 'N'])
        with open(self.file + r'\ini.txt', 'w') as ini:
            ini.write('O(r, azimuth, polar): %.2f %d %d\n' % (r1, azi1, ele1))
            ini.write('N(r, azimuth, polar): %.2f %d %d\n' % (r2, azi2, ele2))
            ini.write(
                'O(x, y, z): %.2f %.2f %.2f\n' % (self.coordinate0[1][0], self.coordinate0[1][1], self.coordinate0[1][2]))
            ini.write(
                'N(x, y, z): %.2f %.2f %.2f\n' % (self.coordinate0[2][0], self.coordinate0[2][1], self.coordinate0[2][2]))

    def create_Cu_TiO2(self, init, init_ele, spherical=True):
        from ase.spacegroup import crystal
        from ase.build import cut, add_adsorbate
        if init.size == 0:
            init = np.array([[1.84, 255, 137]])
            init_ele = np.array(['Cu'])
        a = 4.6
        c = 2.95
        tio2 = crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                       spacegroup=136, cellpar=[a, a, c, 90, 90, 90], size=(4, 4, 4))
        self.atoms = cut(tio2, (0, 0, -1), (1, -1, 0), (1, 1, 0), nlayers=2)
        self.atoms.rotate(45, 'z')
        self.atoms.rotate(90, 'x')
        self.atoms.rotate(90, 'z')
        root = randrange(self.atoms.positions.shape[0])
        if spherical:
            absorbing = np.array([self.atoms.positions[root][0] + init[0][0] * sin((180 - init[0][2]) / 180 * pi)
                                  * cos((360 - init[0][1]) / 180 * pi),
                                  self.atoms.positions[root][1] + init[0][0] * sin((180 - init[0][2]) / 180 * pi)
                                  * sin((360 - init[0][1]) / 180 * pi),
                                  self.atoms.positions[root][2] + init[0][0] * cos((180 - init[0][2]) / 180 * pi)])
        else:
            absorbing = init[0]
        add_adsorbate(self.atoms, init_ele[0], 2, (absorbing[0], absorbing[1]))
        self.atoms.positions[-1][2] = absorbing[2]
        for i in range(1, init_ele.size):
            if spherical:
                satellite = np.array([self.atoms.positions[root][0] + init[i][0] * sin(init[i][2] / 180 * pi)
                                      * cos(init[i][1] / 180 * pi),
                                      self.atoms.positions[root][1] + init[i][0] * sin(init[i][2] / 180 * pi)
                                      * sin(init[i][1] / 180 * pi),
                                      self.atoms.positions[root][2] + init[i][0] * cos(init[i][2] / 180 * pi)])
            else:
                satellite = init[i]
            add_adsorbate(self.atoms, init_ele[i], 2, (satellite[0], satellite[1]))
            self.atoms.positions[-1][2] = satellite[2]
        satellite = init_ele.size
        self.atoms.translate(-self.atoms.positions[49])
        self.coordinate0 = np.round(self.atoms.get_positions(), 3)
        self.coordinate_temp = self.coordinate0.copy()
        self.element0 = np.asarray(self.atoms.get_chemical_symbols())
        distance = np.array([sqrt(((self.coordinate0[_] - self.coordinate0[-satellite]) ** 2).sum())
                             for _ in range(self.coordinate0.shape[0])])
        dangling = []
        for i in range(self.element0.size):
            if self.element0[i] == 'Ti' and (self.coordinate0[i][2] == -13.011 or
                                             abs(self.coordinate0[i][2]) == 6.505 or self.coordinate0[i][2] == 0):
                dangling.append(self.coordinate0[i])
        self.dangling = np.asarray(dangling)

        self.coordinate = self.coordinate0[-satellite:].copy()
        self.element = self.element0[-satellite:].copy()
        for i in range(self.element0.size - satellite):
            if self.element0[i] == 'O' and distance[i] < 3:
                self.coordinate = np.vstack((self.coordinate, self.coordinate0[i]))
                self.element = np.append(self.element, 'O')
        for i in range(self.element0.size - satellite):
            if self.element0[i] == 'Ti' and distance[i] < 3:
                self.coordinate = np.vstack((self.coordinate, self.coordinate0[i]))
                self.element = np.append(self.element, 'Ti')
        self.coordinate -= self.coordinate[0]
        self.distance = np.array([sqrt((_ ** 2).sum()) for _ in self.coordinate])
        self.coordinate_new = self.coordinate.copy()
        self.element_new = self.element.copy()
        with open(self.file + r'\ini.txt', 'w') as ini:
            ini.write('Cu(x,y,z): %.5f %.5f %.5f\n' % (self.coordinate0[-satellite][0], self.coordinate0[-satellite][1],
                                                       self.coordinate0[-satellite][2]))
            for i in range(1, satellite):
                ini.write('%s(x,y,z): %.5f %.5f %.5f\n' % (init_ele[i], self.coordinate0[1 - satellite][0],
                                                           self.coordinate0[1 - satellite][1],
                                                           self.coordinate0[1 - satellite][2]))

    def create_CuAlO(self, init):
        ele = np.array([])
        coor = np.array([])
        with open(self.file.split(self.file.split('\\')[-1])[0] + r'CuAlO_large.xyz', 'r') as f:
            f.readline()
            f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                temp = line.split()
                ele = np.append(ele, temp[0])
                coor = np.append(coor, np.array([float(temp[1]), float(temp[2]), float(temp[3])]))
        if init.size == 0:
            init = np.array([0, 0, 2])
        coor = np.append(coor, init)
        ele = np.append(ele, 'Cu')
        self.coordinate0 = np.round(np.reshape(coor, (int(coor.size / 3), 3)), 3)
        adsorb = np.array([])
        for i in range(ele.size):
            if ele[i] == 'O':
                if abs(self.coordinate0[i][0]) < 3.5 and abs(self.coordinate0[i][1]) < 5 \
                        and abs(self.coordinate0[i][2]) > 0.8:
                    adsorb = np.append(adsorb, self.coordinate0[i])
        adsorb = np.reshape(adsorb, (int(adsorb.size/3), 3))
        self.coordinate0[-1] = adsorb[randrange(adsorb.shape[0])]
        self.coordinate0[-1][2] += init[2]
        self.coordinate_temp = self.coordinate0.copy()
        self.element0 = ele.copy()
        distance = np.array([sqrt(((self.coordinate0[_] - self.coordinate0[-1]) ** 2).sum())
                             for _ in range(self.coordinate0.shape[0])])
        self.coordinate = self.coordinate0[-1:].copy()
        self.element = self.element0[-1:].copy()
        for i in range(self.element0.size - 2):
            if self.element0[i] == 'O' and distance[i] < 3:
                self.coordinate = np.vstack((self.coordinate, self.coordinate0[i]))
                self.element = np.append(self.element, 'O')
        for i in range(self.element0.size - 2):
            if self.element0[i] == 'Al' and distance[i] < 3:
                self.coordinate = np.vstack((self.coordinate, self.coordinate0[i]))
                self.element = np.append(self.element, 'Al')
        self.coordinate -= self.coordinate[0]
        self.distance = np.array([sqrt((_ ** 2).sum()) for _ in self.coordinate])
        self.coordinate_new = self.coordinate.copy()
        self.element_new = self.element.copy()
        with open(self.file + r'\ini.txt', 'w') as ini:
            ini.write('Cu(x,y,z): %.5f %.5f %.5f\n' % (self.coordinate0[-1][0], self.coordinate0[-1][1],
                                                       self.coordinate0[-1][2]))'''

    def write(self, coor, ele):
        distance = np.array([sqrt((_ ** 2).sum()) for _ in coor])
        with open(self.file + r'\feff.inp', 'r+') as file_feff:
            file_feff.seek(0)
            while True:
                lines = file_feff.readline(6)
                if not lines or not lines.find('ATOMS') == -1:
                    break
            file_feff.seek(file_feff.tell())
            file_feff.write('\n')
            file_feff.seek(file_feff.tell())
            file_feff.write('   %.5f     %.5f     %.5f    %d  %s1              %.5f\n'
                            % (coor[0][0], coor[0][1], coor[0][2], 0, ele[0], distance[0]))
            for i in range(1, coor.shape[0]):
                file_feff.seek(file_feff.tell())
                file_feff.write('   %.5f     %.5f     %.5f    %d  %s1              %.5f\n'
                                % (coor[i][0], coor[i][1],
                                   coor[i][2], np.where(self.symbol == ele[i])[0][0]+1, ele[i], distance[i]))
            file_feff.seek(file_feff.tell())
            file_feff.write('END\n')
            file_feff.truncate()

    def read(self, index):
        os.popen('copy "%s" "%s"' % (self.file + r'\result\result%d.txt' % index,
                                     self.file + r'\backup\result%d.txt' % index))
        with open(self.file + r'\result\result%d.txt' % index, 'r') as f:
            coordinate_best = np.array([])
            element_best = np.array([])
            distance_best = np.array([])
            coordinate = np.array([])
            element = np.array([])
            distance = np.array([])
            while True:
                lines = f.readline()
                if not lines.find('Best') == -1:
                    break
            while True:
                data = f.readline()
                if data.isspace() or not data:
                    break
                temp = data.split()
                coordinate_best = np.append(coordinate_best, np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
                element_best = np.append(element_best, temp[4][:-1])
                distance_best = np.append(distance_best, float(temp[5]))
            while True:
                lines = f.readline()
                if not lines.find('final') == -1:
                    break
            while True:
                data = f.readline()
                if data.isspace() or not data:
                    break
                temp = data.split()
                coordinate = np.append(coordinate, np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
                element = np.append(element, temp[4][:-1])
                distance = np.append(distance, float(temp[5]))
        print('data read')
        if not self.surface == '':
            self.coordinate_whole = coordinate.reshape(distance.size, 3)
            self.cw_temp = self.coordinate_whole.copy()
            self.element_whole = element.copy()
            self.surface_c = self.coordinate_whole[:-self.local_size]
            self.surface_e = self.element_whole[:-self.local_size]
            distance = get_distance(self.coordinate_whole - self.coordinate_whole[-self.local_size])
            self.coordinate = self.coordinate_whole[-self.local_size:].copy()
            self.element = self.element_whole[-self.local_size:].copy()
            surface_symbol = np.unique(self.surface_e)
            for j in range(surface_symbol.size):
                for i in range(self.surface_e.size):
                    if self.surface_e[i] == surface_symbol[j] and distance[i] < self.local_range:
                        self.coordinate = np.vstack((self.coordinate, self.surface_c[i]))
                        self.element = np.append(self.element, self.surface_e[i])
            self.coordinate -= self.coordinate[0]
            self.distance = get_distance(self.coordinate)
            self.cw_temp = self.coordinate_whole.copy()
            self.ew_temp = self.element_whole.copy()
            self.center_c = self.coordinate_whole[-self.local_size].copy()
            self.center_e = self.element_whole[-self.local_size].copy()
            self.satellite_c = self.coordinate_whole[-self.local_size - 1:].copy()
            self.satellite_e = self.element_whole[-self.local_size - 1:].copy()
        else:
            self.coordinate = coordinate.reshape(distance.size, 3)
            self.distance = distance.copy()
            self.element = element.copy()
            self.center_c = self.coordinate[0].copy()
            self.center_e = self.element[0].copy()
            self.satellite_c = self.coordinate[1:].copy()
            self.satellite_e = self.element[1:].copy()
        self.c_best = coordinate_best.reshape(distance_best.size, 3)
        self.e_best = element_best.copy()
        self.c_temp = self.coordinate.copy()
        self.e_temp = self.element.copy()
        print('parameter set up')

    def moving(self, target_i=None):
        trials = 50
        while trials > 0:
            self.c_temp = self.coordinate.copy()
            target = randrange(1, self.distance.size) if target_i is None else target_i
            self.c_temp[target][randrange(3)] += round(randrange(-100, 101) * 0.001, 3)
            distance = get_distance(self.c_temp - self.c_temp[target])
            error_flag = False
            for _ in range(self.distance.size):
                if not (_ == target) and not (self.min_distance < distance[_]):
                    error_flag = True
                    break
            if error_flag:
                trials -= 1
                continue
            self.distance[target] = sqrt((self.c_temp[target] ** 2).sum())
            break
        flag = True if trials == 0 else False
        return target, flag

    def moving_spherical(self, target_i=None):
        trials = 20
        while trials > 0:
            self.c_temp = self.coordinate.copy()
            target = randrange(1, self.distance.size) if target_i is None else target_i
            ri = sqrt(((self.c_temp[target]) ** 2).sum())
            if not self.c_temp[target][0] == 0:
                azimuth = atan(self.c_temp[target][1] / self.c_temp[target][0])
                if self.c_temp[target][0] < 0:
                    azimuth += pi
            else:
                azimuth = pi / 2
                if self.c_temp[target][1] < 0:
                    azimuth += pi
            elevation = acos(self.c_temp[target][2] / ri)
            axis = randrange(3)
            if axis == 0:
                ri += round(randrange(-100, 101) * 0.001, 2)
            elif axis == 1:
                azimuth += randrange(-100, 101) / 1800 * pi
            else:
                elevation += randrange(-100, 101) / 1800 * pi
                if elevation > pi:
                    elevation = 2 * pi - elevation
                    azimuth += pi
                elif elevation < 0:
                    elevation = -elevation
                    azimuth += pi
            if azimuth > 2 * pi:
                azimuth -= 2 * pi
            elif azimuth < 0:
                azimuth += 2 * pi
            self.c_temp[target][0] = round(ri * sin(elevation) * cos(azimuth), 3)
            self.c_temp[target][1] = round(ri * sin(elevation) * sin(azimuth), 3)
            self.c_temp[target][2] = round(ri * cos(elevation), 3)
            distance = get_distance(self.c_temp - self.c_temp[target])
            error_flag = False
            for _ in range(distance.size):
                if not (_ == target) and not (self.min_distance < distance[_]):
                    error_flag = True
                    break
            if error_flag:
                trials -= 1
                continue
            self.distance[target] = sqrt((self.c_temp[target] ** 2).sum())
            break
        flag = True if trials == 0 else False
        return target, flag

    def moving_center(self):
        trials = 20
        while trials > 0:
            self.cw_temp = self.coordinate_whole.copy()
            step = round(randrange(-100, 101) * 0.01, 2)
            axis = randrange(3)
            for i in range(self.local_size):
                self.cw_temp[self.surface_e.size+i][axis] += step
            if self.surface == 'TiO2':
                if not (-4.7 < self.cw_temp[-self.local_size][0] < 6.2) or \
                        not (-13.3 < self.cw_temp[-self.local_size][1] < 12) or \
                        not (1.5 < self.cw_temp[-self.local_size][2] < 5):
                    trials -= 1
                    continue
            elif self.surface == 'Al2O3':
                if not (-8 < self.cw_temp[-self.local_size][0] < 5.6) or \
                        not (-10 < self.cw_temp[-self.local_size][1] < 10) or \
                        not (0 < self.cw_temp[-self.local_size][2] < 3):
                    trials -= 1
                    continue
            distance = get_distance(self.cw_temp - self.cw_temp[-self.local_size])
            self.c_temp = self.cw_temp[-self.local_size:].copy()
            self.e_temp = self.ew_temp[-self.local_size:].copy()
            surface_symbol = np.unique(self.surface_e)
            for j in range(surface_symbol.size):
                for i in range(self.surface_e.size):
                    if self.surface_e[i] == surface_symbol[j] and distance[i] < self.local_range:
                        self.c_temp = np.vstack((self.c_temp, self.surface_c[i]))
                        self.e_temp = np.append(self.e_temp, self.surface_e[i])
            self.c_temp -= self.c_temp[0]
            self.distance = get_distance(self.c_temp)
            error_flag = False
            for _ in range(self.local_size, self.distance.size):
                if not self.min_distance < self.distance[_]:
                    error_flag = True
                    break
            if error_flag:
                trials -= 1
                continue
            if self.distance.size > self.local_size and np.min(self.distance[self.local_size:]) > 2.2:
                trials -= 1
                continue
            break
        return True if trials == 0 else False

    '''def tca_filter(self, coor_s):
        coor_rot = np.zeros(2)
        temp = self.dangling - coor_s
        qualify = False
        for j in range(self.dangling.shape[0]):
            if not temp[j][0] == 0:
                azi = atan(temp[j][1] / temp[j][0])
                if temp[j][0] < 0:
                    azi += pi
            else:
                azi = pi/2
            coor_rot[0] = temp[j][0] * cos(azi) + temp[j][1] * sin(azi)
            coor_rot[1] = temp[j][2]
            if self.y2l < coor_rot[1] < self.y1:
                vect = coor_rot - self.eccenter
                if self.rx < sqrt((vect ** 2).sum()) < self.rxl:
                    if not vect[0] == 0 and 0 < -atan(vect[1] / vect[0]) < self.elevation:
                        qualify = True
                        break
        return qualify'''

    def add_atom(self, rx, ry, rz):
        #center = np.array([self.coordinate[:, _].sum() for _ in range(3)]) / self.distance.size
        i_vector = np.array([rx, ry, rz])
        surface = np.array([], dtype=int)
        for i in range(self.coordinate.shape[0]):
            if np.where(np.array(([sqrt(((self.coordinate[i] - self.coordinate[_]) ** 2).sum()
                                        ) for _ in range(self.coordinate.shape[0])])) < 3)[0].size < 12:
                surface = np.append(surface, i)
        off_center = self.coordinate[surface] #- center[0]
        d_off_center = np.array([sqrt((off_center[_] ** 2).sum()) for _ in range(off_center.shape[0])])
        vect_angle = np.array([acos(((np.abs(off_center[_]) * i_vector).sum()) / sqrt((i_vector ** 2).sum())
                                    / d_off_center[_]) for _ in range(d_off_center.size)])
        d_off_center = d_off_center ** 2
        print('distance_index:', np.argsort(d_off_center))
        print('distance:', np.sort(d_off_center))
        print('angle_index:', np.argsort(vect_angle))
        print('angle:', np.sort(vect_angle))
        vect_final = vect_angle * d_off_center
        ad_index = np.argmin(vect_final)#randrange(surface.size)
        print('final_index:', np.argsort(vect_final))
        print('final:', np.sort(vect_final))

        #neighbor = np.where(np.array(([sqrt(((self.coordinate[surface[ad_index]] - self.coordinate[_]) ** 2).sum()
        #                                    ) for _ in range(self.coordinate.shape[0])])) < 3)[0]
        possible_pos = self.coordinate0[1:13]

        # single random add
        '''while True:
            overlap = False
            ad_pos = self.coordinate[surface[ad_index]] + possible_pos[randrange(1, 13)]
            for i in range(self.distance.size):
                if not i == surface[ad_index]:
                    if sqrt(((ad_pos - self.coordinate[i]) ** 2).sum()) < self.min_distance[1]:
                        overlap = True
            if overlap:
                continue
            break
        self.coordinate = np.vstack((self.coordinate, ad_pos))
        self.distance = np.append(self.distance, sqrt((self.coordinate[-1] ** 2).sum()))
        self.element = np.append(self.element, self.element[-1])'''
        # cluster add
        '''for i in range(possible_pos.shape[0]):
            overlap = False
            ad_pos = self.coordinate[surface[ad_index]] + possible_pos[i]
            for j in range(self.distance.size):
                if not j == surface[ad_index]:
                    if sqrt(((ad_pos - self.coordinate[j]) ** 2).sum()) < self.min_distance[1]:
                        overlap = True
            if not overlap:
                self.coordinate = np.vstack((self.coordinate, ad_pos))
                self.distance = np.append(self.distance, sqrt((self.coordinate[-1] ** 2).sum()))
                self.element = np.append(self.element, self.element[-1])'''

        # single vector add
        off_center = possible_pos #- center[0]
        d_off_center = np.array([sqrt((off_center[_] ** 2).sum()) for _ in range(off_center.shape[0])])
        vect_angle = np.array([acos(((np.abs(off_center[_]) * i_vector).sum()) / sqrt((i_vector ** 2).sum())
                                    / d_off_center[_]) for _ in range(d_off_center.size)])
        for i in np.argsort(vect_angle):
            overlap = False
            ad_pos = self.coordinate[surface[ad_index]] + possible_pos[i]
            for i in range(self.distance.size):
                if not i == surface[ad_index]:
                    if sqrt(((ad_pos - self.coordinate[i]) ** 2).sum()) < self.min_distance[1]:
                        overlap = True
            if not overlap:
                self.coordinate = np.vstack((self.coordinate, ad_pos))
                self.distance = np.append(self.distance, sqrt((self.coordinate[-1] ** 2).sum()))
                self.element = np.append(self.element, self.element[-1])
                break

        '''distance_sum = np.zeros(self.coordinate0.shape[0])
        for i in range(1, self.coordinate0.shape[0]):
            ad_pos = self.coordinate[surface[ad_index]] + self.coordinate0[i]
            neighbor = np.where(np.array(([sqrt(((self.coordinate[surface[ad_index]] - self.coordinate[_]) ** 2).sum()
                                                ) for _ in range(self.coordinate.shape[0])])) < 3)[0]
            for j in neighbor:
                distance_sum[i] += sqrt(((ad_pos - self.coordinate[j]) ** 2).sum())'''


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.fftpack import fft

    '''expe = np.array([EXP(r'D:\simulation\Cu202_sum.rex', 3, 9, 1, 2.7),
                     EXP(r'D:\simulation\Cu207_sum.rex', 3, 9, 1, 2.7),
                     EXP(r'D:\simulation\Cu211_sum_d.rex', 3, 9, 1, 2.7)])'''
    expe = EXP(r'D:\simulation\Cu202_sum.rex', 3, 9, 1, 2.7)
    ax = plt.subplot(projection='3d')
    x, y = np.meshgrid(expe.k, expe.r_cut)
    ax.plot_surface(x, y, expe.cross, cmap='magma')
    ax.plot(expe.k, expe.chi*np.max(expe.ft_cut), zs=3, zdir='y', c='k')
    ax.plot(expe.r_cut, expe.ft_cut * np.max(expe.chi), zs=2.5, zdir='x', c='k')

    #chi = CHI(r'D:\tri_angle\chi.dat', 0, 20, 0, 6)
    #plt.plot(chi.r, chi.ft, c='black')
    '''rep = RMC4(0, expe, 0, [16, -3], 1, os.getcwd() + r'\cuos\substrate\test', 'CuTiO2', ini_flag=True)
    plt.plot(rep.table[0].k, rep.table[0].chi, c='red')
    rep.walk(tau=4e-2)
    plt.plot(rep.table[0].k, rep.table[0].chi, c='blue')
    rep.accept()
    plt.plot(rep.table[0].k, rep.table[0].chi, c='green', linestyle='--')'''
    # plt.title(r'invert Fourier Transform Spectrum')
    plt.show()
