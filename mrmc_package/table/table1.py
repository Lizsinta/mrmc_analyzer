from mrmc_package import intp1D
import numpy as np
from math import sqrt, acos, pi


def used_file(file_name):
    with open(file_name + r'\chi.dat', 'r') as f_chi:
        file_index = np.array([], dtype=int)
        while True:
            line = f_chi.readline()
            if not line.find('file') == -1:
                break
        while True:
            temp = f_chi.readline()
            if not temp.find('paths used') == -1:
                break
            file_index = np.append(file_index, int(temp.split()[1]))
        return file_index


def read_path(file_name, index):
    with open(file_name + r'\paths.dat', 'r') as f_path:
        path_r = []
        path_angle = []
        f_path.readline()
        f_path.readline()
        i = 0
        while True:
            line = f_path.readline()
            if not line:
                break
            if not line.find('index') == -1:
                if index[i] == int(line.split()[0]):
                    f_path.readline()
                    r = ()
                    angle = ()
                    i += 1
                    for step in range(int(line.split()[1])):
                        data = f_path.readline().split()
                        r += (float(data[6]),)
                        angle += (round(float(data[7]), 0),)
                    path_r.append(r)
                    path_angle.append(angle)
                    del r
                    del angle
        return path_r, path_angle


def cal_angle(coor1, coor2, coor3):
    vect1 = coor2 - coor1
    vect2 = coor3 - coor2
    transvection = round(np.sum(vect1 * vect2), 6)
    module = round(sqrt(np.sum(vect1**2) * np.sum(vect2**2)), 6)
    return round(acos(transvection / module) / pi * 180, 2)

class FEFF:
    def __init__(self, file_name='feff0001.dat', k=np.array([]), reduce=1.0):
        temp = file_name.split('\\')[-1].split('.')[0]
        self.index = int(temp[-4]) if len(temp) == 8 else int(temp.split('feff')[1])
        #self.phase_c = np.array([])
        #self.phase_s = np.array([])
        self.phase = np.array([])
        self.amp = np.array([])
        self.lamb = np.array([])
        self.atom = np.array([])
        self.chi = np.array([])
        self.distance = 0
        #self.realp = np.array([])
        self.k = k

        self.read(file_name, reduce)

    def read(self, file_name, reduce):
        with open(file_name, 'r') as target:
            k0 = np.array([])
            amp = np.array([])
            #phase_c = np.array([])
            #phase_s = np.array([])
            phase = np.array([])
            lamb = np.array([])
            #realp = np.array([])
            while True:
                lines = target.readline()
                if not lines or not lines.find('------') == -1:
                    break
            info = target.readline().split()
            leg = int(info[0]) - 1
            self.distance = float(info[2])
            target.readline()
            target.readline()
            for i in range(leg):
                temp = target.readline().split()
                self.atom = np.append(self.atom, np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
            self.atom = self.atom.reshape(leg, 3)
            target.readline()
            while True:
                data = target.readline()
                if not data:
                    break
                temp = [i for i in data.split()]
                k0 = np.append(k0, float(temp[0]))
                amp = np.append(amp, float(temp[2]) * float(temp[4]))
                phase = np.append(phase, float(temp[1]) + float(temp[3]))
                lamb = np.append(lamb, float(temp[5]))
                #realp = np.append(realp, float(temp[6]))

            if self.k.size == 0:
                self.k = k0.copy()
            #self.phase_c = intp1D(k0, self.k, phase_c)
            #self.phase_s = intp1D(k0, self.k, phase_s)
            self.phase = intp1D(k0, self.k, phase)
            self.amp = np.multiply(intp1D(k0, self.k, amp), self.k ** 2) * reduce
            self.lamb = np.reciprocal(intp1D(k0, self.k, lamb))
            self.chi = np.multiply(np.multiply(np.sin(2 * self.distance * self.k + self.phase),
                                               np.exp(-2 * self.distance * self.lamb)), self.amp) / (self.distance ** 2)
            #self.chi = np.multiply(np.multiply(np.sin(2 * self.distance * self.k + intp1D(k0, self.k, phase)),
            #                                   np.exp(-2 * self.distance * np.reciprocal(intp1D(k0, self.k, lamb)))),
            #                       np.multiply(intp1D(k0, self.k, amp), self.k ** 2) * reduce) / (self.distance ** 2)
            #self.realp = intp1D(k0, self.k, realp)

