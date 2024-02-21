import os
from math import sqrt, pi, atan, cos, sin, acos

import matplotlib.pyplot as plt
import numpy as np

from mrmc_package import cal_angle, scatter, cylinder, line, bar


def read_shell_data():
    distance = np.array([])
    rdf = []
    for index in range(1):
        with open(r'D:\Monte Carlo\cu\temperature-3\result\result%d.txt' % index, 'r') as f:
            while True:
                lines = f.readline()
                if not lines.find('Best') == -1:
                    break
            f.readline()
            while True:
                data = f.readline()
                if data.isspace():
                    break
                distance = np.append(distance, float(data.split()[5]))
    rdf = list(np.unique(np.trunc(distance * 10) / 10, return_counts=True))
    r1_aver = distance[:12].mean()
    c21 = distance[:12].var()
    c31 = ((distance[:12] - r1_aver) ** 3).mean()
    plt.bar(rdf[0], rdf[1], width=0.075, align='center')
    print(r1_aver, c21, c31)
    return r1_aver, c21, c31, rdf


def read_coordinate():
    group = np.array([])
    for index in range(50):
        with open(r'D:\Monte Carlo\cu\total-6 each-3\result\result%d.txt' % index, 'r') as f:
            coordinate = np.array([])
            while True:
                lines = f.readline()
                if not lines.find('final') == -1:
                    break
            while True:
                data = f.readline()
                if data.isspace():
                    break
                temp = data.split()
                coordinate = np.append(coordinate, [float(temp[0]), float(temp[1]), float(temp[2])])
            coordinate = coordinate.reshape((55, 3))
        '''if group.size == 0:
            group = coordinate.copy()
        else:
            group = np.concatenate((group, coordinate))'''
        if len(group) == 0:
            group = [coordinate]
        else:
            group.append(coordinate)
    return group


def collect_pov():
    coordinate = np.array([])
    for index in range(50):
        with open(r'image\image%d.pov' % index, 'r') as f:
            while True:
                lines = f.readline()
                if not lines.find('no cell vertices') == -1:
                    break
            while True:
                data = f.readline()
                if not data.find('cylinder') == -1:
                    break
                coordinate = np.append(coordinate, data)
    f = open(r'../image/image_origin.pov', 'r')
    w = open(r'../image/image_write.pov', 'w')
    while True:
        lines = f.readline()
        if not lines.find('Rbond') == -1:
            w.write(lines.split('=')[0] + '= 0.010;')
        else:
            w.write(lines)
        if not lines.find('no cell vertices') == -1:
            break
    while True:
        data = f.readline()
        if not data.find('cylinder') == -1:
            temp = data
            break
        alpha = data.split(', 0.0,')
        w.write(alpha[0] + ', 0.7,' + alpha[1])
    for a_line in coordinate:
        w.write(a_line)
    rgb = temp.split('color rgb <0.78, 0.50, 0.20>')
    w.write(rgb[0] + 'color rgb <0, 0, 0>' + rgb[1])
    while True:
        data = f.readline()
        if data.find('cylinder') == -1:
            w.write(data)
            break
        rgb = data.split('color rgb <0.78, 0.50, 0.20>')
        w.write(rgb[0] + 'color rgb <0, 0, 0>' + rgb[1])
    f.close()
    w.close()


def read_chi(f_rep, rep_size, pol_size):
    chi_x = []
    chi_y = []
    chi_z = []
    with open(f_rep + r'\log.txt', 'r') as f:
        for i in range(rep_size):
            f.readline()
            for j in range(pol_size):
                f.readline()
                chi_i = np.array([])
                while True:
                    data = f.readline()
                    if data.isspace() or not data:
                        break
                    temp = data.split()
                    chi_i = np.append(chi_i, float(temp[1]))
                if j == 0:
                    chi_x.append(chi_i)
                elif j == 1:
                    chi_y.append(chi_i)
                else:
                    chi_z.append(chi_i)

    if j == 0:
        return np.array([np.asarray(chi_x)])
    elif j == 1:
        return np.array([np.asarray(chi_x), np.asarray(chi_y)])
    else:
        return np.array([np.asarray(chi_x), np.asarray(chi_y), np.asarray(chi_z)])


def read_rep(f_rep, choice, local_size):
    rep_size = len(os.listdir(f_rep + r'\result'))
    coor_rep = []
    print(choice)
    if rep_size <= 0:
        fail = True
    else:
        fail = False
        for i in range(rep_size):
            with open(f_rep + r'\result\result%d.txt' % i, 'r') as f:
                coor = np.array([])
                while True:
                    lines = f.readline()
                    if not lines.find(choice) == -1:
                        break
                while True:
                    data = f.readline()
                    if data.isspace() or not data:
                        break
                    temp = data.split()
                    coor = np.append(coor, np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
            coor_rep.append(coor.reshape((int(coor.size / 3), 3))[-local_size:])
    return np.asarray(coor_rep), fail


def load_substrate(f_rep, local_size):
    if not f_rep.find('.') == -1:
        name = f_rep
        seq = [1, 2, 3, 0]
    else:
        name = f_rep + r'\result\result0.txt'
        seq = [0, 1, 2, 4]
    with open(name, 'r') as f:
        coordinate_substrate = np.array([])
        ele_substrate = np.array([])
        f.readline()
        if name.split('.')[1] == 'xyz':
            f.readline()
        while True:
            temp = f.readline()
            if not temp or temp.isspace():
                break
            temp = temp.split()
            coordinate_substrate = np.append(coordinate_substrate,
                                             np.array([float(temp[seq[0]]), float(temp[seq[1]]), float(temp[seq[2]])]))
            ele_substrate = np.append(ele_substrate, temp[seq[3]][:-1] if seq[3] == 4 else temp[seq[3]])
        coordinate_substrate = coordinate_substrate.reshape(int(coordinate_substrate.size / 3), 3)

    if name.split('.')[1] == 'xyz':
        return coordinate_substrate, ele_substrate
    else:
        surface_size = coordinate_substrate.shape[0] - local_size
        return coordinate_substrate[:surface_size], ele_substrate[:surface_size]


def tca_filter(surface_c, surface_e, shift_c):
    y1 = -1.93
    y2l = -6.45
    rx = 4.95
    rxl = 6.09
    elevation = 65 / 180 * pi
    center = np.array([-0.52, y1])
    filtered = np.array([False] * shift_c.shape[0])
    index = np.array([], dtype=int)
    coor_rot = np.zeros(2)
    '''for i in range(len(coor)):
        coordinate = coor[i][-2].copy()
        ne = np.array([], dtype=int)
        for j in range(ele.size - 2):
            if ele[j] == 'O' and sqrt(((coor[i][-2] - substrate[j]) ** 2).sum()) < 2.4:
                coordinate = np.vstack((coordinate, substrate[j]))
                ne = np.append(ne, j)
        k = np.argmin(
            np.array([sqrt(((coordinate[_] - coor[i][-2]) ** 2).sum()) for _ in range(1, coordinate.shape[0])]))
        coor[i][-2] -= (substrate[ne[k]] - substrate[45])
        coor[i][-1] -= (substrate[ne[k]] - substrate[45])'''
    for i in range(surface_e.size):
        if surface_e[i] == 'Ti' and (surface_c[i][1] == -13.011 or abs(surface_c[i][1]) == 6.505 or surface_c[i][1] == 0):
            index = np.append(index, i)
    for i in range(shift_c.shape[0]):
        temp = shift_c[i][1].copy()
        for j in index:
            # if sqrt(((coor[0][j] - temp[0]) ** 2).sum()) < 7 and sqrt((coor[0][j][1] - temp[0][1]) ** 2) < 3.253:
            temp = np.vstack((temp, surface_c[j]))
            temp = np.vstack((temp, surface_c[j] + np.array([11.8, 0, 0])))
        temp -= temp[0]
        for j in range(1, temp.shape[0]):
            azi = atan(temp[j][1] / temp[j][0]) if not temp[j][0] == 0 else 0
            if temp[j][0] < 0:
                azi += pi
            coor_rot[0] = temp[j][0] * cos(azi) + temp[j][1] * sin(azi)
            coor_rot[1] = temp[j][2]
            if y2l < coor_rot[1] < y1:
                vect = coor_rot - center
                if rx < sqrt((vect ** 2).sum()) < rxl:
                    if 0 < -atan(vect[1] / vect[0]) < elevation:
                        filtered[i] = True
                        break
    return filtered


def select_atom(surface_c, surface_e, local_c, local_e, rpath):
    coordinate = []
    dist = []
    element = []
    adsorb = np.array([], dtype=int)
    for i in range(local_c.shape[0]):
        temp = local_c[i].copy()
        temp_ele = local_e.copy()
        temp_index = np.array([], dtype=int)
        for j in range(surface_e.size):
            if sqrt(((surface_c[j] - temp[0]) ** 2).sum()) < rpath[surface_e[j]]:
                temp = np.vstack((temp, surface_c[j]))
                temp_ele = np.append(temp_ele, surface_e[j])
                temp_index = np.append(temp_index, j)
        temp -= temp[0]
        temp_d = np.array([sqrt((temp[_] ** 2).sum()) for _ in range(local_e.size, temp.shape[0])])
        align = np.append(np.arange(local_e.size), np.argsort(temp_d) + local_e.size)
        adsorb = np.append(adsorb, temp_index[np.argsort(temp_d)[0]])
        coordinate.append(temp[align])
        element.append(temp_ele[align])
        print(temp_ele[align])
        dist.append(np.array([sqrt((coordinate[-1][_] ** 2).sum()) for _ in range(coordinate[-1].shape[0])]))
    return coordinate, dist, element, adsorb


def plot_bondanlge(coor, ele):
    angle = np.array([])
    for i in range(coor.shape[0]):
        if ele[i][1] == 'O':
            angle = np.append(angle, 180 - cal_angle(coor[i][1], np.zeros(3), coor[i][2]))
        else:
            for j in range(2, coor[i].shape[0]):
                angle = np.append(angle, 180 - cal_angle(coor[i][j], np.zeros(3), coor[i][1]))
    return bar(angle[0], angle[1], width=0.01)


def plot_TiO2(coor, ele, graph, local=False):
    index = np.arange(ele.size, dtype=int) if not local else np.array([], dtype=int)
    for i in range(ele.size):
        if ele[i] == 'Ti':
            color = 'grey'
        else:
            if coor[i][2] > 0:
                color = 'purple'
            else:
                color = 'red'
        size = 0.6 if ele[i] == 'O' else 0.4
        if local:
            if abs(coor[i][0] - coor[42][0]) <= 2.95 and abs(coor[i][1] - coor[42][1]) <= 4 and abs(
                    coor[i][2] - coor[42][2]) <= 2.95:
                index = np.append(index, i)
                graph.addItem(scatter(coor[i][0], coor[i][1], coor[i][2], c=color, scale=size))
        else:
            graph.addItem(scatter(coor[i][0], coor[i][1], coor[i][2], c=color, scale=size))

    for i in range(index.size):  # substrate bonds
        for j in range(i):
            if not (ele[index[i]] == ele[index[j]]) and sqrt(((coor[index[i]] - coor[index[j]]) ** 2).sum()) < 3.5:
                graph.addItem(cylinder([coor[index[i]][0], coor[index[j]][0]], [coor[index[i]][1], coor[index[j]][1]],
                                       [coor[index[i]][2], coor[index[j]][2]], c='black', width=0.1))
    graph.addItem(line([-5, 5], [0, 0], [0, 0], c='red', width=3))
    graph.addItem(line([0, 0], [-5, 5], [0, 0], c='green', width=3))
    graph.addItem(line([0, 0], [0, 0], [-5, 5], c='blue', width=3))


def plot_Al2O3(coor, ele, graph):
    for i in range(ele.size):
        color = 'red'
        if ele[i] == 'Al':
            if coor[i][2] < 0:
                color = 'dimgray'
            elif coor[i][2] == 0:
                color = 'darkgray'
            else:
                color = 'lightgray'
        size = 0.4 if ele[i] == 'O' else 0.6
        graph.addItem(scatter(coor[i][0], coor[i][1], coor[i][2], c=color, scale=size))
    '''for i in range(ele.size - 2):  # substrate bonds
        for j in range(i):
            if not (ele[i] == ele[j]) and sqrt(((coor[i] - coor[j]) ** 2).sum()) < 3.5:
                graph.addItem(cylinder([coor[i][0], coor[j][0]], [coor[i][1], coor[j][1]],
                                       [coor[i][2], coor[j][2]], c='black', width=0.1))'''
    graph.addItem(line([-5, 5], [0, 0], [0, 0], c='red', width=3))
    graph.addItem(line([0, 0], [-5, 5], [0, 0], c='green', width=3))
    graph.addItem(line([0, 0], [0, 0], [-5, 5], c='blue', width=3))


def plot_on_substrate(surface_c, surface_e, local_c, rpath, graph):
    item_scatter = np.array([])
    item_cylinder = []
    color = ['brown', 'yellow', 'green', 'orange', 'cyan', 'red', 'purple']
    for i in range(local_c.shape[0]):
        for j in range(local_c.shape[1]):
            item_scatter = np.append(item_scatter, scatter(local_c[i][j][0], local_c[i][j][1],
                                                           local_c[i][j][2], c=color[j], alpha=1, scale=0.3))
            graph.addItem(item_scatter[-1])

        item_cylinder_rep = np.array([])
        dist = np.array([])
        add_list = np.array([], dtype=int)
        for j in range(surface_e.size):
            if sqrt(((local_c[i][0] - surface_c[j]) ** 2).sum()) < rpath[surface_e[j]]:
                dist = np.append(dist, sqrt(((local_c[i][0] - surface_c[j]) ** 2).sum()))
                add_list = np.append(add_list, j)
        nearest = add_list[np.argmin(dist)]
        item_cylinder_rep = np.append(item_cylinder_rep, cylinder([local_c[i][0][0], surface_c[nearest][0]],
                                                                  [local_c[i][0][1], surface_c[nearest][1]],
                                                                  [local_c[i][0][2], surface_c[nearest][2]],
                                                                  c='black', alpha=1, width=0.05))
        #graph.addItem(item_cylinder_rep[-1])
        for j in range(1, local_c.shape[1]):
            item_cylinder_rep = np.append(item_cylinder_rep, cylinder([local_c[i][0][0], local_c[i][j][0]],
                                                                      [local_c[i][0][1], local_c[i][j][1]],
                                                                      [local_c[i][0][2], local_c[i][j][2]],
                                                                      c='black', alpha=1, width=0.05))
            #graph.addItem(item_cylinder_rep[-1])
        item_cylinder.append(item_cylinder_rep)
    # item_scatter = item_scatter.reshape(local_c.shape[0], local_c.shape[1])
    # graph.addItem(cylinder([rep[0][0], rep[1][0]], [rep[0][1], rep[1][1]],
    # [rep[0][2], rep[1][2]], c='red', width=0.05))
    return item_scatter, item_cylinder


def rdf_polarization(coor, dist, ele, sym, select=np.array([])):
    if select.size == 0:
        select = np.arange(len(ele))
    sym = np.append(sym, 'Ol')
    size = sym.size if sym.size <= 3 else 3
    sym_dict = {sym[_]: _ for _ in range(size)}
    stat_d = [np.array([])] * (sym.size + 2)
    stat_a = [np.array([])] * (sym.size + 2)
    stat_e = [np.array([])] * (sym.size + 2)
    label = [np.array([], dtype=int)] * sym.size

    for rep in select:
        for i in range(1, ele[rep].size):
            if ele[rep][i] == 'O' and dist[rep][i] > 2:
                locate = 2
            else:
                locate = sym_dict[ele[rep][i]]
            label[locate] = np.append(label[locate], rep)
            if not coor[rep][i][0] == 0:
                azi = abs(round(atan(coor[rep][i][1] / coor[rep][i][0]) / pi * 180, 0))
                azi = azi if azi < 90 else 180 - azi
                #azi = round(atan(coor[rep][i][1] / coor[rep][i][0]) / pi * 180, 0)
            else:
                azi = 90
            polar = acos(coor[rep][i][2] / dist[rep][i]) / pi * 180
            stat_a[locate] = np.append(stat_a[locate], azi)
            stat_e[locate] = np.append(stat_e[locate], round(polar if polar < 90 else (180 - polar), 0))
            stat_d[locate] = np.append(stat_d[locate], round(dist[rep][i], 2))
            if ele[rep][i] == 'O':
                if dist[rep][i] > 2:
                    stat_a[-1] = np.append(stat_a[-1], azi)
                    stat_e[-1] = np.append(stat_e[-1], round(polar if polar < 90 else (180 - polar), 0))
                    stat_d[-1] = np.append(stat_d[-1], round(dist[rep][i], 2))
                elif dist[rep][i] <= 2:
                    stat_a[-2] = np.append(stat_a[-2], azi)
                    stat_e[-2] = np.append(stat_e[-2], round(polar if polar < 90 else (180 - polar), 0))
                    stat_d[-2] = np.append(stat_d[-2], round(dist[rep][i], 2))
    print('O1:%.3f-%.3f, %.1f-%.1f, %.1f-%.1f' % (stat_d[-2].mean(), stat_d[-2].std(),
                                                  stat_a[-2].mean(), stat_a[-2].std(),
                                                  stat_e[-2].mean(), stat_e[-2].std()))
    print('O2:%.3f-%.3f, %.1f-%.1f, %.1f-%.1f' % (stat_d[-1].mean(), stat_d[-1].std(),
                                                  stat_a[-1].mean(), stat_a[-1].std(),
                                                  stat_e[-1].mean(), stat_e[-1].std()))

        # bond_angle = np.append(bond_angle, 180 - cal_angle(coor[i][1], coor[i][0], coor[i][2]))
    rdf = [stat_d[0], stat_a[0], stat_e[0]]
    for i in range(1, size):
        rdf.append(stat_d[i])
        rdf.append(stat_a[i])
        rdf.append(stat_e[i])
    '''ax = np.array([plt.subplot(2, 3, _ + 1) for _ in range(6)])
    name = ['r-O', 'φ-O', 'θ-O', 'r-S', 'φ-S', 'θ-S']
    urdf = [np.unique(rdf[_], return_counts=True) for _ in range(7)]
    for i in range(6):
        if i == 0 or i == 3:
            w = 0.005
            ax[i].set_xlabel('distance [Å]')
            ax[i].set_ylabel('Frequency')
        else:
            w = 0.5
            ax[i].set_xlabel('angle [°]')
        ax[i].bar(urdf[i][0], urdf[i][1], width=w, align='center', color='k')
        ax[i].set_title(name[i])
    plt.show()'''

    # print('%.1f %.1f %.1f %.1f' % (elevation_o.mean(), sqrt(elevation_o.var()), elevation_s.mean(), sqrt(elevation_s.var())))

    return rdf, label


def plot_rotate(surface_c, surface_e, local_c, local_e, rpath, surface, graph):
    color = ['brown', 'yellow', 'purple', 'green', 'orange', 'cyan', 'red']
    # color = ['blue', 'yellow', 'red', 'red', 'red', 'red', 'red']
    graph.addItem(scatter(0, 0, 0, c=color[0], scale=0.3))
    graph.addItem(line([-3, 3], [0, 0], [0, 0], c='r'))
    graph.addItem(line([0, 0], [-3, 3], [0, 0], c='g'))
    graph.addItem(line([0, 0], [0, 0], [-3, 3], c='b'))
    item_list = []
    print(local_e)
    for i in range(local_c.shape[0]):
        item_rep = np.array([])
        for j in range(1, local_e.size):
            temp = local_c[i][j] - local_c[i][0]
            if j == 1:
                c = color[j]
            if j > 1:
                c = color[2] if np.sqrt((temp**2).sum()) < 2 else color[3]
            item_rep = np.append(item_rep, scatter(temp[0], temp[1], temp[2], c=c, scale=0.3))
            graph.addItem(item_rep[-1])
        if not surface == '':
            for j in range(surface_e.size):
                if sqrt(((surface_c[j] - local_c[i][0]) ** 2).sum()) < rpath[surface_e[j]]:
                    temp = surface_c[j] - local_c[i][0]
                    if surface_e[j] == 'O':
                        c = 'purple' if surface == 'TiO2' and surface_c[j][2] > 0 else 'red'
                    else:
                        if surface_e[j] == 'Al':
                            if surface_c[j][2] < 0:
                                c = 'dimgray'
                            elif surface_c[j][2] == 0:
                                c = 'darkgray'
                            else:
                                c = 'lightgray'
                        else:
                            c = 'silver'
                    item_rep = np.append(item_rep, scatter(temp[0], temp[1], temp[2], c=c, scale=0.3))
                    graph.addItem(item_rep[-1])
        item_list.append(item_rep)
    return item_list
