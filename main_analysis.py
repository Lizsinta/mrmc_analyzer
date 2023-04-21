from mrmc_package import EXP, init_3dplot, TABLE_POL, norm_fft
from mrmc_package.analysis import *
from mrmc_package.ui_analysis import Ui_MainWindow
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from PyQt5.QtGui import QIcon, QFont
from pyqtgraph.Qt import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QMenu
from PyQt5.Qt import QCursor
from PyQt5.QtCore import Qt
import numpy as np
import os
import matplotlib.pyplot as plt
import icon_rc


class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self.setWindowTitle('mRMC')
        self.setWindowIcon(QIcon(':\\mRMC.ico'))
        self.setMinimumSize(0, 0)

        self.subs = gl.GLViewWidget()
        self.subs.setObjectName('subs')
        init_3dplot(self.subs, grid=False, view=30, title='substrate')
        self.g3dLayout.addWidget(self.subs)
        self.rota = gl.GLViewWidget()
        self.rota.setObjectName('rota')
        init_3dplot(self.rota, grid=True, view=10, title='angle distribution')
        self.g3dLayout.addWidget(self.rota)
        self.subs.setContextMenuPolicy(Qt.CustomContextMenu)
        self.subs.customContextMenuRequested.connect(self.save_3d_menu)
        self.rota.setContextMenuPolicy(Qt.CustomContextMenu)
        self.rota.customContextMenuRequested.connect(self.save_3d_menu)


        self.plotx = pg.PlotWidget(background=(255, 255, 255, 255))
        self.plotx.setTitle('χ(k) [001]', color='#000000', size='20pt')
        self.plotx.setLabel('left', 'χ(k)')
        self.plotx.setLabel('bottom', 'k', unit='Å-1')
        self.plotx.addLegend((0, 0))
        self.chixLayout.addWidget(self.plotx)
        self.ploty = pg.PlotWidget(background=(255, 255, 255, 255))
        self.ploty.setTitle('χ(k) [1-10]', color='#000000', size='20pt')
        self.ploty.setLabel('left', 'χ(k)')
        self.ploty.setLabel('bottom', 'k', unit='Å-1')
        self.ploty.addLegend((0, 0))
        self.chiyLayout.addWidget(self.ploty)
        self.plotz = pg.PlotWidget(background=(255, 255, 255, 255))
        self.plotz.setTitle('χ(k) [110]', color='#000000', size='20pt')
        self.plotz.setLabel('left', 'χ(k)')
        self.plotz.setLabel('bottom', 'k', unit='Å-1')
        self.plotz.addLegend((0, 0))
        self.chizLayout.addWidget(self.plotz)
        self.widget_plot = {0: self.plotx, 1: self.ploty, 2: self.plotz}

        fn = QFont()
        fn.setPointSize(15)
        # xy label
        self.do = pg.PlotWidget(background=(255, 255, 255, 255))
        self.Ao = pg.PlotWidget(background=(255, 255, 255, 255))
        self.eo = pg.PlotWidget(background=(255, 255, 255, 255))
        self.ds = pg.PlotWidget(background=(255, 255, 255, 255))
        self.As = pg.PlotWidget(background=(255, 255, 255, 255))
        self.es = pg.PlotWidget(background=(255, 255, 255, 255))
        self.dti = pg.PlotWidget(background=(255, 255, 255, 255))
        self.Ati = pg.PlotWidget(background=(255, 255, 255, 255))
        self.eti = pg.PlotWidget(background=(255, 255, 255, 255))


        self.widget_dict = {0: self.do, 1: self.Ao, 2: self.eo,
                            3: self.ds, 4: self.As, 5: self.es,
                            6: self.dti, 7: self.Ati, 8: self.eti}
        self.layout_dict = {0: self.doLayout, 1: self.aoLayout, 2: self.eoLayout,
                            3: self.dsLayout, 4: self.asLayout, 5: self.esLayout,
                            6: self.dtiLayout, 7: self.atiLayout, 8: self.etiLayout}
        self.action_dict = {0: self.d_oAction, 1: self.a_oAction, 2: self.e_oAction,
                            3: self.d_sAction, 4: self.a_sAction, 5: self.e_sAction, 6: self.b_aAction}
        self.actioname_dict = {'d_oAction': 0, 'a_oAction': 1, 'e_oAction': 2,
                               'd_sAction': 3, 'a_sAction': 4, 'e_sAction': 5, 'b_aAction': 6}
        self.info_dict = {0: self.b1rLine, 1: self.b1aLine, 2: self.b1eLine,
                          3: self.b2rLine, 4: self.b2aLine, 5: self.b2eLine, 6: self.baLine}
        self.rdf_name = ['r-O', 'φ-O', 'θ-O', 'r-S', 'φ-S', 'θ-S', 'r-Ti', 'φ-Ti', 'θ-Ti']

        for i in range(9):
            self.widget_dict[i].getAxis('bottom').setTickFont(fn)
            self.widget_dict[i].getAxis('bottom').setTextPen('black')
            self.widget_dict[i].getAxis('left').setTickFont(fn)
            self.widget_dict[i].getAxis('left').setTextPen('black')
            self.widget_dict[i].setLabel('left', 'count')
            if i == 0 or i == 3 or i == 6:
                self.widget_dict[i].setLabel('bottom', 'distance', unit='Å')
            else:
                self.widget_dict[i].setLabel('bottom', 'angle', unit='°')
            self.layout_dict[i].addWidget(self.widget_dict[i])

        self.do.setTitle('R', color='#000000', size='20pt')
        self.Ao.setTitle('φ', color='#000000', size='20pt')
        self.eo.setTitle('θ', color='#000000', size='20pt')
        self.ds.setTitle('R', color='#000000', size='20pt')
        self.As.setTitle('φ', color='#000000', size='20pt')
        self.es.setTitle('θ', color='#000000', size='20pt')
        self.dti.setTitle('R', color='#000000', size='20pt')
        self.Ati.setTitle('φ', color='#000000', size='20pt')
        self.eti.setTitle('θ', color='#000000', size='20pt')
        #self.ba.setTitle('Bond angle', color='#000000', size='20pt')

        self.region_line = pg.LinearRegionItem()
        self.region_enable = True
        self.regionAction.triggered.connect(self.region_switch)
        self.hmin, self.hmax = 0, 0

        for i in range(7):
            self.action_dict[i].triggered.connect(self.select)
        self.heightAction.triggered.connect(self.select)
        self.expAction.triggered.connect(self.load_exp)
        self.readAction.triggered.connect(self.read)
        self.saveAction.triggered.connect(self.save)
        self.folder = ''
        self.rangeMenu.setEnabled(False)
        self.d_oAction.setChecked(True)
        self.range_target = 0
        self.substrate = True
        self.material = ''

        self.exp_name = np.array([os.getcwd() + r'\source\Cu202_sum.rex',
                                  os.getcwd() + r'\source\Cu207_sum.rex',
                                  os.getcwd() + r'\source\Cu211_sum_d.rex'])

        self.exp = np.array([], dtype=EXP)

    def load_exp(self):
        file_list, _ = QFileDialog.getOpenFileNames(self, 'select experimental data...', r'C:/Monte Carlo/',
                                              'REX Files(*.rex *.ex3)')
        if len(file_list) == 0:
            return
        if len(file_list) > 3:
            self.statusbar.showMessage('too many files', 3000)
            return
        self.statusbar.showMessage('experimental data loaded', 3000)
        self.exp_name = np.asarray(file_list)
        if self.exp_name.size == 1:
            self.plotx.setTitle('χ(k)', color='#000000', size='20pt')
            self.ploty.setTitle('', color='#000000', size='20pt')
            self.plotz.setTitle('', color='#000000', size='20pt')
        elif self.exp_name.size == 2:
            self.plotx.setTitle('χ(k)-s', color='#000000', size='20pt')
            self.ploty.setTitle('χ(k)-p', color='#000000', size='20pt')
            self.plotz.setTitle('', color='#000000', size='20pt')
        elif self.exp_name.size == 3:
            self.plotx.setTitle('χ(k) [001]', color='#000000', size='20pt')
            self.ploty.setTitle('χ(k) [1-10]', color='#000000', size='20pt')
            self.plotz.setTitle('χ(k) [110]', color='#000000', size='20pt')
        print(self.exp_name)

    def read(self):
        address = self.folder.split(self.folder.split('/')[-1])[0] if not self.folder == '' \
            else r'D/CuAlO'#r'C:/Monte Carlo/cuos/substrate'
        folder = QFileDialog.getExistingDirectory(self, 'select folder...', address)
        if folder == '':
            return
        if not self.folder == '':
            self.subs.clear()
            self.rota.clear()
            init_3dplot(self.rota, view=10, title='angle distribution')
            self.plotx.clear()
            self.ploty.clear()
            self.plotz.clear()
            for i in range(7):
                self.widget_dict[i].clear()
        self.folder = folder
        with open(self.folder + r'\info.txt', 'r') as f:
            self.material = f.readline().split()[1]
        coordinate_cu, distance_cu, self.element_cu, flag = read_rep_substrate(self.folder, self.material, filter=False)
        if not flag:
            self.statusbar.showMessage('No result file found', 3000)
            self.folder = ''
            return
        self.substrate = True if coordinate_cu[0].shape[0] > 3 else False
        absorb = 1 if self.material == 'CuAlO' else 2
        self.coordinate_substrate = coordinate_cu[0][:-absorb].copy() if self.substrate else self.load_substrate()


        if self.material == 'CuAlO':
            plot_Al2O3(self.coordinate_substrate, self.element_cu, self.subs)
        else:
            plot_TiO2(self.coordinate_substrate, self.element_cu, self.subs, local=False)
        '''self.subs.addItem(line([self.coordinate_substrate[45][0]-5, self.coordinate_substrate[45][0]+5], [self.coordinate_substrate[45][1], self.coordinate_substrate[45][1]], [self.coordinate_substrate[45][2], self.coordinate_substrate[45][2]], c='red', width=3))
        self.subs.addItem(line([self.coordinate_substrate[45][0], self.coordinate_substrate[45][0]], [self.coordinate_substrate[45][1]-5, self.coordinate_substrate[45][1]+5], [self.coordinate_substrate[45][2], self.coordinate_substrate[45][2]], c='green', width=3))
        self.subs.addItem(line([self.coordinate_substrate[45][0], self.coordinate_substrate[45][0]], [self.coordinate_substrate[45][1], self.coordinate_substrate[45][1]], [self.coordinate_substrate[45][2]-5, self.coordinate_substrate[45][2]+5], c='blue', width=3))
        self.subs.addItem(scatter(self.coordinate_substrate[45][0], self.coordinate_substrate[45][1], self.coordinate_substrate[45][2], c='purple', scale=0.6))'''
        self.rep = len(coordinate_cu)
        self.amountLine.setText('%d / %d' % (self.rep, len(os.listdir(self.folder + r'/result'))))
        if self.substrate:
            self.coordinate_np = np.array([coordinate_cu[0][-absorb:].copy()])
            for i in range(1, self.rep):
                self.coordinate_np = np.vstack((self.coordinate_np, np.array([coordinate_cu[i]])))
        else:
            self.coordinate_np = np.asarray(coordinate_cu)
            for i in range(self.rep):
                self.coordinate_np[i] = self.coordinate_np[i] - self.coordinate_np[i][2] + self.coordinate_substrate[45]
            self.coordinate_np = self.coordinate_np[:, :2, :]

        self.coordinate_select, self.distance_select, self.element_select = select_atom(self.coordinate_substrate,
                                                                   self.coordinate_np, self.element_cu, nearest=True)
        print(self.coordinate_np.shape)



        surface_symbol = np.array(['O', 'Ti'])
        dist = np.array([])
        for rep in range(self.rep):
            distance = np.array([np.sqrt(((_ - self.coordinate_np[rep][0]) ** 2).sum()) for _ in self.coordinate_substrate])
            temp_c = self.coordinate_np[rep].copy()
            temp_e = self.element_cu[-2:].copy()
            for j in range(surface_symbol.size):
                for i in range(self.element_cu.size - 2):
                    if self.element_cu[i] == surface_symbol[j] and distance[i] < 2.7:
                        temp_c = np.vstack((temp_c, self.coordinate_substrate[i]))
                        temp_e = np.append(temp_e, self.element_cu[i])
                        if self.element_cu[i] == 'O':
                            dist = np.append(dist, sqrt(((self.coordinate_substrate[i] - temp_c[0]) ** 2).sum()))
            if rep == 0:
                self.coor_cal = [temp_c - temp_c[0]]
                self.element_cal = [temp_e]
            else:
                self.coor_cal.append(temp_c - temp_c[0])
                self.element_cal.append(temp_e)

        if not self.material == 'CuAlO':
            self.rdf = rdf_polarization(self.coordinate_select, self.distance_select, self.element_select)
            urdf = [np.unique(self.rdf[_], return_counts=True) for _ in range(9)]
            self.bar_spher = [bar(urdf[_][0], urdf[_][1], width=0.5) for _ in range(9)]
            for i in range(9):
                if i < 6:
                    self.info_dict[i].setText('%.2f-%.2f' % (self.rdf[i].mean(), self.rdf[i].std()))
                self.widget_dict[i].addItem(self.bar_spher[i])
                self.widget_dict[i].setLabel('left', 'count')
                if i == 0 or i == 3 or i == 6:
                    self.widget_dict[i].setLabel('bottom', 'distance', unit='Å')
                    self.bar_spher[i].setOpts(width=0.005)
                    self.widget_dict[i].setXRange(np.min(urdf[i][0]) - 0.01, np.max(urdf[i][0]) + 0.01)
                else:
                    self.widget_dict[i].setLabel('bottom', 'angle', unit='°')
            self.do.plotItem.enableAutoRange(axis='x', enable=False)
            with open(self.folder + r'/info.txt', 'r') as info:
                for _ in range(4):
                    info.readline()
                k = info.readline().split()[1:]
                r = info.readline().split()[1:]
                print(k, r)
                for _ in range(5):
                    info.readline()
                b1name = 'Cu-O'#info.readline().split()[0]
                b2name = 'Cu-S'#info.readline().split()[0]
                print(b1name, b2name)
                self.do.setTitle('R (%s)' % b1name, color='#000000', size='20pt')
                self.Ao.setTitle('φ (%s)' % b1name, color='#000000', size='20pt')
                self.eo.setTitle('θ (%s)' % b1name, color='#000000', size='20pt')
                self.ds.setTitle('R (%s)' % b2name, color='#000000', size='20pt')
                self.As.setTitle('φ (%s)' % b2name, color='#000000', size='20pt')
                self.es.setTitle('θ (%s)' % b2name, color='#000000', size='20pt')
                self.b1rLabel.setText('Bond distance (%s)' % b1name)
                self.b1aLabel.setText('Azimuth angle (%s)' % b1name)
                self.b1eLabel.setText('Polar angle (%s)' % b1name)
                self.b2rLabel.setText('Bond distance (%s)' % b2name)
                self.b2aLabel.setText('Azimuth angle (%s)' % b2name)
                self.b2eLabel.setText('Polar angle (%s)' % b2name)
                self.exp = np.array([EXP(self.exp_name[i], float(k[0]), float(k[1]), float(r[0]), float(r[1]))
                                     for i in range(self.exp_name.size)])

        table = np.zeros((3, self.rep), dtype=TABLE_POL)
        for pol in range(3):
            for i in range(self.rep):
                table[pol][i] = TABLE_POL(self.exp[pol].k_start, self.exp[pol].k_end, self.exp[pol].r_start,
                                          self.exp[pol].r_end, 0, [13, -3, 13], 1, self.exp[pol].k0,
                                          self.coor_cal[i].copy(), self.element_cal[i].copy(),
                                          r'J:\Monte Carlo\cutio2', pol, ms_en=False)
        chi_sum = np.zeros((3, self.exp[0].k.size))
        ax = np.array([plt.subplot(3, 1, _+1) for _ in range(3)])
        direct = ['[001]', '[1-10]', '[110]']
        for pol in range(3):
            for index in range(self.rep):
                chi_sum[pol] += table[pol][index].chi
            chi_sum[pol] /= self.rep
            print(self.exp[pol].r_factor(chi_sum[pol]))
            ax[pol].plot(self.exp[pol].r, self.exp[pol].ft, c='k', label='expe %s'%direct[pol])
            ax[pol].plot(self.exp[pol].r, np.abs(norm_fft(chi_sum[pol], self.exp[pol].r.size)), c='r', label='average %s'%direct[pol])
            ax[pol].legend(loc='lower right')
            ax[pol].set_ylabel('χ %s'%direct[pol])
        ax[2].set_xlabel(r'$k[Å^{-1}]$')
        plt.show()

        self.item_r = plot_rotate(self.coordinate_substrate, self.coordinate_np, self.element_cu, self.material,
                                  self.rota, nearest=True, color_assign='' if self.substrate else 'purple')

        self.item_s, self.item_c = plot_on_substrate(self.coordinate_np, self.coordinate_substrate, self.element_cu,
                                                     self.material, self.subs)

        self.line_exp = np.array([line(x=self.exp[_].k, y=self.exp[_].chi, c='black', width=3)
                                  for _ in range(self.exp.size)])
        self.chi = read_chi(self.folder, self.rep, self.exp.size)
        self.line_chi = np.zeros((self.exp.size, self.rep), dtype=pg.PlotDataItem)
        for i in range(self.exp.size):
            self.widget_plot[i].addItem(self.line_exp[i])
            for j in range(self.rep):
                self.line_chi[i][j] = line(x=self.exp[0].k, y=table[i][j].chi, c='red', alpha=0.3)
                self.widget_plot[i].addItem(self.line_chi[i][j])

        #ax = np.array([plt.subplot(1, 3, _ + 1) for _ in range(3)])
        #direct = ['[001]', '[1-10]', '[110]']

        self.line_aver = np.zeros(3, dtype=pg.PlotDataItem)
        for i in range(self.exp.size):
            chi = np.sum(self.chi[i], axis=0) / self.rep
            #ax[i].plot(self.exp[i].k, self.exp[i].chi, c='k', label='expe %s'%direct[i])
            #ax[i].plot(self.exp[i].k, chi, c='r', label='average %s' % direct[i])
            #ax[i].legend(loc='lower right')
            #ax[i].set_ylabel('χ %s'%direct[i])
            self.line_aver[i] = line(x=self.exp[i].k, y=chi_sum[i], c='blue', width=3)
            self.line_aver[i].opts['name'] = '%.3f' % self.exp[i].r_factor(chi_sum[i])
            self.widget_plot[i].addItem(self.line_aver[i])
            rfactor = np.array([self.exp[i].r_factor(self.chi[i][_]) for _ in range(self.rep)])
            rdf = np.unique(np.round(rfactor, 1), return_counts=True)
            #ax[i].bar(rdf[0], rdf[1], width=0.05, align='center', color='k')
            #ax[i].set_xlabel('R-factor')
            #ax[i].set_title(direct[i])
        #ax[2].set_xlabel(r'$k[Å^{-1}]$')
        #ax[0].set_ylabel(r'frequency')
        #plt.show()

        self.hmin = np.min(self.rdf[self.range_target]) - 0.01
        self.hmax = np.max(self.rdf[self.range_target]) + 0.01
        self.region_line.setRegion([self.hmin, self.hmax])
        self.do.addItem(self.region_line)

        self.region_line.sigRegionChanged.connect(lambda: self.update(False))
        self.rangeMenu.setEnabled(True)

    def save(self):
        address = self.folder if not self.folder == '' else os.getcwd()
        name = QFileDialog.getSaveFileName(self, 'select path...', address + '/distribution.png', 'png(*.png)')
        if name[0] == '':
            return
        self.widget_dict[self.range_target].removeItem(self.region_line)
        fig = self.g2dWidget.grab()
        self.widget_dict[self.range_target].addItem(self.region_line)
        fig.save(name[0], 'PNG')

    def select(self):
        target = self.sender()
        self.hmin = np.min(self.rdf[self.range_target]) - 0.01
        self.hmax = np.max(self.rdf[self.range_target]) + 0.01
        self.update(True)
        self.region_line.setRegion([self.hmin, self.hmax])
        self.widget_dict[self.range_target].removeItem(self.region_line)
        self.action_dict[self.range_target].setChecked(False)
        self.widget_dict[self.range_target].plotItem.enableAutoRange(axis='x')

        self.range_target = self.actioname_dict[target.objectName()]

        self.hmin = np.min(self.rdf[self.range_target]) - 0.01
        self.hmax = np.max(self.rdf[self.range_target]) + 0.01
        self.region_line.setRegion([self.hmin, self.hmax])

        self.widget_dict[self.range_target].addItem(self.region_line)
        self.widget_dict[self.range_target].plotItem.disableAutoRange(axis='x')

        self.region_enable = True
        self.regionAction.setChecked(True)

    def update(self, manual=False):
        data = self.rdf[self.range_target]
        if manual:
            hmin, hmax = self.region_line.getRegion()
            index0 = np.where((hmin < data) & (data < hmax), True, False)
            index = np.array([True] * self.rep)
        else:
            index0 = np.where((self.hmin < data) & (data < self.hmax), True, False)
            self.hmin, self.hmax = self.region_line.getRegion()
            index = np.where((self.hmin < data) & (data < self.hmax), True, False)
        for i in range(self.rep):
            if index0[i] and not index[i]:
                self.subs.removeItem(self.item_s[i * 2])
                self.subs.removeItem(self.item_s[i * 2 + 1])
                for j in range(self.item_c[i].size):
                    self.subs.removeItem(self.item_c[i][j])
                for j in range(self.item_r[i].size):
                    self.rota.removeItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].removeItem(self.line_chi[pol][i])
            elif not index0[i] and index[i]:
                self.subs.addItem(self.item_s[i * 2])
                self.subs.addItem(self.item_s[i * 2 + 1])
                for j in range(self.item_c[i].size):
                    self.subs.addItem(self.item_c[i][j])
                for j in range(self.item_r[i].size):
                    self.rota.addItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].addItem(self.line_chi[pol][i])

        self.amountLine.setText('%d / %d' % (np.where(index == True)[0].size,
                                             len(os.listdir(self.folder + r'/result'))))
        rdf = rdf_polarization(self.coordinate_select[index], self.distance_select[index])
        urdf = [np.unique(rdf[_], return_counts=True) for _ in range(7)]
        urdf[self.range_target] = np.unique(self.rdf[self.range_target], return_counts=True)
        for i in range(7):
            self.info_dict[i].setText('%.2f-%.2f' % (rdf[i].mean(), rdf[i].std()))
            self.bar_spher[i].setOpts(x=urdf[i][0], height=urdf[i][1])
        if not urdf[0][0].size == 0:
            self.do.setXRange(np.min(urdf[0][0] - 0.01), np.max(urdf[0][0] + 0.01))
        if not urdf[3][0].size == 0:
            self.ds.setXRange(np.min(urdf[3][0] - 0.01), np.max(urdf[3][0] + 0.01))
        if index.any():
            count = np.where(index == True)[0].size
            for i in range(self.exp.size):
                chi = np.sum(self.chi[i][index], axis=0) / count
                self.line_aver[i].setData(x=self.exp[0].k, y=chi)
                self.widget_plot[i].plotItem.legend.items[0][1].setText('%.3f' % self.exp[i].r_factor(chi))

    def load_substrate(self):
        with open(os.getcwd() + r'\tio2_coor.txt', 'r') as f:
            coordinate_substrate = np.array([])
            ele_substrate = np.array([])
            while True:
                temp = f.readline()
                if not temp or temp.isspace():
                    break
                temp = temp.split()
                coordinate_substrate = np.append(coordinate_substrate,
                                                 np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
                ele_substrate = np.append(ele_substrate, temp[4][:-1])
            coordinate_substrate = coordinate_substrate.reshape(int(coordinate_substrate.size/3), 3)
            self.element_cu = np.append(ele_substrate, np.array(['Cu', 'S']))
        return coordinate_substrate

    def save_3d_menu(self, pos):
        target = self.sender()
        menu = QMenu()
        action = menu.addAction('save')
        action.triggered.connect(lambda: self.save_3d_action(target))
        menu.exec_(QCursor.pos())

    def save_3d_action(self, target):
        address = self.folder if not self.folder == '' else os.getcwd()
        name = '/substrate.jpg' if target.objectName() == 'subs' else '/central.jpg'
        size = (1920, 1080) if target.objectName() == 'subs' else (1080, 1080)
        name = QFileDialog.getSaveFileName(self, 'select path...', address + name, 'jpg(*.jpg)')
        if name[0] == '':
            return
        image = target.renderToArray(size)
        pg.makeQImage(image.transpose(1, 0, 2)).save(name[0])
        # image.save(name[0], 'JPG', 100)

    def region_switch(self):
        if self.regionAction.isChecked() and not self.region_enable:
            self.widget_dict[self.range_target].addItem(self.region_line)
            self.region_enable = True
        elif not self.regionAction.isChecked() and self.region_enable:
            self.widget_dict[self.range_target].removeItem(self.region_line)
            self.region_enable = False



if __name__ == '__main__':
    from sys import flags
    from sys import argv, exit

    app = QApplication(argv)
    main = MainWindow()
    main.show()
    exit(app.exec_())

    '''app = QApplication([])
    w = gl.GLViewWidget()
    w.show()
    init_3dplot(w, grid=False)
    with open(r'C:\Monte Carlo\cuos\substrate\from_0_point3\feff.inp', 'r') as f:
        coor = np.array([])
        ele = np.array([])
        while True:
            lines = f.readline(6)
            if not lines or not lines.find('ATOMS') == -1:
                break
        while True:
            data = f.readline()
            if not data.find('END') == -1:
                break
            temp = data.split()
            coor = np.append(coor, np.array([float(temp[0]), float(temp[1]), float(temp[2])]))
            ele = np.append(ele, temp[4][:-1])
        coor = [coor.reshape((int(coor.size / 3), 3))]
        ele = [ele]
    plot_with_substrate(coor, ele, w)'''

    if (flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QApplication.instance().exec_()
