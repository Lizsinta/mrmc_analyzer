from mrmc_package import EXP, init_3dplot, TABLE_POL, norm_fft, cal_angle
from mrmc_package.analysis import *
from mrmc_package.ui_analysis import Ui_MainWindow
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from PyQt5.QtGui import QIcon, QFont, QCursor
from pyqtgraph.Qt import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QMenu, QMessageBox
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


        fn = QFont()
        fn.setPointSize(12)
        fn.setBold(True)
        labelStyle = {'color': '#000000', 'font-size': '25t'}
        self.plotx = pg.PlotWidget(background=(255, 255, 255, 255))
        self.plotx.setTitle('χ(k) [001]', color='#000000', size='20pt')
        #self.plotx.setLabel('left', 'χ(k)')
        self.plotx.getAxis('left').setLabel(text='χ(k)', unit='10nm$^-1$', **labelStyle)
        self.plotx.getAxis('bottom').setLabel(text='k', unit='Å-1', **labelStyle)
        self.plotx.addLegend((0, 0), labelTextSize='12pt', labelTextColor='k')
        self.chixLayout.addWidget(self.plotx)
        self.ploty = pg.PlotWidget(background=(255, 255, 255, 255))
        self.ploty.setTitle('χ(k) [1-10]', color='#000000', size='20pt')
        self.ploty.getAxis('left').setLabel(text='χ(k)', **labelStyle)
        self.ploty.getAxis('bottom').setLabel(text='k', unit='Å-1', **labelStyle)
        self.ploty.addLegend((0, 0), labelTextSize='12pt', labelTextColor='k')
        self.chiyLayout.addWidget(self.ploty)
        self.plotz = pg.PlotWidget(background=(255, 255, 255, 255))
        self.plotz.setTitle('χ(k) [110]', color='#000000', size='20pt')
        self.plotz.getAxis('left').setLabel(text='χ(k)', **labelStyle)
        self.plotz.getAxis('bottom').setLabel(text='k', unit='Å-1', **labelStyle)
        self.plotz.addLegend((0, 0), labelTextSize='12pt', labelTextColor='k')
        self.chizLayout.addWidget(self.plotz)
        self.widget_plot = {0: self.plotx, 1: self.ploty, 2: self.plotz}
        self.rota.setSizePolicy(self.plotz.sizePolicy())
        self.subs.setSizePolicy(self.plotz.sizePolicy())
        for i in range(3):
            self.widget_plot[i].getAxis('bottom').setTickFont(fn)
            self.widget_plot[i].getAxis('bottom').setTextPen('black')
            self.widget_plot[i].getAxis('left').setTickFont(fn)
            self.widget_plot[i].getAxis('left').setTextPen('black')


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
                            3: self.d_sAction, 4: self.a_sAction, 5: self.e_sAction,
                            6: self.d_tiAction, 7: self.a_tiAction, 8: self.e_tiAction}
        self.actioname_dict = {'d_oAction': 0, 'a_oAction': 1, 'e_oAction': 2,
                               'd_sAction': 3, 'a_sAction': 4, 'e_sAction': 5,
                               'd_tiAction': 6, 'a_tiAction': 7, 'e_tiction': 8}
        self.info_dict = {0: self.b1rLine, 1: self.b1aLine, 2: self.b1eLine,
                          3: self.b2rLine, 4: self.b2aLine, 5: self.b2eLine}
        self.rdf_name = ['r-O', 'φ-O', 'θ-O', 'r-S', 'φ-S', 'θ-S', 'r-Ti', 'φ-Ti', 'θ-Ti']

        for i in range(9):
            self.widget_dict[i].getAxis('bottom').setTickFont(fn)
            self.widget_dict[i].getAxis('bottom').setTextPen('black')
            self.widget_dict[i].getAxis('left').setTickFont(fn)
            self.widget_dict[i].getAxis('left').setTextPen('black')
            self.widget_dict[i].getAxis('left').setLabel(text='count', **labelStyle)
            if i == 0 or i == 3 or i == 6:
                self.widget_dict[i].getAxis('bottom').setLabel(text='distance', unit='Å', **labelStyle)
            else:
                self.widget_dict[i].getAxis('bottom').setLabel(text='angle', unit='°', **labelStyle)
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

        for i in range(9):
            self.action_dict[i].triggered.connect(self.select)
        self.heightAction.triggered.connect(self.select)
        self.readAction.triggered.connect(self.read)
        self.saveAction.triggered.connect(self.save)
        self.outputAction.triggered.connect(self.output)
        self.plot_rdfAction.triggered.connect(self.save_matplotlib_action)
        self.plot_chiAction.triggered.connect(self.save_matplotlib_action)
        self.plot_baAction.triggered.connect(self.save_matplotlib_action)
        self.rfSlider.valueChanged.connect(self.rf_filter)
        self.folder = ''
        self.rangeMenu.setEnabled(False)
        self.d_oAction.setChecked(True)
        self.range_target = 0
        self.substrate = True
        self.material = ''

        '''self.exp_name = np.array([os.getcwd() + r'\source\Cu202_sum.rex',
                                  os.getcwd() + r'\source\Cu207_sum.rex',
                                  os.getcwd() + r'\source\Cu211_sum_d.rex'])'''

        self.exp = np.array([], dtype=EXP)

        self.surface = ''
        self.rpath = {}
        self.fit_space = 'k'
        self.local_c, self.local_e = np.array([]), np.array([])
        self.surface_c, self.surface_e = np.array([]), np.array([])
        self.select_c, self.select_e = np.array([]), np.array([])
        self.rep = self.local_c.shape[0]
        self.current_rep = np.arange(self.rep)

    def read(self):
        address = self.folder.split(self.folder.split('/')[-1])[0] if not self.folder == '' \
            else r'D:/CuAlO'#r'C:/Monte Carlo/cuos/substrate'
        folder = QFileDialog.getExistingDirectory(self, 'select folder...', address)
        if folder == '':
            return

        self.clear()
        self.folder = folder
        k, r, dw, dE = self.load_info()
        weight, s0, ms, f_material, surface_file, surface_range = self.load_inp(k, r)
        self.load_rep()
        '''for i in range(self.rep):
            if self.local_c[i][1][0] < 0:
                self.local_c[i][1][0] = -self.local_c[i][1][0]
            if self.local_c[i][2][0] > 0:
                self.local_c[i][2][0] = -self.local_c[i][2][0]'''
        self.load_atom()
        chi_sum = self.cal_spectrum(f_material, dw, dE, s0, ms, weight)
        self.cal_rdf()

        '''shift_c = self.local_c.copy()
        if not self.surface == '':
            for rep in range(self.rep):
                shift_c[rep] = self.local_c[rep] - self.surface_c[self.adsorb[rep]] + self.surface_c[41]
        else:
            for rep in range(self.rep):
                self.local_c[rep] = self.local_c[rep][[0, 2, 1], :]'''

        self.item_r = plot_rotate(self.surface_c, self.surface_e, self.local_c, self.local_e, self.rpath, self.surface,
                                  self.rota)

        '''if self.surface == '':
            for rep in range(self.rep):
                shift_c[rep] = self.local_c[rep] - self.local_c[rep][2] + self.surface_c[41]
            shift_c = shift_c[:, :2, :]
            self.local_e = self.local_e[[0, 2]]
            self.local_size = self.local_e.size'''

        if not self.surface == '':
            self.item_s, self.item_c = plot_on_substrate(self.surface_c, self.surface_e, self.local_c, self.rpath, self.subs)

        self.set_plot(chi_sum)
        self.setWindowTitle('mRMC (%s)' % self.folder)

    def load_exp(self, file_list, k, r):
        exp_name = np.asarray(file_list)
        self.exp = np.array([EXP(exp_name[i], k[0], k[1], r[0], r[1]) for i in range(exp_name.size)])

        if exp_name.size == 1:
            self.plotx.setTitle('χ(k)', color='#000000', size='20pt')
            self.ploty.setTitle('', color='#000000', size='20pt')
            self.plotz.setTitle('', color='#000000', size='20pt')
        elif exp_name.size == 2:
            self.plotx.setTitle('χ(k)-s', color='#000000', size='20pt')
            self.ploty.setTitle('χ(k)-p', color='#000000', size='20pt')
            self.plotz.setTitle('', color='#000000', size='20pt')
        elif exp_name.size == 3:
            self.plotx.setTitle('χ(k) [001]', color='#000000', size='20pt')
            self.ploty.setTitle('χ(k) [1-10]', color='#000000', size='20pt')
            self.plotz.setTitle('χ(k) [110]', color='#000000', size='20pt')

    def clear(self):
        if not self.folder == '':
            self.subs.clear()
            self.rota.clear()
            init_3dplot(self.rota, view=10, title='angle distribution')
            self.plotx.clear()
            self.ploty.clear()
            self.plotz.clear()
            for i in range(9):
                self.widget_dict[i].clear()

    def load_info(self):
        with open(self.folder + r'\info.txt', 'r') as f:
            self.local_e = np.array([f.readline().split()[1]])
            self.local_e = np.append(self.local_e, f.readline().split()[1:])
            self.local_size = self.local_e.size
            temp = f.readline().split()
            self.surface = temp[1] if len(temp) > 1 else ''
            sym = np.unique(self.local_e[1:])
            if not self.surface == '':
                if self.surface == 'TiO2':
                    sym = np.append(sym, np.array(['O', 'Ti']))
                elif self.surface == 'Al2O3':
                    sym = np.append(sym, np.array(['O', 'Al']))
            temp = f.readline().split()
            if not temp[0].find('dw') == -1 or not temp[0].find('SIG2') == -1:
                temp = temp[1:]
                if temp[0].find('=') == -1:
                    dw = {i: float(temp[0]) for i in sym}
                else:
                    dw = {_[0]: float(_[1]) for _ in np.char.split(np.asarray(temp, dtype=str), '=')}
                temp = f.readline().split()[1:]
            else:
                dw = {i: 0.0 for i in sym}
                temp = temp[1:]
            if temp[0].find('=') == -1:
                dE = {i:float(j) for i, j in zip(sym, temp)}
            else:
                dE = {_[0]:float(_[1]) for _ in np.char.split(np.asarray(temp, dtype=str), '=')}
            k = [float(_) for _ in f.readline().split()[1:]]
            r = [float(_) for _ in f.readline().split()[1:]]
            print(k, r)
        return k, r, dw, dE

    def load_inp(self, k, r):
        with (open(self.folder + r'/mrmc.inp', 'r') as f):
            f.readline()
            exp_path = np.array([])
            for i in range(3):
                temp = f.readline()
                if not temp or temp.find('exp') == -1:
                    self.sig_warning.emit('exp formal error')
                    return False
                temp = temp.split(temp.split(':')[0] + ':')[1].strip()
                if not temp or len(temp.split()) == 0:
                    break
                exp_path = np.append(exp_path, temp.split('\n')[0])
            self.load_exp(exp_path, k, r)
            while True:
                lines = f.readline()
                if not lines.find('weight') == -1:
                    break
            weight = int(lines.split(':')[1])

            while True:
                lines = f.readline()
                if not lines.find('surface:') == -1:
                    break
            temp = lines.split('surface:')[1].strip()
            self.surface = temp.split(' ')[0]
            surface_file = temp[len(self.surface):].strip()
            if self.surface == 'TiO2':
                self.symbol = np.array(['O', 'Ti'])
            elif self.surface == 'Al2O3':
                self.symbol = np.array(['O', 'Al'])

            while True:
                lines = f.readline()
                if not lines.find('surface_range') == -1:
                    break
            temp = lines.split(':')[1].strip()
            surface_range = np.array([float(_) for _ in temp.split()]) if not len(temp) == 0 else np.array([])

            while True:
                lines = f.readline()
                if not lines.find('S0') == -1:
                    break
            s0 = float(lines.split(':')[1])
            for i in range(4):
                f.readline()
            if not len(self.surface) == 0:
                temp = f.readline().split(':')[1]
                if temp.find('=') == -1:
                    rpath = np.array([float(_) for _ in temp.strip().split()])
                    self.rpath ={self.symbol[i]:rpath[i] for i in range(rpath.size)}
                else:
                    self.rpath = {_[0]:float(_[1]) for _ in np.char.split(np.asarray(temp.strip().split(), dtype=str), '=')}
            else:
                f.readline()
                self.rpath = {}
            temp = f.readline().split(':')[1].strip()
            if temp == 'True' or temp == 'true' or temp == '1':
                ms = True
            elif temp == 'False' or temp == 'false' or temp == '0':
                ms = False
            fspace = f.readline().split(':')
            if not fspace[0].find('fitting_space') == -1:
                temp = fspace[1].strip()
                if temp == 'K' or temp == 'k':
                    self.fit_space = 'k'
                elif temp == 'R' or temp == 'r':
                    self.fit_space = 'r'
                elif temp == 'X' or temp == 'x':
                    self.fit_space = 'x'

            while True:
                lines = f.readline()
                if not lines.find('material_folder') == -1:
                    break
            f_material = lines.split(lines.split(':')[0] + ':')[1].split('\n')[0].strip()
            return weight, s0, ms, f_material, surface_file, surface_range

    def load_rep(self):
        self.local_c, flag = read_rep(self.folder, self.choice_window(), self.local_size)
        if flag:
            self.statusbar.showMessage('No result file found', 3000)
            self.folder = ''
            return
        self.rep = self.local_c.shape[0]
        self.current_rep = np.arange(self.rep)
        # filter = tca_filter(self.surface_c, self.surface_e, self.local_c)
        # self.local_c = self.local_c[filter, :, :]
        self.amountLine.setText('%d / %d' % (self.rep, len(os.listdir(self.folder + r'/result'))))

    def load_atom(self):
        if not self.surface == '':
            self.surface_c, self.surface_e = load_substrate(self.folder, self.local_size)
            if self.surface == 'Al2O3':
                plot_Al2O3(self.surface_c, self.surface_e, self.subs)
            elif self.surface == 'TiO2':
                plot_TiO2(self.surface_c, self.surface_e, self.subs, local=False)
            self.select_c, self.select_d, self.select_e, self.adsorb = select_atom(
                self.surface_c, self.surface_e, self.local_c, self.local_e, self.rpath)
            self.symbol = np.unique(np.append(self.local_e[1:], self.symbol))
        else:
            self.symbol = np.unique(self.local_e[1:])#np.unique(np.append(self.local_e[1:], ['O', 'Ti']))
            self.select_c = self.local_c.copy()
            self.select_d = np.array([])
            for i in range(self.local_c.shape[0]):
                self.select_d = np.append(self.select_d, np.array([sqrt((self.local_c[i][_] ** 2).sum())
                                                                   for _ in range(self.local_c.shape[1])]))
            self.select_d = np.reshape(self.select_d, (self.rep, self.local_size))
            self.select_e = np.tile(self.local_e, (self.rep, 1))
            '''for rep in range(self.rep):
                self.local_c[rep] = self.local_c[rep][[0, 2, 1], :] - self.local_c[rep][1] + self.surface_c[41]
            self.local_c = self.local_c[:, :2, :]
            self.local_e = self.local_e[[0, 2]]
            self.local_size = self.local_e.size
            self.select_c, self.select_d, self.select_e, self.adsorb = select_atom(
                self.surface_c, self.surface_e, self.local_c, self.local_e, self.rpath)'''
        '''short = np.array([])
        long = np.array([])
        sa = np.array([])
        la = np.array([])
        se = np.array([])
        le = np.array([])
        si = []
        li = []
        sulf = []
        sud = np.array([])
        sua = np.array([])
        sue = np.array([])
        for i in range(len(self.select_c)):
            if self.select_d[i][1] < 2.3 and self.select_d[i][2] < 2 and self.select_d[i][3] < 2:
                sulf.append(self.select_c[i][1])
                if not self.select_c[i][1][0] == 0:
                    azi1 = round(atan(self.select_c[i][1][1] / self.select_c[i][1][0]) / pi * 180, 0)
                    azi1 = azi1 if azi1 < 90 else 180 - azi1
                else:
                    azi1 = 90
                sud = np.append(sud, self.select_d[i][1])
                sua = np.append(sua, azi1)
                sue = np.append(sue, round(acos(self.select_c[i][1][2] / self.select_d[i][1]) / pi * 180, 0))
                if not self.select_c[i][2][0] == 0:
                    azi2 = round(atan(self.select_c[i][2][1] / self.select_c[i][2][0]) / pi * 180, 0)
                    azi2 = azi2 if azi2 < 90 else 180 - azi2
                else:
                    azi2 = 90
                if not self.select_c[i][3][0] == 0:
                    azi3 = round(atan(self.select_c[i][3][1] / self.select_c[i][3][0]) / pi * 180, 0)
                    azi3 = azi3 if azi3 < 90 else 180 - azi3
                else:
                    azi3 = 90
                if self.select_d[i][2] > self.select_d[i][3]:
                    si.append(self.select_c[i][3])
                    li.append(self.select_c[i][2])
                    short = np.append(short, self.select_d[i][3])
                    long = np.append(long, self.select_d[i][2])
                    sa = np.append(sa, azi3)
                    la = np.append(la, azi2)
                    se = np.append(se, round(acos(self.select_c[i][3][2] / self.select_d[i][3]) / pi * 180, 0))
                    le = np.append(le, round(acos(self.select_c[i][2][2] / self.select_d[i][2]) / pi * 180, 0))
                else:
                    si.append(self.select_c[i][2])
                    li.append(self.select_c[i][3])
                    short = np.append(short, self.select_d[i][2])
                    long = np.append(long, self.select_d[i][3])
                    sa = np.append(sa, azi2)
                    la = np.append(la, azi3)
                    se = np.append(se, round(acos(self.select_c[i][2][2] / self.select_d[i][3]) / pi * 180, 0))
                    le = np.append(le, round(acos(self.select_c[i][3][2] / self.select_d[i][2]) / pi * 180, 0))
        si = np.asarray(si)
        li = np.asarray(li)
        sulf = np.asarray(sulf)
        print(short.size, short.mean(), short.var(), long.mean(), long.var(), sud.mean(), sud.var())
        print(short.size, sa.mean(), sa.var(), la.mean(), la.var(), sua.mean(), sua.var())
        print(short.size, se.mean(), se.var(), le.mean(), le.var(), sue.mean(), sue.var())
        print(si.mean(0))
        print(li.mean(0))
        print(sulf.mean(0))
        dshort = np.unique(np.round(short, 2), return_counts=True)
        dlong = np.unique(np.round(long, 2), return_counts=True)
        dsa = np.unique(np.round(sa, 0), return_counts=True)
        dla = np.unique(np.round(la, 0), return_counts=True)
        dse = np.unique(np.round(se, 0), return_counts=True)
        dle = np.unique(np.round(le, 0), return_counts=True)
        dsud = np.unique(np.round(sud, 2), return_counts=True)
        dsua = np.unique(np.round(sua, 0), return_counts=True)
        dsue = np.unique(np.round(sue, 0), return_counts=True)
        plt.subplot(331)
        plt.bar(dshort[0], dshort[1], width=0.005)
        plt.title('%.3f-%.5f' % (short.mean(), short.std()))
        plt.subplot(332)
        plt.bar(dsa[0], dsa[1], width=0.5)
        plt.title('%.1f-%.3f' % (sa.mean(), sa.std()))
        plt.subplot(333)
        plt.bar(dse[0], dse[1], width=0.5)
        plt.title('%.1f-%.3f' % (se.mean(), se.std()))
        plt.subplot(334)
        plt.bar(dlong[0], dlong[1], width=0.005)
        plt.title('%.3f-%.5f' % (long.mean(), long.std()))
        plt.subplot(335)
        plt.bar(dla[0], dla[1], width=0.5)
        plt.title('%.1f-%.3f' % (la.mean(), la.std()))
        plt.subplot(336)
        plt.bar(dle[0], dle[1], width=0.5)
        plt.title('%.1f-%.3f' % (le.mean(), le.std()))
        plt.subplot(337)
        plt.bar(dsud[0], dsud[1], width=0.005)
        plt.title('%.3f-%.5f' % (sud.mean(), sud.std()))
        plt.subplot(338)
        plt.bar(dsua[0], dsua[1], width=0.5)
        plt.title('%.1f-%.3f' % (sua.mean(), sua.std()))
        plt.subplot(339)
        plt.bar(dsue[0], dsue[1], width=0.5)
        plt.title('%.1f-%.3f' % (sue.mean(), sue.std()))
        plt.show()'''

    def cal_spectrum(self, folder_material, dw, dE, s0, ms, weight):
        self.table = np.zeros((self.exp.size, self.rep), dtype=TABLE_POL)
        pol = np.arange(3) if self.exp.size == 3 else (np.arange(2) + 2)
        for i in range(self.exp.size):
            for j in range(self.rep):
                self.table[i][j] = TABLE_POL(self.exp[i].k_start, self.exp[i].k_end, self.exp[i].r_start,
                                             self.exp[i].r_end, dw, dE, s0, self.exp[i].k0, self.select_c[j],
                                             self.select_e[j], folder_material, pol[i], ms_en=ms, weight=weight)
        chi_sum = np.zeros((self.exp.size, self.exp[0].k.size))
        self.chi = np.zeros((self.exp.size, self.rep, self.exp[0].k.size))
        for pol in range(self.exp.size):
            for i in range(self.rep):
                self.chi[pol][i] = self.table[pol][i].chi
                if self.fit_space == 'k':
                    chi_sum[pol] += self.table[pol][i].chi
                elif self.fit_space == 'r':
                    chi_sum[pol] += self.table[pol][i].ft
                else:
                    chi_sum[pol] += self.table[pol][i].chi * np.transpose([self.table[pol][i].ft[self.exp[pol].r_range]])
            chi_sum[pol] /= self.rep
            print(self.exp[pol].r_factor_chi(chi_sum[pol]))
        return chi_sum

    def cal_rdf(self):
        self.rdf, self.rdf_label = rdf_polarization(self.select_c, self.select_d, self.select_e, self.symbol)
        graph_size = 0
        for i in self.rdf:
            if i.size > 0:
                graph_size += 1
        urdf = [np.unique(self.rdf[_], return_counts=True) for _ in range(graph_size)]
        self.bar_spher = [bar(urdf[_][0], urdf[_][1], width=0.5) for _ in range(graph_size)]

        '''ba = np.array([180 - cal_angle(self.select_c[_, 1, :], np.zeros(3), self.select_c[_, 2, :]) for _ in range(self.rep)])
        if self.select_c.shape[1] > 3:
            ba = np.append(ba, np.array([
                180 - cal_angle(self.select_c[_, 1, :], np.zeros(3), self.select_c[_, 3, :]) for _ in range(self.rep)]))
        ba = ba[np.where(ba > 130)[0]]
        urdf_ba = np.unique(np.round(ba, 0), return_counts=True)

        self.g3dLayout.removeWidget(self.subs)
        self.ba = pg.PlotWidget(background=(255, 255, 255, 255))
        fn = QFont()
        fn.setPointSize(20)
        self.ba.getAxis('bottom').setTickFont(fn)
        self.ba.getAxis('bottom').setTextPen('black')
        self.ba.getAxis('left').setTickFont(fn)
        self.ba.getAxis('left').setTextPen('black')
        self.g3dLayout.addWidget(self.ba)
        self.g3dLayout.setStretch(0, 1)
        self.g3dLayout.setStretch(1, 1)
        self.ba_bar = bar(urdf_ba[0], urdf_ba[1], width=0.5)
        self.ba.addItem(self.ba_bar)
        self.ba.setLabel('left', 'count')
        self.ba.setLabel('bottom', 'angle', unit='Å')
        self.ba.setTitle('S-Cu-O angle', color='#000000', size='30pt')
        print('%.1f-%.1f' % (ba.mean(), ba.std()))'''


        '''third_ele = True if graph_size > 6 else False
        for i in range(6, 9):
            self.action_dict[i].setEnabled(third_ele)'''
        for i in range(graph_size):
            if i < 6:
                self.info_dict[i].setText('%.3f-%.3f' % (self.rdf[i].mean(), self.rdf[i].std()))
            self.widget_dict[i].addItem(self.bar_spher[i])
            self.widget_dict[i].setLabel('left', 'count')
            if i == 0 or i == 3 or i == 6:
                self.widget_dict[i].setLabel('bottom', 'distance', unit='Å')
                self.bar_spher[i].setOpts(width=0.005)
                self.widget_dict[i].setXRange(np.min(urdf[i][0]) - 0.01, np.max(urdf[i][0]) + 0.01)
            else:
                self.widget_dict[i].setLabel('bottom', 'angle', unit='°')
        symbol = np.append(self.symbol, 'Ol')
        print(symbol)
        for i in range(graph_size // 3):
            self.widget_dict[i * 3].setTitle('R (%s-%s)' % (self.local_e[0], symbol[i]),
                                             color='#000000', size='20pt')
            self.widget_dict[i * 3 + 1].setTitle('φ (%s-%s)' % (self.local_e[0], symbol[i]),
                                                 color='#000000', size='20pt')
            self.widget_dict[i * 3 + 2].setTitle('θ (%s-%s)' % (self.local_e[0], symbol[i]),
                                                 color='#000000', size='20pt')
        self.do.plotItem.enableAutoRange(axis='x', enable=False)

    def set_plot(self, chi_sum):
        self.line_exp = np.array([line(x=self.exp[_].k, y=self.exp[_].chi, c='black', width=3)
                                  for _ in range(self.exp.size)])
        self.line_chi = np.zeros((self.exp.size, self.rep), dtype=pg.PlotDataItem)
        for i in range(self.exp.size):
            self.widget_plot[i].addItem(self.line_exp[i])
            for j in range(self.rep):
                self.line_chi[i][j] = line(x=self.exp[0].k, y=self.table[i][j].chi, c='red', alpha=0.3)
                self.widget_plot[i].addItem(self.line_chi[i][j])

        self.r_factor = np.zeros(self.rep)
        for i in range(self.rep):
            for j in range(self.exp.size):
                self.r_factor[i] += self.exp[j].r_factor_chi(self.table[j][i].chi)
        self.rfSlider.setValue(99)
        self.rf_range = np.linspace(np.min(self.r_factor), np.max(self.r_factor), 100)
        self.rfLabel.setText('r-factor filter (<=%.3f)' % self.rf_range[-1])
        self.rf_rep = np.arange(self.rep)

        self.line_aver = np.zeros(self.exp.size, dtype=pg.PlotDataItem)
        self.r_factor_temp = np.zeros(self.exp.size)
        for i in range(self.exp.size):
            self.line_aver[i] = line(x=self.exp[i].k, y=chi_sum[i], c='blue', width=3)
            self.r_factor_temp[i] = self.exp[i].r_factor_chi(chi_sum[i])
            self.line_aver[i].opts['name'] = '%.3f' % self.r_factor_temp[i]
            self.widget_plot[i].addItem(self.line_aver[i])

        self.hmin = np.min(self.rdf[self.range_target]) - 0.01
        self.hmax = np.max(self.rdf[self.range_target]) + 0.01
        self.region_line.setRegion([self.hmin, self.hmax])
        self.do.addItem(self.region_line)
        # self.flush(np.arange(self.rep)[filter], np.arange(self.rep)[np.invert(filter)])

        self.region_line.sigRegionChanged.connect(lambda: self.update(False))
        self.rangeMenu.setEnabled(True)
        self.plot_rdfAction.setEnabled(True)
        self.plot_chiAction.setEnabled(True)

    def save(self):
        address = self.folder if not self.folder == '' else os.getcwd()
        name = QFileDialog.getSaveFileName(self, 'select path...', address + r'\distribution.png', 'png(*.png)')
        if name[0] == '':
            return
        self.widget_dict[self.range_target].removeItem(self.region_line)
        fig = self.g2dWidget.grab()
        self.widget_dict[self.range_target].addItem(self.region_line)
        fig.save(name[0], 'PNG')
        fig = self.chiWidget.grab()
        fig.save(name[0].split('.')[0] + '_chi.png', 'PNG')
        with open(self.folder + r'/rdf.dat', 'w') as f:
            f.write('Simulation folder: %s\n' % self.folder)
            f.write('Selected replicas:\n')
            for i in range(self.current_rep.size):
                f.write('%d ' % self.current_rep[i])
            f.write('\n')
            f.write('R-factor: ')
            for i in range(self.exp.size):
                f.write('%.3f ' % self.r_factor_temp[i])
            f.write('\n')
            for i in range(self.symbol.size if self.symbol.size <= 2 else 2):
                temp = self.info_dict[i * 3].text().split('-')
                f.write('Bond distance[%s-%s]: %s±%s\n' % (self.local_e[0], self.symbol[i], temp[0], temp[1]))
                temp = self.info_dict[i * 3 + 1].text().split('-')
                f.write('Azimuth angle[%s-%s]: %s±%s\n' % (self.local_e[0], self.symbol[i], temp[0], temp[1]))
                temp = self.info_dict[i * 3 + 2].text().split('-')
                f.write('Polar angle[%s-%s]: %s±%s\n' % (self.local_e[0], self.symbol[i], temp[0], temp[1]))

    def output(self):
        address = self.folder + r'/output'
        index = int(1)
        while True:
            if os.path.isdir(address):
                name_new = address + str(index)
                if os.path.isdir(name_new):
                    index += 1
                else:
                    address = name_new
                    break
            else:
                break
        os.makedirs(address)
        os.makedirs(address + r'\result')
        os.popen('copy "%s" "%s"' % (self.folder + r'\mrmc.inp', address + r'\mrmc.inp'))
        os.popen('copy "%s" "%s"' % (self.folder + r'\info.txt', address + r'\info.txt'))
        os.popen('copy "%s" "%s"' % (self.folder + r'\model.dat', address + r'\model.dat'))
        for i in range(self.current_rep.size):
            os.popen('copy "%s" "%s"' % (self.folder + r'\result\result%d.txt' % self.current_rep[i],
                                         address + r'\result\result%d.txt' % i))
        while True:
            try:
                with open(address + r'\mrmc.inp', 'r+') as f:
                    while True:
                        posi = f.tell()
                        line = f.readline()
                        if not line.find('replicas_size') == -1:
                            f.seek(posi)
                            f.write(('replicas_size:%d' % self.current_rep.size).ljust(len(line) - 1) + '\n')
                            break
            except PermissionError:
                continue
            else:
                break
        chi_sum = np.zeros((self.exp.size, self.exp[0].k.size))
        ft_sum = np.zeros((self.exp.size, self.exp[0].r.size))
        for pol in range(self.exp.size):
            for i in self.current_rep:
                self.chi[pol][i] = self.table[pol][i].chi
                chi_sum[pol] += self.table[pol][i].chi
                ft_sum[pol] += self.table[pol][i].ft
            chi_sum[pol] /= self.current_rep.size
            ft_sum[pol] /= self.current_rep.size
        cross_sum = np.array(
                [chi_sum[_] * np.transpose([ft_sum[_][self.exp[_].r_range]]) for _ in range(self.exp.size)])
        while True:
            try:
                with open(address + r'/info.txt', 'r+') as f:
                    while True:
                        posi = f.tell()
                        line = f.readline()
                        if not line.find('best_R_factor') == -1:
                            f.seek(posi)
                            if self.fit_space == 'k':
                                best_r = np.array([self.exp[_].r_factor_chi(chi_sum[_]) for _ in range(self.exp.size)])
                            elif self.fit_space == 'r':
                                best_r = np.array([self.exp[_].r_factor_ft(ft_sum[_]) for _ in range(self.exp.size)])
                            else:
                                best_r = np.array([self.exp[_].r_factor_cross(cross_sum[_]) for _ in range(self.exp.size)])
                            data = 'best_R_factor([001], [1-10], [110]): %.6f (' % best_r.sum()
                            for i in best_r:
                                data += '%.6f ' % i
                            data = (data[:-1] + ')').ljust(len(line) - 1) + '\n'
                            f.write(data)
                            break
            except PermissionError:
                continue
            else:
                break
        with open(self.folder + r'\model.dat', 'r') as fi:
            name = ['initial', 'final', 'best']
            with open(address + r'\model.dat', 'w') as fo:
                for _ in range(3):
                    if not self.surface == '':
                        while True:
                            line = fi.readline()
                            if not line.find('[%s surface model]' % name[_]) == -1:
                                break
                            if not line.find('[best model]') == -1:
                                break
                        fo.write(line)
                        temp_c = np.array([])
                        for i in range(self.surface_e.size):
                            fo.write(fi.readline())
                        for i in range(self.local_e.size):
                            for rep in range(self.rep):
                                line = fi.readline()
                                if np.where(self.current_rep == rep)[0].size > 0:
                                    fo.write(line)
                                    temp = line.split()
                                    temp_c = np.append(temp_c, np.array([float(temp[1]), float(temp[2]), float(3)]))
                        temp_c = temp_c.reshape((self.local_e.size, self.current_rep.size, 3))
                        temp_c = temp_c.transpose(1, 0, 2)
                        fo.write('\n')
                    else:
                        temp_c = self.local_c
                    fo.write('[%s local model]\n' % name[_])
                    fo.write('%s %.6f %.6f %.6f\n' % (self.local_e[0], 0, 0, 0))
                    for rep in range(temp_c.shape[0]):
                        for i in range(1, temp_c.shape[1]):
                            temp = temp_c[rep][i] - temp_c[rep][0]
                            fo.write('%s %.6f %.6f %.6f\n' % (self.local_e[i], temp[0], temp[1], temp[2]))
                        for i in range(self.surface_e.size):
                            temp = self.surface_c[i] - temp_c[rep][0]
                            if sqrt((temp**2).sum()) < self.rpath[self.surface_e[i]]:
                                fo.write('%s %.6f %.6f %.6f\n' % (self.surface_e[i], temp[0], temp[1], temp[2]))
                    fo.write('\n')
            QMessageBox.information(self, 'Output succeed', 'Output folder:\n%s' % address)

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
        graph_size = self.symbol.size if self.symbol.size <= 3 else 3
        data = self.rdf[self.range_target]
        label = self.rdf_label[self.range_target//graph_size]

        if manual:
            select_remove = np.array([])
            select_add = np.arange(self.rep)[self.rf_rep]
        else:
            self.hmin, self.hmax = self.region_line.getRegion()
            index1 = np.where((self.hmin < data) & (data < self.hmax), False, True)
            select0 = np.unique(label[index1])
            for i in select0:
                index1[np.where(label == i)[0]] = True
            select_remove = np.unique(label[index1])
            select_add = np.unique(label[np.invert(index1)])
        self.flush(select_add, select_remove)

    def rf_filter(self):
        rf_max = self.rfSlider.value()
        index0 = np.where(self.r_factor <= self.rf_range[rf_max], False, True)
        rf_remove = np.arange(self.rep)[index0]
        rf_add = np.arange(self.rep)[np.invert(index0)]
        self.rfLabel.setText('r-factor filter (<=%.3f)' % self.rf_range[rf_max])
        if self.rf_rep.size < rf_add.size:
            self.rf_rep = rf_add.copy()
        self.flush(rf_add, rf_remove)
        if self.rf_rep.size > rf_add.size:
            self.rf_rep = rf_add.copy()

    def flush(self, rep_add, rep_remove):
        flag = False
        for i in range(self.rep):
            if np.where(self.current_rep == i)[0].size > 0 and np.where(rep_remove == i)[0].size > 0\
                    and np.where(self.rf_rep == i)[0].size > 0:
                flag = True
                if not self.surface == '':
                    for j in range(self.local_size):
                        self.subs.removeItem(self.item_s[i * self.local_size + j])
                        self.subs.removeItem(self.item_c[i * 2 + j])
                for j in range(self.item_r[i].size):
                    self.rota.removeItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].removeItem(self.line_chi[pol][i])
                self.current_rep = np.delete(self.current_rep, np.where(self.current_rep == i)[0][0])
            elif np.where(self.current_rep == i)[0].size == 0 and np.where(rep_add == i)[0].size > 0\
                    and np.where(self.rf_rep == i)[0].size > 0:
                flag = True
                if not self.surface == '':
                    for j in range(self.local_size):
                        self.subs.addItem(self.item_s[i * self.local_size + j])
                        self.subs.addItem(self.item_c[i * 2 + j])
                for j in range(self.item_r[i].size):
                    self.rota.addItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].addItem(self.line_chi[pol][i])
                self.current_rep = np.insert(self.current_rep, 0, i)
        self.current_rep = np.sort(self.current_rep)

        if flag:
            self.amountLine.setText('%d / %d' % (self.current_rep.size, len(os.listdir(self.folder + r'/result'))))
            graph_size = 0
            if self.current_rep.size > 0:
                rdf, label = rdf_polarization(self.select_c, self.select_d, self.select_e, self.symbol,
                                              select=self.current_rep)
                for i in rdf:
                    if i.size > 0:
                        graph_size += 1
            urdf = [np.unique(rdf[_], return_counts=True) for _ in range(graph_size)]
            # urdf[self.range_target] = np.unique(self.rdf[self.range_target], return_counts=True)

            '''ba = np.array(
                [180 - cal_angle(self.select_c[_, 1, :], np.zeros(3), self.select_c[_, 2, :]) for _ in self.current_rep])
            if self.select_c.shape[1] > 3:
                ba = np.append(ba, np.array([
                    180 - cal_angle(self.select_c[_, 1, :], np.zeros(3), self.select_c[_, 3, :]) for _ in
                    self.current_rep]))
            ba = ba[np.where(ba > 130)[0]]
            urdf_ba = np.unique(np.round(ba, 0), return_counts=True)
            self.ba_bar.setOpts(x=urdf_ba[0], height=urdf_ba[1])
            print('%.1f-%.1f' % (ba.mean(), ba.std()))'''

            for i in range(graph_size):
                if i < 6 and self.current_rep.size > 0:
                    if rdf[i].size > 0:
                        self.info_dict[i].setText('%.3f-%.3f' % (rdf[i].mean(), rdf[i].std()))
                    else:
                        self.info_dict[i].setText('%.3f-%.3f' % (0, 0))
                self.bar_spher[i].setOpts(x=urdf[i][0], height=urdf[i][1])
                '''if (i == 0 or i == 3 or i == 6) and not urdf[i][0].size == 0:
                    self.widget_dict[i].setXRange(np.min(urdf[i][0] - 0.01), np.max(urdf[i][0] + 0.01))'''
            if self.current_rep.size > 0:
                for i in range(self.exp.size):
                    chi = np.sum(self.chi[i][self.current_rep], axis=0) / self.current_rep.size
                    self.line_aver[i].setData(x=self.exp[0].k, y=chi)
                    self.r_factor_temp[i] = self.exp[i].r_factor_chi(chi)
                    self.widget_plot[i].plotItem.legend.items[0][1].setText('%.3f' % self.r_factor_temp[i])

    def save_matplotlib_action(self):
        target = self.sender()
        if target.objectName() == 'plot_rdfAction':
            graph_size = 0
            if self.current_rep.size > 0:
                rdf, label = rdf_polarization(self.select_c, self.select_d, self.select_e, self.symbol,
                                              select=self.current_rep)
                for i in rdf:
                    if i.size > 0:
                        graph_size += 1
            urdf = [np.unique(rdf[_], return_counts=True) for _ in range(graph_size)]
            ax = np.array([plt.subplot(graph_size // 3, 3, _ + 1) for _ in range(graph_size)])
            urdf0 = [np.unique(self.rdf[_], return_counts=True) for _ in range(graph_size)]
            for i in range(graph_size):
                ax[i].bar(urdf[i][0], urdf[i][1], align='center', width=0.005 if i % 3 == 0 else 0.5, color='k')
                ax[i].set_xlabel('distance/Å' if i % 3 == 0 else 'angle/°')
                ax[i].set_ylabel('count')
                #ax[i].set_xlim(np.min(urdf0[i][0])-(0.01 if i % 3 == 0 else 1), np.max(urdf0[i][0])+(0.01 if i % 3 == 0 else 1))
                ax[i].set_ylim(0, np.max(urdf0[i][1])*1.1)
                ax[i].text(0.05, 0.8, self.rdf_name[i], fontsize=14, transform=ax[i].transAxes)
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.42, wspace=0.3)
        elif target.objectName() == 'plot_chiAction':
            plt.figure(figsize=(5, 10))
            ax = np.array([plt.subplot(3, 1, _ + 1) for _ in range(3)])
            title = ['E//[001]', 'E//[1-10]', 'E//[110]']
            if self.current_rep.size > 0:
                for i in range(self.exp.size):
                    chi = np.sum(self.chi[i][self.current_rep], axis=0) / self.current_rep.size
                    self.line_aver[i].setData(x=self.exp[0].k, y=chi)
                    self.r_factor_temp[i] = self.exp[i].r_factor_chi(chi)
                    self.widget_plot[i].plotItem.legend.items[0][1].setText('%.3f' % self.r_factor_temp[i])
                    ax[i].plot(self.exp[i].k[0:1], self.chi[i][0][0:1], c='b', linestyle='--', label='replica χ(k)')
                    for j in self.current_rep:
                        ax[i].plot(self.exp[i].k, self.chi[i][j], c='b', linestyle='--', alpha=0.05)
                    ax[i].plot(self.exp[i].k, chi, c='r', label='average χ(k)', linewidth=2.5)
                    ax[i].plot(self.exp[i].k, self.exp[i].chi, c='k', label='experiment χ(k)', linewidth=2.5)
                    ax[i].legend(loc='upper right')
                    ax[i].text(0.05, 0.9, title[i], fontsize=18, transform=ax[i].transAxes)
                    ax[i].set_xlabel('k/$Å^{-1}$', fontsize=18)
                    ax[i].set_ylabel('χ(k)', fontsize=18)
                    ax[i].tick_params(axis='both', labelsize=18)
                    ax[i].set_xlim(self.exp[i].k[0], self.exp[i].k[-1])
                    ax[i].set_ylim(np.min(self.chi[:, self.current_rep]) * 1.1, np.max(self.chi[:, self.current_rep]) * 1.1)
                    for j in range(int(self.exp[i].k[0])+1, int(self.exp[i].k[-1])+1):
                        ax[i].plot([j, j], [np.min(self.chi[:, self.current_rep])*1.1, np.max(self.chi[:, self.current_rep])*1.1], c='k', alpha=0.4, linestyle='--')
                plt.tight_layout()
                #plt.subplots_adjust(hspace=0.4, wspace=0.3)
        elif target.objectName() == 'plot_baAction':
            ba = np.array([])
            for rep in self.current_rep:
                angle = cal_angle(self.select_c[rep][1], self.select_c[rep][0], self.select_c[rep][2])
                ba = np.append(ba, angle if angle > 90 else (180 - angle))
            urdf = np.unique(np.round(ba, 0), return_counts=True)
            ax = plt.subplot()
            ax.bar(urdf[0], urdf[1], align='center', width=0.5, color='k')
            ax.set_xlabel('angle/°')
            ax.set_ylabel('count')
            ax.text(0.05, 0.8, 'S-Cu-O', fontsize=14, transform=ax.transAxes)
            ax.set_ylim(0, np.max(urdf[1]) * 1.1)
            plt.tight_layout()
        plt.show()

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

    def choice_window(self):
        msg = QMessageBox(self)
        msg.setWindowTitle('resulted model')
        msg.setText('Please select a set of model')
        msg.setIcon(QMessageBox.Question)
        msg.addButton('Best', QMessageBox.YesRole)
        msg.addButton('Final', QMessageBox.NoRole)
        msg.exec_()
        if msg.clickedButton().text() == 'Best':
            return 'Best'
        else:
            return 'final'


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
