from mrmc_package import EXP, init_3dplot, TABLE_POL, norm_fft
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

        for i in range(9):
            self.action_dict[i].triggered.connect(self.select)
        self.heightAction.triggered.connect(self.select)
        self.readAction.triggered.connect(self.read)
        self.saveAction.triggered.connect(self.save)
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

    def load_exp(self, file_list):
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
            for i in range(9):
                self.widget_dict[i].clear()
        self.folder = folder
        with open(self.folder + r'\info.txt', 'r') as f:
            self.local_e = np.array([f.readline().split()[1]])
            self.local_e = np.append(self.local_e, f.readline().split()[1:])
            self.local_size = self.local_e.size
            self.surface = f.readline().split()[1]
            dE = [float(_) for _ in f.readline().split()[1:]]
            k = [float(_) for _ in f.readline().split()[1:]]
            r = [float(_) for _ in f.readline().split()[1:]]
            print(k, r)
        self.local_c, flag = read_rep(self.folder, self.choice_window(), self.local_size)
        if flag:
            self.statusbar.showMessage('No result file found', 3000)
            self.folder = ''
            return
        self.surface_c, self.surface_e = load_substrate(self.folder, self.local_size)

        if self.surface == 'Al2O3':
            plot_Al2O3(self.surface_c, self.surface_e, self.subs)
        elif self.surface == 'TiO2':
            plot_TiO2(self.surface_c, self.surface_e, self.subs, local=False)

        self.rep = self.local_c.shape[0]
        self.amountLine.setText('%d / %d' % (self.rep, len(os.listdir(self.folder + r'/result'))))

        with open(self.folder + r'/mrmc.inp', 'r') as f:
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
            self.load_exp(exp_path)
            while True:
                lines = f.readline()
                if not lines.find('weight') == -1:
                    break
            weight = int(lines.split(':')[1])
            while True:
                lines = f.readline()
                if not lines.find('S0') == -1:
                    break
            s0 = float(lines.split(':')[1])
            sig2 = float(f.readline().split(':')[1])
            for i in range(3):
                f.readline()
            rpath = float(f.readline().split(':')[1])
            temp = f.readline().split(':')[1].strip()
            if temp == 'True' or temp == 'true' or temp == '1':
                ms = True
            elif temp == 'False' or temp == 'false' or temp == '0':
                ms = False
            while True:
                lines = f.readline()
                if not lines.find('material_folder') == -1:
                    break
            f_material = lines.split(lines.split(':')[0] + ':')[1].split('\n')[0].strip()
        if not self.surface == '':
            self.select_c, self.select_d, self.select_e = select_atom(self.surface_c, self.surface_e,
                                                                      self.local_c, self.local_e, rpath)
        print(self.local_c.shape)

        if self.surface == 'Al2O3':
            self.symbol = np.append(self.local_e[1:], ['O', 'Al'])
        elif self.surface == 'TiO2':
            self.symbol = np.append(self.local_e[1:], ['O', 'Ti'])
        else:
            self.symbol = self.local_e[1:].copy()
        self.symbol = np.unique(self.symbol)
        if not self.material == 'CuAlO':
            graph_size = self.symbol.size if self.symbol.size <= 3 else 3
            self.rdf, self.rdf_label = rdf_polarization(self.select_c, self.select_d, self.select_e, self.symbol)
            urdf = [np.unique(self.rdf[_], return_counts=True) for _ in range(graph_size * 3)]
            self.bar_spher = [bar(urdf[_][0], urdf[_][1], width=0.5) for _ in range(graph_size * 3)]
            for i in range(graph_size * 3):
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
            for i in range(graph_size):
                self.widget_dict[i * 3].setTitle('R (%s-%s)' % (self.local_e[0], self.symbol[i]),
                                                 color='#000000', size='20pt')
                self.widget_dict[i * 3 + 1].setTitle('φ (%s-%s)' % (self.local_e[0], self.symbol[i]),
                                                     color='#000000', size='20pt')
                self.widget_dict[i * 3 + 2].setTitle('θ (%s-%s)' % (self.local_e[0], self.symbol[i]),
                                                     color='#000000', size='20pt')
            self.do.plotItem.enableAutoRange(axis='x', enable=False)
            self.exp = np.array([EXP(self.exp_name[i], k[0], k[1], r[0], r[1]) for i in range(self.exp_name.size)])

        table = np.zeros((3, self.rep), dtype=TABLE_POL)
        for pol in range(3):
            for i in range(self.rep):
                table[pol][i] = TABLE_POL(self.exp[pol].k_start, self.exp[pol].k_end, self.exp[pol].r_start,
                                          self.exp[pol].r_end, sig2, dE, s0, self.exp[pol].k0,
                                          self.select_c[i], self.select_e[i], f_material, pol, ms_en=ms, weight=weight)
        chi_sum = np.zeros((3, self.exp[0].k.size))
        #ax = np.array([plt.subplot(3, 1, _ + 1) for _ in range(3)])
        #direct = ['[001]', '[1-10]', '[110]']
        for pol in range(3):
            for i in range(self.rep):
                chi_sum[pol] += table[pol][i].chi
            chi_sum[pol] /= self.rep
            print(self.exp[pol].r_factor(chi_sum[pol]))
            '''ax[pol].plot(self.exp[pol].k, self.exp[pol].chi, c='k', label='expe %s'%direct[pol])
            #ax[pol].plot(self.exp[pol].r, np.abs(norm_fft(chi_sum[pol], self.exp[pol].r.size)), c='r', label='average %s'%direct[pol])
            ax[pol].plot(self.exp[pol].k, chi_sum[pol], c='r', label='average %s' % direct[pol])
            ax[pol].legend(loc='lower right')
            ax[pol].set_ylabel('χ %s'%direct[pol])'''
        #ax[2].set_xlabel(r'$k[Å^{-1}]$')
        #plt.show()

        self.item_r = plot_rotate(self.surface_c, self.surface_e, self.local_c, self.local_e, rpath, self.surface,
                                  self.rota)

        self.item_s, self.item_c = plot_on_substrate(self.surface_c, self.surface_e, self.local_c, rpath, self.subs)

        self.line_exp = np.array([line(x=self.exp[_].k, y=self.exp[_].chi, c='black', width=3)
                                  for _ in range(self.exp.size)])
        self.chi = read_chi(self.folder, self.rep, self.exp.size)
        self.line_chi = np.zeros((self.exp.size, self.rep), dtype=pg.PlotDataItem)
        for i in range(self.exp.size):
            self.widget_plot[i].addItem(self.line_exp[i])
            for j in range(self.rep):
                self.line_chi[i][j] = line(x=self.exp[0].k, y=table[i][j].chi, c='red', alpha=0.3)
                self.widget_plot[i].addItem(self.line_chi[i][j])

        self.line_aver = np.zeros(3, dtype=pg.PlotDataItem)
        for i in range(self.exp.size):
            self.line_aver[i] = line(x=self.exp[i].k, y=chi_sum[i], c='blue', width=3)
            self.line_aver[i].opts['name'] = '%.3f' % self.exp[i].r_factor(chi_sum[i])
            self.widget_plot[i].addItem(self.line_aver[i])

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
        graph_size = self.symbol.size if self.symbol.size <= 3 else 3
        data = self.rdf[self.range_target]
        label = self.rdf_label[self.range_target//graph_size]
        if manual:
            hmin, hmax = self.region_line.getRegion()
            index0 = np.where((hmin < data) & (data < hmax), False, True)
            select0 = np.unique(label[index0])
            for i in select0:
                index0[np.where(label == i)[0]] = True
            index1 = np.array([False] * self.rdf[self.range_target].size)
            select1 = np.arange(self.rep)
        else:
            index0 = np.where((self.hmin < data) & (data < self.hmax), False, True)
            select0 = np.unique(label[index0])
            for i in select0:
                index0[np.where(label == i)[0]] = True
            self.hmin, self.hmax = self.region_line.getRegion()
            index1 = np.where((self.hmin < data) & (data < self.hmax), False, True)
            select1 = np.unique(label[index1])
            for i in select1:
                index1[np.where(label == i)[0]] = True

        for i in range(self.rep):
            if np.where(select0 == i)[0].size == 0 and np.where(select1 == i)[0].size > 0:
                self.subs.removeItem(self.item_s[i * self.local_size])
                self.subs.removeItem(self.item_s[i * self.local_size + 1])
                for j in range(self.item_c[i].size):
                    self.subs.removeItem(self.item_c[i][j])
                for j in range(self.item_r[i].size):
                    self.rota.removeItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].removeItem(self.line_chi[pol][i])
            elif np.where(select1 == i)[0].size == 0 and np.where(select0 == i)[0].size > 0:
                self.subs.addItem(self.item_s[i * self.local_size])
                self.subs.addItem(self.item_s[i * self.local_size + 1])
                for j in range(self.item_c[i].size):
                    self.subs.addItem(self.item_c[i][j])
                for j in range(self.item_r[i].size):
                    self.rota.addItem(self.item_r[i][j])
                for pol in range(self.exp.size):
                    self.widget_plot[pol].addItem(self.line_chi[pol][i])

        if (index0 == index1).all():
            index1 = np.invert(index1)
            select = np.unique(label[index1])
            self.amountLine.setText('%d / %d' % (select.size, len(os.listdir(self.folder + r'/result'))))
            rdf, label = rdf_polarization(self.select_c, self.select_d, self.select_e, self.symbol, select=select)
            urdf = [np.unique(rdf[_], return_counts=True) for _ in range(graph_size * 3)]
            urdf[self.range_target] = np.unique(self.rdf[self.range_target], return_counts=True)
            for i in range(graph_size * 3):
                if i < 6:
                    self.info_dict[i].setText('%.2f-%.2f' % (rdf[i].mean(), rdf[i].std()))
                self.bar_spher[i].setOpts(x=urdf[i][0], height=urdf[i][1])
                if (i == 0 or i == 3 or i == 6) and not urdf[i][0].size == 0:
                    self.widget_dict[i].setXRange(np.min(urdf[i][0] - 0.01), np.max(urdf[i][0] + 0.01))
            if index1.any():
                for i in range(self.exp.size):
                    chi = np.sum(self.chi[i][select], axis=0) / select.size
                    self.line_aver[i].setData(x=self.exp[0].k, y=chi)
                    self.widget_plot[i].plotItem.legend.items[0][1].setText('%.3f' % self.exp[i].r_factor(chi))

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
