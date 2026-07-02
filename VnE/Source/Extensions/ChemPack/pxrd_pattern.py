from pymatgen.io.cif import CifParser
from powdiffrac.simulation.core import Powder
from pymatgen.analysis.diffraction.xrd import XRDCalculator, WAVELENGTHS
import os

import numpy as np
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget,
    QListWidgetItem, QGroupBox, QFormLayout, QDoubleSpinBox, QSpinBox,
    QComboBox, QCheckBox, QColorDialog, QFileDialog, QSplitter, QLabel
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QPen, QColor
import pyqtgraph as pg
from ..ChemPack.ui.select_mol_dialog import SelectMolDialog
from . import MOLECULE_SYSTEMS

WID = None

class Pxrd:
    def __init__(self, file_name):
        self.pipe = [self.loadPxrd,
                     self.loadCif,
                     self.loadStr,
                     self.loadPattern,
                     self.loadPowder,
                     self.loadPlotData,
                     ]
        self.pipe_state = {self.setFile: False,
                           self.loadPxrd: False,
                           self.loadCif: False,
                           self.loadStr: False,
                           self.loadPattern: False,
                           self.loadPowder: False,
                           self.loadPlotData: False,
                           }
        self.file_name = None
        self.setFile(file_name)
        self.waves = WAVELENGTHS
        self.wave_len = self.waves['CuKa1']
        self.min_2th = 3
        self.max_2th = 120
        self.pymatgen_str = None
        self.pxrd = None
        self.pxrd_pattern = None
        self.powder = None
        self.min_dom, self.max_dom = None, None
        self.min_2th, self.max_2th = None, None
        self.step = None
        self.x, self.y = None, None

    def setFile(self, file_name):
        if os.path.splitext(file_name)[-1] == '.cif':
            self.file_name = file_name
            self.pipe_state[self.setFile] = True

    def loadCif(self):
        if self.pipe_state[self.setFile]:
            self.cif = CifParser.from_str(open(self.file_name, 'r').read())
            self.pymatgen_str = None
            self.pxrd_pattern = None
            self.pipe_state[self.loadCif] = True
            self.pipe_state[self.loadStr] = False
            self.pipe_state[self.loadPattern] = False

    def loadStr(self):
        if self.pipe_state[self.loadCif]:
            self.pymatgen_str = self.cif.parse_structures(symmetrized=True, primitive=False)[0]
            self.pxrd_pattern = None
            self.pipe_state[self.loadStr] = True
            self.pipe_state[self.loadPattern] = False

    def loadPxrd(self):
        self.pxrd = XRDCalculator(wavelength=self.wave_len, symprec=0.02)
        self.pxrd_pattern = None
        self.pipe_state[self.loadPxrd] = True
        self.pipe_state[self.loadPattern] = False

    def loadPattern(self):
        if self.pipe_state[self.loadPxrd] and self.pipe_state[self.loadStr]:
            self.pxrd_pattern = self.pxrd.get_pattern(self.pymatgen_str, scaled=True)
            self.pipe_state[self.loadPattern] = True

    def loadPowder(self):
        if self.pipe_state[self.loadStr]:
            self.powder = Powder(
                self.pymatgen_str, two_theta=(self.min_2th, self.max_2th),
                step_size=self.step, max_strain=0.04, max_texture=0.0,
                min_domain_size=self.min_dom, max_domain_size=self.max_dom,
                peak_shape='pseudo_voigt', vary_domain=True, vary_strain=True
            )
            self.pipe_state[self.loadPowder] = True

    def loadPlotData(self):
        if self.pipe_state[self.loadPowder]:
            self.y = self.powder.get_signal(vary=True, scaling=10000)
            self.x = self.powder.steps
            self.pipe_state[self.loadPlotData] = True

    def change2th(self, min=3, max=120):
        self.min_2th = min
        self.max_2th = max
        self.pipe_state[self.loadPowder] = False

    def changeWave(self, wave):
        if isinstance(wave, str):
            self.wave_len = self.waves[wave]
            self.pipe_state[self.loadPowder] = False
        elif isinstance(wave, float):
            self.wave_len = wave
            self.pipe_state[self.loadPowder] = False

    def changeDom(self, min=47, max=53):
        self.min_dom = min
        self.max_dom = max
        self.pipe_state[self.loadPowder] = False

    def changeStep(self, step=0.02):
        self.step = step
        self.pipe_state[self.loadPowder] = False

    def process(self):
        for st in self.pipe:
            if not self.pipe_state[st]:
                st()

    def default(self):
        self.change2th()
        self.changeWave('CuKa1')
        self.changeDom()
        self.changeStep()


class PowderPattern:
    """
    Модель данных одной порошкограммы.
    Хранит исходные/текущие данные, параметры и графический элемент.
    """
    def __init__(self, name, mol_sys=None, wavelength=1.54184, size_params=(43, 54),
                 two_theta_range=(10.0, 90.0), num_points=5000, x=None, y=None):
        self.mol_sys = mol_sys
        self.name = name
        self.x = x
        self.y = y
        self.wave_len = wavelength
        self.size_params = size_params
        self.two_theta_range = two_theta_range
        self.num_points = num_points
        self.visible = True
        self.pen = pg.mkPen(color='b', width=1.5)
        self.plot_item: pg.PlotDataItem = None

    def recalculate(self):
        pxrd = XRDCalculator(wavelength=self.wave_len, symprec=0.02)
        cif = CifParser.from_str(open(self.mol_sys.file_name, 'r').read())
        pymatgen_str = cif.parse_structures(symmetrized=True, primitive=False)[0]
        powder = Powder(pymatgen_str, two_theta=self.two_theta_range,
            step_size=(self.two_theta_range[1] - self.two_theta_range[0])/self.num_points, max_strain=0.00, max_texture=0.0,
            min_domain_size=self.size_params[0], max_domain_size=self.size_params[1],
            peak_shape='gaussian', vary_domain=True, vary_strain=True, seed=0
        )
        powder._calculator = pxrd
        powder.radiation = self.wave_len
        self.y = powder.get_signal(vary=True, scaling=10000)
        self.x = powder.steps
        ...

    def update_plot(self):
        if self.plot_item is not None:
            self.plot_item.setData(self.x, self.y)


class PowderDiffractionWindow(QWidget):
    """
    Окно для отображения, редактирования и сравнения порошкограмм.
    """
    pattern_selected = Signal(object)   # выбранный объект PowderPattern

    def __init__(self, parent=None):
        super().__init__(parent)
        #self.setWindowTitle("Порошковый дифракционный анализ")
        self.patterns: list[PowderPattern] = []
        self.selected_pattern: PowderPattern = None
        self.cursor_line = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen('r', width=1))
        self.cursor_point = pg.ScatterPlotItem(size=10, pen=pg.mkPen('r'), brush=pg.mkBrush('r'))
        self.cursor_label = pg.TextItem(color='r', anchor=(0.5, 1.5))
        self.cursor_visible = False

        self.dialog = None

        self._init_ui()
        self._connect_signals()

    def _init_ui(self):
        main_layout = QHBoxLayout(self)

        # ---------- Левая часть: график ----------
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setLabel('bottom', '2θ', units='°')
        self.plot_widget.setLabel('left', 'Intensity')
        self.plot_widget.showGrid(x=True, y=True)
        self.plot_widget.setBackground('w')
        # Добавляем элементы курсора (изначально скрыты)
        self.plot_widget.addItem(self.cursor_line)
        self.plot_widget.addItem(self.cursor_point)
        self.plot_widget.addItem(self.cursor_label)
        self.cursor_line.hide()
        self.cursor_point.hide()
        self.cursor_label.hide()

        # ---------- Правая часть: панель управления ----------
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)

        # Кнопка добавления новой порошкограммы
        btn_add = QPushButton("Add a powder pattern")
        control_layout.addWidget(btn_add)

        # Список всех порошкограмм
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QListWidget.SingleSelection)
        control_layout.addWidget(QLabel("List of powder patterns:"))
        control_layout.addWidget(self.list_widget)

        # Удаление выбранной
        btn_remove = QPushButton("Delete selected")
        control_layout.addWidget(btn_remove)

        # Группа параметров выбранной порошкограммы
        params_group = QGroupBox("Parameters of the selected pattern")
        form = QFormLayout()

        # Длина волны
        self.wavelength_combo = QComboBox()
        waves = [f'{x} ({WAVELENGTHS[x]} Å)' for x in WAVELENGTHS]
        waves.append("Other")
        self.wavelength_combo.addItems(waves)
        self.wavelength_spin = QDoubleSpinBox()
        self.wavelength_spin.setDecimals(5)
        self.wavelength_spin.setRange(0.1, 10.0)
        self.wavelength_spin.setValue(1.54184)
        self.wavelength_spin.setSingleStep(0.00001)
        self.wavelength_spin.setEnabled(False)
        wl_layout = QHBoxLayout()
        wl_layout.addWidget(self.wavelength_combo)
        wl_layout.addWidget(self.wavelength_spin)
        form.addRow("Wavelength:", wl_layout)

        # Параметры уширения (размер кристаллита)
        self.size1_spin = QDoubleSpinBox()
        self.size1_spin.setRange(1.0, 250.0)
        self.size1_spin.setSingleStep(1)
        self.size1_spin.setValue(43)
        self.size2_spin = QDoubleSpinBox()
        self.size2_spin.setRange(1, 250.0)
        self.size2_spin.setSingleStep(1)
        self.size2_spin.setValue(54)
        size_layout = QHBoxLayout()
        size_layout.addWidget(QLabel("Min:"))
        size_layout.addWidget(self.size1_spin)
        size_layout.addWidget(QLabel("Max:"))
        size_layout.addWidget(self.size2_spin)
        form.addRow("Crystallite:", size_layout)

        # Область расчёта 2θ
        self.tth_min_spin = QDoubleSpinBox()
        self.tth_min_spin.setRange(0.0, 180.0)
        self.tth_min_spin.setValue(10.0)
        self.tth_max_spin = QDoubleSpinBox()
        self.tth_max_spin.setRange(0.0, 180.0)
        self.tth_max_spin.setValue(90.0)
        tth_layout = QHBoxLayout()
        tth_layout.addWidget(QLabel("from"))
        tth_layout.addWidget(self.tth_min_spin)
        tth_layout.addWidget(QLabel("to"))
        tth_layout.addWidget(self.tth_max_spin)
        form.addRow("2θ range:", tth_layout)

        # Число точек
        self.num_points_spin = QSpinBox()
        self.num_points_spin.setRange(100, 50000)
        self.num_points_spin.setSingleStep(100)
        self.num_points_spin.setValue(5000)
        form.addRow("Number of points:", self.num_points_spin)

        # Видимость
        self.visible_check = QCheckBox("Visible")
        self.visible_check.setChecked(True)
        form.addRow(self.visible_check)

        params_group.setLayout(form)
        control_layout.addWidget(params_group)

        # Кнопка применения параметров
        btn_apply_params = QPushButton("Apply")
        control_layout.addWidget(btn_apply_params)

        # Кнопка изменения стиля линии
        btn_style = QPushButton("Line style")
        control_layout.addWidget(btn_style)

        # Группа для разностного/суммарного графика
        combo_group = QGroupBox("Combinations")
        combo_form = QFormLayout()

        self.combo_pattern1 = QComboBox()
        self.combo_pattern2 = QComboBox()
        self.weight1_spin = QDoubleSpinBox()
        self.weight1_spin.setRange(-10.0, 10.0)
        self.weight1_spin.setValue(1.0)
        self.weight1_spin.setSingleStep(0.1)
        self.weight2_spin = QDoubleSpinBox()
        self.weight2_spin.setRange(-10.0, 10.0)
        self.weight2_spin.setValue(1.0)
        self.weight2_spin.setSingleStep(0.1)

        combo_form.addRow("Pattern 1:", self.combo_pattern1)
        combo_form.addRow("Weight 1:", self.weight1_spin)
        combo_form.addRow("Pattern 2:", self.combo_pattern2)
        combo_form.addRow("Weight 2:", self.weight2_spin)

        btn_sum = QPushButton("Sum (Weight1*I1 + Weight2*I2)")
        btn_diff = QPushButton("Difference (Weight1*I1 - Weight2*I2)")
        combo_form.addRow(btn_sum)
        combo_form.addRow(btn_diff)

        combo_group.setLayout(combo_form)
        control_layout.addWidget(combo_group)

        # Цвет фона
        btn_bg = QPushButton("Background color")
        control_layout.addWidget(btn_bg)

        control_layout.addStretch()

        # Размещаем график и панель в сплиттере
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.plot_widget)
        splitter.addWidget(control_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter)

        # Сохраняем ссылки на кнопки для подключения сигналов
        self.btn_add = btn_add
        self.btn_remove = btn_remove
        self.btn_apply_params = btn_apply_params
        self.btn_style = btn_style
        self.btn_sum = btn_sum
        self.btn_diff = btn_diff
        self.btn_bg = btn_bg

    def _connect_signals(self):
        self.btn_add.clicked.connect(self.on_add_pattern)
        self.btn_remove.clicked.connect(self.on_remove_pattern)
        self.btn_apply_params.clicked.connect(self.on_apply_parameters)
        self.btn_style.clicked.connect(self.on_change_style)
        self.btn_sum.clicked.connect(self.on_build_sum)
        self.btn_diff.clicked.connect(self.on_build_difference)
        self.btn_bg.clicked.connect(self.on_change_background)
        self.list_widget.currentItemChanged.connect(self.on_selection_changed)
        self.wavelength_combo.currentIndexChanged.connect(self.on_wavelength_combo_changed)
        # Курсор мыши
        self.plot_widget.scene().sigMouseMoved.connect(self.on_mouse_moved)

    # ===================== ЛОГИКА =====================

    @Slot()
    def on_add_pattern(self):
        """
        Открывает диалог выбора файла с данными (2 колонки: 2θ, интенсивность),
        создаёт объект PowderPattern, добавляет в список и отрисовывает.
        При необходимости можно использовать генерацию модельных данных.
        """

        def process(args):
            point_list, mol_sys = args
            pxrd = PowderPattern(point_list.name, mol_sys)
            pxrd.recalculate()
            self.patterns.append(pxrd)
            self._add_plot(pxrd)
            self._update_list()
            self._update_combo_boxes()
            self.on_apply_parameters(pxrd)


        mol_sys = {k: v for k, v in MOLECULE_SYSTEMS.items() if v.file_name and os.path.splitext(v.file_name)[-1] == '.cif'}
        self.dialog = SelectMolDialog(mol_sys, process)
        self.dialog.show()

    def _add_plot(self, pattern: PowderPattern):
        """Добавляет PlotDataItem на график и сохраняет ссылку."""
        item = self.plot_widget.plot(pattern.x, pattern.y, pen=pattern.pen, name=pattern.name)
        pattern.plot_item = item
        item.setVisible(pattern.visible)

    @Slot()
    def on_remove_pattern(self):
        """Удаляет выбранную порошкограмму из списка и с графика."""
        if self.selected_pattern is None:
            return
        self.plot_widget.removeItem(self.selected_pattern.plot_item)
        self.patterns.remove(self.selected_pattern)
        self.selected_pattern = None
        self._update_list()
        self._update_combo_boxes()

    @Slot()
    def on_apply_parameters(self, patt=None):
        """
        Считывает значения из полей параметров и применяет их к выбранной порошкограмме.
        Вызывает пересчёт данных и обновление графика.
        """
        if patt:
            p = patt
        elif self.selected_pattern is None:
            return
        else:
            p = self.selected_pattern
        # Длина волны
        if self.wavelength_combo.currentText() == "Другая...":
            p.wave_len = self.wavelength_spin.value()
        else:
            # В зависимости от выбора можно задать словарь
            wl_dict = {f'{k} ({v} Å)': v for k,v in WAVELENGTHS.items()}
            p.wave_len = wl_dict.get(self.wavelength_combo.currentText(), 1.54184)
        # Уширение
        p.size_params = (self.size1_spin.value(), self.size2_spin.value())
        # Диапазон
        p.two_theta_range = (self.tth_min_spin.value(), self.tth_max_spin.value())
        # Число точек
        p.num_points = self.num_points_spin.value()
        # Видимость
        p.visible = self.visible_check.isChecked()
        p.plot_item.setVisible(p.visible)

        p.recalculate()
        p.update_plot()

    @Slot()
    def on_change_style(self):
        """Открывает диалог выбора цвета, толщины и стиля линии для выделенной порошкограммы."""
        if self.selected_pattern is None:
            return
        color = QColorDialog.getColor(initial=self.selected_pattern.pen.color(), parent=self)
        if color.isValid():
            pen = pg.mkPen(color=color, width=2, style=Qt.SolidLine)  # стиль можно менять дополнительно
            self.selected_pattern.pen = pen
            self.selected_pattern.plot_item.setPen(pen)

    @Slot()
    def on_change_background(self):
        """Меняет цвет фона графика."""
        color = QColorDialog.getColor(initial=self.plot_widget.backgroundBrush().color(), parent=self)
        if color.isValid():
            self.plot_widget.setBackground(color)

    @Slot()
    def on_wavelength_combo_changed(self, idx):
        """Включает/отключает поле ввода своей длины волны."""
        if self.wavelength_combo.currentText() == "Другая...":
            self.wavelength_spin.setEnabled(True)
        else:
            self.wavelength_spin.setEnabled(False)

    def _update_list(self):
        """Обновляет QListWidget именами порошкограмм."""
        self.list_widget.blockSignals(True)
        self.list_widget.clear()
        for p in self.patterns:
            item = QListWidgetItem(p.name)
            item.setData(Qt.UserRole, p)
            self.list_widget.addItem(item)
        self.list_widget.blockSignals(False)

    def _update_combo_boxes(self):
        """Заполняет комбобоксы для выбора порошкограмм в секции комбинаций."""
        self.combo_pattern1.clear()
        self.combo_pattern2.clear()
        for p in self.patterns:
            self.combo_pattern1.addItem(p.name, p)
            self.combo_pattern2.addItem(p.name, p)

    @Slot(QListWidgetItem, QListWidgetItem)
    def on_selection_changed(self, current, previous):
        """Обрабатывает выбор порошкограммы в списке: обновляет selected_pattern и поля параметров."""
        if current is None:
            self.selected_pattern = None
            return
        pattern = current.data(Qt.UserRole)
        self.selected_pattern = pattern
        # Заполняем поля параметров
        self._populate_parameters(pattern)
        self.pattern_selected.emit(pattern)

    def _populate_parameters(self, pattern: PowderPattern):
        """Записывает параметры выбранной порошкограммы в элементы управления."""
        # Длина волны
        wl = pattern.wave_len
        for k, v in WAVELENGTHS.items():
            if abs(wl - v) < 1e-4:
                self.wavelength_combo.setCurrentText(f'{k} ({v} Å)')
                break
        else:
            self.wavelength_combo.setCurrentText("Другая...")
            self.wavelength_spin.setValue(wl)
        self.size1_spin.setValue(pattern.size_params[0])
        self.size2_spin.setValue(pattern.size_params[1])
        self.tth_min_spin.setValue(pattern.two_theta_range[0])
        self.tth_max_spin.setValue(pattern.two_theta_range[1])
        self.num_points_spin.setValue(pattern.num_points)
        self.visible_check.setChecked(pattern.visible)

    @Slot(object)
    def on_mouse_moved(self, pos):
        """
        Реализует курсор: вертикальная линия, точка и подпись значения Y
        для выделенного графика.
        """
        if self.selected_pattern is None or not self.selected_pattern.visible:
            self.cursor_line.hide()
            self.cursor_point.hide()
            self.cursor_label.hide()
            return

        # Преобразование координат сцены в координаты графика
        mouse_point = self.plot_widget.plotItem.vb.mapSceneToView(pos)
        x = mouse_point.x()

        # Находим ближайшую точку данных выделенной порошкограммы
        pattern = self.selected_pattern
        x_data = pattern.x
        y_data = pattern.y
        if len(x_data) == 0:
            return

        idx = np.searchsorted(x_data, x)
        idx = np.clip(idx, 0, len(x_data) - 1)
        # Уточнение ближайшей точки
        if idx > 0 and abs(x_data[idx - 1] - x) < abs(x_data[idx] - x):
            idx -= 1
        nearest_x = x_data[idx]
        nearest_y = y_data[idx]

        # Обновляем элементы курсора
        self.cursor_line.setPos(nearest_x)
        self.cursor_point.setData([nearest_x], [nearest_y])
        self.cursor_label.setText(f"{nearest_y:.3f}")
        self.cursor_label.setPos(nearest_x, nearest_y)

        if not self.cursor_visible:
            self.cursor_line.show()
            self.cursor_point.show()
            self.cursor_label.show()
            self.cursor_visible = True

    @Slot()
    def on_build_sum(self):
        """Строит суммарный график I = w1*I1 + w2*I2."""
        self._build_combination(is_difference=False)

    @Slot()
    def on_build_difference(self):
        """Строит разностный график I = w1*I1 - w2*I2."""
        self._build_combination(is_difference=True)

    def _build_combination(self, is_difference):
        """Общая логика построения комбинации двух выбранных порошкограмм."""
        p1 = self.combo_pattern1.currentData()
        p2 = self.combo_pattern2.currentData()
        if p1 is None or p2 is None:
            return
        w1 = self.weight1_spin.value()
        w2 = self.weight2_spin.value()

        # Интерполяция на общую ось x (используем ось первого, либо общую)
        x_common = np.linspace(
            max(p1.x[0], p2.x[0]),
            min(p1.x[-1], p2.x[-1]),
            max(len(p1.x), len(p2.x))
        )
        y1_interp = np.interp(x_common, p1.x, p1.y)
        y2_interp = np.interp(x_common, p2.x, p2.y)

        if is_difference:
            y_comb = w1 * y1_interp - w2 * y2_interp
            name = f"{p1.name} - {p2.name}"
        else:
            y_comb = w1 * y1_interp + w2 * y2_interp
            name = f"{p1.name} + {p2.name}"

        # Создаём новый объект PowderPattern для результирующего графика
        combo_pattern = PowderPattern(name, x=x_common, y=y_comb)
        combo_pattern.recalculate = lambda: None  # пересчёт не меняет данные
        self.patterns.append(combo_pattern)
        self._add_plot(combo_pattern)
        self._update_list()
        self._update_combo_boxes()

    def leaveEvent(self, event):
        """Скрывает курсор при уходе мыши с виджета (дополнительная мера)."""
        self.cursor_line.hide()
        self.cursor_point.hide()
        self.cursor_label.hide()
        self.cursor_visible = False
        super().leaveEvent(event)


def process(args):
    point_list, mol_sys = args
    pxrd = Pxrd(mol_sys.file_name)
    pxrd.default()
    pxrd.process()
    a = 0
    pw = pg.plot(pxrd.x, pxrd.y)
    ...

def execute():
    global WID
    WID = PowderDiffractionWindow()
    WID.show()

def setup(menu, model, *args, **kwargs):
    from PySide6.QtGui import QAction

    global TREE_MODEL
    TREE_MODEL = model

    action = QAction('Powder')
    action.triggered.connect(execute)
    menu.addAction(action)

    actions = [action]
    return actions