# Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
# This file is part of ASID - Atomistic Simulation Instruments and Database
# For more information see <https://github.com/ASID-Production/ASID>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# *****************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# *****************************************************************************************

import os.path

from django.db import models
from django.contrib.auth import get_user_model

User = get_user_model()

CENTRINGS = [
    (1, 'P'),
    (2, 'I'),
    (3, 'A'),
    (4, 'B'),
    (5, 'C'),
    (6, 'F'),
    (7, 'R'),
]
SYSTEMS = [
    (1, 'triclinic'),
    (2, 'monoclinic'),
    (3, 'orthorombic'),
    (4, 'tetragonal'),
    (5, 'trigonal'),
    (6, 'hexagonal'),
    (7, 'cubic'),
    (8, 'rhombohedron'),
]
SHAPES = [
    ('block', 'block'),
    ('plate', 'plate'),
    ('needle', 'needle'),
    ('prism', 'prism'),
    ('irregular', 'irregular'),
    ('cube', 'cube'),
    ('trapezoid', 'trapezoid'),
    ('rect. prism', 'rect. prism'),
    ('rhombohedral', 'rhombohedral'),
    ('hexagonal', 'hexagonal'),
    ('octahedral', 'octahedral'),
    ('plank', 'plank')
]


class Author(models.Model):
    '''Table with authors.'''
    family_name = models.CharField(blank=True, null=True, verbose_name='author', max_length=250)
    initials = models.CharField(blank=True, null=True, verbose_name='author', max_length=250)

    def __str__(self):
        if self.initials:
            return f'{self.family_name} {self.initials}'
        return self.family_name


class StructureCode(models.Model):
    '''Table with codes.'''
    refcode = models.CharField(verbose_name='Refcode', max_length=17, unique=True)
    CCDC_number = models.CharField(verbose_name='CCDC number', max_length=10, blank=True, null=True)
    ICSD = models.BooleanField(verbose_name='ICSD structure', blank=True, null=True)
    COD = models.BooleanField(verbose_name='Crystallography Open Database structure', blank=True, null=True)
    user = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        related_name='refcode',
        null=True,
        blank=True,
        verbose_name='Owner'
    )
    authors = models.ManyToManyField(
        Author,
        related_name='refcodes',
        blank=True,
        verbose_name='Authors'
    )

    def __str__(self):
        return self.refcode


class CalculationProgram(models.Model):
    '''Table with caclulation types.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='calc_prog',
        on_delete=models.CASCADE
    )
    orca = models.BooleanField(verbose_name='Orca calculation', blank=True, null=True)
    vasp = models.BooleanField(verbose_name='VASP calculation', blank=True, null=True)
    qchem = models.BooleanField(verbose_name='QChem calculation', blank=True, null=True)
    gaussian = models.BooleanField(verbose_name='Gaussian calculation', blank=True, null=True)
    nwchem = models.BooleanField(verbose_name='NWChem calculation', blank=True, null=True)
    abinit = models.BooleanField(verbose_name='ABINIT calculation', blank=True, null=True)
    crystal = models.BooleanField(verbose_name='CRYSTAL calculation', blank=True, null=True)
    gamess = models.BooleanField(verbose_name='GAMESS calculation', blank=True, null=True)
    mopac = models.BooleanField(verbose_name='MOPAC calculation', blank=True, null=True)
    quantum_espresso = models.BooleanField(verbose_name='Quantum ESPRESSO calculation', blank=True, null=True)


def user_directory_path(instance, filename):
    return os.path.join('cifs', f'user_{instance.refcode.user.id}', f'{instance.refcode}.cif')


class CifFile(models.Model):
    '''User cif files.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='cif_file',
        on_delete=models.CASCADE
    )
    file = models.FileField(
        upload_to=user_directory_path,
        max_length=200,
        verbose_name='Cif file'
    )
    old_file_name = models.CharField(
        verbose_name='File name when upload',
        max_length=200,
        blank=True,
        null=True
    )


class Spacegroup(models.Model):
    '''Space group list.'''
    number = models.IntegerField(db_column='Number')
    system = models.IntegerField(choices=SYSTEMS, verbose_name='lattice system')
    name = models.CharField(verbose_name='space group', max_length=20)
    hall_name = models.CharField(verbose_name='Hall space group', max_length=100, blank=True, null=True)
    symops = models.TextField(verbose_name='Symmetry operations', blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        unique_together = (('name', 'hall_name'),)


class Cell(models.Model):
    '''Table with refcode-cell correspondence.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='cell',
        on_delete=models.CASCADE
    )
    centring = models.SmallIntegerField(choices=CENTRINGS, default=1)
    a = models.FloatField(db_column='LengthA', verbose_name='a')
    b = models.FloatField(db_column='LengthB', verbose_name='b')
    c = models.FloatField(db_column='LengthC', verbose_name='c')
    al = models.FloatField(db_column='AngleAlpha', verbose_name='alpha')
    be = models.FloatField(db_column='AngleBeta', verbose_name='betta')
    ga = models.FloatField(db_column='AngleGamma', verbose_name='gamma')
    a_err = models.CharField(db_column='LengthAErr', verbose_name='a', max_length=50, null=True)
    b_err = models.CharField(db_column='LengthBErr', verbose_name='b', max_length=50, null=True)
    c_err = models.CharField(db_column='LengthCErr', verbose_name='c', max_length=50, null=True)
    al_err = models.CharField(db_column='AngleAlphaErr', verbose_name='alpha', max_length=50, null=True)
    be_err = models.CharField(db_column='AngleBetaErr', verbose_name='betta', max_length=50, null=True)
    ga_err = models.CharField(db_column='AngleGammaErr', verbose_name='gamma', max_length=50, null=True)

    spacegroup = models.ForeignKey(
        Spacegroup,
        related_name='cell',
        on_delete=models.PROTECT,
        verbose_name='space group'
    )
    zvalue = models.FloatField(verbose_name='Z')
    zprime = models.FloatField(verbose_name="Z'", null=True, blank=True)


class NormalisedReducedCell(models.Model):
    '''Normalized Niggli cell.'''
    refcode = models.ForeignKey(
        StructureCode,
        related_name='normalised_reduced_cells',
        on_delete=models.CASCADE,
        verbose_name='Refcode'
    )
    a = models.FloatField(db_column='LengthA', verbose_name='a')
    b = models.FloatField(db_column='LengthB', verbose_name='b')
    c = models.FloatField(db_column='LengthC', verbose_name='c')
    al = models.FloatField(db_column='AngleAlpha', verbose_name='alpha')
    be = models.FloatField(db_column='AngleBeta', verbose_name='betta')
    ga = models.FloatField(db_column='AngleGamma', verbose_name='gamma')
    volume = models.FloatField(db_column='Volume', verbose_name='Volume')

    class Meta:
        unique_together = (('refcode', 'a', 'b', 'c', 'al', 'be', 'ga'),)


class ReducedCell(models.Model):
    '''Niggli cell.'''
    refcode = models.ForeignKey(
        StructureCode,
        related_name='reduced_cells',
        on_delete=models.CASCADE,
        verbose_name='Refcode'
    )
    a = models.FloatField(db_column='LengthA', verbose_name='a')
    b = models.FloatField(db_column='LengthB', verbose_name='b')
    c = models.FloatField(db_column='LengthC', verbose_name='c')
    al = models.FloatField(db_column='AngleAlpha', verbose_name='alpha')
    be = models.FloatField(db_column='AngleBeta', verbose_name='betta')
    ga = models.FloatField(db_column='AngleGamma', verbose_name='gamma')
    volume = models.FloatField(db_column='Volume', verbose_name='Volume')

    class Meta:
        unique_together = (('refcode', 'a', 'b', 'c', 'al', 'be', 'ga'),)


class CompoundName(models.Model):
    '''Compound names.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='name',
        on_delete=models.CASCADE)
    systematic_name = models.CharField(
        blank=True,
        null=True,
        max_length=1500,
        verbose_name='Compound name'
    )
    trivial_name = models.CharField(
        blank=True,
        null=True,
        max_length=1000,
        verbose_name='Trivial name'
    )

    def __str__(self):
        if self.trivial_name:
            return self.trivial_name
        return self.systematic_name


class ExperimentalInfo(models.Model):
    '''Information about experiment.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='experimental_info',
        on_delete=models.CASCADE)
    structure_determination_temperature = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Temperature',
        help_text='Temperature in Kelvin scale: 120.0'
    )
    measurement_reflns = models.PositiveIntegerField(
        blank=True,
        null=True,
        verbose_name='Number of measurement reflections'
    )
    measurement_theta_max = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Maximum value of measurement theta'
    )
    measurement_theta_min = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Minimum value of measurement theta'
    )
    calculated_density_value = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Calculated density'
    )
    absorpt_correction_type = models.CharField(
        blank=True,
        null=True,
        verbose_name='Absorption correction type',
        max_length=500,
        help_text='multi-scan'
    )
    measurement_device = models.CharField(
        blank=True,
        null=True,
        verbose_name='Measurement device',
        max_length=250,
        help_text='Bruker APEX-II CCD'
    )
    measurement_method = models.CharField(
        blank=True,
        null=True,
        verbose_name='Measurement method',
        max_length=250,
        help_text=r'\w and \f scans'
    )
    monochromator = models.CharField(
        blank=True,
        null=True,
        verbose_name='Monochromator',
        max_length=100,
        help_text='graphite'
    )
    radiation_type = models.CharField(
        blank=True,
        null=True,
        verbose_name='Radiation type',
        max_length=500,
        help_text=r'CuK\a'
    )
    wavelength = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Wavelength',
        help_text='1.54178'
    )
    radiation_source = models.CharField(
        blank=True,
        null=True,
        verbose_name='Radiation source',
        max_length=250,
        help_text='sealed tube'
    )

    class Meta:
        verbose_name_plural = 'Experimental info'


class RefinementInfo(models.Model):
    '''Information about structure refinment details.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='refinement_info',
        on_delete=models.CASCADE
    )
    r_factor = models.FloatField(
        blank=True,
        null=True,
        verbose_name='R1'
    )  # _refine_ls_R_factor_gt

    wR_factor = models.FloatField(
        blank=True,
        null=True,
        verbose_name='wR2'
    )  # _refine_ls_wR_factor_ref
    gof = models.FloatField(
        blank=True,
        null=True,
        verbose_name='GooF'
    )  # _refine_ls_goodness_of_fit_ref
    diff_density_max = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Maximum residual electron density'
    )
    diff_density_min = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Minimum residual electron density'
    )
    extinction_coef = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Extinction coefficient'
    )  # _refine_ls_extinction_coef
    flack = models.FloatField(
        blank=True,
        null=True,
        verbose_name='Flack'
    )  # _refine_ls_abs_structure_Flack

    class Meta:
        verbose_name_plural = 'Refinement info'


class CoordinatesBlock(models.Model):
    '''Block of coordinates and other structural information.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='coordinates',
        on_delete=models.CASCADE
    )
    coordinates = models.TextField(
        verbose_name='Coordinates'
    )
    smiles = models.TextField(
        verbose_name='Smiles',
        blank=True,
        null=True,
    )
    graph = models.TextField(
        verbose_name='Graph',
        help_text='string graph for structure search',
        blank=True,
        null=True,
        editable=False
    )


class Substructure1(models.Model):
    '''Substructure information.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='substructure1',
        on_delete=models.CASCADE
    )
    only_CHNO = models.BooleanField(verbose_name='Only CHNO atoms', default=False)
    no_C = models.BooleanField(verbose_name='No C atoms', default=False)
    hetero = models.BooleanField(verbose_name='Heteroatoms (B, Si, P, S, As, Se, Ge, Sb, Te)', default=False)
    d_Me4 = models.BooleanField(verbose_name='d-Metalls 4 group', default=False)
    d_Me5 = models.BooleanField(verbose_name='d-Metalls 5 group', default=False)
    d_Me67 = models.BooleanField(verbose_name='d-Metalls 6 group', default=False)
    f_Me = models.BooleanField(verbose_name='f-Metalls', default=False)
    AEMet = models.BooleanField(verbose_name='Alkali earth metalls', default=False)
    AMet = models.BooleanField(verbose_name='Alkali metalls', default=False)
    halogens = models.BooleanField(verbose_name='Halogen atoms', default=False)
    other_met = models.BooleanField(verbose_name='Other metalls (Al, Ga, In, Sn, Tl, Pb, Bi)', default=False)
    NO2 = models.BooleanField(verbose_name='N-(O)2 group', default=False)
    SO2 = models.BooleanField(verbose_name='S-(O)2 group', default=False)
    CS = models.BooleanField(verbose_name='C-S bond', default=False)
    # C_Met = models.BooleanField(verbose_name='C-Metall bond', default=False)
    # N_Met = models.BooleanField(verbose_name='N-Metall bond', default=False)

    class Meta:
        verbose_name_plural = 'Substructures'


class Substructure2(models.Model):
    '''Substructure information.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='substructure2',
        on_delete=models.CASCADE
    )
    Ph = models.BooleanField(verbose_name='Ph group', default=False)
    ring6 = models.BooleanField(verbose_name='6 atoms ring', default=False)
    ring5 = models.BooleanField(verbose_name='5 atoms ring', default=False)
    ring6N1 = models.BooleanField(verbose_name='6 atoms ring with one N', default=False)
    ring6N2 = models.BooleanField(verbose_name='6 atoms ring with two N', default=False)
    ring5N1 = models.BooleanField(verbose_name='5 atoms ring with one N', default=False)
    ring5N2 = models.BooleanField(verbose_name='5 atoms ring with two N', default=False)
    iPr = models.BooleanField(verbose_name='i-Pr group', default=False)
    tBu = models.BooleanField(verbose_name='t-Bu group', default=False)
    C3N = models.BooleanField(verbose_name='(C)3-N group', default=False)
    C2O = models.BooleanField(verbose_name='(C)2-O group', default=False)
    CNO = models.BooleanField(verbose_name='C-N-O group', default=False)
    C3P = models.BooleanField(verbose_name='(C)3-P group', default=False)
    CO2 = models.BooleanField(verbose_name='C-(O)2 group', default=False)
    # C_Hal = models.BooleanField(verbose_name='C-Halogen bond', default=False)

    class Meta:
        verbose_name_plural = 'Substructures'


class CrystalAndStructureInfo(models.Model):
    '''Information about crystal and structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='crystal_and_structure_info',
        on_delete=models.CASCADE
    )
    color = models.CharField(
        blank=True,
        null=True,
        verbose_name='Color',
        max_length=250
    )
    bioactivity = models.CharField(
        blank=True,
        null=True,
        verbose_name='Bioactivity',
        max_length=1000
    )
    crystal_shape = models.CharField(
        blank=True,
        null=True,
        verbose_name='Crystal shape',
        choices=SHAPES,
        max_length=200
    )
    phase_transitions = models.CharField(
        blank=True,
        null=True,
        verbose_name='Phase transitions',
        max_length=1000
    )
    polymorph = models.CharField(
        blank=True,
        null=True,
        verbose_name='Polymorph',
        max_length=500
    )
    melting_point = models.CharField(
        blank=True,
        null=True,
        verbose_name='Melting point',
        max_length=100
    )
    sensitivity = models.CharField(
        blank=True,
        null=True,
        verbose_name='Sensitivity',
        max_length=1000
    )
    pressure = models.CharField(
        blank=True,
        null=True,
        verbose_name='Pressure',
        max_length=100
    )
    disorder = models.TextField(db_column='Disorder', blank=True, null=True, verbose_name='Disorder')
    recrystallisation_solvent = models.CharField(
        blank=True,
        null=True,
        verbose_name='Recrystallisation solvent',
        max_length=500
    )
    size_max = models.FloatField(blank=True, null=True, verbose_name='Maximum size')
    size_min = models.FloatField(blank=True, null=True, verbose_name='Minimum size')
    size_mid = models.FloatField(blank=True, null=True, verbose_name='Middle size')

    class Meta:
        verbose_name_plural = 'Crystal and structure info'


class Formula(models.Model):
    '''Formula.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='formula',
        on_delete=models.CASCADE
    )
    formula_moiety = models.CharField(
        max_length=500,
        verbose_name='Moiety Formula',
        blank=True,
        null=True,
    )
    formula_sum = models.CharField(
        max_length=500,
        verbose_name='Sum Formula',
        blank=True,
        null=True,
    )


class ElementsSet1(models.Model):
    '''First set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_1',
        on_delete=models.CASCADE
    )
    H = models.FloatField(blank=True, null=True)
    C = models.FloatField(blank=True, null=True)
    N = models.FloatField(blank=True, null=True)
    O = models.FloatField(blank=True, null=True)
    Cl = models.FloatField(blank=True, null=True)
    F = models.FloatField(blank=True, null=True)
    S = models.FloatField(blank=True, null=True)
    P = models.FloatField(blank=True, null=True)
    D = models.FloatField(blank=True, null=True)
    T = models.FloatField(blank=True, null=True)


class ElementsSet2(models.Model):
    '''Second set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_2',
        on_delete=models.CASCADE
    )
    Pt = models.FloatField(blank=True, null=True)
    Zn = models.FloatField(blank=True, null=True)
    Mo = models.FloatField(blank=True, null=True)
    I = models.FloatField(blank=True, null=True)
    Ru = models.FloatField(blank=True, null=True)
    Ni = models.FloatField(blank=True, null=True)
    Co = models.FloatField(blank=True, null=True)
    Br = models.FloatField(blank=True, null=True)
    Si = models.FloatField(blank=True, null=True)
    Fe = models.FloatField(blank=True, null=True)
    B = models.FloatField(blank=True, null=True)
    Cu = models.FloatField(blank=True, null=True)
    Sn = models.FloatField(blank=True, null=True)
    W = models.FloatField(blank=True, null=True)
    Mn = models.FloatField(blank=True, null=True)
    Pd = models.FloatField(blank=True, null=True)


class ElementsSet3(models.Model):
    '''Third set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_3',
        on_delete=models.CASCADE
    )
    Os = models.FloatField(blank=True, null=True)
    Zr = models.FloatField(blank=True, null=True)
    V = models.FloatField(blank=True, null=True)
    Ir = models.FloatField(blank=True, null=True)
    Al = models.FloatField(blank=True, null=True)
    Au = models.FloatField(blank=True, null=True)
    K = models.FloatField(blank=True, null=True)
    Cd = models.FloatField(blank=True, null=True)
    Ti = models.FloatField(blank=True, null=True)
    Na = models.FloatField(blank=True, null=True)
    Li = models.FloatField(blank=True, null=True)
    Cr = models.FloatField(blank=True, null=True)
    Se = models.FloatField(blank=True, null=True)
    Re = models.FloatField(blank=True, null=True)
    Ag = models.FloatField(blank=True, null=True)
    Rh = models.FloatField(blank=True, null=True)


class ElementsSet4(models.Model):
    '''Fourth set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_4',
        on_delete=models.CASCADE
    )
    Nb = models.FloatField(blank=True, null=True)
    Bi = models.FloatField(blank=True, null=True)
    Nd = models.FloatField(blank=True, null=True)
    Yb = models.FloatField(blank=True, null=True)
    Sm = models.FloatField(blank=True, null=True)
    In = models.FloatField(blank=True, null=True)
    Mg = models.FloatField(blank=True, null=True)
    Pb = models.FloatField(blank=True, null=True)
    U = models.FloatField(blank=True, null=True)
    Ge = models.FloatField(blank=True, null=True)
    Te = models.FloatField(blank=True, null=True)
    Ga = models.FloatField(blank=True, null=True)
    Hg = models.FloatField(blank=True, null=True)
    As = models.FloatField(blank=True, null=True)
    Sb = models.FloatField(blank=True, null=True)
    La = models.FloatField(blank=True, null=True)


class ElementsSet5(models.Model):
    '''Fifth set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_5',
        on_delete=models.CASCADE
    )
    Dy = models.FloatField(blank=True, null=True)
    Sr = models.FloatField(blank=True, null=True)
    Rb = models.FloatField(blank=True, null=True)
    Hf = models.FloatField(blank=True, null=True)
    Tb = models.FloatField(blank=True, null=True)
    Pr = models.FloatField(blank=True, null=True)
    Ce = models.FloatField(blank=True, null=True)
    Er = models.FloatField(blank=True, null=True)
    Tl = models.FloatField(blank=True, null=True)
    Cs = models.FloatField(blank=True, null=True)
    Ba = models.FloatField(blank=True, null=True)
    Gd = models.FloatField(blank=True, null=True)
    Y = models.FloatField(blank=True, null=True)
    Ca = models.FloatField(blank=True, null=True)
    Eu = models.FloatField(blank=True, null=True)
    Ta = models.FloatField(blank=True, null=True)


class ElementsSet6(models.Model):
    '''Sixth set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_6',
        on_delete=models.CASCADE
    )
    Pm = models.FloatField(blank=True, null=True)
    Cm = models.FloatField(blank=True, null=True)
    Pa = models.FloatField(blank=True, null=True)
    Kr = models.FloatField(blank=True, null=True)
    Am = models.FloatField(blank=True, null=True)
    Ar = models.FloatField(blank=True, null=True)
    Xe = models.FloatField(blank=True, null=True)
    Pu = models.FloatField(blank=True, null=True)
    Np = models.FloatField(blank=True, null=True)
    Tm = models.FloatField(blank=True, null=True)
    Be = models.FloatField(blank=True, null=True)
    Th = models.FloatField(blank=True, null=True)
    Ho = models.FloatField(blank=True, null=True)
    Sc = models.FloatField(blank=True, null=True)
    Lu = models.FloatField(blank=True, null=True)
    Tc = models.FloatField(blank=True, null=True)


class ElementsSet7(models.Model):
    '''Seventh set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_7',
        on_delete=models.CASCADE
    )
    At = models.FloatField(blank=True, null=True)
    Es = models.FloatField(blank=True, null=True)
    Fm = models.FloatField(blank=True, null=True)
    Fr = models.FloatField(blank=True, null=True)
    He = models.FloatField(blank=True, null=True)
    Lr = models.FloatField(blank=True, null=True)
    Md = models.FloatField(blank=True, null=True)
    Ne = models.FloatField(blank=True, null=True)
    No = models.FloatField(blank=True, null=True)
    Po = models.FloatField(blank=True, null=True)
    Ra = models.FloatField(blank=True, null=True)
    Rn = models.FloatField(blank=True, null=True)
    Ac = models.FloatField(blank=True, null=True)
    Bk = models.FloatField(blank=True, null=True)
    Cf = models.FloatField(blank=True, null=True)
    Rf = models.FloatField(blank=True, null=True)


class ElementsSet8(models.Model):
    '''Eighth set of elements with quantity of each atom in structure.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='element_set_8',
        on_delete=models.CASCADE
    )
    Db = models.FloatField(blank=True, null=True)
    Sg = models.FloatField(blank=True, null=True)
    Bh = models.FloatField(blank=True, null=True)
    Hs = models.FloatField(blank=True, null=True)
    Mt = models.FloatField(blank=True, null=True)
    Ds = models.FloatField(blank=True, null=True)
    Rg = models.FloatField(blank=True, null=True)
    Cn = models.FloatField(blank=True, null=True)
    Nh = models.FloatField(blank=True, null=True)
    Fl = models.FloatField(blank=True, null=True)
    Mc = models.FloatField(blank=True, null=True)
    Lv = models.FloatField(blank=True, null=True)
    Ts = models.FloatField(blank=True, null=True)
    Og = models.FloatField(blank=True, null=True)


class Journal(models.Model):
    '''Information about journal.'''
    international_coden = models.CharField(
        max_length=10,
        verbose_name='Internations journal code',
        blank=True,
        null=True,
    )
    name = models.CharField(max_length=500, verbose_name='Journal')
    discontinued = models.BooleanField(verbose_name='Discontinued', blank=True, null=True,)
    fullname = models.CharField(max_length=1000, verbose_name='Journal')
    translated_name = models.CharField(max_length=500, blank=True, null=True, verbose_name='Journal')
    abbreviated_translated_name = models.CharField(
        max_length=500,
        blank=True,
        null=True,
        verbose_name='Journal'
    )

    def __str__(self):
        return self.name

    class Meta:
        unique_together = (('name', 'fullname'),)


class Publication(models.Model):
    '''Information about publication.'''
    journal = models.ForeignKey(
        Journal,
        related_name='publication',
        on_delete=models.CASCADE,
        verbose_name='Journal',
        blank=True,
        null=True
    )
    authors = models.ManyToManyField(
        Author,
        related_name='publications',
        verbose_name='Author',
        blank=True
    )
    page = models.CharField(blank=True, null=True, max_length=100, verbose_name='Page')
    volume = models.CharField(blank=True, null=True, max_length=100, verbose_name='Volume')
    year = models.PositiveIntegerField(verbose_name='Year')
    doi = models.CharField(blank=True, null=True, max_length=100, verbose_name='DOI', unique=True)

    class Meta:
        unique_together = (('journal', 'page', 'volume', 'year', 'doi'),)

    def __str__(self):
        return self.doi


class RefcodePublicationConnection(models.Model):
    '''Table of refcode-publication correspondence.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='publication',
        on_delete=models.CASCADE,
        verbose_name='Refcode'
    )
    publication = models.ForeignKey(
        Publication,
        related_name='refcode',
        on_delete=models.CASCADE,
        verbose_name='Publication',
        blank=True,
        null=True
    )


class Other(models.Model):
    '''Other informatrion about structures.'''
    refcode = models.OneToOneField(
        StructureCode,
        related_name='characteristics',
        on_delete=models.CASCADE
    )
    has_3d_structure = models.BooleanField(default=True, verbose_name='Has 3d structure')
    number_atoms_with_sites = models.IntegerField(verbose_name='Number of atoms with sites', null=True, blank=True)
    maximum_atomic_number = models.IntegerField(verbose_name='Maximum atomic number', null=True, blank=True)

    class Meta:
        verbose_name_plural = 'Other'
