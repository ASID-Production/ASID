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

from structure.models import (Spacegroup, AbstractCell, AbstractReducedCell,
                              AbstractCompoundName, AbstractCoordinatesBlock,
                              AbstractSubstructure1, AbstractSubstructure2,
                              AbstractFormula, AbstractElementsManager,
                              ElementsSet1, ElementsSet2,
                              ElementsSet3, ElementsSet4, ElementsSet5,
                              ElementsSet6, ElementsSet7, ElementsSet8)
from django.db import models
from django.contrib.auth import get_user_model
import os

User = get_user_model()


class QCStructureCode(models.Model):
    '''Table with codes.'''
    refcode = models.CharField(verbose_name='Refcode', max_length=17, unique=True)
    user = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        related_name='qc_refcode',
        null=True,
        blank=True,
        verbose_name='Owner'
    )

    def __str__(self):
        return self.refcode


class QCProgram(models.Model):
    '''Calculation program name.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_prog',
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


def vasp_user_directory_path(instance, filename):
    return os.path.join('vasp', f'user_{instance.refcode.user.id}', f'{instance.refcode}_vasprun.xml')


class VaspFile(models.Model):
    '''VASP .xml files.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='vasp_file',
        on_delete=models.CASCADE
    )
    file = models.FileField(
        upload_to=vasp_user_directory_path,
        max_length=200,
        verbose_name='VASP .xml file'
    )


class QCCell(AbstractCell):
    '''Table with refcode-cell correspondence.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_cell',
        on_delete=models.CASCADE
    )
    spacegroup = models.ForeignKey(
        Spacegroup,
        related_name='qc_cell',
        on_delete=models.PROTECT,
        verbose_name='space group'
    )


class QCReducedCell(AbstractReducedCell):
    '''Niggli cell.'''
    refcode = models.ForeignKey(
        QCStructureCode,
        related_name='qc_reduced_cells',
        on_delete=models.CASCADE,
        verbose_name='Refcode'
    )

    class Meta:
        unique_together = (('refcode', 'a', 'b', 'c', 'al', 'be', 'ga'),)


class QCCompoundName(AbstractCompoundName):
    '''Compound names.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_name',
        on_delete=models.CASCADE)


class QCCoordinatesBlock(AbstractCoordinatesBlock):
    '''Block of coordinates and other structural information.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_coordinates',
        on_delete=models.CASCADE
    )

    class Meta:
        verbose_name_plural = 'QCCoordinates'


class QCSubstructure1(AbstractSubstructure1):
    '''Substructure information.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_substructure1',
        on_delete=models.CASCADE
    )


class QCSubstructure2(AbstractSubstructure2):
    '''Substructure information.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_substructure2',
        on_delete=models.CASCADE
    )


class QCFormula(AbstractFormula):
    '''Formula.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_formula',
        on_delete=models.CASCADE
    )


class QCElementsManager(AbstractElementsManager):
    '''Refcode to ElementsSets relations.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='elements',
        on_delete=models.PROTECT
    )
    element_set_1 = models.ForeignKey(
        ElementsSet1,
        related_name='qc_element_set_1',
        on_delete=models.PROTECT,
        verbose_name='element set 1',
        null=True
    )
    element_set_2 = models.ForeignKey(
        ElementsSet2,
        related_name='qc_element_set_2',
        on_delete=models.PROTECT,
        verbose_name='element set 2',
        null=True
    )
    element_set_3 = models.ForeignKey(
        ElementsSet3,
        related_name='qc_element_set_3',
        on_delete=models.PROTECT,
        verbose_name='element set 3',
        null=True
    )
    element_set_4 = models.ForeignKey(
        ElementsSet4,
        related_name='qc_element_set_4',
        on_delete=models.PROTECT,
        verbose_name='element set 4',
        null=True
    )
    element_set_5 = models.ForeignKey(
        ElementsSet5,
        related_name='qc_element_set_5',
        on_delete=models.PROTECT,
        verbose_name='element set 5',
        null=True
    )
    element_set_6 = models.ForeignKey(
        ElementsSet6,
        related_name='qc_element_set_6',
        on_delete=models.PROTECT,
        verbose_name='element set 6',
        null=True
    )
    element_set_7 = models.ForeignKey(
        ElementsSet7,
        related_name='qc_element_set_7',
        on_delete=models.PROTECT,
        verbose_name='element set 7',
        null=True
    )
    element_set_8 = models.ForeignKey(
        ElementsSet8,
        related_name='qc_element_set_8',
        on_delete=models.PROTECT,
        verbose_name='element set 8',
        null=True
    )


class QCProperties(models.Model):
    '''Quantum chemistry calculation properties.'''
    refcode = models.OneToOneField(
        QCStructureCode,
        related_name='qc_properties',
        on_delete=models.CASCADE
    )
    energy = models.FloatField(verbose_name='Final energy in eV')
