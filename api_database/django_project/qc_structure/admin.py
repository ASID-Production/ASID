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

from .models import (QCStructureCode, QCCell, QCProgram,
                     QCCompoundName, QCFormula, QCReducedCell,
                     QCCoordinatesBlock, QCProperties,
                     QCSubstructure1, QCSubstructure2)
from django.contrib import admin


@admin.register(QCStructureCode)
class QCStructureCodeAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'user')
    search_fields = ('refcode',)
    empty_value_display = '-empty-'


@admin.register(QCCell)
class QCCellAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'spacegroup', 'a', 'b', 'c', 'al', 'be', 'ga')
    search_fields = ('refcode__qc_refcode__startswith', 'spacegroup__name')
    list_filter = ('spacegroup',)
    empty_value_display = '-empty-'


@admin.register(QCReducedCell)
class QCReducedCellAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'a', 'b', 'c', 'al', 'be', 'ga', 'volume')
    search_fields = ('refcode__qc_refcode__startswith',)
    empty_value_display = '-empty-'


@admin.register(QCCompoundName)
class QCCompoundNamesAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'systematic_name', 'trivial_name')
    search_fields = ('refcode__qc_refcode__startswith', 'systematic_name', 'trivial_name')
    empty_value_display = '-empty-'


@admin.register(QCCoordinatesBlock)
class QCCoordinatesBlockAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'smiles', 'graph')
    search_fields = ('refcode__qc_refcode__startswith', 'smiles')
    empty_value_display = '-empty-'


@admin.register(QCSubstructure1)
class QCSubstructure1Admin(admin.ModelAdmin):
    list_display = (
        'id', 'refcode', 'only_CHNO', 'no_C',
        'hetero', 'd_Me4', 'd_Me5', 'd_Me67',
        'f_Me', 'AEMet', 'AMet', 'halogens',
        'other_met', 'NO2', 'SO2', 'CS'
    )
    search_fields = ('refcode__qc_refcode__startswith', )
    empty_value_display = '-empty-'


@admin.register(QCSubstructure2)
class QCSubstructure2Admin(admin.ModelAdmin):
    list_display = (
        'id', 'refcode', 'Ph', 'ring6',
        'ring5', 'ring6N1', 'ring6N2', 'ring5N1',
        'ring5N2', 'iPr', 'tBu', 'C3N',
        'C2O', 'CNO', 'C3P', 'CO2'
    )
    search_fields = ('refcode__qc_refcode__startswith', )
    empty_value_display = '-empty-'


@admin.register(QCFormula)
class QCFormulaAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'formula_moiety', 'formula_sum')
    search_fields = ('refcode__qc_refcode__startswith', 'formula_moiety', 'formula_sum')
    empty_value_display = '-empty-'


@admin.register(QCProgram)
class QCProgramAdmin(admin.ModelAdmin):
    list_display = (
        'id', 'refcode', 'orca', 'vasp', 'qchem', 'gaussian',
        'nwchem', 'abinit', 'crystal', 'gamess', 'mopac', 'quantum_espresso'
    )
    search_fields = ('refcode__qc_refcode__startswith',)
    list_filter = (
        'orca', 'vasp', 'qchem', 'gaussian', 'nwchem',
        'abinit', 'crystal', 'gamess', 'mopac', 'quantum_espresso'
    )
    empty_value_display = '-empty-'


@admin.register(QCProperties)
class QCPropertiesAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'energy', 'calculated_density')
    search_fields = ('refcode__qc_refcode__startswith',)
    empty_value_display = '-empty-'
