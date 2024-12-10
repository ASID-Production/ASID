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

from .models import (StructureCode, Author, Spacegroup, Cell,
                     CompoundName, Formula, Publication,
                     NormalisedReducedCell, ReducedCell,
                     ExperimentalInfo, RefinementInfo,
                     CoordinatesBlock, CrystalAndStructureInfo,
                     Journal, RefcodePublicationConnection, Other, InChI)
from django.contrib import admin


@admin.register(Author)
class AuthorAdmin(admin.ModelAdmin):
    list_display = ('id', 'family_name', 'initials')
    search_fields = ('family_name',)
    empty_value_display = '-empty-'


@admin.register(StructureCode)
class StructureCodeAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'CCDC_number', 'user')
    search_fields = ('refcode', 'CCDC_number')
    empty_value_display = '-empty-'


@admin.register(Spacegroup)
class SpacegroupAdmin(admin.ModelAdmin):
    list_display = ('id', 'number', 'system', 'name', 'hall_name')
    search_fields = ('name', 'number')
    list_filter = ('system',)
    empty_value_display = '-empty-'


@admin.register(Cell)
class CellAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'spacegroup', 'a', 'b', 'c', 'al', 'be', 'ga')
    search_fields = ('refcode__refcode__startswith', 'spacegroup__name')
    list_filter = ('spacegroup',)
    empty_value_display = '-empty-'


@admin.register(NormalisedReducedCell)
class NormalisedReducedCellAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'a', 'b', 'c', 'al', 'be', 'ga', 'volume')
    search_fields = ('refcode__refcode__startswith',)
    empty_value_display = '-empty-'


@admin.register(ReducedCell)
class ReducedCellAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'a', 'b', 'c', 'al', 'be', 'ga', 'volume')
    search_fields = ('refcode__refcode__startswith',)
    empty_value_display = '-empty-'


@admin.register(CompoundName)
class CompoundNamesAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'systematic_name', 'trivial_name')
    search_fields = ('refcode__refcode__startswith', 'systematic_name', 'trivial_name')
    empty_value_display = '-empty-'


@admin.register(ExperimentalInfo)
class ExperimentalInfoAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'structure_determination_temperature',
                    'calculated_density_value', 'monochromator', 'radiation_type', 'wavelength')
    search_fields = ('refcode__refcode__startswith', 'structure_determination_temperature', 'monochromator', 'radiation_type', 'wavelength')
    empty_value_display = '-empty-'


@admin.register(RefinementInfo)
class RefinementInfoAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'r_factor', 'gof', 'diff_density_max', 'diff_density_min')
    search_fields = ('refcode__refcode__startswith',)
    empty_value_display = '-empty-'


@admin.register(CoordinatesBlock)
class CoordinatesBlockAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'smiles')
    search_fields = ('refcode__refcode__startswith', 'smiles')
    empty_value_display = '-empty-'


@admin.register(CrystalAndStructureInfo)
class CrystalAndStructureInfoAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'color', 'crystal_shape', 'disorder',)
    search_fields = ('refcode__refcode__startswith', 'color', 'crystal_shape')
    empty_value_display = '-empty-'


@admin.register(Formula)
class FormulaAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'formula_moiety', 'formula_sum')
    search_fields = ('refcode__refcode__startswith', 'formula_moiety', 'formula_sum')
    empty_value_display = '-empty-'


@admin.register(Journal)
class JournalAdmin(admin.ModelAdmin):
    list_display = ('id', 'international_coden', 'name', 'translated_name', 'discontinued')
    search_fields = ('international_coden', 'name', 'translated_name')
    list_filter = ('discontinued',)
    empty_value_display = '-empty-'


@admin.register(Publication)
class PublicationsAdmin(admin.ModelAdmin):
    list_display = ('id', 'journal', 'page', 'volume', 'year', 'doi')
    search_fields = ('doi', 'refcode__refcode__refcode__startswith', 'journal__name__startswith')
    empty_value_display = '-empty-'


@admin.register(RefcodePublicationConnection)
class RefcodePublicationConnectionAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode')
    search_fields = ('publication__doi',)
    empty_value_display = '-empty-'


@admin.register(Other)
class OtherAdmin(admin.ModelAdmin):
    list_display = ('id', 'refcode', 'has_3d_structure', 'number_atoms_with_sites', 'maximum_atomic_number')
    search_fields = ('refcode__refcode__startswith',)
    list_filter = ('has_3d_structure',)
    empty_value_display = '-empty-'


@admin.register(InChI)
class InChIAdmin(admin.ModelAdmin):
    list_display = (
        'id', 'refcode', 'version', 'formula', 'connectivity',
        'hydrogens', 'q_charge', 'p_charge', 'b_stereo',
        't_stereo', 'm_stereo', 's_stereo', 'i_isotopic'
    )
    search_fields = ('refcode__refcode__startswith', 'formula')
    empty_value_display = '-empty-'
