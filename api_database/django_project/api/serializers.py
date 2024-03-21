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

from rest_framework import serializers
from rest_framework.fields import SerializerMethodField
from structure.models import (StructureCode, Author, Spacegroup, Cell,
                              CompoundName, Formula, Publication,
                              RefcodePublicationConnection, ExperimentalInfo,
                              ReducedCell, ExperimentalInfo, RefinementInfo,
                              CoordinatesBlock, CrystalAndStructureInfo,
                              CifFile)
from qc_structure.models import (QCStructureCode, QCCell, QCCompoundName, QCFormula,
                                 QCReducedCell, QCCoordinatesBlock, QCProgram,
                                 QCProperties, VaspFile)
from djoser.serializers import UserSerializer, UserCreateSerializer
from django.contrib.auth import get_user_model
from .fields import NodesListField, EdgesListField

User = get_user_model()


#########################################################################
#                    Structure Serializers                            #
#########################################################################


class AuthorSerializer(serializers.ModelSerializer):
    class Meta:
        model = Author
        fields = ('id', 'family_name', 'initials')


class SpacegroupSerializer(serializers.ModelSerializer):
    system = SerializerMethodField(read_only=True)

    class Meta:
        model = Spacegroup
        fields = ('id', 'number', 'system', 'name', 'hall_name', 'symops')

    def get_system(self, obj):
        return obj.get_system_display()


class CellSerializer(serializers.ModelSerializer):
    spacegroup = SpacegroupSerializer(read_only=True)
    centring = SerializerMethodField(read_only=True)

    class Meta:
        model = Cell
        fields = ('id', 'spacegroup', 'centring', 'a', 'b', 'c', 'al', 'be', 'ga', 'zvalue')

    def get_centring(self, obj):
        return obj.get_centring_display()


class ReducedCellSerializer(serializers.ModelSerializer):
    class Meta:
        model = ReducedCell
        fields = ('id', 'a', 'b', 'c', 'al', 'be', 'ga', 'volume')


class CompoundNameSerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundName
        fields = ('id', 'systematic_name', 'trivial_name')


class FormulaSerializer(serializers.ModelSerializer):
    class Meta:
        model = Formula
        fields = ('id', 'formula_moiety', 'formula_sum')


class PublicationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Publication
        fields = ('id', 'journal', 'page', 'volume', 'year', 'doi')


class ExperimentalInfoSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExperimentalInfo
        fields = ('id', 'structure_determination_temperature', 'calculated_density_value',
                  'measurement_device', 'radiation_type', 'wavelength')


class RefinementInfoSerializer(serializers.ModelSerializer):
    class Meta:
        model = RefinementInfo
        fields = ('id', 'r_factor', 'wR_factor', 'gof')


class CoordinatesSerializer(serializers.ModelSerializer):
    class Meta:
        model = CoordinatesBlock
        fields = ('id', 'coordinates', 'smiles')


class CrystalInfoSerializer(serializers.ModelSerializer):
    class Meta:
        model = CrystalAndStructureInfo
        fields = ('id', 'color', 'bioactivity', 'crystal_shape', 'phase_transitions',
                  'polymorph', 'melting_point', 'sensitivity', 'pressure',
                  'disorder', 'recrystallisation_solvent', 'size_max', 'size_min',
                  'size_mid')


class RefcodeShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = StructureCode
        fields = ('id', 'refcode')


class RefcodeFullSerializer(serializers.ModelSerializer):
    formula = FormulaSerializer(read_only=True)
    authors = AuthorSerializer(many=True, read_only=True)
    publication = SerializerMethodField(read_only=True)
    cell = CellSerializer(read_only=True)
    reduced_cells = ReducedCellSerializer(many=True, read_only=True)
    compound_name = CompoundNameSerializer(read_only=True, source='name')
    experiment = ExperimentalInfoSerializer(read_only=True, source='experimental_info')
    refinement_info = RefinementInfoSerializer(read_only=True)
    coordinates = CoordinatesSerializer(read_only=True)
    crystal_info = CrystalInfoSerializer(read_only=True, source='crystal_and_structure_info')

    class Meta:
        model = StructureCode
        fields = (
            'id', 'refcode', 'CCDC_number', 'cell',
            'reduced_cells', 'compound_name', 'formula', 'authors',
            'publication', 'experiment', 'refinement_info', 'coordinates',
            'crystal_info'
        )

    def get_publication(self, obj):
        try:
            publication = RefcodePublicationConnection.objects.get(refcode=obj).publication
            return PublicationSerializer(publication).data
        except Exception:
            pass


class CifUploadSerializer(serializers.ModelSerializer):
    class Meta:
        model = CifFile
        fields = ('file', 'refcode')


class CustomUserCreateSerializer(UserCreateSerializer):
    class Meta(UserCreateSerializer.Meta):
        fields = (
            'email', 'id', 'username', 'first_name', 'last_name', 'password'
        )


class CustomUserSerializer(UserSerializer):
    class Meta:
        model = User
        fields = (
            'email',
            'id',
            'username',
            'first_name',
            'last_name',
        )


class SearchSerializer(serializers.Serializer):
    search_type = serializers.ChoiceField(
        choices=['exact', 'substructure'],
        required=True,
        error_messages={"invalid_choice": "Unsupported value: use 'exact' or 'substructure' keywords"}
    )
    # nodes input: ['(1, {"type": "C", "Hnum": "2",})', '(2, {"type": "N", "Hnum": "1", })', ...]
    nodes = NodesListField(child=serializers.CharField(max_length=100), allow_empty=False, required=True)
    # edges input: ['(1, 2, {"distance": "1.250",})', '(2, 3, {"distance": "1.340",})', ...]
    edges = EdgesListField(child=serializers.CharField(max_length=100), allow_empty=True)
    chunk_size = serializers.IntegerField(required=False, default=0, min_value=0)
    iter_num = serializers.IntegerField(required=False, default=0, min_value=0)


#########################################################################
#                    QCStructure Serializers                            #
#########################################################################


def get_fields_list(model):
    return [field.name for field in model._meta.get_fields()]


class QCCellSerializer(serializers.ModelSerializer):
    spacegroup = SpacegroupSerializer(read_only=True)
    centring = SerializerMethodField(read_only=True)

    class Meta:
        model = QCCell
        fields = ('id', 'spacegroup', 'centring', 'a', 'b', 'c', 'al', 'be', 'ga', 'zvalue')

    def get_centring(self, obj):
        return obj.get_centring_display()


class QCReducedCellSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCReducedCell
        fields = ('id', 'a', 'b', 'c', 'al', 'be', 'ga', 'volume')


class QCCompoundNameSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCCompoundName
        fields = ('id', 'systematic_name', 'trivial_name')


class QCFormulaSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCFormula
        fields = ('id', 'formula_moiety', 'formula_sum')


class QCCoordinatesSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCCoordinatesBlock
        fields = ('id', 'coordinates', 'smiles', 'graph')


class QCProgramSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCProgram
        fields = ('id', 'orca', 'vasp', 'qchem', 'gaussian',
                  'nwchem', 'abinit', 'crystal', 'gamess',
                  'mopac', 'quantum_espresso')


class QCPropertiesSerializer(serializers.ModelSerializer):
    energy = SerializerMethodField(read_only=True)

    class Meta:
        model = QCProperties
        fields = ('id', 'energy')

    def get_energy(self, obj):
        try:
            energy = obj.energy
            if energy:
                return str(energy) + ' eV'
        except Exception:
            pass


class QCRefcodeFullSerializer(serializers.ModelSerializer):
    formula = QCFormulaSerializer(read_only=True, source='qc_formula')
    cell = QCCellSerializer(read_only=True, source='qc_cell')
    reduced_cells = QCReducedCellSerializer(many=True, read_only=True, source='qc_reduced_cells')
    compound_name = QCCompoundNameSerializer(read_only=True, source='qc_name')
    coordinates = QCCoordinatesSerializer(read_only=True, source='qc_coordinates')
    properties = QCPropertiesSerializer(read_only=True, source='qc_properties')
    programs = QCProgramSerializer(read_only=True, source='qc_prog')

    class Meta:
        model = QCStructureCode
        fields = (
            'id', 'refcode', 'cell', 'reduced_cells', 'compound_name',
            'formula', 'coordinates', 'properties', 'programs'
        )


class QCRefcodeShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = QCStructureCode
        fields = ('id', 'refcode')


class VaspUploadSerializer(serializers.ModelSerializer):
    systematic_name = serializers.CharField(max_length=200, required=False, default='')
    trivial_name = serializers.CharField(max_length=200, required=False, default='')

    class Meta:
        model = VaspFile
        fields = ('file', 'systematic_name', 'trivial_name')
