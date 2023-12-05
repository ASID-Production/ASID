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
from django.shortcuts import get_object_or_404
from structure.models import (StructureCode, Author, Spacegroup, Cell,
                              CompoundName, Formula, Publication,
                              RefcodePublicationConnection, ExperimentalInfo,
                              ReducedCell, ExperimentalInfo, RefinementInfo,
                              CoordinatesBlock, CrystalAndStructureInfo,
                              CifFile)
from djoser.serializers import UserSerializer, UserCreateSerializer
from django.contrib.auth import get_user_model
from .fields import NodesListField, EdgesListField

User = get_user_model()


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
    formula = SerializerMethodField(read_only=True)
    temperature = SerializerMethodField(read_only=True)

    class Meta:
        model = StructureCode
        fields = (
            'id', 'refcode', 'CCDC_number', 'formula', 'temperature'
        )

    def get_formula(self, obj):
        formula = get_object_or_404(Formula, refcode=obj)
        if formula.formula_moiety:
            return formula.formula_moiety
        return formula.formula_sum

    def get_temperature(self, obj):
        exp_info = get_object_or_404(ExperimentalInfo, refcode=obj)
        return exp_info.structure_determination_temperature


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
        publication = get_object_or_404(RefcodePublicationConnection, refcode=obj).publication
        return PublicationSerializer(publication).data


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
