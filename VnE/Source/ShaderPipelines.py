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
# ******************************************************************************************
#  Author:      Alexander A. Korlyukov (head)
#  ORCID:       0000-0002-5600-9886
#  Author:      Alexander D. Volodin (author of cpplib)
#  ORCID:       0000-0002-3522-9193
#  Author:      Petr A. Buikin (author of api_database)
#  ORCID:       0000-0001-9243-9915
#  Author:      Alexander R. Romanenko (author of VnE)
#  ORCID:       0009-0003-5298-6836
#
# ******************************************************************************************


from OpenGL.GL import *
import numpy as np
from abc import ABC, abstractmethod
from typing import List
import ctypes
import os.path as opath

from .ShaderDataObjects import ShaderData, ShaderDataCreator, ShaderDataText

import debug


class Uniform(ABC):

    setter_func = {np.float32: glUniform1f}

    def __init__(self, dtype, name, program):
        self.data = None
        self.dtype = dtype
        self.program = program
        self.setter = type(self).setter_func.get(dtype, None)
        self.location = glGetUniformLocation(program, name)

    def set(self, value):
        if self.setter is not None:
            glUseProgram(self.program)
            self.setter(*(self.location, value))
            self.data = value
            return True
        else:
            return False


class ShaderProgramsCreator:

    def __init__(self):
        self.programs_list = {}

    def createProgram(self, source, shader):
        if type(source) is str:
            source = source.encode('utf-8')
        elif type(source) is bytes:
            pass
        else:
            raise TypeError
        program = self.programs_list.get(source, None)
        if program:
            return (program, shader)
        else:
            '''source_line_p = ctypes.c_char_p(source)
            source_line_p_p = ctypes.cast(ctypes.addressof(source_line_p), ctypes.POINTER((ctypes.POINTER(ctypes.c_char))))
            try:
                program = glCreateShaderProgramv(shader, 1, source_line_p_p)
            except Exception as e:
                print(str(e).encode('utf-8').decode('unicode_escape'))
                raise SystemExit()
            if glGetProgramiv(program, GL_LINK_STATUS) == 0:
                log = glGetProgramInfoLog(program).decode('utf-8')
                print('{:-^30}'.format('START'))
                print(log, end='')
                print('{:-^30}'.format('END'))
                raise SystemExit()
            self.programs_list[source] = (program, shader)
            return (program, shader)'''
            shader_id = glCreateShader(shader)
            if shader_id:
                glShaderSource(shader_id, source)
                glCompileShader(shader_id)
                program = glCreateProgram()
                if program:
                    compiled = glGetShaderiv(shader_id, GL_COMPILE_STATUS)
                    glProgramParameteri(program, GL_PROGRAM_SEPARABLE, GL_TRUE)
                    if compiled:
                        glAttachShader(program, shader_id)
                        glLinkProgram(program)
                        glDetachShader(program, shader_id)
                    else:
                        log = glGetShaderInfoLog(shader_id)
                        glDeleteShader(shader_id)
                        raise SystemExit('Failed to compile shader')
                else:
                    raise SystemExit('Failed to create program')
                glDeleteShader(shader_id)
                return (program, shader)
            else:
                raise SystemExit('Failed to create shader')

    def getProgramById(self, id):
        prog_list = self.programs_list.values()
        for prog in prog_list:
            if prog[1] == id:
                return prog
        return None


class aShaderPipeline(ABC):
    shader_bit = {GL_VERTEX_SHADER: GL_VERTEX_SHADER_BIT,
                  GL_TESS_CONTROL_SHADER: GL_TESS_CONTROL_SHADER_BIT,
                  GL_TESS_EVALUATION_SHADER: GL_TESS_EVALUATION_SHADER_BIT,
                  GL_FRAGMENT_SHADER: GL_FRAGMENT_SHADER_BIT,
                  GL_COMPUTE_SHADER: GL_COMPUTE_SHADER_BIT}

    @abstractmethod
    def __init__(self):
        self.shader_data: List[ShaderData]
        self.VAOFormat: List[List[int, str]]
        self.vertex_source: str
        self.fragment_source: str
        self.draw_mode: int
        self.uniforms: dict = {}

    def add_shader_data(self, shader_data_inst=None):
        if shader_data_inst:
            if isinstance(shader_data_inst, ShaderData):
                self.shader_data.append(shader_data_inst)
        shader_data_inst = ShaderDataCreator.createShaderData(self.VAOFormat)
        self.shader_data.append(shader_data_inst)
        return shader_data_inst

    @abstractmethod
    def draw(self):
        pass

    def getInfo(self):
        return self.VAOFormat

    def changeUniform(self, property, value):
        if property in self.uniforms:
            self.uniforms[property].set(value)

    @abstractmethod
    def changeShaderProgram(self):
        return


class BallsShaderPipeline(aShaderPipeline):

    shader_bit = {GL_VERTEX_SHADER: GL_VERTEX_SHADER_BIT,
                  GL_TESS_CONTROL_SHADER: GL_TESS_CONTROL_SHADER_BIT,
                  GL_TESS_EVALUATION_SHADER: GL_TESS_EVALUATION_SHADER_BIT,
                  GL_FRAGMENT_SHADER: GL_FRAGMENT_SHADER_BIT,
                  GL_COMPUTE_SHADER: GL_COMPUTE_SHADER_BIT}

    def __init__(self):
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (4, np.float32), (1, np.float32), (1, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/base.vert', 'r').read()],
                         [GL_TESS_CONTROL_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/tesc/balls.tesc', 'r').read()],
                         [GL_TESS_EVALUATION_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/tese/balls.tese', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/base.frag', 'r').read()]]

        self.pipeline = glGenProgramPipelines(1)
        glBindProgramPipeline(self.pipeline)

        for i, program in enumerate(self.programs):
            self.programs[i][1] = SHADER_PROGRAM_CREATOR.createProgram(program[2], program[0])[0]
            glUseProgramStages(self.pipeline, BallsShaderPipeline.shader_bit[self.programs[i][0]], self.programs[i][1])

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        glPatchParameteri(GL_PATCH_VERTICES, 1)
        for shader_data in self.shader_data:
            shader_data.draw(GL_PATCHES)

    def changeShaderProgram(self, shader_bit=None, source=None, id=None):
        if id:
            prog = SHADER_PROGRAM_CREATOR.getProgramById(id)
            if prog:
                glUseProgramStages(self.pipeline, prog[1], prog[0])

        if source and shader_bit:
            prog = SHADER_PROGRAM_CREATOR.createProgram(source, shader_bit)
            glUseProgramStages(self.pipeline, prog[1], prog[0])


class BondShaderPipeline(BallsShaderPipeline):

    def __init__(self):
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (4, np.float32), (1, np.float32), (1, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/base.vert', 'r').read()],
                         [GL_TESS_CONTROL_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/tesc/bonds.tesc', 'r').read()],
                         [GL_TESS_EVALUATION_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/tese/bonds.tese', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/base.frag', 'r').read()]]

        self.pipeline = glGenProgramPipelines(1)
        glBindProgramPipeline(self.pipeline)

        for i, program in enumerate(self.programs):
            self.programs[i][1] = SHADER_PROGRAM_CREATOR.createProgram(program[2], program[0])[0]
            glUseProgramStages(self.pipeline, BallsShaderPipeline.shader_bit[self.programs[i][0]], self.programs[i][1])

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        glPatchParameteri(GL_PATCH_VERTICES, 2)
        for shader_data in self.shader_data:
            shader_data.draw(GL_PATCHES)


class TextShaderPipeline(aShaderPipeline):

    def __init__(self):
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (2, np.float32), (2, np.float32), (2, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/text.vert', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/text.frag', 'r').read()]]
        self.textures = []

        self.pipeline = glGenProgramPipelines(1)

        glBindProgramPipeline(self.pipeline)

        for i, program in enumerate(self.programs):
            self.programs[i][1] = SHADER_PROGRAM_CREATOR.createProgram(program[2], program[0])[0]
            glUseProgramStages(self.pipeline, BallsShaderPipeline.shader_bit[self.programs[i][0]], self.programs[i][1])

        self.uniforms = {'const_scale': Uniform(np.float32, 'const_scale', self.programs[0][1])}
        self.uniforms['const_scale'].set(ctypes.c_float(150.0))
        glUseProgram(0)

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        for shader_data in self.shader_data:
            shader_data.draw(GL_TRIANGLES, self.textures)

    def add_shader_data(self, shader_data_inst=None):
        if shader_data_inst:
            if isinstance(shader_data_inst, ShaderData):
                self.shader_data.append(shader_data_inst)
        shader_data_inst = ShaderDataCreator.createShaderData(self.VAOFormat, shader_data=ShaderDataText)
        self.shader_data.append(shader_data_inst)
        return shader_data_inst

    def changeShaderProgram(self):
        super().changeShaderProgram()


class LinesShaderPipeline(aShaderPipeline):
    def __init__(self):
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (4, np.float32), (1, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/lines.vert', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/lines.frag', 'r').read()]]

        self.pipeline = glGenProgramPipelines(1)
        glBindProgramPipeline(self.pipeline)

        for i, program in enumerate(self.programs):
            self.programs[i][1] = SHADER_PROGRAM_CREATOR.createProgram(program[2], program[0])[0]
            glUseProgramStages(self.pipeline, BallsShaderPipeline.shader_bit[self.programs[i][0]], self.programs[i][1])

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        glDisable(GL_LINE_STIPPLE)
        for shader_data in self.shader_data:
            shader_data.draw(GL_LINES)

    def changeShaderProgram(self):
        pass


class DashedLineShaderPipeline(LinesShaderPipeline):

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        glEnable(GL_LINE_STIPPLE)
        glLineStipple(10, 0xAAAA)
        for shader_data in self.shader_data:
            shader_data.draw(GL_LINES)
        pass


class PlaneShaderPipeline(aShaderPipeline):
    def __init__(self):
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (4, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/plane.vert', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/base.frag', 'r').read()]]

        self.pipeline = glGenProgramPipelines(1)
        glBindProgramPipeline(self.pipeline)

        for i, program in enumerate(self.programs):
            self.programs[i][1] = SHADER_PROGRAM_CREATOR.createProgram(program[2], program[0])[0]
            glUseProgramStages(self.pipeline, super().shader_bit[self.programs[i][0]], self.programs[i][1])

    def draw(self):
        glBindProgramPipeline(self.pipeline)
        for shader_data in self.shader_data:
            shader_data.draw(GL_TRIANGLES)

    def changeShaderProgram(self, shader_bit=None, source=None, id=None):
        if id:
            prog = SHADER_PROGRAM_CREATOR.getProgramById(id)
            if prog:
                glUseProgramStages(self.pipeline, prog[1], prog[0])

        if source and shader_bit:
            prog = SHADER_PROGRAM_CREATOR.createProgram(source, shader_bit)
            glUseProgramStages(self.pipeline, prog[1], prog[0])


class TestShader(PlaneShaderPipeline):

    def __init__(self):
        from OpenGL.GL.shaders import compileProgram, compileShader
        self.shader_data = []
        self.VAOFormat = [(3, np.float32), (4, np.float32)]
        self.programs = [[GL_VERTEX_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/vert/plane.vert', 'r').read()],
                         [GL_FRAGMENT_SHADER, None, open(f'{opath.dirname(__file__)}/shaders/frag/base.frag', 'r').read()]]
        try:
            vert = compileShader(self.programs[0][2], GL_VERTEX_SHADER)
            frag = compileShader(self.programs[1][2], GL_FRAGMENT_SHADER)
            self.id = compileProgram(vert, frag, validate=False)
        except Exception as e:
            print(str(e).encode('utf-8').decode('unicode_escape'))
            raise SystemExit()

    def draw(self):
        glUseProgram(self.id)
        for shader_data in self.shader_data:
            shader_data.draw(GL_LINES)


SHADER_PROGRAM_CREATORS = {}
SHADER_PROGRAM_CREATOR = ShaderProgramsCreator()


def setCreator(context):
    global SHADER_PROGRAM_CREATORS
    creator = SHADER_PROGRAM_CREATORS.get(context, None)
    if creator is None:
        creator = ShaderProgramsCreator()
    global SHADER_PROGRAM_CREATOR
    SHADER_PROGRAM_CREATOR = creator
