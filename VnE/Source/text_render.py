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
import freetype
import numpy as np

openGL_context = None


class Font:
    """Fonts class"""

    # line simbols
    line = 'qwertyuiop[]asdfghjkl;zxcvbnm,./"QWERTYUIOP{}ASDFGHJKL:ZXCVBNM<>?1234567890*-+!@#$%^&()_ '

    def __init__(self):
        self.char_map = {}
        self.create_chars()
        return

    def create_chars(self):
        """
        Creates chars textures and assign texture options
        :return:
        """

        for chr in self.line:
            face = freetype.Face(f'./Source/fonts/arial.ttf')
            face.set_char_size(width=0, height=int((bin(32)+'0'*6), 2), hres=576, vres=0)
            face.load_char(chr)
            bitmap = face.glyph.bitmap
            char = CharTex(
                buffer=bitmap.buffer,
                width=bitmap.width,
                height=bitmap.rows,
                bearingX=face.glyph.bitmap_left,
                bearingY=face.glyph.bitmap_top,
                advance=face.glyph.advance.x,
                texture=glGenTextures(1)
            )
            self.char_map[chr] = char

            glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
            glBindTexture(GL_TEXTURE_2D, char.texture)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
            glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                char.width,
                char.height,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                char.buffer)
        return


class CharTex:
    """Char texture class, just dict"""

    def __init__(self, **kwargs):
        self.buffer = kwargs['buffer']
        self.width = kwargs['width']
        self.height = kwargs['height']
        self.bearingX = kwargs['bearingX']
        self.bearingY = kwargs['bearingY']
        self.advance = kwargs['advance']
        self.texture = kwargs['texture']


class Letter:
    def __init__(self, chr, vert, font: Font):
        self.chr = chr
        self.texture = font.char_map[self.chr].texture
        self.vert = vert


class Word:

    def __init__(self, xyz, wh, font, word=None):
        self._xyz = xyz
        self._word = word
        self._font = font
        self._buffer_data = np.array([], dtype=np.float32)
        self._wh = wh
        if word is not None:
            self.gen_letters(self._word, self._xyz, self._wh, self._font)

    def gen_letters(self, word, xyz, wh, font):
        self._letters = []
        x,y,z = xyz
        x_shift = 0
        w_s,h_s = wh
        for chr in word:
            xpos = x
            ypos = y
            y_shift = -(font.char_map[chr].height - font.char_map[chr].bearingY)

            w = font.char_map[chr].width/w_s
            h = font.char_map[chr].height/h_s
            x_shift_p_bearing = x_shift + font.char_map[chr].bearingX/w_s
            text_vert = np.array([[xpos, ypos, z, 0.0, 1.0, 0.0, 0.0, x_shift_p_bearing, y_shift],
                                  [xpos, ypos, z, 0.0, 0.0, 0.0, h, x_shift_p_bearing, y_shift],
                                  [xpos, ypos, z, 1.0, 1.0, w, 0.0, x_shift_p_bearing, y_shift],
                                  [xpos, ypos, z, 0.0, 0.0, 0.0, h, x_shift_p_bearing, y_shift],
                                  [xpos, ypos, z, 1.0, 0.0, w, h, x_shift_p_bearing, y_shift],
                                  [xpos, ypos, z, 1.0, 1.0, w, 0.0, x_shift_p_bearing, y_shift]], dtype=np.float32)
            x_shift += int(bin(font.char_map[chr].advance)[:-6], 2)/w_s
            self._letters.append(Letter(chr, text_vert, font))
        self.gen_buffer()

    def delete_letter(self, pos_lett):
        if type(pos_lett) == int:
            for lett in self._letters[pos_lett+1:]:
                lett.vert[:,0] -= self._font.char_map[self._letters[pos_lett].chr].advance/self._wh[0]
        elif type(pos_lett) == Letter:
            for lett in self._letters[self._letters.index(pos_lett)+1:]:
                lett.vert[:,0] -= self._font.char_map[pos_lett.chr].advance/self._wh[0]
            self._letters.remove(pos_lett)

    def add_letter(self, chr, pos):
        w_s,h_s = self._wh
        font = self._font
        try:
            x,y,z = self._letters[pos].vert[0,0] - font.char_map[self._letters[pos].chr].bearingX/w_s,\
                    self._xyz[1],\
                    self._xyz[2]
        except IndexError:
            x,y,z = self._letters[pos - 1].vert[0, 0] - font.char_map[self._letters[pos].chr].bearingX/w_s +\
                                                       font.char_map[self._letters[pos-1].chr].advance/w_s, \
                    self._xyz[1], \
                    self._xyz[2]

        xpos = x
        ypos = y
        y_shift = -(font.char_map[chr].height - font.char_map[chr].bearingY)

        w = font.char_map[chr].width/w_s
        h = font.char_map[chr].heighth_s
        x_shift_p_bearing = font.char_map[chr].bearingX / w_s
        text_vert = np.array([[xpos, ypos, z, 0.0, 1.0, 0.0, 0.0, x_shift_p_bearing, y_shift],
                              [xpos, ypos, z, 0.0, 0.0, 0.0, h, x_shift_p_bearing, y_shift],
                              [xpos, ypos, z, 1.0, 1.0, w, 0.0, x_shift_p_bearing, y_shift],
                              [xpos, ypos, z, 0.0, 0.0, 0.0, h, x_shift_p_bearing, y_shift],
                              [xpos, ypos, z, 1.0, 0.0, w, h, x_shift_p_bearing, y_shift],
                              [xpos, ypos, z, 1.0, 1.0, w, 0.0, x_shift_p_bearing, y_shift]], dtype=np.float32)
        for lett in self._letters[pos:]:
            lett.vert[:,0] += font.char_map[chr].advance/w_s
        letter = Letter(chr, text_vert, font)
        self._letters.insert(pos, letter)

    def move_to(self, xyz):
        vec = xyz - self._xyz
        for lett in self._letters:
            lett.vert[:,0] += vec[0]
            lett.vert[:,1] += vec[1]
            lett.vert[:,2] += vec[2]
        if self._word is not None:
            self._buffer_data[:,0] += vec[0]
            self._buffer_data[:,1] += vec[1]
            self._buffer_data[:,2] += vec[2]

    def gen_buffer(self):
        self._buffer_data = None
        for lett in self._letters:
            if self._buffer_data is None:
                self._buffer_data = lett.vert
            else:
                self._buffer_data = np.append(self._buffer_data, lett.vert, axis=0)
        return self._buffer_data

    def get_textures(self):
        return [x.texture for x in self._letters]

    @property
    def word(self):
        return self._word

    @word.setter
    def word(self, word):
        self._word = word
        self._letters = self.gen_letters(self._word, self._xyz, self._wh, self._font)

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self.move_to(xyz)
        self._xyz = xyz


class TextBuffer:

    def __init__(self, font):
        self.char_vert = {}
        self.words = []
        self.font = font
        self.vert_format = (6,8)

    def gen_char_vert(self):
        self.char_vert = {}
        for word in self.words:
            for letter in word.letters:
                if self.char_vert.get(letter.chr, None) is not None:
                    self.char_vert[letter.chr] = np.append(self.char_vert[letter.chr], letter.vert.reshape((1,6,letter.vert.shape[-1])), axis=0)
                else:
                    self.char_vert[letter.chr] = np.array([letter.vert])

    def add_word(self, word, xyz, wh, font=None, index=None):
        if font is None:
            font = self.font
        word = Word(word, xyz, wh, font, self)
        if index is None:
            self.words.append(word)
        else:
            self.words.insert(index, word)
        return word

    def delete_word(self, word):
        self.words.remove(word)
        for letter in word.letters:
            ind = [i for i, x in enumerate(np.where(self.char_vert[letter.chr] == letter.vert, True, False)) if x.all()]
            np.delete(self.char_vert[letter.chr], ind[0], None)

    def add_letter(self, letter):
        if self.char_vert.get(letter.chr, None):
            self.char_vert[letter.chr].append(letter.vert)
        else:
            self.char_vert[letter] = [letter.vert]

    def add_var_to_vert(self):
        self.vert_format = (self.vert_format[0],self.vert_format[1]+1)
        for word in self.words:
            for letter in word.letters:
                np.append(letter.vert, np.zeros((6,1)), axis=1)

    def change_vert_par(self, pos: int, data: int, *args, word=None, letter=None, stay=None):
        """

        :param pos:
        :param data: byte/int
        :param args:
        :param word:
        :param letter:
        :param stay: stay mask (0b1111)
        :return:
        0.0 - nothing
        1.0 - highlighted by mouse
        2.0 - picked
        """
        if stay is None:
            stay = 0b0000

        if letter is not None:
            letter.vert[:,pos] = [float((int(x) & stay) | data) for x in letter.vert[:,pos]]
        elif word is not None:
            for letter in word.letters:
                letter.vert[:,pos] = [float((int(x) & stay) | data) for x in letter.vert[:,pos]]
        else:
            for word in self.words:
                for letter in word.letters:
                    letter.vert[:,pos] = [float((int(x) & stay) | data) for x in letter.vert[:,pos]]

    def remove_letter(self, letter):
        self.char_vert[letter.chr].remove(letter.vert)