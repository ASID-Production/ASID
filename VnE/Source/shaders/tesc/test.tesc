// Copyright 2023 Alexander A. Korlyukov, Alexander D. Volodin, Petr A. Buikin, Alexander R. Romanenko
// This file is part of ASID - Atomistic Simulation Instruments and Database
// For more information see <https://github.com/ASID-Production/ASID>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// ******************************************************************************************
//  Author:      Alexander A. Korlyukov (head)
//  ORCID:       0000-0002-5600-9886
//  Author:      Alexander D. Volodin (author of cpplib)
//  ORCID:       0000-0002-3522-9193
//  Author:      Petr A. Buikin (author of api_database)
//  ORCID:       0000-0001-9243-9915
//  Author:      Alexander R. Romanenko (author of VnE)
//  ORCID:       0009-0003-5298-6836
//
// ******************************************************************************************


#version 460

layout(vertices = 1) out;

in float rad_tcs[];
in vec4 color_tcs[];
in float pick_tcs[];

layout(std140, binding = 0) uniform Matrices
    {
        mat4 rotation;
        mat4 scale;
        mat4 scaledef;
        mat4 translation;
        mat4 persp;
    };

uniform float shift;

patch out float rad_tes;
patch out vec4 color_tes;
patch out float pick_tes;

void main()
    {
        //Multiples of eight can cause crash on some systems
        gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
        rad_tes = rad_tcs[gl_InvocationID];
        color_tes = color_tcs[gl_InvocationID];
        pick_tes = pick_tcs[gl_InvocationID];
        float grade;
        mat4 transl = mat4(1.0);
        transl[3][2] = shift;
        grade = ((scale * persp * transl * scaledef * rotation * vec4(1.0, 0.0, 0.0, 1.0))/gl_in[gl_InvocationID].gl_Position.w).x;
        int tess = int(round(48 * pow(4 * grade, 2.0)) + 16);

        /*if (grade < 0.5) {
            tess = int(round(56 * pow(2 * grade, 3.0)) + 8);
        }*/

        if (grade > 1.0) {
            tess = 64;
        }

        gl_TessLevelOuter[0] = 8;
        gl_TessLevelOuter[1] = 8;
        gl_TessLevelOuter[2] = 8;
        gl_TessLevelOuter[3] = 8;
        gl_TessLevelInner[0] = 8;
        gl_TessLevelInner[1] = 8;

        /*if (grade < 0.01) {
            gl_TessLevelOuter[0] = 12;
            gl_TessLevelOuter[1] = 12;
            gl_TessLevelOuter[2] = 12;
            gl_TessLevelOuter[3] = 12;
            gl_TessLevelInner[0] = 12;
            gl_TessLevelInner[1] = 12;
        }

        if (grade < 0.05) {
            gl_TessLevelOuter[0] = 16;
            gl_TessLevelOuter[1] = 16;
            gl_TessLevelOuter[2] = 16;
            gl_TessLevelOuter[3] = 16;
            gl_TessLevelInner[0] = 16;
            gl_TessLevelInner[1] = 16;
        }

        if (grade < 0.1) {
            gl_TessLevelOuter[0] = 24;
            gl_TessLevelOuter[1] = 24;
            gl_TessLevelOuter[2] = 24;
            gl_TessLevelOuter[3] = 24;
            gl_TessLevelInner[0] = 24;
            gl_TessLevelInner[1] = 24;
        }

        if (grade < 0.25) {
            gl_TessLevelOuter[0] = 32;
            gl_TessLevelOuter[1] = 32;
            gl_TessLevelOuter[2] = 32;
            gl_TessLevelOuter[3] = 32;
            gl_TessLevelInner[0] = 32;
            gl_TessLevelInner[1] = 32;
        }

        if (grade < 0.5) {
            gl_TessLevelOuter[0] = 64;
            gl_TessLevelOuter[1] = 64;
            gl_TessLevelOuter[2] = 64;
            gl_TessLevelOuter[3] = 64;
            gl_TessLevelInner[0] = 64;
            gl_TessLevelInner[1] = 64;
        }*/
    }