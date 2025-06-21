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

out gl_PerVertex
{
  vec4 gl_Position;
  float gl_PointSize;
  float gl_ClipDistance[];
};

layout(location = 0) in vec3 pos_vert;
layout(location = 1) in vec4 color;
layout(location = 2) in float rad;
layout(location = 3) in float freq;
layout(location = 4) in float hfreq;
layout(location = 5) in float pick;
layout(std140, binding = 0) uniform Matrices
    {
        mat4 scale;
        mat4 translation;
        mat4 rotation;
        mat4 aspect_ratio;
        mat4 clip_distance;
        mat4 perspective;
        mat4 scene_shift;
    };

out vec4 color_geom;
out float rad_geom;
out float freq_geom;
out float hfreq_geom;

void main()
    {
        if (pick == 1.0) {
            color_geom = vec4(0.0,0.0,0.0,1.0);
        }
        else {
            color_geom = color;
        }
        rad_geom = rad;
        freq_geom = freq;
        hfreq_geom = hfreq;
        gl_Position = translation * perspective * aspect_ratio * scale * rotation * scene_shift * vec4(pos_vert, 1.0);
    }