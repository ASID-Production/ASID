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

out gl_PerVertex { vec4 gl_Position;};

layout (location = 0) in vec3 vertex; // <vec3 pos>
layout (location = 1) in vec2 tex_cords; // <vec2 tex_coord>
layout (location = 2) in vec2 size;
layout (location = 3) in vec2 shifts;

out vec2 TexCoords;

layout(std140, binding = 0) uniform Matrices
{
    mat4 scale;
    mat4 translation;
    mat4 rotation;
    mat4 aspect_ratio;
    mat4 clip_distance;
    mat4 perspective;
    mat4 scene_shift;
    vec2 wh;
};

uniform float const_scale;

void main()
    {
        vec4 pos = translation * perspective * aspect_ratio * scale * rotation * scene_shift * vec4(vertex, 1.0);
        pos.x = pos.x + (size.x + shifts.x) * pos.w;
        pos.y = pos.y + shifts.y + size.y * pos.w;
        pos.z = -pos.w;
        gl_Position = pos;
        TexCoords = tex_cords;
    }