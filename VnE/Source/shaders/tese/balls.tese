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
#define PI 3.1415926538

in gl_PerVertex
{
  vec4 gl_Position;
  float gl_PointSize;
  float gl_ClipDistance[];
} gl_in[gl_MaxPatchVertices];

out gl_PerVertex
{
    vec4 gl_Position;
    float gl_PointSize;
    float gl_ClipDistance[];
};

layout(quads, equal_spacing, ccw) in;

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

patch in vec4 color_tes;
patch in float rad_tes;
patch in float pick_tes;

out vec4 color_frag;
out vec3 normals_frag;
out vec4 frag_pos;

vec3 sphere_int(in float cord_x, in float cord_y)
{
    float cos_a = cos(cord_x * 2 * PI);
    float sin_a = sin(cord_x * 2 * PI);
    float cos_b = cos((cord_y - 0.5) * PI);
    float sin_b = sin((cord_y - 0.5) * PI);
    return vec3(cos_a, cos_b * sin_a, sin_b * sin_a);
}


void main()
{
    if (bool(pick_tes)){
        color_frag = vec4(0.1,0.1,0.1,color_tes.w);
    }
    else{
        color_frag = color_tes;
    }
    normals_frag = mat3(rotation) * sphere_int(gl_TessCoord[0], gl_TessCoord[1]);

    frag_pos = translation * perspective * aspect_ratio * scale * (vec4(normals_frag * rad_tes, 0.0) + rotation * gl_in[0].gl_Position);
    gl_Position = frag_pos;
}