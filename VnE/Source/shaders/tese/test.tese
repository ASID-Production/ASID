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

layout(quads, equal_spacing, ccw) in;

layout(std140, binding = 0) uniform Matrices
{
    mat4 rotation;
    mat4 scale;
    mat4 scaledef;
    mat4 translation;
    mat4 persp;
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

vec3 tube_int(in float cord_x, in vec3 dest)
{
    vec3 point = vec3(cos(cord_x * 2 * PI), sin(cord_x * 2 * PI), 0.0);

    float cos_b = dot(vec3(0.0, 0.0, 1.0), normalize(dest));
    float sin_b = sqrt(1-pow(cos_b, 2.0));

    vec3 a = normalize(cross(vec3(0.0,0.0,1.0), dest));
    float cos_bm = 1-cos_b;

    mat3 rot;
    rot[0] = vec3(cos_b + pow(a.x, 2.0)*cos_bm, a.x*a.y*cos_bm - a.z*sin_b, a.x*a.z*cos_bm + a.y*sin_b);
    rot[1] = vec3(a.y*a.x*cos_bm + a.z*sin_b, cos_b + pow(a.y, 2.0)*cos_bm, a.y*a.z*cos_bm - a.x*sin_b);
    rot[2] = vec3(a.z*a.x*cos_bm - a.y*sin_b, a.z*a.y*cos_bm + a.x*sin_b, cos_b + pow(a.z, 2.0)*cos_bm);

    point = rot * point;
    return point;
}

void main()
{
    if (bool(pick_tes)){
        color_frag = vec4(0.1,0.1,0.1,1.0);
    }
    else{
        color_frag = color_tes;
    }
    normals_frag = mat3(scaledef * rotation) * sphere_int(gl_TessCoord[0], gl_TessCoord[1]);

    frag_pos.xyz = rad_tes * normals_frag + gl_in[0].gl_Position.xyz;
    frag_pos.w = gl_in[0].gl_Position.w;
    gl_Position =  translation * scale * frag_pos;
}