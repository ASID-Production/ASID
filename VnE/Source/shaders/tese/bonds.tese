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
    mat4 scale;
    mat4 translation;
    mat4 rotation;
    mat4 aspect_ratio;
    mat4 clip_distance;
    mat4 perspective;
    mat4 scene_shift;
};

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

in vec4 color_tes[];
patch in float rad_tes;

out vec4 color_frag;
out vec3 normals_frag;
out vec4 frag_pos;

const vec3 s = normalize(vec3(1,1,1));

vec3 tube_int(in float cord_x, in float cord_y, in vec3 dest)
{
    //vec3 point = (sin(cord_x * 2 * PI) * normalize(vec3(1.0, 0.0, -dest.x/dest.z))) + (cos(cord_x * 2 * PI) * normalize(cross(dest,vec3(1.0, 0.0, -dest.x/dest.z))));

    vec3 ndest = normalize(dest);
    vec3 point = sin(cord_x * 2 * PI) * normalize(s - dot(ndest, s) * ndest) + cos(cord_x * 2 * PI) * normalize(cross((s - dot(ndest, s) * ndest), dest));

    return point;
}

void main()
{
    color_frag = color_tes[0] * (1 - gl_TessCoord[1]) + color_tes[1] * gl_TessCoord[1];
    vec3 dest = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
    normals_frag = tube_int(gl_TessCoord[0], gl_TessCoord[1], dest);
    //frag_pos.xyz = mat3(aspect_ratio * scale) * rad_tes * normals_frag + gl_in[0].gl_Position.xyz + dest * gl_TessCoord[1];
    frag_pos = translation * perspective * aspect_ratio * scale * rotation * vec4(((normals_frag * rad_tes) + gl_in[0].gl_Position.xyz + (dest * gl_TessCoord[1])), 1.0);
    normals_frag = mat3(rotation) * normals_frag;
    //frag_pos.w = gl_in[0].gl_Position.w * (1 - gl_TessCoord[1]) + gl_in[1].gl_Position.w * gl_TessCoord[1];
    gl_Position = frag_pos;
}