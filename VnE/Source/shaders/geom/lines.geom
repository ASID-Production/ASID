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

#version 460 core

in gl_PerVertex
{
  vec4 gl_Position;
  float gl_PointSize;
  float gl_ClipDistance[];
} gl_in[];

out gl_PerVertex
{
  vec4 gl_Position;
  float gl_PointSize;
  float gl_ClipDistance[];
};

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

layout(lines) in;
in vec4 color_geom[];
in float rad_geom[];
in float freq_geom[];
in float hfreq_geom[];
layout(triangle_strip, max_vertices = 4) out;
flat out vec4 origin;
out vec4 color_frag;
out vec4 pos;
out vec4 hpos;
out float freq_frag;
out float hfreq_frag;
out float tlen;
out float htlen;

void main() {
    vec3 n = normalize(cross(vec3(gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w - gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w, 0.0), vec3(0.0,0.0,1.0)));
    n.x *= aspect_ratio[0][0];
    n.y *= aspect_ratio[1][1];
    origin = gl_in[0].gl_Position - vec4(n.xy, 0.0, 0.0)*rad_geom[0]*scale[0][0];
    tlen = length(gl_in[1].gl_Position/gl_in[1].gl_Position.w-gl_in[0].gl_Position/gl_in[0].gl_Position.w);
    htlen = 2 * length(n) * (rad_geom[0]/gl_in[0].gl_Position.w)*scale[0][0];

    pos = gl_in[0].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    hpos = origin;
    hpos.xy = origin.xy + 2*(n.xy * rad_geom[0]*scale[0][0]);
    color_frag = color_geom[0];
    gl_Position = gl_in[0].gl_Position;
    gl_Position.xy = gl_Position.xy + (n.xy * rad_geom[0]*scale[0][0]);
    freq_frag = freq_geom[0];
    hfreq_frag = hfreq_geom[0];

    EmitVertex();

    pos = gl_in[1].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    hpos = origin;
    hpos.xy = origin.xy + 2*(n.xy * rad_geom[0]*scale[0][0]);
    color_frag = color_geom[1];
    gl_Position = gl_in[1].gl_Position;
    gl_Position.xy = gl_Position.xy + (n.xy * rad_geom[1]*scale[0][0]);
    freq_frag = freq_geom[1];
    hfreq_frag = hfreq_geom[1];

    EmitVertex();

    pos = gl_in[0].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    hpos = origin;
    color_frag = color_geom[0];
    gl_Position = gl_in[0].gl_Position;
    gl_Position.xy = gl_Position.xy - (n.xy * rad_geom[0]*scale[0][0]);
    freq_frag = freq_geom[0];
    hfreq_frag = hfreq_geom[0];

    EmitVertex();

    pos = gl_in[1].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    hpos = origin;
    color_frag = color_geom[1];
    gl_Position = gl_in[1].gl_Position;
    gl_Position.xy = gl_Position.xy - (n.xy * rad_geom[1]*scale[0][0]);
    freq_frag = freq_geom[1];
    hfreq_frag = hfreq_geom[1];

    EmitVertex();
    EndPrimitive();
}