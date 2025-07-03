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
out PerVertex {
 vec4 origin;
 vec4 color_frag;
 vec4 pos;
 vec4 hpos;
 float freq_frag;
 float hfreq_frag;
 float tlen;
 float htlen;
} out_v;

void main() {
    vec3 n = normalize(cross(vec3(gl_in[1].gl_Position.xy/gl_in[1].gl_Position.w - gl_in[0].gl_Position.xy/gl_in[0].gl_Position.w, 0.0), vec3(0.0,0.0,1.0)));
    n.x *= aspect_ratio[0][0];
    n.y *= aspect_ratio[1][1];
    out_v.origin = gl_in[0].gl_Position - vec4(n.xy, 0.0, 0.0)*rad_geom[0]*scale[0][0];
    out_v.tlen = length(gl_in[1].gl_Position/gl_in[1].gl_Position.w-gl_in[0].gl_Position/gl_in[0].gl_Position.w);
    out_v.htlen = 2 * length(n) * (rad_geom[0]/gl_in[0].gl_Position.w)*scale[0][0];

    out_v.pos = gl_in[0].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    out_v.hpos = out_v.origin;
    out_v.hpos.xy = out_v.origin.xy + 2*(n.xy * rad_geom[0]*scale[0][0]);
    out_v.color_frag = color_geom[0];
    gl_Position = gl_in[0].gl_Position;
    gl_Position.xy = gl_Position.xy + (n.xy * rad_geom[0]*scale[0][0]);
    out_v.freq_frag = freq_geom[0];
    out_v.hfreq_frag = hfreq_geom[0];

    EmitVertex();

    out_v.pos = gl_in[1].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    out_v.hpos = out_v.origin;
    out_v.hpos.xy = out_v.origin.xy + 2*(n.xy * rad_geom[0]*scale[0][0]);
    out_v.color_frag = color_geom[1];
    gl_Position = gl_in[1].gl_Position;
    gl_Position.xy = gl_Position.xy + (n.xy * rad_geom[1]*scale[0][0]);
    out_v.freq_frag = freq_geom[1];
    out_v.hfreq_frag = hfreq_geom[1];

    EmitVertex();

    out_v.pos = gl_in[0].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    out_v.hpos = out_v.origin;
    out_v.color_frag = color_geom[0];
    gl_Position = gl_in[0].gl_Position;
    gl_Position.xy = gl_Position.xy - (n.xy * rad_geom[0]*scale[0][0]);
    out_v.freq_frag = freq_geom[0];
    out_v.hfreq_frag = hfreq_geom[0];

    EmitVertex();

    out_v.pos = gl_in[1].gl_Position - vec4(n * rad_geom[0]*scale[0][0], 0.0);
    out_v.hpos = out_v.origin;
    out_v.color_frag = color_geom[1];
    gl_Position = gl_in[1].gl_Position;
    gl_Position.xy = gl_Position.xy - (n.xy * rad_geom[1]*scale[0][0]);
    out_v.freq_frag = freq_geom[1];
    out_v.hfreq_frag = hfreq_geom[1];

    EmitVertex();
    EndPrimitive();
}