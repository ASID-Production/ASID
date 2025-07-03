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

in PerVertex {
 in vec4 origin;
 vec4 color_frag;
 vec4 pos;
 vec4 hpos;
 float freq_frag;
 float hfreq_frag;
 float tlen;
 float htlen;
} in_f;
out vec4 color_out;

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

void main()
    {
        vec3 x = in_f.pos.xyz/in_f.pos.w - in_f.origin.xyz/in_f.origin.w;
        //x.x /= aspect_ratio[0][0];
        //x.y /= aspect_ratio[1][1];
        float len = length(x);

        vec3 hx = in_f.hpos.xyz/in_f.hpos.w - in_f.origin.xyz/in_f.origin.w;
        //hx.x /= aspect_ratio[0][0];
        //hx.y /= aspect_ratio[1][1];
        float hlen = length(hx);

        color_out = in_f.color_frag;

        if (int(in_f.freq_frag+0.0001) == 0) {
        }
        else {
            int n = int(floor(in_f.freq_frag+0.0001)*2-1);
            float p = in_f.tlen/n;
            if (mod(int(len/p+0.0001), 2) == 1){
                discard;
            }
        }

        if (int(in_f.hfreq_frag+0.0001) == 0) {
        }
        else {
            int hn = int(floor(in_f.hfreq_frag+0.0001)*2-1);
            float hp = in_f.htlen/hn;
            if (mod(int(hlen/hp+0.0001), 2) == 1){
                discard;
            }
        }
    }