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

out vec4 out_color;

in vec4 color_frag;
in vec3 normals_frag;
in vec4 frag_pos;

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

uniform vec3 light_pos = vec3(-10, 10, -5);
vec3 viewPos = vec3(0.0,0.0,0.0);

uniform vec3 lightColor = vec3(1.0, 1.0, 1.0);
//vec3 lightColor = vec3(1.0,0.0,0.0);
uniform float specularStrength = 0.5;
uniform float ambientStrength = 0.1;

vec3 pal(in float t) {
    vec3 a = vec3(0.5,0.5,0.5);
    vec3 b = vec3(0.5,0.5,0.5);
    vec3 c = vec3(1.0,0.7,0.4);
    vec3 d = vec3(-0.392,-0.242,-0.192);
    return a + b*cos(6.28318*(c*t + d));
}

void main()
{
    vec3 ambient = ambientStrength * lightColor;

    vec3 normals_frag = normalize(normals_frag);
    vec3 lightDir = normalize(light_pos - frag_pos.xyz);
    float diff = max(dot(normals_frag, lightDir), 0.5);
    vec3 diffuse = diff * lightColor;

    vec3 viewDir = normalize(viewPos - frag_pos.xyz);
    vec3 reflectDir = reflect(-lightDir, normals_frag);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 4.0);
    vec3 specular = specularStrength * spec * lightColor;

    //vec3 color_1 = vec3(0.0,0.0,1.0*((frag_pos.z/frag_pos.w)));
    vec3 out_color3 = color_frag.xyz * (ambient + diffuse + specular);
    //out_color3 = color_frag.xyz * pal(gl_FragCoord.x/wh.x) * (ambient + diffuse + specular);
    out_color3 = color_frag.xyz * (ambient + diffuse + specular);
    out_color = vec4(out_color3, color_frag.w);
    //out_color = vec4(abs(normals_frag), 1.0);
    //out_color = vec4(1,0,0,1);
}