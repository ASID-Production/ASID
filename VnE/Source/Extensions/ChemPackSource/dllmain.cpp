#if defined(_WIN32)
#define EXPORT __declspec(dllexport)
#elif defined(__linux__)
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

using namespace std;

decltype(auto) genMap()
{
    
    cout << "GenMap\n";
    ifstream stream;
    stream.open("./Source/Extensions/ChemPackSource/BondLength.ini");
    if (!stream.is_open()) {
        cout << "File not open\n";
    }
    unsigned int first;
    unsigned int second;
    float min;
    float max;
    static array<array<array<float, 2>, 120>, 120> ret;

    while (stream >> first) {
        if ((((first == 0) || !(stream >> second)) || !(stream >> min)) || !(stream >> max)) {
            break;
        }
        if (first > second) {
            swap(first, second);
        }
        ret[first][second][0] = min;
        ret[first][second][1] = max;
        }
    return ret;
}

array<array<array<float, 2>, 120>, 120> bond_map = genMap();

inline float dist(float coord1[4], float coord2[4]) {
    float distance = sqrt(pow(coord2[1] - coord1[1], 2) + pow(coord2[2] - coord1[2], 2) + pow(coord2[3] - coord1[3], 2));
    return distance;
}

extern "C" EXPORT const char* genBonds(const int size, float** atoms) {
    static string line = "";
    line.clear();
    for (int a1 = 0; a1 < size - 1; a1++) {
        unsigned int a1_type = (unsigned int)atoms[a1][0];
        for (int a2 = a1 + 1; a2 < size; a2++) {
            unsigned int a2_type = (unsigned int)atoms[a2][0];
            float distance = dist(atoms[a1], atoms[a2]);
            float min;
            float max;
            if (a1_type <= a2_type) {
                min = bond_map[a1_type][a2_type][0];
                max = bond_map[a1_type][a2_type][1];
            }
            else {
                min = bond_map[a2_type][a1_type][0];
                max = bond_map[a2_type][a1_type][1];
            }
            if ((min < distance) && (distance < max)) {
                line += to_string(a1) + ":" + to_string(a2) + "\n";
            }
        }
    }
    return line.c_str();
}