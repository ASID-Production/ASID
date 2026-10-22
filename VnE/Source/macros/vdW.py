'''
env:
observers: list[str] = ['Sphere', 'Bond', 'Label', 'Plane', 'Line', 'Ellipsoid']
Point
PointsList
attach(point: Point, observer: str)
detach(point: Point, observer: str)
addProp(key: str, value, point: Point)
selected: list[Point | PointList]
root: PointList'''

#VdW rad:
#Charry J, Tkatchenko A. van der Waals Radii of Free and Bonded Atoms from Hydrogen (Z = 1) to Oganesson (Z = 118). J Chem Theory Comput. 2024 Sep 10;20(17):7469-7478. doi: 10.1021/acs.jctc.4c00784. Epub 2024 Aug 29. PMID: 39208255; PMCID: PMC11391583.

import json
import os
f = open(os.path.normpath(os.path.join(os.path.dirname(__file__), 'macros/vdW.json')), 'r')
vdW_d = json.load(f)
f.close()
if selected:
    for l in selected:
        if isinstance(l, PointsList):
            for p in l.children:
                if p.atom_type:
                    p.rad = vdW_d.get(str(p.atom_type), 0)
        else:
            p = l
            if p.atom_type:
                p.rad = vdW_d.get(str(p.atom_type), 0)
else:
    selection = self.selection_model.selection()
    for i in selection:
        inds = i.indexes()
        for ind in inds:
            l = ind.internalPointer()
            if isinstance(l, PointsList):
                for p in l.children:
                    if p.atom_type:
                        p.rad = vdW_d.get(str(p.atom_type), 0)
            else:
                p = l
                if p.atom_type:
                    p.rad = vdW_d.get(str(p.atom_type), 0)