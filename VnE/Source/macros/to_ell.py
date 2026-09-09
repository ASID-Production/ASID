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

isot_types = [1,6]
for l in selected:
    detach(l, 'Sphere')
    if isinstance(l, PointsList):
        for p in l.children:
            if p.atom_type and p.atom_type in isot_types:
                attach(p, 'Sphere')
            elif p.atom_type:
                attach(p, 'Ellipsoid')
    else:
        p = l
        if p.atom_type and p.atom_type in isot_types:
            attach(p, 'Sphere')
        elif p.atom_type:
            attach(p, 'Ellipsoid')

