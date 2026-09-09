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

for l in selected:
    if isinstance(l, PointsList):
        for p in l.children:
            if p.atom_type and p.atom_type:
                for obs in p.observers:
                    detach(p, obs.NAME)
        attach(l, 'Sphere')
    else:
        p = l
        if p.atom_type and p.atom_type in isot_types:
            for obs in p.observers:
                detach(p, obs.NAME)
