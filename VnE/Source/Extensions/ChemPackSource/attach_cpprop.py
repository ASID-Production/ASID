DIALOG = None

def execute():
    from .ui import select_mol_dialog
    from ..ChemPack import MOLECULE_SYSTEMS
    global DIALOG
    DIALOG = select_mol_dialog.SelectMolDialog(MOLECULE_SYSTEMS, attachCPprop)
    DIALOG.show()


def attachCPprop(molsys):
    from PySide6.QtWidgets import QFileDialog
    from .parsers import PARSER
    molsys = molsys[1]
    file_path = QFileDialog.getOpenFileName()
    if file_path:
        cp_props = PARSER.parsCpProp(file_path[0])[0]
        cp_props = cp_props.children[0].children
        vis_cps = molsys.children[0].children.copy()
        vis_cps.sort(key=lambda x: x.coord[0])
        cp_props.sort(key=lambda x: x.coord[0])
        if len(vis_cps) != len(cp_props):
            return
        else:
            for i in range(len(cp_props)):
                vis_cp = vis_cps[i]
                cp_prop = cp_props[i]
                for key in cp_prop.__dict__:
                    if key[0] != '_' and key not in vis_cp.__dict__:
                        vis_cp.__dict__[key] = cp_prop.__dict__[key]
                    if vis_cp.point() is not None and key[0] != '_' and key not in vis_cp.point().__dict__:
                        vis_cp.point().addProperty(key, cp_prop.__dict__[key], False)
        a = 0
    else:
        return