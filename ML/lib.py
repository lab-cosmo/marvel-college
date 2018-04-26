import ipywidgets
import nglview as nv
import ase

def view(sdict, tu):
    ase_str=ase.Atoms(symbols=sdict['species'], positions=sdict['positions'], cell=sdict['cell'])
    return ase_view(ase_str*tu)
def ase_view(s):
    view=nv.show_ase(s, gui=True)
    view.clear_representations()
    view.add_ball_and_stick(aspectRatio=4)
    return view
def traj_view(t):
    view=nv.show_asetraj(t, gui=True)
    view.clear_representations()
    view.add_ball_and_stick(aspectRatio=4)
    return view