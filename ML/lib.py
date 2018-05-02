import ipywidgets
import nglview as nv
import ase

def view(sdict, rep=(1,1,1)):
    ase_str=ase.Atoms(symbols=sdict['species'], positions=sdict['positions'], cell=sdict['cell'])
    return ase_view(ase_str*rep)
def ase_view(s, gui=False):
    view=nv.show_ase(s, gui=gui)
    view.clear_representations()
    view.add_ball_and_stick(aspectRatio=4)
    return view
def traj_view(t):
    view=nv.show_asetraj(t)
    view.clear_representations()
    view.add_ball_and_stick(aspectRatio=4)
    return view