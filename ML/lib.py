import ipywidgets
import nglview as nv

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