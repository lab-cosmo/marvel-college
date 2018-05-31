import ipywidgets
import nglview as nv
import ase
import numpy as np

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

def aff_bruit(a,x):
    return a*x+np.array([np.random.normal(scale=0.2) for i in range(len(x))])
def aff_bruit2D(a,x):
    return a[0]*x[0,:]+a[1]*x[1,:]+np.array([np.random.normal(scale=0.2*np.sqrt(2)) for i in range(len(x[0,:]))])
def aff_bruitND(A,X,Ndim):
    return np.dot(X.T,A)+np.array([np.random.normal(scale=np.sqrt(Ndim)*0.2) for i in range(len(X.T))])
def aff_bruit_biaisND(A,X,Ndim):
    Np=len(X.T)
    return np.dot(X.T,A)+np.array([np.random.normal(loc=np.rint(2*(float(i)-Np/2)/Np), scale=np.sqrt(Ndim)*0.2) for i in range(Np)])