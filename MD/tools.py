from nglview import show_asetraj,show_ase
import numpy as np
import ase
from ipywidgets import interactive,FloatSlider,interact
import matplotlib.pyplot as plt

def visualiser_trajectoire(pos,stride):
    from ase.visualize import view
    from nglview import show_asetraj,show_ase
    from ase import Atoms
    Natom,_ = pos[0].shape
    num = [1]*Natom
    xsize,ysize = 400,300
    frames = [Atoms(numbers=num,positions=pp,pbc=False) for pp in pos[::stride]]
    view = show_asetraj(frames, gui=False)
    view.clear_representations()
    view.representations = [
    {"type": "ball+stick", "params": {
        "aspectRatio": 4,'color': "atomindex",
    }},
    {"type": "distance", "params": {
        "atomPair": [[it,it+1] for it in range(Natom-1)],
        'colorScheme': "black",
    'labelSize': 4,
    'labelColor': "white",
    'labelVisible': False,
    'labelUnit': "angstrom",
    'radiusScale': 1,
    'opacity': 1,
    'name': "link",
    'side': "front",
    'useCylinder': True
    }}
    ]
    view._remote_call('setSize', target='Widget',
                               args=['%dpx' % (xsize,), '%dpx' % (ysize,)])
    return view

def faire_une_chaine_circulaire(N,r_m):
    r = N * r_m / (2*np.pi)
    n = np.arange(N)
    pos = np.asarray([r * np.cos(2*np.pi* n/N), r * np.sin(2*np.pi* n/N),np.zeros(N)]).T
    return pos.reshape((-1,3))

def faire_une_chaine_lineaire(N,r_m):
    n = np.arange(N)
    pos = np.asarray([n*r_m,n*0,n*0]).T
    return pos.reshape((-1,3))


def manipulation_hist():
    the_figure, the_plot = plt.subplots(figsize=(5,5))
    the_plot.set_xlabel("vitesse [l/t]")
    the_plot.set_ylabel("Nombre de particules")

    def f(temperature):
        the_plot.axes.clear()
        velocities = np.zeros((500,1))
        velocities = np.random.normal(loc=0.0, scale=np.sqrt(temperature),size=velocities.shape)
        the_plot.hist(velocities);
        the_plot.set_xlim(-10,10)
        the_figure.canvas.draw()
        the_figure.canvas.flush_events()
    return interactive(f, temperature=FloatSlider(min=0.1, max=7, step=0.3,description=r'\(T\)',value=1))

def manipulation_LJ_force(potentiel,force):
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    r_m,epsilon = 1,1
    def f(d):
        fig, ax = plt.subplots(figsize=(5,5))
        circle1 = plt.Circle((0, 2), 0.2, color='b')
        circle2 = plt.Circle((d, 2), 0.2, color='r')
        
        ax.add_artist(circle1)
        ax.add_artist(circle2)
        st = np.array([d,2])
        
        if d != r_m:
            end = np.array([force(d,r_m,epsilon),2])
            nnn = end - st 
            ax.arrow(st[0], st[1],2*nnn[0]/force(r_m-0.3,r_m,epsilon), nnn[1], head_width=0.25, head_length=0.15, fc='k', ec='k',width=0.05)
            plt.text(2, 0.65, r'$F_{LJ}=$'+'{val:.1f}'.format(val=end[0]), fontdict=font)
        else:
            plt.text(2, 0.65, r'$F_{LJ}=$'+'{val:.1f}'.format(val=0), fontdict=font)
        r = np.linspace(0.01, 4, num=1000)
        plt.plot(r, potentiel(r,r_m,epsilon))
        
        plt.xlim(-0.5,4)
        plt.ylim(-1.5,3)
        plt.show()
    return interactive(f, d=FloatSlider(min=r_m-0.3, max=4, step=0.1,description=r'\(d\)',value=r_m))

def manipulation_Harmonique_force(potentiel,force):
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    r_m,k_spring = 2,1
    def f(d):
        fig, ax = plt.subplots(figsize=(5,5))
        circle1 = plt.Circle((0, 2), 0.2, color='b')
        circle2 = plt.Circle((d, 2), 0.2, color='r')
        
        ax.add_artist(circle1)
        ax.add_artist(circle2)
        st = np.array([d,2])
        
        if d != r_m:
            end = np.array([force(d,r_m,k_spring),2])
            nnn = end - np.array([0,2])
            ax.arrow(st[0], st[1],nnn[0]/abs(force(0.5,r_m,k_spring)), nnn[1], head_width=0.25, head_length=0.15, fc='k', ec='k',width=0.05)
            plt.text(2, 3, r'$F_{Harm}=$'+'{val:.1f}'.format(val=end[0]), fontdict=font)
        else:
            plt.text(2, 3, r'$F_{Harm}=$'+'{val:.1f}'.format(val=0), fontdict=font)
        r = np.linspace(-0.5, 4, num=50)
        plt.plot(r, potentiel(r,r_m,k_spring))
        
        plt.xlim(-0.5,4)
        plt.ylim(-0.5,4)
        plt.show()
    return interactive(f, d=FloatSlider(min=0.5, max=3.5, step=0.1,description=r'\(d\)',value=r_m))

def manipulation_Harmonique(func):
    def f(r_m,k):
        fig, ax = plt.subplots(figsize=(5,5))
        circle1 = plt.Circle((0, 2), 0.2, color='b')
        circle2 = plt.Circle((r_m, 2), 0.2, color='r')
        ax.add_artist(circle1)
        ax.add_artist(circle2)
        r = np.linspace(-2, 6, num=100)
        plt.plot(r, func(r,r_m,k))
        
        plt.xlim(-0.5,4)
        plt.ylim(-0.5,4)
        plt.show()
    return interactive(f, r_m=FloatSlider(min=0.5, max=5, step=0.2,description=r'\(r_m\)',value=2), 
                       k=FloatSlider(min=0.5,max=5.,step=0.2,description=r'\(k\)',value=2))

def manipulation_LJ(func):
    def f(r_m,epsilon):
        fig, ax = plt.subplots(figsize=(5,5))
        circle1 = plt.Circle((0, 2), 0.2, color='b')
        circle2 = plt.Circle((r_m, 2), 0.2, color='r')
        ax.add_artist(circle1)
        ax.add_artist(circle2)
        r = np.linspace(0.01, 6, num=1000)
        plt.plot(r, func(r,r_m,epsilon))
        
        plt.xlim(-0.5,5)
        plt.ylim(-2.5,3)
        plt.show()
    return interactive(f, r_m=FloatSlider(min=0.1, max=5, step=0.2,description=r'\(r_m\)',value=2), 
                       epsilon=FloatSlider(min=0.1,max=2.,step=0.2,description=r'\(\epsilon\)',value=2))

def get_numerical_force(pot_func,r,*args):
    e_x,e_y,e_z = np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])
    h = 1e-3
    F_x = - ( pot_func(norm(r+e_x*h),*args) - pot_func(norm(r-e_x*h),*args) )  / (2*h)
    F_y = - ( pot_func(norm(r+e_y*h),*args) - pot_func(norm(r-e_y*h),*args) )  / (2*h)
    F_z = - ( pot_func(norm(r+e_z*h),*args) - pot_func(norm(r-e_z*h),*args) )  / (2*h)
    return np.array([F_x,F_y,F_z])
