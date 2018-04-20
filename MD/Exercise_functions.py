import numpy as np



def get_LJ_forces(positions,r_m,epsilon):
    Npart,_ = positions.shape
    forces = np.zeros(positions.shape)
    energy = 0.0
    r_m6 = r_m**6
    r_m12 = r_m**12
    factor = -12*epsilon
    
    dirVec = np.zeros((int(Npart*(Npart+1)/2-Npart),3))
    stride = np.cumsum([0]+[it for it in range(Npart-1,0,-1)])
    iparts = range(len(stride[:-1]))
    for ipart, st, nd in zip(iparts, stride[:-1], stride[1:]):
        dirVec[st:nd,:] = positions[ipart+1:] - positions[ipart]
    
    d2 = np.power(dirVec,2).sum(axis=1)
    c2 = np.power(d2,-1)
    c6 = r_m6*np.power(c2,3)
    c12 = np.power(c6,2)
    energy = epsilon * (c12 - 2*c6).sum()
    fs = (factor*(c12-c6)*c2 ).reshape((-1,1))* dirVec
    
    for ipart, st, nd in zip(iparts, stride[:-1], stride[1:]):
        forces[ipart] += fs[st:nd].sum(axis=0)
        forces[ipart+1:] -= fs[st:nd]
    energy = energy
    return forces,energy
def get_Harm_forces(positions,r_m,k_spring):
    Npart,_ = positions.shape
    forces = np.zeros(positions.shape)
    energy = 0.0
    r = np.subtract(positions[1:], positions[:-1])
    d = np.linalg.norm(r,axis=1)
    dx = np.subtract(d,r_m)
    
    f = - 2 * k_spring * np.divide(dx,d).reshape((-1,1))
    rf = np.multiply(r, f)
    
    forces[:Npart-1] -= rf
    forces[1:] += rf
    energy = k_spring * np.power(dx,2).sum()
    return forces,energy
def get_forces(positions,r_m,epsilon,k_spring):
    LJ_forces,LJ_pot = get_LJ_forces(positions,r_m,epsilon)
    H_forces,H_pot = get_Harm_forces(positions,r_m,k_spring)
    return LJ_forces+H_forces,LJ_pot,H_pot

def andersen_thermostat(velocities,temperature,freq,dt):
    if temperature > 0:
        vshape = velocities.shape
        mask = np.random.rand(velocities.size) < 1 - np.exp(-freq*dt)
        Nupdate = np.sum(mask)
        velocities = velocities.flatten()
        velocities[mask] = np.sqrt(temperature)*np.random.normal(loc=0.0, scale=1,size=(Nupdate,))
        velocities = velocities.reshape(vshape)
    return velocities


def MD_NVT_simulator_rapide(positions,velocities,mass,temperature,r_m,epsilon,k_spring,Nstep,dt,enregistrement_stride=10):
    from tqdm import tqdm_notebook as tqdm_cs
    Nparticule, _ = positions.shape
    accelerations = np.zeros(positions.shape)
    pos = []
    vel = []
    masses = np.ones((Nparticule))*mass
    thermostat_frequency = np.sqrt(2)*np.sqrt(2*k_spring)
    Nrecord = Nstep//enregistrement_stride 
    diagnostic = dict(E_variation=np.zeros((Nrecord,)),Temperature=np.zeros((Nrecord,)),
                      E_system=np.zeros((Nrecord,)),E_potentiel=np.zeros((Nrecord,)),
                     E_cinetique=np.zeros((Nrecord,)),time=np.zeros((Nrecord,)))
    dt_half = 0.5*dt
    th_en = 0.5*mass*np.power(velocities,2).sum()
    sys_en = 0
    econs = 0.0
    # Calculates the initial potential energy
    forces,ljpot,sppot = get_forces(positions,r_m,epsilon,k_spring)
    sys_en += ljpot + sppot + th_en
    econs += sys_en 
    ii = 0
    for it in tqdm_cs(range(Nstep)):
        
        #Apply thermostat
        econs += 0.5*mass*np.power(velocities,2).sum()
        velocities = andersen_thermostat(velocities,temperature,thermostat_frequency,dt)
        CoM = np.average(velocities,weights=masses,axis=0).reshape((1,3))
        velocities = velocities - CoM
        econs -= 0.5*mass*np.power(velocities,2).sum()
        
        # half update of velocities
        velocities = velocities + dt_half * accelerations
        # update of positions
        positions = positions + dt * velocities
        # update forces from new posittions
        forces,ljpot,sppot = get_forces(positions,r_m,epsilon,k_spring)
        # update acceleration
        accelerations = forces / mass
        # half update of velocities
        velocities = velocities + dt_half * accelerations
        
        Ekin_tot = 0.5*mass*np.power(velocities,2).sum()
        econs +=  ljpot + sppot + Ekin_tot  - sys_en
        sys_en =  ljpot + sppot + Ekin_tot
        
        if it % enregistrement_stride == 0:
            diagnostic['E_variation'][ii] = econs
            diagnostic['E_system'][ii] = sys_en
            diagnostic['E_cinetique'][ii] = Ekin_tot
            diagnostic['time'][ii] = it * dt
            diagnostic['Temperature'][ii] = Ekin_tot/1.5/len(positions)
            diagnostic['E_potentiel'][ii] = ljpot+sppot
            ii += 1
            
            vel.append(velocities)
            pos.append(positions)
            
        
    return pos,vel,diagnostic