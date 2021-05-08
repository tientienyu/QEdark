import numpy as np
from scipy.special import gamma, iv, modstruve, hyp2f1, erf
from scipy import integrate, interpolate
from scipy.integrate import quad, dblquad, nquad
from QEdark_constants import *


def vmin(EE,qin,mX):
    q = qin * alpha *me_eV
    return (EE/q+q/(2*mX))*c_light*1e-3 # to convert to km/s

def etaSHM(vmin, _params):
    """
    Standard Halo Model with sharp cutoff. 
    Fiducial values are v0=220 km/s, vE=232 km/s, vesc= 544 km/s
    params = [v0, vE, vesc]
    """
    v0 = _params[0]
    vE = _params[1]
    vesc = _params[2]
    KK=v0**3*(-2.0*np.exp(-vesc**2/v0**2)*np.pi*vesc/v0+np.pi**1.5*erf(vesc/v0))
#    print('KK=',KK)
    def func(vx2):
        return np.exp(-vx2/(v0**2))

    if vmin <= vesc - vE:
        # eq. B4 from 1509.01598
        def bounds_cosq():
            return [-1,1]
        def bounds_vX(cosq):
            return [vmin, -cosq*vE+np.sqrt((cosq**2-1)*vE**2+vesc**2)]
        def eta(vx,cosq):
            return (2*np.pi/KK)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_vX,bounds_cosq])[0]
        
    elif vesc - vE < vmin <= vesc + vE:
        # eq. B5 from 1509.01598
        def bounds_cosq(vx):
            return [-1, (vesc**2-vE**2-vx**2)/(2*vx*vE)] 
        def bounds_vX():
            return [vmin, vE+vesc]
        def eta(cosq,vx):
            return (2*np.pi/KK)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_cosq,bounds_vX])[0]
    else:
        return 0

def etaTsa(vmin, _params):
    """
    Tsallis Model, q = .773, v0 = 267.2 km/s, and vesc = 560.8 km/s
    give best fits from arXiv:0902.0009. 
    params = [v0, vE, q]
    """
    v0 = _params[0]
    vE = _params[1]
    q = _params[2]
    if q <1:
        vesc = v0/np.sqrt(1-q)
    else:
        vesc = 560.8e5 # cm/s
#    vesc = 544*kmpers_nu ## to test against SHM    
    def func(vx2):
        if q == 1:
            return np.exp(-vx2/v0**2)
        else:
            return (1-(1-q)*vx2/v0**2)**(1/(1-q))
    " calculate normalization constant "
    def inttest(vx):
        if q == 1:
            if vx <= vesc:
                return vx**2*np.exp(-vx**2/v0**2)
            else:
                return 0 
        else:
            if vx <= vesc:
                return vx**2*(1-(1-q)*vx**2/v0**2)**(1/(1-q))
            else:
                return 0            
    def bounds():
        return [0.,vesc]
    K_=4*np.pi*nquad(inttest,[bounds])[0]

#    K_ = 4/3*np.pi*vesc**3*hyp2f1(3/2, 1/(q-1), 5/2, (1-q)*vesc**2/v0**2) # analytic expression, runs faster
    
    if vmin <= vesc - vE:
        def bounds_cosq():
            return [-1,1]
        def bounds_vX(cosq):
            return [vmin, -cosq*vE+np.sqrt((cosq**2-1)*vE**2+vesc**2)]
        def eta(vx,cosq):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_vX,bounds_cosq])[0]
        
    elif vesc - vE < vmin <= vesc + vE:
        def bounds_cosq(vx):
            return [-1, (vesc**2-vE**2-vx**2)/(2*vx*vE)] 
        def bounds_vX():
            return [vmin, vE+vesc]
        def eta(cosq,vx):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_cosq,bounds_vX])[0]
    else:
        return 0   

def etaDPL(vmin, _params):
    """
    Double Power Law Profile, 1.5 <= k <= 3.5 found to give best fit to N-body
    simulations. 
    takes input velocities in km/s
    params = [vmin, v0, vE, vesc, k]
    """
    v0 = _params[0]
    vE = _params[1]    
    vesc = _params[2]

    
    def func(vx2):
        return (np.exp((vesc**2-vx2)/(k*v0**2))-1)**k
    " calculate normalization constant "
    def inttest(vx):
        if vx <= vesc:
            return vx**2*(np.exp((vesc**2-vx**2)/(k*v0**2))-1)**k
        else:
            return 0
    def bounds():
        return [0.,vesc]
    K_=4*np.pi*nquad(inttest,[bounds])[0]   
    
    if vmin <= vesc - vE:
        def bounds_cosq():
            return [-1,1]
        def bounds_vX(cosq):
            return [vmin, -cosq*vE+np.sqrt((cosq**2-1)*vE**2+vesc**2)]
        def eta(vx,cosq):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_vX,bounds_cosq])[0]
        
    elif vesc - vE < vmin <= vesc + vE:
        def bounds_cosq(vx):
            return [-1, (vesc**2-vE**2-vx**2)/(2*vx*vE)] 
        def bounds_vX():
            return [vmin, vE+vesc]
        def eta(cosq,vx):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_cosq,bounds_vX])[0]
    else:
        return 0

def etaMSW(vmin_in, params):
    """
    empirical model by Mao, Strigari, Weschler 
    takes input velocities in km/s
    params = [vmin, v0, vE, vesc, p]
    """
    v0 = _params[0]
    vE = _params[1]    
    vesc = _params[2]
    p = _params[3]
    
    def func(vx2):
        return np.exp(-np.sqrt(vx2)/v0)*(vesc**2-vx2)**p

    " calculate normalization constant "       
    def inttest(vx):
        if vx <= vesc:
            return vx**2*np.exp(-vx/v0)*(vesc**2-vx**2)**p
        else:
            return 0
    def bounds():
        return [0.,vesc]
    K_=4*np.pi*nquad(inttest,[bounds])[0]   
        
    if vmin <= vesc - vE:
        def bounds_cosq():
            return [-1,1]
        def bounds_vX(cosq):
            return [vmin, -cosq*vE+np.sqrt((cosq**2-1)*vE**2+vesc**2)]
        def eta(vx,cosq):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_vX,bounds_cosq])[0]
        
    elif vesc - vE < vmin <= vesc + vE:
        def bounds_cosq(vx):
            return [-1, (vesc**2-vE**2-vx**2)/(2*vx*vE)] 
        def bounds_vX():
            return [vmin, vE+vesc]
        def eta(cosq,vx):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_cosq,bounds_vX])[0]
    else:
        return 0

def etaDebris(vmin, params):
    vE = params[0]
    vflow = params[1]

    if vmin < vflow-vE:
        return 1/vflow
    elif vflow-vE <= vmin < vflow+vE:
        return (vflow+vE-vmin)/(2*vflow*vE)
    else:
        return 0

def etaFile(vmin,infile):
    """
    takes as input the infile which gives {v,eta(v)}
    outputs eta(vmin)
    """
    data = np.loadtxt(infile)

    x = data[:,0]
    y = data[:,1]

    # Define an interpolation function
    func = interpolate.interp1d(x,y,kind='linear')
    
    return func(vmin)

def etaFromF(vmin,infile): 
    """
    takes as input the infile which gives {v,f(v)}
    outputs function eta(vmin)
    
    n.b. takes vE = 232 km/s and vesc = 544 km/s
    """
    vE = 232*kmpers_nu
    vesc = 544*kmpers_nu
    
    data = np.loadtxt(infile)
    x = data[:,0]
    y = data[:,1]
    func = interpolate.interp1d(x,y,kind='linear') # interpolate to get f(v)
    
    " calculate normalization constant "
    def inttest(vx):
        if vx <= vesc:
            return vx**2*func(vx)
        else:
            return 0
    def bounds():
        return [0.,vesc]
    K_=4*np.pi*nquad(inttest,[bounds])[0]   
    
    if vmin <= vesc - vE:
        def bounds_cosq():
            return [-1,1]
        def bounds_vX(cosq):
            return [vmin, -cosq*vE+np.sqrt((cosq**2-1)*vE**2+vesc**2)]
        def eta(vx,cosq):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_vX,bounds_cosq])[0]
        
    elif vesc - vE < vmin <= vesc + vE:
        def bounds_cosq(vx):
            return [-1, (vesc**2-vE**2-vx**2)/(2*vx*vE)] 
        def bounds_vX():
            return [vmin, vE+vesc]
        def eta(cosq,vx):
            return (2*np.pi/K_)*vx*func(vx**2+vE**2+2*vx*vE*cosq)
        return nquad(eta, [bounds_cosq,bounds_vX])[0]
    else:
        return 0