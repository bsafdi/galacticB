import numpy as np 
import numpy.linalg as LA


#B-field model from https://arxiv.org/pdf/1204.3662.pdf

iopen=11.5 #degrees
rmx_array = np.array([5.1,6.3,7.1,8.3,9.8,11.4,12.7,15.5]) #kpc

def return_B(x,y):
    '''
    x,y in Galactic coords in kpc
    Earth at (x,y) = (-8.5,0)
    '''
    r = np.sqrt(x,y)
    phi = np.arctanM(y,x)
    
    r_hat = np.array([np.cos(phi),np.sin(phi)])
    phi_hat = np.array([-np.sin(phi),np.cos(phi)])
    
    if r<5.0:
        B = b_ring*phi_hat
        
    bv_hat = np.sin(iopen)*r_hat+np.cos(iopen)*phi_hat
    
    rs = rmx_array*np.exp(phi*np.tan(np.pi/2.-i_open))
    entry = np.searchsorted(rs,r)
    
# Use paper https://academic.oup.com/mnras/article/431/1/683/1050400

B0 = 1. # \muG
Rscale = 20.0 # kpc
hg = 6.0 #kpc
Rmax = 20.0 # kpc

Rmol = 5.0 #kpc
theta_p = -11.5*np.pi/180. #radians

def Br(r):
    return B0*np.exp(-r**2/Rscale**2)

def Bcoh(z):
    return 1/np.cosh(z/hg)**2

def Bhat(x,y):
    r = np.sqrt(x**2+y**2)
    phi = arctanM(y,x)
    if r < Rmol:
        return np.array([np.cos(phi+np.pi/2.),np.sin(phi+np.pi/2.)])
    else:
        r_hat = np.array([np.cos(phi),np.sin(phi)])
        phi_hat = np.array([-np.sin(phi),np.cos(phi)])
        #print r_hat, phi_hat
        return np.sin(theta_p)*r_hat+np.cos(theta_p)*phi_hat
    
def B_ASS(x,y,z):
    r = np.sqrt(x**2+y**2)
    return Br(r)*Bcoh(z)*Bhat(x,y)


ai = np.array([3, .5, -4, 1.2, -.8]) # spiral arm amplitudes
phi_0i = np.deg2rad(10+90*np.arange(1, 6)) # aximuthal orientation of the rotation of the spiral
Rcomp = 7.1 # kpc scale radius of compressed spiral arms
C0 = 2.5 # compression arm amplitude
rcc = 12 #kpc, region of constant compression
d0 = .3 # kpc, base width of arm enhancement
hc = 2. # kpc, scaleheight of the spiral compression
ThetaP = np.deg2rad(-11.5)
Rscale = 20. 

def spiral_arm(phi, phi0):
    beta = 1 /np.tan(ThetaP)
    radius = Rcomp * np.exp((phi-phi0) / beta)
    
    x = np.cos(phi)*radius
    y = np.sin(phi)*radius
    
    return np.array([x, y])

def min_distance(x, y, phi0):
    phi = np.arctan(float(y)/float(x)) + np.arange(-10, 11) * np.pi
    dists = np.zeros_like(phi)
    
    best_phi = -10*np.pi
    min_dist = 1e10
    
    point_loc = np.array([x, y])
    
    for i in range(len(phi)):
        spiral_arm_loc = spiral_arm(phi[i], phi0)
        dist = LA.norm(spiral_arm_loc - point_loc)
        
        if dist < min_dist:
            best_phi = phi[i]
            min_dist = dist
            
    return min_dist

def Barm(x,y,z):
    r = np.sqrt(x**2 + y**2)
    beta = 1 /np.tan(ThetaP)
    phi = np.arctan2(y, x)
    
    ri = np.zeros_like(ai)
    for i in range(len(ri)):
        ri[i] = min_distance(x, y, phi_0i[i])
    Br = B0 * np.exp(-r**2 / Rscale**2)
    cr = C0 * np.minimum(1, (r / rcc)**(-3))
    d0_r = d0 / cr / Br
    
    Bcomp = 1./np.cosh(z / hc)**2
    
    rhoC = cr * Bcomp * np.exp(-ri**2 / d0_r**2)
    
    if r < Rmol:
        BVec = np.array([np.cos(phi+np.pi/2), np.sin(phi+np.pi/2), 0])
        
    else:
        r_hat = np.array([np.cos(phi),np.sin(phi), 0])
        phi_hat = np.array([-np.sin(phi),np.cos(phi), 0])
        
        BVec = np.sin(ThetaP)*r_hat+np.cos(ThetaP)*phi_hat
        
    return np.sum((Br*ai*rhoC)[:, None]*BVec[None, :], axis = 0)

def arctanM(x,y):
    tmp = np.arctan2(x,y)
    if tmp<0:
        res= 2*np.pi+tmp
    else:
        res = tmp
    return res

######## B_ASS

# ais = np.array([3.0,0.5,-4.0,1.2,-0.8])
# phi0_is = np.array([10+90*1,10+90*2,10+90*3,10+90*4,10+90*5])*np.pi/180.
# Rcomp = 7.1 #kpc
# C0 = 2.5
# rcc =12.0 #kpc
# d0 = 0.3 #kpc
# hc = 2.0 #kpc

# def arctanM(x,y):
#     tmp = np.arctan2(x,y)
#     if tmp<0:
#         res= 2*np.pi+tmp
#     else:
#         res = tmp
#     return res
        

# def Bcomp(z):
#     return 1./np.cosh(z/hc)**2

# def c(r):
#     if r<rcc:
#         return C0
#     else:
#         return C0*(r/rcc)**(-3.)

# def d0f(r):
#     return d0/(c(r)*Br(r))

# def ri(phi,i):
#     return 7.1*np.exp((phi-phi0_is[i-1])*np.tan(theta_p))

# def di(r,phi,i):
#     riA = ri(phi,i)
#     return np.abs(r-riA)

# def rhoc(x,y,z,d,i):
#     r = np.sqrt(x**2+y**2)
#     phi = arctanM(y,x)
#     rI = ri(phi,i)
#     return c(rI)*Bcomp(z)*np.exp(-d**2/d0f(rI)**2)

# def Barmi(x,y,z,i):
#     r = np.sqrt(x**2+y**2)
#     phi = arctanM(y,x)
    
#     Bi = Br(r)
#     ai = ais[i-1]
#     #print ai
#     d = di(r,phi,i)
#     #print d, d0f(ri(phi,i))
#     rhoci = rhoc(x,y,z,d,i)
#     #print d, ri(phi,i), rhoci #, d0f(r)
#     Bh = Bhat(x,y)
#     return Bi*ai*rhoci*Bh

# def Barm(x,y,z):
#     res = np.zeros(2)
#     for i in [1,2,3,4,5]:
#         res += Barmi(x,y,z,i)
#     return res




