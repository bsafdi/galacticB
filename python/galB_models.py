import numpy as np 


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
    