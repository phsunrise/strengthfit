import sys
import numpy as np
import matplotlib.pyplot as plt

def voigt(chi_deg, exx, eyy, ezz, d0, wavelength, \
          phi_array=[], npoints=300, tol=1.e-12):
    # exx, eyy, ezz are the compression strains,
    # d0 is the lattice plane distance (uncompressed) 
    
    ## calculate parameters
    chi = chi_deg / 180. * np.pi
    R = np.matrix([[np.cos(chi), 0, -np.sin(chi)],
                   [     0,      1,     0       ],
                   [np.sin(chi), 0,  np.cos(chi)]])
            # this defines a clockwise rotation about y-axis by chi
            # to transform from X-ray coord to sample coord
    a = np.array(np.diag([1.-exx, 1.-eyy, 1.-ezz]).T.dot(R))
    ## calculate the coefficients
    A1 = np.linalg.norm(a[:,0])**2
    A2 = a[:,0].dot(a[:,1])
    A3 = a[:,0].dot(a[:,2])
    A4 = np.linalg.norm(a[:,1])**2
    A5 = a[:,1].dot(a[:,2])
    A6 = np.linalg.norm(a[:,2])**2

    ## calculate thetaB for each phi
    tthB_array = []
    if phi_array == []:
        phi_array = np.linspace(-np.pi, np.pi, npoints)
    for phi in phi_array: 
        if exx == ezz:
            tthB_array.append(2.*np.arcsin(wavelength/2./d0/(1.-exx)))
            continue
        D1 = ( A1*np.cos(phi)**2 
             + 2*A2*np.cos(phi)*np.sin(phi)
             + A4*np.sin(phi)**2 )
        D2 = 2*A3*np.cos(phi) + 2*A5*np.sin(phi)
        D3 = A6
        ## coefficients for quartic equation
        coeff = [D3-D1-1.j*D2, -(4.*D3-2.j*D2), \
                 2.*D1+6.*D3-4.*wavelength**2/d0**2, \
                 -(2.j*D2+4.*D3), D3+1.j*D2-D1]
        roots = np.roots(coeff)
        ## find the correct roots
        sol_exist = False
        for val in np.log(roots):
            if np.abs(val.real) > tol:
                continue
            tthB = val.imag
            if tthB < 0.:
                tthB += 2.*np.pi 
            if 0. < tthB and tthB < np.pi: 
                # root needs to be between 0 and pi
                sol_exist = True
                break

        if sol_exist:
            tthB_array.append(tthB)
        else:
            raise RuntimeError("no solution for phi=%f" % phi)

    return phi_array, np.array(tthB_array)

if __name__=='__main__':
    ## settings
    chi_deg = -30.
    t_array = [0, 50, 100] # GPa
    ## values in paper
    #exx = [0.0937, 0.0805, 0.0667]
    #ezz = [0.093705, 0.121,  0.148]
    ## another set of values
    exx = [0.08906, 0.07881, 0.06857]
    ezz = [0.0890605, 0.1095,  0.1300]
    eyy = exx
    h, k, l = 2, 2, 0
    a0 = 3.56683 # diamond at 300K, in Angstrom
    E_xray = 10. # keV
    tol = 1.e-12 # tolerance for root finding
    colors = ['r', 'g', 'b']

    ## calculate parameters
    d0 = a0 / np.sqrt(h**2+k**2+l**2)
    wavelength = 12.3984193 / E_xray 

    for i_t, t in enumerate(t_array):
        phi, tthB = voigt(chi_deg, exx[i_t], eyy[i_t], ezz[i_t], \
                d0, wavelength, npoints=300, tol=tol)
        plt.plot(phi/np.pi*180., tthB/np.pi*180.,\
                color=colors[i_t], label=r"$t=%d$ GPa" % t)
        np.savez("voigt_t_%d.npz" % t, \
                phi=phi/np.pi*180., tthB=tthB/np.pi*180.)
        
    plt.xlabel(r"$\phi$ (deg)")
    plt.ylabel(r"$2\theta$ (deg)")
    plt.xlim(-180., 180.)
    plt.ylim(58., 70.)
    plt.legend()
    plt.savefig("voigt.pdf")
    plt.show()
