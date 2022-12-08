import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d, RectBivariateSpline

number_of_times=21 # 21 time snapshots including t=0
number_of_amr_levels=4 # levels of AMR to feed BC
# note: x and y here are meant as coordinates on the 2D boundary
#       which we want to drive in a time-dependent fashion
#       we can allow unequal x-y numbers
Nx=8 # total grid resolution in x at level 1, so domain_nx2
Ny=20 # total grid resolution in y at level 1, so domain_nx3
for t in range(0,number_of_times): 
    dx=1./Nx
    dy=1./Ny
    x=np.zeros(Nx)
    y=np.zeros(Ny)
    rho1=np.zeros((Nx,Ny))
    oneline=' {: 0.6E} \n'
    output=open('bc/bc_t'+"%02d" % (t)+'_AMR1.dat','w+')
    for i in range(Nx):
       x[i]=dx/2.+dx*i
       for j in range(Ny):
          y[j]=dy/2.+dy*j
          rho1[i,j]=random.uniform(0.,2.)
          output.write(oneline.format(rho1[i,j]))
    output.close()

    # add a layer in x and y, for plotting purposes
    xx,yy=np.meshgrid(np.append(x,x[Nx-1]+x[1]-x[0]),np.append(y,y[Ny-1]+y[1]-y[0]),indexing='ij')
    plt.pcolor(xx,yy,rho1)
    plt.savefig('bc/bc_t'+"%02d" % (t)+'_AMR1.png')
    plt.close()

    f=RectBivariateSpline(x,y,rho1) # signal on former mesh
    if number_of_amr_levels> 1:
       for ilev in range(2,number_of_amr_levels+1):
          Nx2=Nx*(2**(ilev-1))
          Ny2=Ny*(2**(ilev-1))
          dx=dx/2.
          dy=dy/2.
          xnew=np.zeros(Nx2)
          ynew=np.zeros(Ny2)
          rho2=np.zeros((Nx2,Ny2))
          for i in range(Nx2):
             xnew[i]=dx/2.+dx*i
             for j in range(Ny2):
                ynew[j]=dy/2.+dy*j
          rho2=f(xnew,ynew)
          output=open('bc/bc_t'+"%02d" % (t)+'_AMR'+"%01d" % (ilev)+'.dat','w+')
          for i in range(Nx2):
             for j in range(Ny2):
                output.write(oneline.format(rho2[i,j]))
          output.close()

          xx,yy=np.meshgrid(np.append(xnew,xnew[Nx2-1]+xnew[1]-xnew[0]),np.append(ynew,ynew[Ny2-1]+ynew[1]-ynew[0]),indexing='ij')
          plt.pcolor(xx,yy,rho2)
          plt.savefig('bc/bc_t'+"%02d" % (t)+'_AMR'+"%01d" % (ilev)+'.png')
          plt.close()
