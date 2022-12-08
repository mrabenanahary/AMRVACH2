#--------------------------------------
import sys
#select the path to the python file
import synchro_v2 as synchro

#size and resolution of the box
reslight=[528,1056,528]
sizelight=[100,200,100]

#select snapshot
fileid=05
print('start opening the file')
#load the file (replace datafile with the name of the file)
pd=synchro.get_pointdata(fileid,filenameout='jet2D_srmhd_choc_',type='vtu')
#select phi and theta angles
phi=0
theta=10
print('calculating emission')
em=synchro.emiss(pd,phi,theta,4,0.6,recipe=2,polarized=0,delta=0)
map = shotgun(em,reslight,)sizelight
show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'theta',str(theta_i)[0:4]]))
