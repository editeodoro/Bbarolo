   import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
import sys


outfolder = str(sys.argv[1])
file1 = outfolder+"ringlog1.txt"
file2 = outfolder+"ringlog2.txt"
outplot=outfolder+"plot_parameters.pdf"
twostage = errors = False
if (len(sys.argv)>2 and sys.argv[2].lower()=='true'): twostage = True
if (len(sys.argv)>3 and sys.argv[3].lower()=='true'): errors = True

rad,vrot,disp,inc,pa,z0,xpos,ypos,vsys= np.loadtxt(file1,skiprows=1,usecols=(1,2,3,4,5,7,9,10,11),unpack=True)
color='#B22222'

if twostage==True: 
	vrot2,disp2,inc2,pa2,z02,xpos2,ypos2,vsys2=np.loadtxt(file2,skiprows=1,usecols=(2,3,4,5,7,9,10,11),unpack=True)
	color='#A0A0A0'
	color2='#B22222'

if errors==True:
	vrot_err1,vrot_err2 =np.loadtxt(file1,skiprows=1,usecols=(12,13),unpack=True)
	disp_err1,disp_err2 =np.loadtxt(file1,skiprows=1,usecols=(14,15),unpack=True)
	
	
max_rad = np.max(rad)

plt.rc('font',family='sans-serif',serif='Helvetica',size=10) 

fig1=plt.figure(figsize=(8.27,11.69), dpi=100) 

grid = gridspec.GridSpec(4,2, height_ratios=[2,1,1,1]) 
grid.update(top=0.95, bottom=0.10, left=0.12, right=0.90, wspace=0.26, hspace=0.1)

ax=plt.subplot(grid[0,0])
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.xlim(0,1.2*max_rad)
plt.ylim(0,1.2*np.max(vrot))
plt.ylabel("v$_\mathrm{rot}$ (km/s)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
if errors==True: plt.errorbar(rad,vrot, yerr=[vrot_err1,-vrot_err2],fmt='o', color=color)
else: plt.plot(rad,vrot,'o', color=color)
if twostage==True: 
	if errors==True: plt.errorbar(rad,vrot2, yerr=[vrot2_err1,-vrot2_err2],fmt='o', color=color2)
	else : plt.plot(rad,vrot2,'o', color=color2)

ax=plt.subplot(grid[1,0])
plt.xlim(0,1.2*max_rad)
plt.ylabel("i (deg)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,inc,'o', color=color)
if twostage==True: plt.plot(rad,inc2,'o', color=color2)

ax=plt.subplot(grid[2,0])
plt.xlim(0,1.2*max_rad)
plt.ylabel("$\phi$ (deg)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,pa,'o', color=color)
if twostage==True: plt.plot(rad,pa2,'o', color=color2)

ax=plt.subplot(grid[3,0])
plt.xlim(0,1.2*max_rad)
plt.ylabel("z$_0$ (arcs)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.xlabel("Radius (arcsec)", fontsize=14, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='on',labelleft='on') 
plt.plot(rad,z0,'o', color=color)
if twostage==True: plt.plot(rad,z02,'o', color=color2)



ax=plt.subplot(grid[0,1])
plt.xlim(0,1.2*max_rad)
plt.ylim(0,1.2*np.max(disp))
plt.ylabel("$\sigma_\mathrm{gas}$  (km/s)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
if errors==True : plt.errorbar(rad,disp,yerr=[disp_err1,-disp_err2],fmt='o',color=color)
else: plt.plot(rad,disp,'o', color=color)
if twostage==True: 
	if errors==True: plt.errorbar(rad,disp2,yerr=[disp2_err1,-disp2_err2],fmt='o',color=color2)
	else: plt.plot(rad,disp2,'o', color=color2)

ax=plt.subplot(grid[1,1])
plt.xlim(0,1.2*max_rad)
plt.ylabel("x$_0$ (pix)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,xpos,'o', color=color)
if twostage==True: plt.plot(rad,xpos2,'o', color=color2)

ax=plt.subplot(grid[2,1])
plt.xlim(0,1.2*max_rad)
plt.ylabel("y$_0$ (pix)", fontsize=14)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,ypos,'o', color=color)
if twostage==True: plt.plot(rad,ypos2,'o', color=color2)

ax=plt.subplot(grid[3,1])
plt.xlim(0,1.2*max_rad)
plt.ylabel("v$_\mathrm{sys}$ (km/s)", fontsize=14)
plt.xlabel("Radius (arcsec)", fontsize=14, labelpad=10)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='on',labelleft='on') 
plt.plot(rad,vsys,'o', color=color)
if twostage==True: plt.plot(rad,vsys2,'o', color=color2)

plt.savefig(outplot, orientation = 'landscape', format = 'pdf') 


"""
fig1=plt.figure(figsize=(8.27, 11.69), dpi=100) 

grid = gridspec.GridSpec(3,1, height_ratios=[2,1,1]) 
grid.update(top=0.95, wspace=0.0, hspace=0.1)

plt.subplot(grid[0,0])
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.ylim(0,1.2*np.max(vrot))
plt.ylabel("$V_\mathrm{rot} \; (km/s)$", fontsize=16, labelpad=10)
plt.plot(rad,vrot,'o', color='#B22222')

plt.subplot(grid[1,0])
plt.ylabel("$i \; (deg)$", fontsize=16, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,inc,'o', color='#B22222')

plt.subplot(grid[2,0])
plt.xlabel("$Radius \; (arcsec)$", fontsize=16,labelpad=10)
plt.ylabel("$\phi \; (deg)$", fontsize=16, labelpad=10)
plt.plot(rad,pa,'o', color='#B22222')

plt.savefig('plot.pdf', orientation = 'portrait', format = 'pdf') 
plt.show()

fig2=plt.figure(figsize=(8.27, 11.69), dpi=100)

grid = gridspec.GridSpec(3,1, height_ratios=[1.5,1,1]) 
grid.update(top=0.98, bottom=0.05,left=0.10,right=0.95, wspace=0.0, hspace=0.1)

plt.subplot(grid[0,0])
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.ylim(0,1.2*np.max(disp))
plt.ylabel("$\sigma_\mathrm{gas} \; (km/s)$", fontsize=16, labelpad=10)
plt.plot(rad,disp,'o', color='#B22222')

plt.subplot(grid[1,0])
plt.ylabel("$z_0 \; (arcs)$", fontsize=16, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,z0,'o', color='#B22222')

plt.subplot(grid[2,0])
plt.xlabel("$Radius \; (arcsec)$", fontsize=16,labelpad=10)
plt.ylabel("$i \; (deg)$", fontsize=16, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,inc,'o', color='#B22222')

plt.savefig('plot2.pdf', orientation = 'portrait', format = 'pdf') 
plt.show()

fig3=plt.figure(figsize=(8.27, 11.69), dpi=100)

grid = gridspec.GridSpec(3,1, height_ratios=[1,1,1]) 
grid.update(top=0.95, wspace=0.0, hspace=0.1)

plt.subplot(grid[0,0])
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.ylabel("$x_0 \; (pix)$", fontsize=16, labelpad=10)
plt.plot(rad,xpos,'o', color='#B22222')

plt.subplot(grid[1,0])
plt.ylabel("$y_0 \; (pix)$", fontsize=16, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,ypos,'o', color='#B22222')

plt.subplot(grid[2,0])
plt.xlabel("$Radius \; (arcsec)$", fontsize=16,labelpad=10)
plt.ylabel("$V_\mathrm{sys} \; (km/s)$", fontsize=16, labelpad=10)
plt.tick_params(axis='both',which='both',bottom='on',top='on',labelbottom='off',labelleft='on') 
plt.plot(rad,vsys,'o', color='#B22222')


plt.savefig('plot3.pdf', orientation = 'portrait', format = 'pdf') 
plt.show()
"""
