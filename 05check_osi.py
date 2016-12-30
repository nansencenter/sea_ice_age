import glob
from nansat import *

from iceagelib import *


res = 62500
h = 60 * 60 * 24
factor=5

ifiles = sorted(glob.glob('/Data/sat/downloads/osi405c_demo_archive/*_201[3,4,5]*'))

n0 = Nansat(ifiles[0])
u0,v0,f0 = read_uv_osi(ifiles, 249)
u1,v1,f1 = read_uv_osi(ifiles, 250)
u2,v2,f2 = read_uv_osi(ifiles, 251)
raise

e0 = f0==3
e1 = f1==3
e2 = f2==3

lm = (f0==1).astype(np.float32)
lm[lm==0] = np.nan

r0,r1 = 50,130
vmin=-0.1
vmax=0.1
plt.subplot(2,3,1);plt.imshow(u0[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')
plt.text(50, 10, '2013-09-14')
plt.subplot(2,3,2);plt.imshow(u1[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')
plt.text(50, 10, '2013-09-15')
plt.subplot(2,3,3);plt.imshow(u2[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')
plt.text(50, 10, '2013-09-16')

plt.subplot(2,3,4);plt.imshow(e0[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')
plt.subplot(2,3,5);plt.imshow(e1[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')
plt.subplot(2,3,6);plt.imshow(e2[r0:r1, :], vmin=vmin, vmax=vmax, interpolation='nearest')

for i in range(1,7):
    plt.subplot(2,3,i);plt.imshow(lm[r0:r1, :], vmin=0, vmax=3, cmap='gray', interpolation='nearest')
plt.tight_layout()
plt.savefig('flag_edge.png', dpi=150, bbox_inches='tight', pad_inches=0)
plt.close('all')
