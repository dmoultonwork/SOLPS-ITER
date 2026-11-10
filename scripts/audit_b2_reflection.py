#!/usr/bin/env python3
"""Read-only: how far is a B2 geometry (b2fgmtry) from exact up-down mirror symmetry for given reflection blocks?
usage: audit_b2_reflection.py b2fgmtry LO:HI [LO:HI ...]
Each block is reflected onto itself, x -> LO+HI-x (same radial index); corners swap left/right (0<->1, 2<->3); R equal, Z opposite."""
import sys
from pathlib import Path
from audit_mesh_symmetry import cf_fields
f=cf_fields(Path(sys.argv[1]),{'nx,ny','crx','cry','vol','hx','hy','leftix'})
nx,ny=f['nx,ny'];W=nx+2;n=W*(ny+2)
blocks=[tuple(map(int,b.split(':'))) for b in sys.argv[2:]]
xmap={}
for lo,hi in blocks:
    for x in range(lo,hi+1):xmap[x]=lo+hi-x
assert set(xmap)==set(range(-1,nx+1)),'blocks must cover -1..nx'
swap=[1,0,3,2];crx,cry,vol,hx,hy=f['crx'],f['cry'],f['vol'],f['hx'],f['hy']
dR=dZ=dv=dhx=dhy=0.;wR=wZ=None
for y in range(ny+2):
    for x in range(nx+2):
        px=xmap[x-1]+1;i=x+W*y;j=px+W*y
        if i>j:continue
        for k in range(4):
            a=crx[i+n*k]-crx[j+n*swap[k]];b=cry[i+n*k]+cry[j+n*swap[k]]
            if abs(a)>dR:dR=abs(a);wR=(x-1,y-1)
            if abs(b)>dZ:dZ=abs(b);wZ=(x-1,y-1)
        if x-1 in range(0,nx) and y-1 in range(0,ny) and (vol[i]>0 and vol[j]>0):
            dv=max(dv,abs(vol[i]-vol[j])/(0.5*(vol[i]+vol[j])))
            dhx=max(dhx,abs(hx[i]-hx[j])/(0.5*(hx[i]+hx[j])));dhy=max(dhy,abs(hy[i]-hy[j])/(0.5*(hy[i]+hy[j])))
print('nx,ny',nx,ny,'blocks',blocks)
print('max corner mismatch: R %.3e m at %s, Z %.3e m at %s'%(dR,wR,dZ,wZ))
print('max relative mismatch: volume %.3e, hx %.3e, hy %.3e'%(dv,dhx,dhy))
