#!/usr/bin/env python3
"""Project this STEP short virtual-neutral background, invalidating rate caches.

Eirene x=R,y=Z; toroidal bulk flow changes sign with the odd B2 parallel
velocity in this no-drift magnetic reflection. Preserve particle number and
thermal+ordered energy of each virtual species for equal-volume mirror cells.
The native short reader accepts the legacy moments-only record and rebuilds
collision data rather than loading caches inconsistent with the new moments.
"""
import json,sys
from pathlib import Path
import numpy as np
from scipy.io import FortranFile

def project(source,geometry,destination):
 mapping=np.array(json.loads((geometry/'mesh_generation.json').read_text())['triangle_reflection_zero_based'])
 with FortranFile(source,'r') as f:a=f.read_reals(np.float64)
 nc=a.size//45;assert nc==len(mapping)+6 and a.size==45*nc
 m=np.r_[mapping,np.arange(len(mapping),nc)];assert np.array_equal(m[m],np.arange(nc))
 old=a.reshape((9,nc,5),order='F');new=old.copy();n=old[:,:,1];den=n+n[:,m];assert np.all(n>=0)
 safe=np.where(den>0,den,1.0)
 new[:,:,1]=.5*den
 for k,sign in [(2,1),(3,-1),(4,-1)]:new[:,:,k]=(n*old[:,:,k]+sign*n[:,m]*old[:,m,k])/safe
 mass=np.array([2,2,2,4,4,4,4,4,4])[:,None]*1.67262192369e-27
 oldK=.5*mass*(old[:,:,2:]**2).sum(2)*1e-4/1.602176634e-19
 newK=.5*mass*(new[:,:,2:]**2).sum(2)*1e-4/1.602176634e-19
 new[:,:,0]=(n*(1.5*old[:,:,0]+oldK)+n[:,m]*(1.5*old[:,m,0]+oldK[:,m]))/safe/1.5-newK/1.5
 new[:,:,0]=np.where(den>0,new[:,:,0],.5*(old[:,:,0]+old[:,m,0]))
 # Compute once per pair and copy to eliminate arithmetic-order differences.
 for i,j in enumerate(m):
  if i<j:
   new[:,j,:]=new[:,i,:];new[:,j,3:]*=-1
 for k,sign in [(0,1),(1,1),(2,1),(3,-1),(4,-1)]:assert np.max(np.abs(new[:,:,k]-sign*new[:,m,k]))==0
 assert np.all(np.isfinite(new)) and np.all(new[:,:,:2]>=0)
 with FortranFile(destination,'w') as f:f.write_record(new.ravel(order='F'))
 result={'source':str(source),'destination':str(destination),'parity_max_error':0,'virtual_species_mass_amu':[2,2,2,4,4,4,4,4,4],'cells':nc,'physical_triangles':len(mapping),'auxiliary_cells':'self mapped; odd velocities zero','collision_cache':'omitted: supported moments-only short restart; rates regenerated','density_sum_relative_error':float(abs(new[:,:,1].sum()-n.sum())/n.sum())}
 destination.with_name('neutral_projection.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':print(json.dumps(project(*map(Path,sys.argv[1:])),indent=2))
