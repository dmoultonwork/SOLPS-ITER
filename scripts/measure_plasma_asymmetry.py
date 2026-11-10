#!/usr/bin/env python3
"""Volume-weighted STEP plasma inventory and mirror residual diagnostics."""
import json
import math
from pathlib import Path
import sys
from symmetrise_state import read_state
from audit_mesh_symmetry import cf_fields
state,geometry,output=map(Path,sys.argv[1:])
_,fields=read_state(state);d={name:v for kind,count,name,v in fields}
nx,ny,ns=d['nx,ny,ns'];n=(nx+2)*(ny+2)
g=cf_fields(geometry,{'vol','cry','leftix'})
cuts=[x for x in range(-1,nx+1) if g['leftix'][x+1+(nx+2)]==-2]
assert len(cuts)==2 and cuts[0]==-1,cuts
ids=[x+1+(nx+2)*(y+1) for y in range(ny) for x in range(nx) if x not in [cuts[1]-1,cuts[1]]]
z={i:math.fsum(g['cry'][i+n*k] for k in range(4))/4 for i in ids}
m={i:((cuts[1]-2-(i%(nx+2)-1)) if (i%(nx+2)-1)<cuts[1] else (cuts[1]+nx-(i%(nx+2)-1)))+1+(nx+2)*(i//(nx+2)) for i in ids}
assert all(m[m[i]]==i for i in ids)
results={}
for name,species in [('D_all',[0,1]),('D_ions',[1]),('Ar_all',list(range(2,21))),('Ar_ions',list(range(3,21))),('He_all',[21,22,23]),('He_ions',[22,23])]:
    if any(x>=ns for x in species): continue  # species set not present in this case (e.g. no He)
    den={i:math.fsum(d['na'][i+n*s] for s in species) for i in ids}
    upper=math.fsum(den[i]*g['vol'][i] for i in ids if z[i]>0)
    lower=math.fsum(den[i]*g['vol'][i] for i in ids if z[i]<0)
    total=math.fsum(den[i]*g['vol'][i] for i in ids)
    residual=math.fsum(abs(den[i]-den[m[i]])*g['vol'][i] for i in ids)/total
    results[name]={'upper_inventory':upper,'lower_inventory':lower,'signed_upper_minus_lower_over_total':(upper-lower)/total,'volume_weighted_mirror_L1_over_inventory':residual}
report={'state':str(state),'geometry':str(geometry),'physical_cells':len(ids),'species':results,
        'caveat':'B2 fluid-state inventory only; does not include Eirene neutral molecules or populations stored in virtual BGK backgrounds. Small geometric volume asymmetry can give a roundoff signed inventory offset even for exactly even densities.'}
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:round(v['signed_upper_minus_lower_over_total'],6) for k,v in results.items()}))
