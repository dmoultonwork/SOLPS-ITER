#!/usr/bin/env python3
"""Read the serialized STEP restart and independently check reflected parities."""
import json
from pathlib import Path
import sys
from symmetrise_state import read_state
from audit_mesh_symmetry import cf_fields

state, geometry, output=map(Path,sys.argv[1:])
_,rows=read_state(state)
d={name:v for kind,count,name,v in rows}
nx,ny,ns=d['nx,ny,ns']; n=(nx+2)*(ny+2)
g=cf_fields(geometry,{'leftix','rightix','rightiy'})
# block starts: -1 and the second column without a west neighbour at iy=0 (block 1 is -1..cuts[1]-1, block 2 is cuts[1]..nx)
cuts=[x for x in range(-1,nx+1) if g['leftix'][x+1+(nx+2)]==-2]
assert len(cuts)==2 and cuts[0]==-1,cuts
m=[]
for y in range(-1,ny+1):
    for x in range(-1,nx+1):
        mx=(cuts[1]-2-x) if x<cuts[1] else (cuts[1]+nx-x)
        m.append(mx+1+(nx+2)*(y+1))
assert all(m[m[i]]==i for i in range(n))
f=[]
for i,j in enumerate(m):
    if g['leftix'][i]==-2:f.append(None);continue
    f.append(int(g['rightix'][j])+1+(nx+2)*(int(g['rightiy'][j])+1))
assert all(j is None or f[j]==i for i,j in enumerate(f))
metadata={'nx,ny,ns','label','zamin','zamax','zn','am','time'}
result={}
for kind,count,name,v in rows:
    if name in metadata:continue
    discrepancies=[]
    for b in range(count//n):
        if name in {'na','ne','te','ti','po','kinrgy'}:mapping,sign=m,1
        elif name=='ua':mapping,sign=m,-1
        elif name=='uadia':mapping,sign=m,(-1 if b%2==0 else 1)
        else:
            pol=name=='fch_p' or b%2==0
            mapping,sign=(f,-1) if pol else (m,1)
        discrepancies.append(max(abs(v[i+n*b]-(0 if j is None else sign*v[j+n*b])) for i,j in enumerate(mapping)))
    error=max(discrepancies); scale=max(map(abs,v))
    result[name]={'max_absolute_parity_error':error,'max_norm_relative_error':error/scale if scale else error}
report={'state':str(state),'all_fields_pass':all(r['max_absolute_parity_error']==0 for r in result.values()),'fields':result,'face_mapping_involutive':True}
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='fields'}))
if not report['all_fields_pass']:raise SystemExit(1)
