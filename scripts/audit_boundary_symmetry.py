#!/usr/bin/env python3
"""Check explicitly supplied STEP boundary values under logical reflection.

Reads this case's flat namelist assignments without guessing omitted values.
Native default and input activation checks remain separate runtime checks.
"""
import json
from pathlib import Path
import re
import sys

path,output=map(Path,sys.argv[1:])
text='\n'.join(line.split('!')[0] for line in path.read_text().splitlines())
pattern=re.compile(r'\b([A-Za-z_]\w*)\s*(\([^=]*?\))?\s*=')
assignments=list(pattern.finditer(text)); values={}
for k,match in enumerate(assignments):
    name=match[1].lower(); index=tuple(map(int,match[2].strip('()').split(','))) if match[2] else ()
    end=assignments[k+1].start() if k+1<len(assignments) else text.rfind('/')
    raw=text[match.end():end].strip().rstrip(',')
    row=[]
    for token in raw.split(','):
        token=token.strip()
        if not token:continue
        repeat,value=token.split('*') if '*' in token else ('1',token)
        value=value.strip()
        parsed=value.strip("'\"") if value.startswith(("'",'"')) else float(value.replace('D','E').replace('d','e'))
        row.extend([parsed]*int(repeat))
    values[name,index]=row
mapping={1:1,2:2,3:4,4:3,5:6,6:5,7:8,8:7,9:10,10:9,11:11,12:12}
issues=[]; checked=[]
for name in ['bcene','bceni','bcpot']:
    row=values[name,()]
    if len(row)!=12:raise ValueError('Wrong boundary shape')
    for i,j in mapping.items():
        if row[i-1]!=row[j-1]:issues.append([name,i,j])
    checked.append(name)
for name in ['enepar','enipar','potpar']:
    for (key,index),row in values.items():
        if key!=name:continue
        if index[0]!=1 or len(row)!=12:raise ValueError('Unsupported parameter shape')
        for i,j in mapping.items():
            if row[i-1]!=row[j-1]:issues.append([name,index,i,j])
        checked.append([name,index])
for name in ['bccon','bcmom','conpar','mompar']:
    records={index[1:]:row for (key,index),row in values.items() if key==name}
    for index,row in records.items():
        reflected=(mapping[index[0]],)+index[1:]
        # An absent explicit record uses native defaults; flag for runtime audit.
        if reflected not in records:issues.append([name,index,'reflected record omitted'])
        elif row!=records[reflected]:issues.append([name,index,reflected])
        if len(row)!=24:raise ValueError('Wrong species shape')
    checked.append(name)
# Check the boundary supports themselves, including self-reflecting radial sides.
def mx(x):return 71-x if x<=72 else 267-x
for i,j in mapping.items():
    side=values['bcchar',()][i-1]; other=values['bcchar',()][j-1]
    pos=int(values['bcpos',()][i-1]);start=int(values['bcstart',()][i-1]);end=int(values['bcend',()][i-1])
    if side in 'WE':expected=('E' if side=='W' else 'W',mx(pos),start,end)
    else:expected=(side,pos,mx(end),mx(start))
    actual=(other,int(values['bcpos',()][j-1]),int(values['bcstart',()][j-1]),int(values['bcend',()][j-1]))
    if expected!=actual:issues.append(['boundary support',i,j,expected,actual])
report={'input':str(path),'boundary_reflection':mapping,'explicit_parameter_checks':checked,'issues':issues,
        'all_explicit_values_pass':not issues,'scope':'Explicit boundary supports, codes and parameters; native defaults, global switches and source profiles require separate audit.'}
output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if issues:raise SystemExit(1)
