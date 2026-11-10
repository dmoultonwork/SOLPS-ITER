#!/usr/bin/env python3
"""Mirror STEP surface/source coordinates, preserving surface IDs and physics.

Label permutations are reported separately; their physical use still requires
native interface validation. Does not copy neutral restart data.
"""
import argparse
import copy
import json
from pathlib import Path


def prepare(reference, pairing, destination, mesh):
    surface_doc=json.loads((reference/'eirene_add_surfaces.input.json').read_text())
    surfaces=surface_doc['ADDITIONAL_SURFACES']['SURFACES']
    original=copy.deepcopy(surfaces)
    pairs={int(k):v for k,v in json.loads(pairing.read_text())['candidate_pairs'].items()}
    if sorted(pairs)!=list(range(1,len(surfaces)+1)) or any(pairs[pairs[i]]!=i for i in pairs):
        raise ValueError('Surface reflection must be a complete involution')
    label_map={-4:-3,-3:-4,5:6,6:5}
    changes=[]
    for i,j in pairs.items():
        if i>j: continue
        a,b=surfaces[i-1],surfaces[j-1]
        for k in a:
            if k in {'ILIM','TXTSFL','P1','P2'}: continue
            expected=label_map.get(a[k],a[k]) if k=='ILPLG' else a[k]
            if expected!=b[k]: raise ValueError(f'Surface physics mismatch {i} {j} {k}')
        if i==j:
            if a['P1'][0]!=a['P2'][0] or a['P1'][1]!=-a['P2'][1]:
                raise ValueError('Unexpected asymmetric self-reflecting surface')
            continue
        lower,upper=(a,b) if a['P1'][1]+a['P2'][1]<0 else (b,a)
        for end,other in [('P1','P2'),('P2','P1')]:
            upper[end][:2]=[lower[other][0],-lower[other][1]]
    for before,after in zip(original,surfaces):
        if before!=after: changes.append({'surface':after['ILIM'],'before':[before['P1'],before['P2']],'after':[after['P1'],after['P2']]})
    doc=json.loads((reference/'eirene.input.json').read_text())
    source_changes=[]
    for stratum in doc['SOURCES']['STRATA']:
        rows=stratum['SUBSTRATA']
        if not stratum['NLSRF']: continue
        by_surface={row['INSOR']:row for row in rows}
        if len(by_surface)!=len(rows): raise ValueError('Duplicate source surface')
        for i,a in by_surface.items():
            j=pairs[i]
            if j not in by_surface: raise ValueError('Missing reflected source')
            if i>j:continue
            b=by_surface[j]
            if i==j:
                if a['SORAD1']!=a['SORAD2'] or a['SORAD3']!=-a['SORAD4']:
                    raise ValueError('Unexpected self source')
                continue
            for key in a:
                if key not in {'INUM','INSOR','SORAD1','SORAD2','SORAD3','SORAD4'} and a[key]!=b[key]:
                    raise ValueError(f'Source parameter mismatch {i} {j} {key}')
            lower,upper=(a,b) if a['SORAD3']+a['SORAD4']<0 else (b,a)
            before=copy.deepcopy(upper)
            for k in [1,2]:upper[f'SORAD{k}']=lower[f'SORAD{k}']
            for k in [3,4]:upper[f'SORAD{k}']=-lower[f'SORAD{k}']
            if before!=upper: source_changes.append({'stratum':stratum['ISTRA'],'surface':upper['INSOR'],'before':before,'after':copy.deepcopy(upper)})
    doc['STANDARD_MESH']['RADIAL_GRID']['NRKNOT']=int((mesh/'fort.33').read_text().split()[0])
    destination.mkdir(parents=True,exist_ok=True)
    for name,value in [('eirene_add_surfaces.input.json',surface_doc),('eirene.input.json',doc)]:
        path=destination/name
        if path.exists():raise FileExistsError(path)
        path.write_text(json.dumps(value,indent=2)+'\n')
    (destination/'fort.1').write_text(json.dumps(doc,indent=2)+'\n')
    physics=json.loads((reference/'eirene_physics_model.input.json').read_text())
    (destination/'eirene_physics_model.input.json').write_text(json.dumps(physics,indent=2)+'\n')
    report={'surface_pairs':pairs,'label_permutation':label_map,'surface_coordinate_changes':changes,
            'source_coordinate_changes':source_changes,'source_weights':'unchanged; paired weights checked equal',
            'neutral_restart':'not copied','remaining_validation':['Native wall/triangle intersection and label consistency','Background properties of reflected reservoir/polygon labels','B2 boundary/numerical/source inputs']}
    (destination/'input_projection.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'surfaces_changed':len(changes),'sources_changed':len(source_changes),'destination':str(destination)}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('pairing',type=Path);p.add_argument('destination',type=Path);p.add_argument('mesh',type=Path)
    a=p.parse_args();prepare(a.reference,a.pairing,a.destination,a.mesh)
