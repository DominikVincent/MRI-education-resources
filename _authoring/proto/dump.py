import json,sys,base64,os
nb=json.load(open(sys.argv[1])); out=sys.argv[2]; os.makedirs(out,exist_ok=True)
fi=0
for ci,c in enumerate(nb['cells']):
    if c['cell_type']!='code': continue
    for o in c.get('outputs',[]):
        if o.get('output_type')=='stream': print(f'--- cell {ci}:\n'+''.join(o['text'])[:3000])
        d=o.get('data',{})
        if 'text/plain' in d and 'image/png' not in d: print(f'--- cell {ci} (value):', ''.join(d['text/plain'])[:500])
        if 'image/png' in d:
            p=f'{out}/c{ci:02d}_{fi}.png'; open(p,'wb').write(base64.b64decode(d['image/png'])); fi+=1; print('IMG',p)
