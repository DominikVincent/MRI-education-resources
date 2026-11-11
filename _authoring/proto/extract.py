import sys, json, base64, os
nb=json.load(open(sys.argv[1])); out=sys.argv[2]; os.makedirs(out,exist_ok=True); k=0
for i,c in enumerate(nb['cells']):
    for o in c.get('outputs',[]):
        d=o.get('data',{})
        if 'image/png' in d:
            k+=1; open(f"{out}/c{i:02d}_{k}.png",'wb').write(base64.b64decode(d['image/png']))
        if o.get('output_type')=='stream': print(f"--- cell {i}:", o['text'][:800])
print(sorted(os.listdir(out)))
