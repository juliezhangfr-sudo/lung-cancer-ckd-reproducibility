from pathlib import Path
import csv, json

# Run from the repository root. Uses only the bundled aggregate source files.
root = Path.cwd()
source = root / 'analysis/source'
desk = source / 'wonder_desktop'
old = source / 'archived_analysis'
urban_files = source / 'urbanicity_wonder'
D = json.loads((old / 'data.json').read_text(encoding='utf-8'))
(root / 'analysis/results').mkdir(parents=True, exist_ok=True)
inp=root/'analysis/input';inp.mkdir(exist_ok=True)
rows=list(csv.DictReader((old/'model_data.csv').open()))
with (inp/'national_age_counts.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['year','age','lung','ckd']);w.writeheader()
 for r in rows:
  age=r['age'].replace(' years','');w.writerow(dict(year=r['year'],age=age,lung=r['n'],ckd=r['k']))
# Record exactly which cached model cells agree with archived source records.
checks=[]
for r in rows:
 y=int(r['year']);a=r['age'];era='1999_2020' if y<=2020 else '2018_2024'
 for k,out in [('all','n'),('CKD','k')]:
  name=f'Lung_{k}_'+('LT55_' if a=='<55' else 'age_')+era+'.csv'
  s=[x for x in D[name] if x.get('Year','').strip()==str(y) and (a=='<55' or x.get('Ten-Year Age Groups','').strip()==a)]
  assert len(s)==1,(name,y,a)
  val=int(s[0]['Deaths']);assert val==int(r[out]),(name,y,a,val,r[out])
  checks.append(dict(year=y,age=a,kind=k,count=val,archive_source=name,match=True,raw_file_present=(desk/name).exists()))
with (root/'analysis/results/national_cell_provenance.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=checks[0]);w.writeheader();w.writerows(checks)
# Native overall annual inputs are tab-separated and preserve the original rate SE.
allrows=[]
for name in ['Lung_CKD_1999_2020.csv','Lung_CKD_2018_2024.csv']:
 for r in D[name]:
  y=int(r['Year']);
  if name.endswith('2018_2024.csv') and y<=2020:continue
  allrows.append(dict(series='overall_CKD',year=y,rate=r['Age Adjusted Rate'],se=r['Age Adjusted Rate Standard Error'],reliable='TRUE'))
with (inp/'rate_series.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
# All reported definitions and original detailed records are retained separately.
tidy=[]
for filename,rr in D.items():
 for r in rr:
  y=r.get('Year','').strip()
  if not y.isdigit():continue
  q=dict(source=filename,year=y,deaths=r.get('Deaths',''),population=r.get('Population',''),rate=r.get('Age Adjusted Rate',r.get('Crude Rate','')),se=r.get('Age Adjusted Rate Standard Error',r.get('Crude Rate Standard Error','')),age=r.get('Ten-Year Age Groups','').strip(),sex=r.get('Sex','').strip(),region=r.get('Census Region','').strip(),state=r.get('State','').strip())
  tidy.append(q)
with (inp/'archived_records.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=tidy[0]);w.writeheader();w.writerows(tidy)
# Preserve all valid urbanicity year-age detail rows, retain Not Available separately.
urban={}
for kind,pattern in [('lung','C34 all lung*csv'),('ckd','C34 plus N18*csv')]:
 p=next(urban_files.glob(pattern))
 rr=list(csv.DictReader(p.open(encoding='utf-8-sig')))
 urban[kind]={(r['Year Code'],r['2023 Urbanization'],r['Ten-Year Age Groups Code']):r['Deaths'] for r in rr if r.get('Notes')=='' and r.get('Year Code') and r.get('Ten-Year Age Groups Code')}
print('Urban unmatched cells',urban['lung'].keys()-urban['ckd'].keys(),urban['ckd'].keys()-urban['lung'].keys())
assert {k for k in urban['lung'] if k[1]!='Not Available'}=={k for k in urban['ckd'] if k[1]!='Not Available'}
ur=[]
for (y,u,a),n in urban['lung'].items():ur.append(dict(year=y,urban=u,age=a,lung=n,ckd=urban['ckd'].get((y,u,a),'NOT_RETURNED')))
with (inp/'urbanicity_age_counts.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=ur[0]);w.writeheader();w.writerows(ur)
print('Prepared national cells',len(rows),'provenance checks',len(checks),'urbanicity cells including Not Available',len(ur))
