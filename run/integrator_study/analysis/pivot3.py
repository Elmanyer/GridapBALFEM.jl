# run from $IS_OUT (reads summary.tsv, metrics.tsv, events2.tsv there)
import csv, sys
S={r['name']:r for r in csv.DictReader(open('summary.tsv'),delimiter='\t')}
E={r['name']:r for r in csv.DictReader(open('events2.tsv'),delimiter='\t')}
M={r['name']:r for r in csv.DictReader(open('metrics.tsv'),delimiter='\t')}
short={'theta0.5':'cn','SDIRK_2_2':'sd22','SDIRK_4_3':'sd43','SDIRK_Crouzeix_3_4':'crz','DIRK_TRBDF_3_2':'trbdf','galpha0.0':'ga0','galpha0.5':'ga50','EXRK_RungeKutta_4_4':'rk4','EXRK_SSP_3_3':'rk3'}
order=['cn','sd22','sd43','crz','trbdf','ga0','ga50','rk3','rk4']
code={'DONE':'ok','GUARD':'G','NEWTON':'N','QUOTA':'q'}
fmt=sys.argv[1]; W=int(sys.argv[2]) if len(sys.argv)>2 else 15; filt=sys.argv[3] if len(sys.argv)>3 else ''
cfg={}
for n,r in S.items():
    if r['tier']=='P': continue
    T=float(r['T']); spp=round(T/float(r['dt']))
    key=(r['model'],r['tier'],r['T'],r['kA'],r['ncell'],r['gamma'],spp)
    cfg.setdefault(key,{})[short[r['integ']]]=n
print('%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%('m','tr','T','kA','n','γ','spp')+' '.join(('%-'+str(W)+'s')%o for o in order))
for key in sorted(cfg,key=lambda k:(k[0],k[1],float(k[2]),float(k[3]),int(k[4]),float(k[5]),k[6])):
    line='%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%key
    if filt and filt not in line: continue
    cells=[]
    for o in order:
        n=cfg[key].get(o)
        if n is None: cells.append(('%-'+str(W)+'s')%'.'); continue
        e=E.get(n); m=M.get(n); r=S[n]; c=code[r['outcome']]
        if e is None: cells.append(('%-'+str(W)+'s')%(c+' -')); continue
        def g(k):
            v=e[k] if k in e else m[k]
            return float('inf') if v=='none' else float(v)
        try: s=c+' '+eval(fmt)
        except Exception as ex: s=c+' err'
        cells.append(('%-'+str(W)+'s')%s)
    print(line+' '.join(cells))
