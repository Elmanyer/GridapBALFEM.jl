# run from $IS_OUT (reads summary.tsv, metrics.tsv, events2.tsv there)
import csv, sys, math
S={r['name']:r for r in csv.DictReader(open('summary.tsv'),delimiter='\t')}
M={r['name']:r for r in csv.DictReader(open('metrics.tsv'),delimiter='\t')}
short={'theta0.5':'cn','SDIRK_2_2':'sd22','SDIRK_4_3':'sd43','SDIRK_Crouzeix_3_4':'crz','DIRK_TRBDF_3_2':'trbdf','galpha0.0':'ga0','galpha0.5':'ga50','EXRK_RungeKutta_4_4':'rk4','EXRK_SSP_3_3':'rk3'}
order=['cn','sd22','sd43','crz','trbdf','ga0','ga50','rk3','rk4']
code={'DONE':'ok','GUARD':'G','NEWTON':'N','QUOTA':'q'}
fmt=sys.argv[1]; filt=sys.argv[2] if len(sys.argv)>2 else ''
cfg={}
for n,r in S.items():
    T=float(r['T']); spp=round(T/float(r['dt']))
    key=(r['model'],r['tier'],r['T'],r['kA'],r['ncell'],r['gamma'],spp)
    cfg.setdefault(key,{})[short[r['integ']]]=(r,M.get(n))
W=int(sys.argv[3]) if len(sys.argv)>3 else 16
print('%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%('m','tr','T','kA','n','γ','spp')+' '.join(('%-'+str(W)+'s')%o for o in order))
for key in sorted(cfg,key=lambda k:(k[0],k[1],float(k[2]),float(k[3]),int(k[4]),float(k[5]),k[6])):
    line='%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%key
    if filt and not line.startswith(filt): continue
    cells=[]
    for o in order:
        x=cfg[key].get(o)
        if x is None or x[1] is None: cells.append(('%-'+str(W)+'s')%'.'); continue
        r,m=x; c=code[r['outcome']]
        f=lambda k: float(m[k])
        try: s=c+' '+eval(fmt)
        except Exception as e: s=c+' err'
        cells.append(('%-'+str(W)+'s')%s)
    print(line+' '.join(cells))
