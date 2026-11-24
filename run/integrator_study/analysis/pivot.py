# run from $IS_OUT (reads summary.tsv, metrics.tsv, events2.tsv there)
import csv, sys
R=list(csv.DictReader(open('summary.tsv'),delimiter='\t'))
short={'theta0.5':'cn','SDIRK_2_2':'sd22','SDIRK_4_3':'sd43','SDIRK_Crouzeix_3_4':'crz','DIRK_TRBDF_3_2':'trbdf','galpha0.0':'ga0','galpha0.5':'ga50','EXRK_RungeKutta_4_4':'rk4','EXRK_SSP_3_3':'rk3'}
order=['cn','sd22','sd43','crz','trbdf','ga0','ga50','rk3','rk4']
code={'DONE':'ok','GUARD':'G','NEWTON':'N','QUOTA':'q'}
cfg={}
for r in R:
    T=float(r['T']); spp=round(T/float(r['dt']))
    key=(r['model'],r['tier'] if r['tier']!='P' else 'P',r['T'],r['kA'],r['ncell'],r['gamma'],spp)
    cfg.setdefault(key,{})[short[r['integ']]]=r
mode=sys.argv[1] if len(sys.argv)>1 else 't'
print('%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%('m','tr','T','kA','n','γ','spp')+' '.join('%-15s'%o for o in order))
for key in sorted(cfg,key=lambda k:(k[0],k[1],float(k[2]),float(k[3]),int(k[4]),float(k[5]),k[6])):
    row=cfg[key]; cells=[]
    for o in order:
        r=row.get(o)
        if r is None: cells.append('%-15s'%'.'); continue
        c=code[r['outcome']]
        if mode=='t': s='%s %.0f/%.0f'%(c,float(r['t_end']),float(r['t_target']))
        elif mode=='E': s='%s E%.3g'%(c,float(r['E_end']))
        elif mode=='u': s='%s u%.2f'%(c,float(r['u_max'])/float(r['u0']))
        elif mode=='c': s='%s %.0fs %.1fit'%(c,float(r['s_step']),float(r['nl_mean']))
        cells.append('%-15s'%s)
    print('%-3s %-2s %-4s %-5s %-3s %-4s %-3s | '%key+' '.join(cells))
