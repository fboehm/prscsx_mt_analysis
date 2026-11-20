import csv,sys,glob
for path in sys.argv[1:]:
    f=glob.glob(path+'/*/seed_1/result.csv')
    if not f: print(path,'not finished'); continue
    d={}
    for r in csv.DictReader(open(f[0])): d[(r['phi_mode'],r['pop'],r['trait'],r['method'])]=float(r['corr_pred'])
    for pm in ('auto','fixed'):
        for pop in ('EUR','EAS'):
            v={m:sum(d[(pm,pop,t,m)] for t in sorted({k[2] for k in d}))/2 for m in ('prscsx','prscsx_mt')}
            print(path,pm,pop,'csx %.4f mt %.4f delta %+.4f'%(v['prscsx'],v['prscsx_mt'],v['prscsx_mt']-v['prscsx']))
