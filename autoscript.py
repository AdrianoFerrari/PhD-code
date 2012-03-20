def de(value):
    return str(value).replace(".","_")

def generate_jobs(Ns,Nls,charges,eps,hs,hts,r0s,Ts,tFs,Rs,Pos,sds):
    i = 0
    fscript = open('job_data','w')
    fscript.writelines("R\tkbT\tfilename\n")
    for N in Ns:
        for Nl in Nls:
            for charge in charges:
                for ep in eps:
                    for h in hs:
                        for ht in hts:
                            for r0 in r0s:
                                for T in Ts:
                                    for tF in tFs:
                                        for Po in Pos:
                                            for R in Rs:
                                                for sd in sds:
                                                    f = open('ompjob' + str(i) + '.pbs', 'w')
                                                    filename = "%dq%se%stF%sR%sh%sht%s" % (sd,de(charge),de(ep),de(tF),de(R),de(h),de(ht))
                                                    s = "#!/bin/bash\n#PBS -q debug\nPBS -N %s\n#PBS -l nodes=2:ppn=8,walltime=0:02:30\n\ncd $PBS_O_WORKDIR\nexport OMP_NUM_THREADS=16\n" % filename
                                                    s = s + "./main.exe " + str(N) + " " + str(Nl) + " " + str(charge) + " " + str(ep) + " " + str(h) + " "  + str(ht) + " " + str(r0) + " " + str(T) + " " + str(tF) + " " + str(R) + " " + filename + " " + str(Po) + " 1600 " + str(i+1)
                                                    f.write(s)
                                                    f.close
                                                    i=i+1
                                                    fscript.writelines("\t".join([str(R),str(tF),filename])+"\n")
    fscript.close

generate_jobs([16],[32],[32.0],[0.02],[1.0],[1.0],[0.4],[3000000],[0.1,1.0], [0.2, 0.6, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0], [10000], [1])
              #N    Nl   qci   eps        h    hth   Lmax     T        kT                  R                                                                pos        seeds
