def de(value):
    return str(value).replace(".","_")

def generate_jobs(Ns,Nls,charges,eps,hs,hts,r0s,Ts,tFs,Rs,Pos):
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
                                                f = open('ompjob' + str(i) + '.pbs', 'w')
                                                filename = "N%dNl%dtF%sR%sh%sht%s" % (N,Nl,de(tF),de(R),de(h),de(ht))
                                                s = "#!/bin/bash\n#PBS -N %s\n#PBS -l nodes=2:ppn=8,walltime=0:05:00\n\ncd $PBS_O_WORKDIR\nexport OMP_NUM_THREADS=16\n" % filename
                                                s = s + "./main.exe " + str(N) + " " + str(Nl) + " " + str(charge) + " " + str(ep) + " " + str(h) + " "  + str(ht) + " " + str(r0) + " " + str(T) + " " + str(tF) + " " + str(R) + " " + filename + " " + str(Po)
                                                f.write(s)
                                                f.close
                                                i=i+1
                                                fscript.writelines("\t".join([str(R),str(tF),filename])+"\n")
    fscript.close

generate_jobs([32],[64],[2.0],[1.0],[1.0],[0.0,1.0,10.0],[0.375],[3000000],[0.01],[3.0],[1000])
