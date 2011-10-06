def de(value):
    return str(value).replace(".","_")

def generate_jobs(Ns,charges,eps,Ts,t0s,tFs,Rs,As,wlengths,Gos,Fos,Pos):
    i = 0
    fscript = open('job_data','w')
    fscript.writelines("R\tA\tfilename\n")
    for N in Ns:
        for charge in charges:
            for ep in eps:
                for T in Ts:
                    for t0 in t0s:
                        for tF in tFs:
                            for A in As:
                                for wlength in wlengths:
                                    for Go in Gos:
                                        for Fo in Fos:
                                            for Po in Pos:
                                                for R in Rs:
                                                    f = open('mpijob' + str(i) + '.pbs', 'w')
                                                    filename = "N%dtF%sR%sA%sw%s" % (N,de(tF),de(R),de(A),de(wlength))
                                                    s = "#!/bin/bash\n#PBS -N %s\n#PBS -l nodes=4:ppn=8,walltime=0:00:40\n\ncd $PBS_O_WORKDIR\n" % filename
                                                    s = s + "mpirun -np 32 ./main.exe " + str(N) + " " + str(charge) + " " + str(ep) + " " + str(T) + " " + str(t0) + " " + str(tF) + " " + str(R) + " " + str(A) + " " + str(wlength) + " " + filename + " " + str(Go)+ " " + str(Fo)+ " " + str(Po)
                                                    f.write(s)
                                                    f.close
                                                    i=i+1
                                                    fscript.writelines("\t".join([str(R),str(A),filename])+"\n")
    fscript.close

generate_jobs([32],[2.3],[3.356],[400000],[5.0],[0.001],[0.1,0.5,1.0,2.0,7.0],[0.01,0.0],[4.0],[0],[5000],[0])
