def de(value):
    return str(value).replace(".","_")

def generate_jobs(Ns,charges,eps,Ts,t0s,tFs,Rs,As,wlengths,Gos,Fos,Pos):
    i = 0
    for N in Ns:
        for charge in charges:
            for ep in eps:
                for T in Ts:
                    for t0 in t0s:
                        for tF in tFs:
                            for R in Rs:
                                for A in As:
                                    for wlength in wlengths:
                                        for Go in Gos:
                                            for Fo in Fos:
                                                for Po in Pos:
                                                    f = open('mpijob' + str(i) + '.pbs', 'w')
                                            filename = "N%d%se%sti%stF%sR%sA%sw%s" % (N,de(charge),de(ep),de(t0),de(tF),de(R),de(A),de(wlength))
                                            s = "#!/bin/bash\n#PBS -N %s\n#PBS -l nodes=4:ppn=8,walltime=0:20:00\n\ncd $PBS_O_WORKDIR\n" % filename
                                            s = s + "mpirun -np 32 ./main.exe" + str(N) + " " + str(charge) + " " + str(ep) + " " + str(T) + " " + str(t0) + " " + str(tF) + " " + str(R) + " " + str(A) + " " + str(wlength) + " " + filename + " " + str(Go)+ " " + str(Fo)+ " " + str(Po)
                                            f.write(s)
                                            f.close
                                            i=i+1
generate_jobs([64],[2.3],[3.356],[100000],[5.0],[0.001],[0.1,0.5,1.0,2.0,7.0],[0.01,0.0],[4.0],[10000],[10000],[0])
