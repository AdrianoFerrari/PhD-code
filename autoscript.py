import random

def de(value):
    return str(value).replace(".","_")

def create_script(i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx):
    f = open('ompjob' + str(i) + '.pbs', 'w')

    args = ""
    if N != 16:
        args += "-N " + str(N) + " "
    if Nl != 32:
        args += "-Nl " + str(Nl) + " "
    if charge != 1.0:
        args += "-q " + str(charge) + " "
    if ep != 1.0:
        args += "-e " + str(ep) + " "
    if h != 1.0:
        args += "-h " + str(h) + " "
    if ht != 1.0:
        args += "-hth " + str(ht) + " "
    if r0 != 0.8:
        args += "-L " + str(r0) + " "
    if T != 200000:
        args += "-T " + str(T) + " "
    if tF != 1.0:
        args += "-kf " + str(tF) + " "

    args += "-R " + str(R) + " "

    if Po != 1000:
        args += "-po " + str(Po) + " "
    if Fo != 100:
        args += "-fo " + str(Fo) + " "
    if dx != 0.03:
        args += "-D " + str(dx) + " "
    
    random.seed(args)
    random.jumpahead(sd)
    args += "-s " + str(random.randint(0,9999999999))

    filename = de(args.replace("-","").replace(" ",""))
    s =  "#!/bin/bash\n#PBS -N %s\n" % filename
    #s += "#PBS -q debug\n"
    s += "#PBS -l nodes=2:ppn=8,walltime=0:03:30\n\ncd $PBS_O_WORKDIR\nexport OMP_NUM_THREADS=16\n"
    s = s + "./main.exe " + args + " -f " + filename
    f.write(s)
    f.close

def generate_jobs(Ns,Nls,charges,eps,hs,hts,r0s,Ts,tFs,Rs,Pos,Fos,sds,dxs):
    i = 0
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
                                            for Fo in Fos:
                                              for R in Rs:
                                                  for sd in sds:
                                                      for dx in dxs:
                                                          create_script(i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx)
                                                          i=i+1

generate_jobs([16],[32],[1.0],[1.0],[1.0],[1.0],[0.8],[200000],[0.05,0.08,0.1], [0.2, 0.6, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0], [1000], [100], [1,2,3],[0.03])
              #N    Nl   qci   eps    h    hth   Lmax     T     kT                R               pos   force seeds  steps
