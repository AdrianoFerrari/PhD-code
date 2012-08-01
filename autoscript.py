import random

def de(value):
    return str(value).replace(".","_")

def create_script(i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx,dxc,tr,A,wv,Ly,kc):
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
    if Fo != 1000:
        args += "-fo " + str(Fo) + " "
    if dx != 0.03:
        args += "-D " + str(dx) + " "
    if dxc != 0.003:
        args += "-Dc " + str(dxc) + " "
    if tr != 0.0:
        args += "-tr " + str(tr) + " "
    if A != 0.0:
        args += "-A " + str(A) + " "
    if wv != 0.0:
        args += "-wv " + str(wv) + " "
    if Ly != 24.0:
        args += "-Ly " + str(Ly) + " "
    if kc != 0:
        args += "-kc " + str(kc) + " "
    

    random.seed(args)
    random.jumpahead(sd)

    args += "-s " + str(random.randint(0,9999999999))

    filename = "energyR2_" + str(i)
    s =  "#!/bin/bash\n#PBS -N %s\n" % filename
    #s += "#PBS -q debug\n"
    s += "#PBS -l nodes=1:ppn=8,walltime=00:15:00\n\ncd $PBS_O_WORKDIR\nexport OMP_NUM_THREADS=8\n"
    s = s + "./main.exe " + args + " -f " + filename
    f.write(s)
    f.close

def generate_jobs(Ns,Nls,charges,eps,hs,hts,r0s,Ts,tFs,Rs,Pos,Fos,sds,dxs,dxcs,trs,As,wvs,Lys,kcs):
    i = 0
    fscript = open('energyR2_job_data','w')
    fscript.writelines("N, Nl, qci, eps, h, hth, Lmax, T, kT, R, pos, for, seeds,step, stepc, trs, Amp, wv, Ly, kc\n")
    fscript.writelines( str(Ns)+', '+str(Nls)+', '+str(charges)+', '+str(eps)+', '+str(hs)+', '+str(hts)+', '+str(r0s)+', '+str(Ts)+', '+str(tFs)+', '+str(Rs)+', '+str(Pos)+', '+str(Fos)+', '+str(sds)+', '+str(dxs)+', '+str(dxcs)+', '+str(trs)+', '+str(As)+', '+str(wvs)+', '+str(Lys)+', '+str(kcs))
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
                                                        for dxc in dxcs:
                                                          for tr in trs:
                                                            for A in As:
                                                                for wv in wvs:
                                                                    for Ly in Lys:
                                                                        for kc in kcs:
                                                                          create_script(i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx,dxc,tr,A,wv,Ly,kc)
                                                                        i=i+1
    fscript.close

generate_jobs([32], [32], [1.0], [0.375], [0.2], [0.0], [2.0*12/32.0], [28000], [0.02], [0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0], [1000], [100], range(1,16), [0.03], [0.0], [0.0], [0.0], [12.34], [12.0], [0])
# generate_jobs([32], [32], [1.0], [0.375], [0.2], [0.0], [1.1*12.0/32.0], [21000], [0.002], [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0], [1000], [10], range(1,31), [0.03], [0.0], [0.0], [12.34], [12.0], [0])
              #N     Nl    qci     eps      h     hth     Lmax             T      kT       R                                                                                                          pos    for    seeds       step  stepC  turns Amp    wv Ly    kc
