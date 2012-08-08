import random

def de(value):
    return str(value).replace(".","_")

def create_script(name,i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx,dxc,tr,A,wv,Ly,kc,sig,sigc):
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
    if sig != 0.3:
        args += "-sig " + str(sig) + " "
    if sigc != 0.15:
        args += "-sigc " + str(sigc) + " "
    

    random.seed(args)
    random.jumpahead(sd)

    args += "-s " + str(random.randint(0,9999999999))

    filename = name +"_"+ str(i)
    s =  "#!/bin/bash\n#PBS -N %s\n" % filename
    s += "#PBS -q debug\n"
    s += "#PBS -l nodes=1:ppn=8,walltime=00:05:00\n\ncd $PBS_O_WORKDIR\nexport OMP_NUM_THREADS=8\n"
    s = s + "./main.exe " + args + " -f " + filename + "\n"
    f.write(s)
    f.close

def generate_jobs(name,conditions,Ns,Nls,charges,eps,hs,hts,r0s,Ts,tFs,Rs,Pos,Fos,sds,dxs,dxcs,trs,As,wvs,Lys,kcs,sigs,sigcs):
    i = 0
    fscript = open(name+'_job_data','w')
    fscript.writelines("N, Nl, qci, eps, h, hth, Lmax, T, kT, R, pos, for, seeds,step, stepc, trs, Amp, wv, Ly, kc, speed, speedC\n")
    fscript.writelines( str(Ns)+', '+str(Nls)+', '+str(charges)+', '+str(eps)+', '+str(hs)+', '+str(hts)+', '+str(r0s)+', '+str(Ts)+', '+str(tFs)+', '+str(Rs)+', '+str(Pos)+', '+str(Fos)+', '+str(sds)+', '+str(dxs)+', '+str(dxcs)+', '+str(trs)+', '+str(As)+', '+str(wvs)+', '+str(Lys)+', '+str(kcs)+', '+str(sigs)+', '+str(sigcs))
    fscript.writelines(conditions)

    T = Ts[0]
    Po = Pos[0]
    Fo = Fos[0]
    for N in Ns:
      for Nl in Nls:
        for charge in charges:
          for ep in eps:
            for h in hs:
              for ht in hts:
                for r0 in r0s:
                  for tF in tFs:
                    for R in Rs:
                      for sd in sds:
                        for dx in dxs:
                          for dxc in dxcs:
                            for tr in trs:
                              for A in As:
                                for wv in wvs:
                                  for Ly in Lys:
                                    for kc in kcs:
                                      for sig in sigs:
                                        for sigc in sigcs:
                                          if(eval(conditions)):
                                            create_script(name,i,N,Nl,charge,ep,h,ht,r0,T,tF,Po,Fo,R,sd,dx,dxc,tr,A,wv,Ly,kc,sig,sigc)
                                            i=i+1
    fscript.close
