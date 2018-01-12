import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab, math
from math import sqrt
from math import cos
from math import acos
from math import pi
from math import log

def Scurve():
    """
    Inputs: none
    Return: none
    
    Takes in Sigma and Tgas data, outpus and save figure
    """
    G=6.67e-8
    Msun=1.989e33
    c=2.99792485e10
    kb=1.3807e-16
    sigma=5.67e-5
    mp=1.67e-24
    aR=4.0*sigma/c
    
    alpha=0.02
    M=6.62*Msun
    bhspin = 0.
    mdot = [0.01,3.,10.]
    rg=G*M/(c**2)
    R=10.0*rg
    mu=0.615
    Kr=0.34
    
    S01nx=256
    if R==10.0*rg:
        S01istart=123
    else:
        S01istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S01nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Tgas','r')]
    
    S01SigmaAvg = [0, 0, 0, 0, 0, 0]
    S01TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S01istart,S01istart+10):
        S01SigmaAvg[0] = S01SigmaAvg[0] + f1[i][2]/10
        S01TgasAvg[0]  = S01TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S01SigmaAvg)):
        for n in range(S01nstart,S01nstart+10):
            for i in range(S01istart,S01istart+10):
                S01SigmaAvg[a] = S01SigmaAvg[a] + f1[n*S01nx+i][2]/100
                S01TgasAvg[a]  = S01TgasAvg[a]  + f2[n*S01nx+i][2]/100
        S01nstart = S01nstart + 161
    
    S1nx=256
    if R==10.0*rg:
        S1istart=123
    else:
        S1istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S1nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Tgas','r')]
    
    S1SigmaAvg = [0, 0, 0, 0, 0, 0]
    S1TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S1istart,S1istart+10):
        S1SigmaAvg[0] = S1SigmaAvg[0] + f1[i][2]/10
        S1TgasAvg[0]  = S1TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S1SigmaAvg)):
        for n in range(S1nstart,S1nstart+10):
            for i in range(S1istart,S1istart+10):
                S1SigmaAvg[a] = S1SigmaAvg[a] + f1[n*S1nx+i][2]/100
                S1TgasAvg[a]  = S1TgasAvg[a]  + f2[n*S1nx+i][2]/100
        S1nstart = S1nstart + 161
    
    S3nx=256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') ffor line in open('../shakura_3E_PP/Tgas','r')]
    
    S3SigmaAvg = [0, 0, 0, 0, 0, 0]
    S3TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3SigmaAvg[0] = S3SigmaAvg[0] + f1[i][2]/10
        S3TgasAvg[0]  = S3TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3SigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3SigmaAvg[a] = S3SigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3TgasAvg[a]  = S3TgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S3nx = 256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/Tgas','r')]
    
    S3PSigmaAvg = [0, 0, 0, 0, 0, 0]
    S3PTgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3PSigmaAvg[0] = S3PSigmaAvg[0] + f1[i][2]/10
        S3PTgasAvg[0]  = S3PTgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3PSigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3PSigmaAvg[a] = S3PSigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3PTgasAvg[a]  = S3PTgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S10nx=256
    if R==10.0*rg:
        S10istart=98
    else:
        S10istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S10nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Tgas','r')]
    
    S10SigmaAvg = [0, 0, 0, 0, 0, 0]
    S10TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S10istart,S10istart+10):
        S10SigmaAvg[0] = S10SigmaAvg[0] + f1[i][2]/10
        S10TgasAvg[0]  = S10TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S10SigmaAvg)):
        for n in range(S10nstart,S10nstart+10):
            for i in range(S10istart,S10istart+10):
                S10SigmaAvg[a] = S10SigmaAvg[a] + f1[n*S10nx+i][2]/100
                S10TgasAvg[a]  = S10TgasAvg[a]  + f2[n*S10nx+i][2]/100
        S10nstart = S10nstart + 161
    
    H = []
    SigmaDisk=[[3.25600e+03,3.25545e+03,3.25400e+03,3.27978e+03,4.33925e+03],
               [2.89744e+04,3.46920e+04,2.98701e+04,2.96050e+04],
               [8.63527e+03,4.85676e+03,8.48039e+03,1.16699e+04,1.37747e+04,6.66683e+03]]
    Tgas=[[1.62770e+06,2.86969e+06,2.92055e+06,4.16175e+06,4.26575e+06],
          [2.79669e+07,1.54046e+07,1.41224e+07,1.55943e+07],
          [2.86917e+07,1.06522e+07,9.38209e+06,9.22887e+06,1.12289e+07,1.11367e+07]]
    
    z1  = 1.0 + pow((1. - bhspin*bhspin), (1./3.))*(pow((1. + bhspin), (1./3.))+pow((1. - bhspin), (1./3.)))
    z2  = sqrt(3.*bhspin*bhspin + z1*z1)
    rms = (3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
    y = sqrt(R/rg)
    y0 = sqrt(rms)
    y1 = 2.*cos((acos(bhspin)-pi)/3.)
    y2 = 2.*cos((acos(bhspin)+pi)/3.)
    y3 = -2.*cos((acos(bhspin))/3.)
    m1 = M/Msun
    Ac = 1.+bhspin*bhspin*pow(y,-4)+2.*bhspin*bhspin*pow(y,-6)
    Bc = 1.+bhspin*pow(y,-3)
    Cc = 1.-3./y/y+2.*bhspin*pow(y,-3)
    Dc = 1.-2./y/y+bhspin*bhspin*pow(y,-4)
    Fc = 1.-2.*bhspin/y/y/y+bhspin*bhspin/y/y/y/y
    Gc = 1.-2./y/y+bhspin/y/y/y
    Rc = Fc*Fc/Cc-bhspin*bhspin/y/y*(Gc/sqrt(Cc)-1.)
    Ec = Ac*Ac/Bc/Bc*Cc/Dc*Rc
    Q0 = Bc/sqrt(Cc)/y
    Qc = Q0*(y-y0-1.5*bhspin*log(y/y0)-3.*pow(y1-bhspin,2)/y1/(y1-y2)/(y1-y3)*log((y-y1)/(y0-y1))-3.*pow(y2-bhspin,2)/y2/(y2-y1)/(y2-y3)*log((y-y2)/(y0-y2))-3.*pow(y3-bhspin,2)/y3/(y3-y1)/(y3-y2)*log((y-y3)/(y0-y3)))
    # for i in range(len(mdot)):
    #    if mdot[i] < 0.02:
    ## Middle region
    #        SigmaDisk.append(9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6))
    #        Tgas.append(7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    #        H.append(1.e3*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.2)*pow(y,2.1)*Ac*pow(Bc,-1.2)*sqrt(Cc)*pow(Dc,-0.6)*pow(Ec,-0.5)*pow(Qc,0.2))
    ## Outer region
    #        SigmaDisk.append(4.e5*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.7)*pow(y,-1.5)*pow(Ac,0.1)*pow(Bc,-0.8)*sqrt(Cc)*pow(Dc,-0.85)*pow(Ec,-0.05)*pow(Qc,0.7))
    #        Tgas.append(2.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.3)*pow(y,-1.5)*pow(Ac,-0.1)*pow(Bc,-0.2)*pow(Dc,-0.15)*pow(Ec,0.05)*pow(Qc,0.3))
    #        H.append(4.e2*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.15)*pow(y,2.25)*pow(Ac,0.95)*pow(Bc,-1.1)*sqrt(Cc)*pow(Dc,-0.575)*pow(Ec,-0.475)*pow(Qc,0.15))
    #    else:
    #        SigmaDisk.append(5./alpha/mdot[i]*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc)
    #        Tgas.append(5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-1.5)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    #        H.append(1.e5*m1*mdot[i]*Ac*Ac/Bc/Bc/Bc*sqrt(Cc)/Dc/Ec*Qc)
    #    print H[i]/rg
    print ("0.01E", 9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(0.01,0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6), 7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(0.01,0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    print ("1E", 5./alpha/1.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("3E", 5./alpha/3.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("10E", 5./alpha/10.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    
    Tc=np.arange(3.0e5,4.19e7,1.0e2)
    omega=sqrt(G*M/(R**3))
    #T=Tc - (0.5*10**7)
    Sigma=np.sqrt(1.6*mp/kb*(32.*sigma/(27.*Kr*alpha*omega)*Tc**3-np.sqrt(32.*sigma/(27.*Kr*alpha*omega))/omega*4.*sigma/3./c*Tc**5))
    #S=Sigma +(1*10**4)
    
    ######################################
    # Novikov Thorne 1973 solutions
    ######################################
    Sigma_NT = np.arange(10.,1.e5,1.e2)
    Tin  = 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25)*Sigma_NT/Sigma_NT
    Tmid = 3.5e5*pow(alpha,1./3.)*pow(m1,-1./3.)/y*pow(Cc,-1./3.)*pow(Dc,1./3.)*pow(Sigma_NT,2./3.)
    ######################################
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel('$T_c\,\mathrm{[K]}$',**font)
    plt.xlabel('$\Sigma\,\mathrm{[g/cm^2]}$',**font)
    pylab.xlim([1200,1e5])
    pylab.ylim([1e6,1e8])
    plt.xscale('log')
    plt.yscale('log')
    plot2 = plt.plot(S01SigmaAvg[0], S01TgasAvg[0], "ko", color = 'red', mec = 'red', markersize = 6)
    plot2 = plt.plot(S1SigmaAvg[0], S1TgasAvg[0], "ks", color = 'green', mec = 'green', markersize = 6)
    plot2 = plt.plot(S3SigmaAvg[0], S3TgasAvg[0], "kv", color = 'blue', mec = 'blue', markersize = 6)
    plot2 = plt.plot(S3PSigmaAvg[0], S3PTgasAvg[0], "k<", color = 'cyan', mec = 'cyan', markersize = 6)
    plot2 = plt.plot(S10SigmaAvg[0], S10TgasAvg[0], "kD", color = 'orange', mec = 'orange', markersize = 6)
    plot2 = plt.plot(S01SigmaAvg[1], S01TgasAvg[1], "ko", color = 'red', mec = 'red', markersize = 8)
    plot2 = plt.plot(S01SigmaAvg[2], S01TgasAvg[2], "ko", color = 'red', mec = 'red', markersize = 10)
    plot2 = plt.plot(S01SigmaAvg[3], S01TgasAvg[3], "ko", color = 'red', mec = 'red', markersize = 12)
    plot2 = plt.plot(S01SigmaAvg[4], S01TgasAvg[4], "ko", color = 'red', mec = 'red', markersize = 14)
    plot2 = plt.plot(S01SigmaAvg[5], S01TgasAvg[5], "ko", color = 'red', mec = 'red', markersize = 16)
    plot2 = plt.plot(S1SigmaAvg[1], S1TgasAvg[1], "ks", color = 'green', mec = 'green', markersize = 8)
    plot2 = plt.plot(S1SigmaAvg[2], S1TgasAvg[2], "ks", color = 'green', mec = 'green', markersize = 10)
    plot2 = plt.plot(S1SigmaAvg[3], S1TgasAvg[3], "ks", color = 'green', mec = 'green', markersize = 12)
    plot2 = plt.plot(S1SigmaAvg[4], S1TgasAvg[4], "ks", color = 'green', mec = 'green', markersize = 14)
    plot2 = plt.plot(S1SigmaAvg[5], S1TgasAvg[5], "ks", color = 'green', mec = 'green', markersize = 16)
    plot2 = plt.plot(S3SigmaAvg[1], S3TgasAvg[1], "kv", color = 'blue', mec = 'blue', markersize = 8)
    plot2 = plt.plot(S3SigmaAvg[2], S3TgasAvg[2], "kv", color = 'blue', mec = 'blue', markersize = 10)
    plot2 = plt.plot(S3SigmaAvg[3], S3TgasAvg[3], "kv", color = 'blue', mec = 'blue', markersize = 12)
    plot2 = plt.plot(S3SigmaAvg[4], S3TgasAvg[4], "kv", color = 'blue', mec = 'blue', markersize = 14)
    plot2 = plt.plot(S3SigmaAvg[5], S3TgasAvg[5], "kv", color = 'blue', mec = 'blue', markersize = 16)
    plot2 = plt.plot(S3PSigmaAvg[1], S3PTgasAvg[1], "k<", color = 'cyan', mec = 'cyan', markersize = 8)
    plot2 = plt.plot(S3PSigmaAvg[2], S3PTgasAvg[2], "k<", color = 'cyan', mec = 'cyan', markersize = 10)
    plot2 = plt.plot(S3PSigmaAvg[3], S3PTgasAvg[3], "k<", color = 'cyan', mec = 'cyan', markersize = 12)
    plot2 = plt.plot(S3PSigmaAvg[4], S3PTgasAvg[4], "k<", color = 'cyan', mec = 'cyan', markersize = 14)
    plot2 = plt.plot(S3PSigmaAvg[5], S3PTgasAvg[5], "k<", color = 'cyan', mec = 'cyan', markersize = 16)
    plot2 = plt.plot(S10SigmaAvg[1], S10TgasAvg[1], "kD", color = 'orange', mec = 'orange', markersize = 8)
    plot2 = plt.plot(S10SigmaAvg[2], S10TgasAvg[2], "kD", color = 'orange', mec = 'orange', markersize = 10)
    plot2 = plt.plot(S10SigmaAvg[3], S10TgasAvg[3], "kD", color = 'orange', mec = 'orange', markersize = 12)
    plot2 = plt.plot(S10SigmaAvg[4], S10TgasAvg[4], "kD", color = 'orange', mec = 'orange', markersize = 14)
    plot2 = plt.plot(S10SigmaAvg[5], S10TgasAvg[5], "kD", color = 'orange', mec = 'orange', markersize = 16)
    plot1 = plt.plot(Sigma, Tc, 'k', lw = 2)
    plot3 = plt.plot(Sigma_NT,Tin,'m--',lw = 2)
    plot4 = plt.plot(Sigma_NT,Tmid,'m--',lw = 2)
    
    plt.text(6.e3, 5.e7, r'$Q^+ > Q^-$',**font)
    plt.text(2.e3, 1.e7, r'$Q^+ < Q^-$',**font)
    plt.text(7.e3, 2.5e6, r'$Q^+ > Q^-$',**font)
    
    #plt.ylim(0,1e5)
    plt.legend((r'S01E',r'S1E',r'S3E',r'S3Ep',r'S10E'), 'lower right', shadow=False, numpoints = 1)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    if R==10.0*rg:
        plt.savefig('../figures/S-curve_10rg.eps', format='eps', transparent='True')
    else:
        plt.savefig('../figures/S-curve_15rg.eps', format='eps', transparent='True')
    plt.show()
    return
