#from matplotlib import cm
import streamlit as st
from scipy.integrate import quad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
mpl.use("agg")

#from matplotlib.backends.backend_agg import RendererAgg
#_lock = RendererAgg.lock


# -- Set page config
apptitle = 'TB-ZFC/FC'

st.set_page_config(page_title=apptitle, page_icon=":magnet:")
st.title(r'$T_{B-ZFC/FC}$'+ ' calculation according to: "Determining the Zero-Field Cooling/ Field Cooling'+
' Blocking Temperature from AC-Susceptibility data for Single-Molecule Magnets"'+
' (Inorg. Chem. Front. 2025 DOI:10.1039/D4QI03259D)')
st.write("with "+r"$\tau^{-1}=exp(-U_{eff}/k_{B}T)/\tau_{0}+CT^n+1/\tau_{QT}$")

def simula_zfcfc_TRO(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_TRO(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def simula_zfcfc_RO(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_RO(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def simula_zfcfc_TO(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_TO(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def simula_zfcfc_TR(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_TR(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def simula_zfcfc_R(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_R(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def simula_zfcfc_O(tini,tmax,dt,c1p,tor0,ueff0,c0,n0,tqt0,rh,fijo0,tb,ymax):
    fig, ax1 = plt.subplots()
    xt_Xzfc_T = []
    xt_Xzfc_X = []
    xt_Xfc_X = []
    for v01 in np.arange(tini,tmax,dt):
        if v01 > tini:
            Xzfc = fijo0/v01*(1. - np.exp(-1./rh*integra_rate_func_O(tini,v01,tor0,ueff0,c0,n0,tqt0)))
            xt_Xfc_X.append(fijo0/v01)
            xt_Xzfc_T.append(v01)
            xt_Xzfc_X.append(Xzfc)
                

    ax1.plot(xt_Xzfc_T, xt_Xzfc_X,color="tab:blue",label=str(rh*60)+" (K/min)")
    ax1.plot(xt_Xzfc_T, xt_Xfc_X,color="k")
    plt.legend(fontsize=14,loc="upper right")
    plt.ylabel(r"$\chi$"+" "+r"$(cm^{3} mol^{-1})$",fontsize=18)
    plt.xlabel(r"$T$"+" "+r"$(K)$",fontsize=18)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if tb > 10.:
        plt.xlim(2,tb*3/2)
    else:
        plt.xlim(2,15)
    plt.ylim(-0.1,ymax)
    st.pyplot(fig)
    return(xt_Xzfc_T,xt_Xfc_X,xt_Xzfc_X)

def find_Tb_TRO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_TRO(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_TRO(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb_RO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_RO(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_RO(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb_TO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_TO(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_TO(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb_TR(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_TR(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_TR(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb_R(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_R(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_R(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb_O(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    for i in np.arange(tini+dt,tmax,dt):
        # l.h.s
        lhs = 1.-np.exp(-1./rh*integra_rate_func_O(tini,i,tor,ueff,c,n,tqt))
        rhs = i/(i+rh/rate_func_O(i,tor,ueff,c,n,tqt))
        if (lhs >= rhs)and(rhs > 0.5):
            if i == tini+dt:
                return 0
            else:
                return i
    return -1

def find_Tb100_TRO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_TRO(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_TRO(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def find_Tb100_RO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_RO(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_RO(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def find_Tb100_TO(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_TO(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_TO(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def find_Tb100_TR(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_TR(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_TR(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def find_Tb100_R(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_R(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_R(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def find_Tb100_O(tini,tmax,dt,tor,ueff,c,n,tqt,rh):
    if 1./rate_func_O(tini+dt,tor,ueff,c,n,tqt) < 100.:
        return 0
    else:
        for i in np.arange(tini+dt,tmax,dt):
            taui = 1./rate_func_O(i,tor,ueff,c,n,tqt)
            if taui < 100.:
                return i-dt
        return -1

def rate_func_TRO(temp,tor,ueff,c,n,tqt):
    return 1./tor*np.exp(-1*ueff/temp)+c*np.power(temp,n)+1./tqt

def integra_rate_func_TRO(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_TRO,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def rate_func_RO(temp,tor,ueff,c,n,tqt):
    return 1./tor*np.exp(-1*ueff/temp)+c*np.power(temp,n)

def integra_rate_func_RO(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_RO,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def rate_func_TO(temp,tor,ueff,c,n,tqt):
    return 1./tor*np.exp(-1*ueff/temp)+1./tqt

def integra_rate_func_TO(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_TO,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def rate_func_TR(temp,tor,ueff,c,n,tqt):
    return c*np.power(temp,n)+1./tqt

def integra_rate_func_TR(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_TR,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def rate_func_R(temp,tor,ueff,c,n,tqt):
    return c*np.power(temp,n)

def integra_rate_func_R(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_R,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def rate_func_O(temp,tor,ueff,c,n,tqt):
    return 1./tor*np.exp(-1*ueff/temp)

def integra_rate_func_O(temp0,tempf,tor,ueff,c,n,tqt):
    I = quad(rate_func_O,temp0,tempf,args=(tor,ueff,c,n,tqt))
    return(I[0])

def config_sidebar() -> None:
    with st.sidebar:
        with st.form("my_form",border=False):
            st.info("## Insert "+r"$\tau(T)$"+" parameters")
            tau0 = st.text_input(r"$\tau_{0}$"+" (s) [typical values 1.e-7 - 1.e-12]",value="")
            Ueff = st.text_input(r"$U{eff}$"+" (K) [typical values 10 - 2500]",value="")
            C = st.text_input("C "+r"$(K^{-n}/s)$"+" [typical values 1.e-8 - 1.e2]",value="")
            n = st.text_input("n"+" [typical values 2 - 9]",value="")
            tauqt = st.text_input(r"$\tau_{QT}$"+" (s) [typical values 1.e-4 - 1.e4]",value="")
            rh = st.text_input(r"$R_{H}$"+" (K/min)",value=2.)
            ln = st.text_input("Ln",value="Dy")
            ymax = st.text_input("ymax (maximum value for the y axis)",value=15.)
            submitted = st.form_submit_button("Submit",type="primary",use_container_width=True)
        return submitted,tau0,Ueff,C,n,tauqt,rh,ln,ymax

def main_page(submitted: bool, tau0: float, Ueff: float, C: float, n: float, tauqt: float, rh: float, ln: str, ymax: float) -> None:
    if submitted:
        tau0OK = False
        UeffOK = False
        COK    = False
        nOK    = False
        tauqtOK = False
        rhOK   = False
        lnOK   = False
        ymaxOK   = False
        try:
            if float(tau0) > 0.:
                tau0OK = True
                st.write(r"$\tau_{0}$"+" = "+str(tau0)+ " s")
            else:
                st.write(r"$\tau_{0}$"+" = ----")
        except:
            st.write(r"$\tau_{0}$"+" = ----")
        try:
            if float(Ueff) > 0.:
                UeffOK = True
                st.write(r"$U_{eff}$"+" = "+str(Ueff)+ " K")
            else:
                st.write(r"$U_{eff}$"+" = ----")
        except:
            st.write(r"$U_{eff}$"+" = ----")
        try:
            if float(C) > 0.:
                COK = True
                st.write("C = "+str(C)+ " "+r"$K^{-n}/s$")
            else:
                st.write("C = ----")
        except:
            st.write("C = ----")
        try:
            if float(n) >= 1.:
                nOK = True
                st.write("n = "+str(n))
            else:
                st.write("n = ----")
        except:
            st.write("n = ----")
        try:
            if float(tauqt) > 0.:
                tauqtOK = True
                st.write(r"$\tau_{QT}$"+" = "+str(tauqt)+ " s")
            else:
                st.write(r"$\tau_{QT}$"+" = ----")
        except:
            st.write(r"$\tau_{QT}$"+" = ----")
        try:
            if float(rh) > 0.:
                rhOK = True
            st.write(r"$R_{H}$"+" = "+str(rh)+" (K/min)")
        except:
            st.write("Invalid Heating Rate")
        try:
            lnlist = ["ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb"]
            magmom = [0.806693621,1.602533965,1.638544841,0.898071218,0.088226644,0,7.882830615,11.81336764,14.12887942,14.04924287,11.47551561,7.1463582,2.577228323]
            if ln.lower() in lnlist:
                lnOK = True
                indx = lnlist.index(ln.lower())
                fixXT = magmom[indx]
        except:
            pass
        try:
            if float(ymax) > 0.:
                ymaxOK = True
                ymax = float(ymax)
            st.write("ymax = "+str(ymax))
        except:
            ymax = 15.
        if lnOK == False:
            st.write("Ln not recognized")
        elif rhOK == False:
            st.write("Invalid Heating Rate (<0)")
        else:
            st.write("Ln recognized as "+ln)
            mech = "none"
            if tau0OK and UeffOK and COK and nOK and tauqtOK:
                mech = "TRO"
            elif tau0OK and UeffOK and COK and nOK:
                mech = "RO"
            elif tau0OK and UeffOK and tauqtOK:
                mech = "TO"
            elif COK and nOK and tauqtOK:
                mech = "TR"
            elif COK and nOK:
                mech = "R"
            elif tau0OK and UeffOK:
                mech = "O"
            if mech == "none":
                st.write("not enough parameters for O+R+T, O+R, O+T, R+T, O or R")
                tb = -2
            elif mech == "TRO":
                tb = find_Tb_TRO(1.,30000.,0.1,float(tau0),float(Ueff),float(C),float(n),float(tauqt),float(rh)/60.)
                tb100 = find_Tb100_TRO(1.,30000.,0.1,float(tau0),float(Ueff),float(C),float(n),float(tauqt),float(rh)/60.)
                data = simula_zfcfc_TRO(1.,300.,0.1,1.,float(tau0),float(Ueff),float(C),float(n),float(tauqt),float(rh)/60.,fixXT,tb,ymax)
            elif mech == "RO":
                tb = find_Tb_RO(1.,30000.,0.1,float(tau0),float(Ueff),float(C),float(n),tauqt,float(rh)/60.)
                tb100 = find_Tb100_RO(1.,30000.,0.1,float(tau0),float(Ueff),float(C),float(n),tauqt,float(rh)/60.)
                data = simula_zfcfc_RO(1.,300.,0.1,1.,float(tau0),float(Ueff),float(C),float(n),tauqt,float(rh)/60.,fixXT,tb,ymax)
            elif mech == "TO":
                tb = find_Tb_TO(1.,30000.,0.1,float(tau0),float(Ueff),C,n,float(tauqt),float(rh)/60.)
                tb100 = find_Tb100_TO(1.,30000.,0.1,float(tau0),float(Ueff),C,n,float(tauqt),float(rh)/60.)
                data = simula_zfcfc_TO(1.,300.,0.1,1.,float(tau0),float(Ueff),C,n,float(tauqt),float(rh)/60.,fixXT,tb,ymax)
            elif mech == "TR":
                tb = find_Tb_TR(1.,30000.,0.1,tau0,Ueff,float(C),float(n),float(tauqt),float(rh)/60.)
                tb100 = find_Tb100_TR(1.,30000.,0.1,tau0,Ueff,float(C),float(n),float(tauqt),float(rh)/60.)
                data = simula_zfcfc_TR(1.,300.,0.1,1.,tau0,Ueff,float(C),float(n),float(tauqt),float(rh)/60.,fixXT,tb,ymax)
            elif mech == "R":
                tb = find_Tb_R(1.,30000.,0.1,tau0,Ueff,float(C),float(n),tauqt,float(rh)/60.)
                tb100 = find_Tb100_R(1.,30000.,0.1,tau0,Ueff,float(C),float(n),tauqt,float(rh)/60.)
                data = simula_zfcfc_R(1.,300.,0.1,1.,tau0,Ueff,float(C),float(n),tauqt,float(rh)/60.,fixXT,tb,ymax)
            elif mech == "O":
                tb = find_Tb_O(1.,30000.,0.1,float(tau0),float(Ueff),C,n,tauqt,float(rh)/60.)
                tb100 = find_Tb100_O(1.,30000.,0.1,float(tau0),float(Ueff),C,n,tauqt,float(rh)/60.)
                data = simula_zfcfc_O(1.,300.,0.1,1.,float(tau0),float(Ueff),C,n,tauqt,float(rh)/60.,fixXT,tb,ymax)
            if tb == 0:
                st.write(r"$T_{B-ZFC/FC}$"+" < 1.1 K")
                data = []
            elif tb == -1:
                st.write(r"$T_{B-ZFC/FC}$"+" > 300 K")
                data = []
            elif tb == -2:
                st.write(r"$T_{B-ZFC/FC}$"+" is not calculated")
                data = []
            else:
                st.write(r"$T_{B-ZFC/FC}$"+" = "+"{:5.1f}".format(tb)+" K")
            if tb100 == 0:
                st.write(r"$T_{B-100s}$"+" < 1.1 K")
            elif tb100 == -1:
                st.write(r"$T_{B-100s}$"+" is not calculated")
            else:
                st.write(r"$T_{B-100s}$"+" = "+"{:5.1f}".format(tb100)+" K")
            data1 = "T X_FC X_ZFC\n"
            if len(data) != 0:
                for i in range(len(data[0])):
                    data1 += str(data[0][i])+" "+str(data[1][i])+" "+str(data[2][i])+"\n"
                st.download_button(label="Download Data",data=data1,file_name="ZFCFCdata.txt")
            data1 = ""
    else:
        pass

def main():
    submitted,tau0,Ueff,C,n,tauqt,rh,ln,ymax = config_sidebar()
    main_page(submitted,tau0,Ueff,C,n,tauqt,rh,ln,ymax)

if __name__ == "__main__":
    main()

