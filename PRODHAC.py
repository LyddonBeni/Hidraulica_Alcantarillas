# -*- coding: utf-8 -*-
__author__ = 'LyddonBeni'
import numpy as np
from matplotlib import pyplot as plt
import scipy as sc
import numpy as np
print ("""
             ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
                  UNIVERSIDAD  NACIONAL  DE  HUANCAVELICA
                    FACULTAD  DE  CIENCIAS  DE  INGENIERÍA
                   ESCUELA  ACADÉMICA  PROFESIONAL  DE  CIVIL
             ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
             ░░░░░░   DISEÑO DE ALCANTARILLA Y CANALES       ░░░░
             ░░░░ PARA EL DISEÑO HIDRAULICO DE ALCANTARILLAS  ░░░
             ░ DE PROYECTOS VIALES EN LA REGION DE HUANCAVELICA ░             
             
                       ===== PROYECTO DE TESIS =====
                       
 * AUTOR      : VARJE ESTEBAN, Lyddon Beni
 * ASESOR     : AYALA BIZARRO, Ivan
""")
################     DATA GENERAL    ######################
def fv(Tipo):
    if Tipo == 1:
        print (" Canal De Tipo Rectangular")
        b = float(input("\nBase del Canal(m): ")  )
        Z1 = Z2 = 0.
        error = 0.000001
        y,va,con, cont = 0.8,1.,0.,60.
        while va > error:
            C1=(Qe*n/Sc**(0.5))**(3./2.)
            C2=pow(1+Z1**2,0.5)+pow(1+Z2**2,0.5)
            Z = (Z1+Z2)/2.
            m1 = (b*y+Z*y**2)**(5./2.)/(C2*y+b)-C1
            m2 = 5./2.*(b+y*2*Z)*(b*y+Z*y**2)**(3./2.)/(C2*y+b)-C2*(b*y+Z*y**2)**(5./2.)/(C2*y+b)**2
            yi = y - m1/m2
            va = abs(y-yi)
            y = yi
            yc = y
            #print y
            con = con + 1
            if con > cont:
                break
        A = b*yc + 0.5*Z1*yc**2 + 0.5*Z2*yc**2
        T = b
        V = Qe/A
        return yc,A,V,T
    elif Tipo == 2:
        print (" Canal De Tipo Triangular")
        b = 0.
        Z1 = float(input("\nTalud del Canal Izquierda(m): ")  )
        Z2 = float(input("\nTalud del Canal Derecha(m): ")  )
        error = 0.000001
        y,va,con, cont = 0.8,1.,0.,60.
        while va > error:
            C1=(Qe*n/Sc**(0.5))**(3./2.)
            C2=pow(1+Z1**2,0.5)+pow(1+Z2**2,0.5)
            Z = (Z1+Z2)/2.
            m1 = (b*y+Z*y**2)**(5./2.)/(C2*y+b)-C1
            m2 = 5./2.*(b+y*2*Z)*(b*y+Z*y**2)**(3./2.)/(C2*y+b)-C2*(b*y+Z*y**2)**(5./2.)/(C2*y+b)**2
            yi = y - m1/m2
            va = abs(y-yi)
            y = yi
            yc = y
            #print y
            con = con + 1
            if con > cont:
                break
        A = b*yc + 0.5*Z1*yc**2 + 0.5*Z2*yc**2
        T = y*(Z1+Z2)
        V = Qe/A
        return yc,A,V,T
    elif Tipo == 3:
        print (" Canal De Tipo Trapezoidal")
        b = float(input("\nBase del Canal(m): ")  )
        Z1 = float(input("\nTalud del Canal Izquierda(m): ")  )
        Z2 = float(input("\nTalud del Canal Derecha(m): ")  )
        error = 0.000001
        y,va,con, cont = 0.8,1.,0.,60.
        while va > error:
            C1=(Qe*n/Sc**(0.5))**(3./2.)
            C2=pow(1+Z1**2,0.5)+pow(1+Z2**2,0.5)
            Z = (Z1+Z2)/2.
            m1 = (b*y+Z*y**2)**(5./2.)/(C2*y+b)-C1
            m2 = 5./2.*(b+y*2*Z)*(b*y+Z*y**2)**(3./2.)/(C2*y+b)-C2*(b*y+Z*y**2)**(5./2.)/(C2*y+b)**2
            yi = y - m1/m2
            va = abs(y-yi)
            y = yi
            yc = y
            #print y
            con = con + 1
            if con > cont:
                break
        A = b*yc + 0.5*Z1*yc**2 + 0.5*Z2*yc**2
        V = Qe/A
        T = b+y*Z1+y*Z2
        return yc,A,V,T
print ("\n1.Tipo Rectangular")
print ("2.Tipo Triangular")
print ("3.Tipo Trapezoidal")
Tipo = int(input(u"Que Tipo de Canal Va a Ingresar: ")  )
Qe = float(input(u"Ingrese Caudal de Diseño(m3/s): ")  )    
n = float(input(u"Ingrese Coeficiente de Manning: ")  ) 
Sc = float(input(u"Ingrese Pendiente del Canal: ")  ) 
y,A,V,T = fv(Tipo)
print ("\n1.Tipo Circular")
print ("2.Tipo Abovedado")
Culv = int(input(u"Que Tipo de Alcantarilla Va a Ingresar: ")  )
# Calculo del Diametro de la Alcantarilla:
g = 9.81   # Gravedad
Diam = []
val = []
for i in range(30):
    Diam = np.append(Diam,[12+3*i])
Diam = Diam/100.
if Culv == 1:
    print (" Alcantarilla Tipo Circular")
    CotaA = float(input(u"Ingrese Cota del canal antes de Transicion: ")  )
    Cober = float(input(u"Ingrese Cobertura de Carretera: ")  )
    Borde = float(input(u"Ingrese Borde de Alcantarilla: ")  )
    Talud = float(input(u"Ingrese Talud de la Carretera: ")  )
    Lon = float(input(u"Ingrese Ancho del camino: ")  )
    Ss = float(input(u"Ingrese Pendiente de Alcantarilla;\: ")  )
    n = float(input(u"Ingrese Rugosidad de Alcantarilla;\: ")  )
    Aa = Qe/2.5
    D =pow((4/np.pi)*Aa,0.5)
    for i in range(30):
        if D >= Diam[i]:
            val = np.append(val,[Diam[i]-D])
    Dc = np.max(val)
    Dc = D + Dc
    Ad = np.pi*Dc**2/4
    Vc = Qe/Ad
    hv = Vc**2/(2*g)
    NAEA = CotaA+y
    CotaB = NAEA -1.5*hv-Dc
    CotaF = CotaB+D+Cober
    CotaE = CotaA+Borde+y
    Lal = 2*Talud*(CotaF-CotaE)+Lon
    VZ = Lal*Ss
    CotaC = CotaB - VZ
    Sen = (Vc*n/1.)
   
#CAUDAL, TALUD, RUGOSIDAD, PENDIENTE
Q,Z,Z1,n,S=0.071357,0.5,0.8,0.014,0.012
#CALCULATE
def f(Q,Z,Z1,n,S,y):
    A=y**2*((Z**2+1.)**0.5+(Z1**2+1.)**0.5-(Z+Z1)*0.5)
    P=2*y*(pow(Z1**2+1.,0.5)+pow(Z**2+1.,0.5))-y*(Z1+Z)
    k=Q*n/pow(S,0.5)
    fy=pow(A,5./3.)*pow(P,-2./3.)-k
    dA=y*2*((Z**2+1.)**0.5+(Z1**2+1.)**0.5-(Z+Z1)*0.5)
    dP=2*(pow(Z1**2+1.,0.5)+pow(Z**2+1.,0.5))-(Z1+Z)
    dfy=5./3.*A**(2./3.)*P**(-2./3.)*dA - 2./3.*pow(A,5./3.)*pow(P,-5./3.)*dP
    return fy,dfy,y
y,Imax=0.5,40
Tol=1E-8       # Tolerancia Para la Iteraciones
E,cont=4,0
print ("\t-------------------------------------------------------------")
print ("\t   N°    y_i       f(y_i)      f'(y_i)     y_(i+1)     Error")
print ("\t-------------------------------------------------------------")
while (E>=Tol):
    fy,dfy,y=f(Q,Z,Z1,n,S,y)
    y1=y-fy/dfy
    cont+=1
    E=np.abs(y-y1)
    print ("\t   %.0f    %.5f    %.5f    %.5f     %.5f    %.5f"%(cont,y,fy,dfy,y1,np.abs(y-y1)))
    y=y1
    if (cont>=Imax):
        break
    
print ("\nTIRANTE (y): ", round(y,4), "m")
A=y**2*((Z**2+1.)**0.5+(Z1**2+1.)**0.5-(Z+Z1)*0.5)  
print ("Area: ",round(A,4),"m^2")  
P=2*y*(pow(Z1**2+1.,0.5)+pow(Z**2+1.,0.5))-y*(Z1+Z)
print ("Perimetro Mojado: ",round(P,4),"m")
print ("Velocidad: ",round(Q/A,4),"m/s")
b=(2*A-y**2*(Z+Z1))/(2*y)
Tt=b+y*Z+y*Z1
print ("\nBase de la Seccion: ",b)
print ("Espejo de Agua: ",Tt)
#print "Base Calculado: ",P-y*(pow(Z1**2+1.,0.5)+pow(Z**2+1.,0.5))
F=(Q/A)/pow((9.8106*A/Tt),0.5)
if F == 1:
    print ("Numero de Froude: ",round(F,3))
    print ("Esta en un Regimen Critico esta en crisis")
elif 0 < F < 1:
    print ("Numero de Froude: ",round(F,3))
    print ("Esta en un Regimen Sub Critico, se trata de un rio")
elif 1 < F:
    print ("Numero de Froude: ",round(F,3))
    print ("Esta en un Regimen Super Critico, se trata de un torrente")
print ("Borde Libre segun Boureau of Reclamation: ",0.30)