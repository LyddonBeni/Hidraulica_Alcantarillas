# -*- coding: utf-8 -*-
# IMPORTANDO LIBRERIAS:
from numpy import *   # Importando libreria de matrices
from math import *    # Importando libreria de operadores matematicas
# from numpy import linalg

print ("""
             ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
                  UNIVERSIDAD  NACIONAL  DE  HUANCAVELICA
                    FACULTAD  DE  CIENCIAS  DE  INGENIERÍA
                   ESCUELA  ACADÉMICA  PROFESIONAL  DE  CIVIL
             ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
             ░░░░░░  DETERMINACION DE LA PENDIENTE OPTIMA    ░░░░
             ░░░░ PARA EL DISEÑO HIDRAULICO DE ALCANTARILLAS   ░░
             ░ DE PROYECTOS VIALES EN LA REGION DE HUANCAVELICA ░             
             
                       ===== PROYECTO DE TESIS =====
                       
 * AUTOR      : VARJE ESTEBAN, Lyddon Beni
 * ASESOR     : AYALA BIZARRO, Ivan
""")

# INGRESAR DATOS:
C=0.40   # Coeficiente de escorrimiento (adimensional)
I=34.82   # intensidad  (mm/h)
S=0.5  # Pendiente ()
n=0.014  # Coeficiente de Manning
P=4   # Perimetro mojado (m)
LT=5  # Longitud tipico de tubos para alcantarilla (m)
et=0.1  # Espesor de tuberia (m).

CSE=30.40 # Cota subrasante en el eje
Rell=1.  # Relleno (m)
   
dato=genfromtxt('I_Datos.txt') # importando datos de foramto txt

L=dato[:,1] # longitud de cada alcantarilla
AC=dato[:,2] # area de cada cuenca (Ha)
PG=dato[:,0]

Qm3=C*I*AC/360. 
Qls=Qm3*1000.
D=(Qm3/(S**0.54*C*0.2785))**0.38
Dp=sum(D)/len(D)

print ('================================')
print ('PROMEDIO DE DIAMETRO DE TUBERIA')
print ('================================')
print (Dp)

R=D/4.
A=(Qm3*n)/(R**(2/3.)*S**0.5)
AA=pi*D**2/4.

# CALCULOS
print ('==========================')
print ('CALCULO DE NUMERO DE TUBOS')
print ('==========================')
print ("Cantidad de tubos sin redondeo")
NT=L/LT # numero de tubos
print (NT)

NT=input('Ingrese la cantida de tubos por criterio del especialista [Ntub1,Ntub2,Ntub3,...Ntubn] : ') #[6,8,7,7,7,5,5,6,6]
print ("Cantidad de tubos Redondeados")
print (NT)

print ('=============================')
print ('CALCULO DE COTA DE EXCAVACION')
print ('=============================')

CE=CSE-Rell-Dp
print (CE)

print ('==============================================')
print ('CALCULO DE COTA INVERTIDA DE ENTRADA Y SALIDA')
print ('==============================================')
Di=Dp-2*et
EI=CSE-(Rell+Di+et)
SI=EI-S*LT
EI1=[EI]
SI1=[SI]

for i in range(max(NT)-1):
    EI=round(SI,4)
    EI1.append(EI)
    SI=round((EI-S*LT),4)
    SI1.append(SI)
print (EI1)
print (SI1)
print ('===================================================')
print ('CALCULO DE COTA DE LOMO DE TUBO DE ENTRADA Y SALIDA')
print ('===================================================')
EL=EI1+Dp-et
SL=SI1+Dp-et
print (EL)
print (SL)
    
#RESULTADOS
Dpp=D*100/2.54 
print ("                                              CALCULOS GENERALES")
print ("=======================================================================================================")
print ("  PROG.","│","LONG. (m)","│","A. CU. (Ha)","│","Q (m3)","│"," R (m) ","│","A. SEC. T.(m2)","│","A. ALC. (m2)","│","Diametro Pulg.","│")
print ("=======================================================================================================")
for i in range(len(L)):
    print (repr(round(PG[i],3)).rjust(8),repr(L[i]).rjust(8),repr(round(AC[i],3)).rjust(10),repr(round(Qm3[i],3)).rjust(13),repr(round(R[i],3)).rjust(8),repr(round(A[i],3)).rjust(12),repr(round(AA[i],3)).rjust(14),repr(round(Dpp[i],2)).rjust(14))

print ("=======================================================================================================")









