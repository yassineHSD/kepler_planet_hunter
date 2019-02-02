from __future__ import division
import sys
import math
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import numpy as np
whenMin=range(0,3)
whenMax=range(0,3)
def x_val(list):
   values=range(0,len(list))
   for i in range(0,len(list)):
       values[i]=i
   return values
def file_len(fname):
    with open(fname) as f:
        global len
        for i, l in enumerate(f):
            pass
        return i + 1
def build_transit(inputfile):
        len=file_len(inputfile)
        global Lum
        Lum=range(len)
        with open(inputfile) as f:
            lines = f.readlines()
            for i in range(0,len):
                Lum[i]=float(lines[i])
        return Lum
def maxLum(begin,end,phase):
    max=Lum[begin]
    for i in range(begin,end):
     if((max-Lum[i+1])<0):
         max=Lum[i+1]
         whenMax[phase]=i+1
    return max
def minLum(begin,end,phase):
    min=Lum[begin]
    for i in range(begin,end):
     if((min-Lum[i+1])>0):
         min=Lum[i+1]
         whenMin[phase]=i
    return min
def starmass(Lummax):
    mass=0.0
    root=2 / 7
    mass=math.pow(Lummax,root)
    return mass
def depth(Lummin,Lummax):
    dp=(Lummax-Lummin)/Lummax
    return dp
def bestLummax(begin,end):
    whenbestLummax=0
    compLum=[]
    compH=[]
    Lumprime=[]
    Lumprime=Lum[:]
    diff=[]
    m=0
    for x in Lumprime+[0]:
            if(x==Lummax):
             compH+=[Lumprime.index(Lummax)+m]
             Lumprime.remove(Lummax)
             m+=1
    for i in range(begin,int(len(compH) / 2)):
        diff+=[Lum[compH[i]+1]-Lum[compH[i]]]
    best=diff[0]
    for i in range(0,len(diff)-1):
     if((best-diff[i+1])>0):
         best=diff[i+1]
     gotit=compH[diff.index(best)]
    return gotit
def wavelength():
 try:
    lamda1=float(raw_input('Please type in the peak wavelength (nm):'))
 except ValueError:
    print "Not a number"
    wavelength()
 return lamda1
def albedo():
 try:
    albedo1=float(raw_input('Please type in how much percentage light the Alien planet does reflect:'))
 except ValueError:
    print "Not a number"
    albedo()
 return albedo1
def parallax():
  try:
     parallax1=float(raw_input('Please type in the parallax angle (arcsecs):'))
  except ValueError:
     print "Not a number"
     parallax()
  return parallax1
Lum=build_transit(str(sys.argv[1]))
print """


  _  __          _
 | |/ /         | |
 | ' / ___ _ __ | | ___ _ __
 |  < / _ \ '_ \| |/ _ \ '__|
 | . \  __/ |_) | |  __/ |
 |_|\_\___| .__/|_|\___|_|
          | |
          |_|

                                Copyright: Yassine Ben Alaya
"""
x_list=x_val(Lum)
x = np.array(x_list)
y = np.array(Lum)
plt.plot(x,y,'r')
Lummax=maxLum(0,int((len(Lum)-1)/2),2)
Lummin=minLum(0,len(Lum)-1,2)
massstar=starmass(Lummax)
localLummin1=minLum(0,int((len(Lum)-1)/2),0)
localLummin2=minLum(int((len(Lum)-1)/2),len(Lum)-1,1)
orbitperiod=whenMin[1] - whenMin[0]
root2=1 / 3
k=((orbitperiod/8766)**2 * massstar * 6.67408e-11 )/(4 * math.pi ** 2 )
a=math.pow(k,root2)
massplanet=(((2 * math.pi )** 2 * (a ** 3) * 1e11) / (((orbitperiod /8766) ** 2) * 6.67408) ) - massstar
s=bestLummax(0,whenMin[0]-1)
tF=sqrt((s - whenMin[0])**2)
tT=2 * tF
thedepth=depth(Lummin,Lummax)
g=(1-sqrt(thedepth))**2
h=(tF/tT)**2 * (1 + sqrt(thedepth)**2)
b=sqrt( (g-h) /(1 - (tF/tT)**2))
e=sqrt(abs(1-((b / a)**2)))
Mediumdistance=( a * (abs(1-(e**2))) / (1 + e) + a * (abs(1-(e**2))))/2
l= np.array([0,s,(s+whenMin[0])/2,(s+whenMin[0])/2 + tF,s+tT,2*s+tT])
f=np.array([Lummax,Lummax,Lummin,Lummin,Lummax,Lummax])
plt.plot(l,f,'g')
parallax=parallax()
print "The Alien star is ", 1/parallax * 3.26 , " Light-year away."
print 'Max Lum: ' , Lummax
print 'Min Lum: ' , Lummin
print 'Star Mass : ' , massstar , '(Solar mass) ' , massstar * 1.9891e30 ,'(Kg)'
print 'Depth: ', depth(Lummin,Lummax)
print 'The Orbit period (The alien year):' , orbitperiod , '(hours) ' , orbitperiod / 24 ,'(days)'
print 'Major semi-axis: ', a ,'(AU) ' , a * 1.496e+8 ,'(Km)'
print 'Senior semi-axis:',b,'(AU) '  , b * 1.496e+8 ,'(Km)'
print 'Planet Mass: ' , massplanet , '(Solar mass) ' , massplanet * 1.9891e30 , '(Kg)'
print 'Medium Distance between the Alien planet and its Star:',Mediumdistance,'(AU)'
xlabel('Time(hours)')
ylabel('Luminosity(Solar Luminosity)')
title('The Luminosity transit curve:')
axis([0,len(Lum) - 1, Lummin - (Lummax-Lummin)/2, Lummax + (Lummax-Lummin)/2])
savefig('fig.pdf')
lamda=600
lamda=wavelength()
Temperature= (2.9 / lamda) * 1e6
print "Stellar Temperature (Kelvin):" , Temperature
Rs=sqrt(Lummax * 3.828e26/(4 * math.pi * Temperature**4 * 5.670367e-8))/1000
print "Stellar radius", Rs ,"(Km)" , Rs/695.508, "(Solar radius)"
Rp=Rs * sqrt(depth(Lummin,Lummax))
print "Planet Radius ", Rp , " (km)", Rp/695.508 , " (Solar radius) "
pg=(massplanet* 1.9891e30)/(Rp**2) * 6.67408e-17
print "The Alien Planet's gravity: " ,pg, "(m/s^2)"
Density=massplanet * 1.9891e30 / (4/3 * math.pi * (Mediumdistance*1.496e11)**3)
print "The Alien Planet's Density ", Density , "(Kg/m^3)"
albedo=albedo()
aux=Lummax*(1-albedo/100)/(16*math.pi*5.670367e-8*(Mediumdistance*1.496e11)**2)
pt=math.pow(aux,0.25)
ptc= pt - 237.15
print "The Planets Temperature considering the Green House effect " , ptc , " (Celsius)"
if (5<=ptc<=40):
    print "The planet is Habitable !!! liquid water may exist on its surface <3"
else:
    print "The planet is not in the habitable zone :'("
plt.show()
