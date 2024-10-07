#Solar Panel Mount Design

import math

#known variables
rcirc = 9.78/100 #m
widthsa = 0.95 #m
lengthsa = 3.5 #m
msa = 30.5 #kg

longload = 10.08*9.81 #m/s^2
latload = 3.36*9.81 #m/s^2
longfreq = 28 #Hz
latfreq = 8.4 #Hz

#force and moment exerted on mount by the solar panel
fbendlong = msa * longload
fbendlat = msa * latload
mtors = msa * latload * lengthsa/2

#materials and properties
#density in kg/m^3, E modulus in GPa, Yield strength in MPa
M1 = [2790, 73.1, 469]
M2 = [2710, 68.9, 255]
M3 = [4430, 120, 924]
M4 = [2810, 71.7, 503]
M5 = [4540, 117, 1210]

materials = [M1, M2, M3, M4, M5]

#define cross sectional area function
def arearing(r, t):
    Aring = math.pi*(r**2-(r-t)**2)
    return(Aring)

def arearect(w, h, t):
    Arect = w*h-(w-2*t)*(h-2*t)
    return(Arect)

#define area moment of inertia calculations
def areamomentring(r, t):
    Iring = math.pi*0.25*(r**4-(r-t)**4)
    return(Iring)

def areamomentrect(w, h, t):
    Irect = (w*h**3-(w-2*t)*(h-2*t)**3)/12
    return(Irect)

#define first moment of area calculations
def qring(r, t, area):
    centroid1 = ((4*r)/(3*math.pi)*math.pi*r**2*0.5-(4*(r-t))/(3*math.pi)*math.pi*(r-t)**2*0.5)/(math.pi*r**2*0.5-math.pi*(r-t)**2*0.5)
    qri = area/2 * centroid1
    return(qri)

def qrect(w, h, t, area):
    centroid = ((h/4)*w*h*0.5-((h-t)*0.5)*(w-2*t)*(h*0.5-t))/(w*h*0.5-(w-2*t)*(h*0.5-t))
    qre = area/2 * centroid
    return(qre)

#define shear stress functions
def shearstress(F, I, Q, t):
    shear = (F*Q)/(I*t*2)
    return(shear)

#define functions for natural frequency
def fnatlateral(E, I, m, l, mb):
    fnatlat = 0.276*math.sqrt((E*I)/(m*l**3+0.236*mb*l**3))
    return(fnatlat)

def fnatlongitudinal(E, A, m, l, mb):
    fnatlong = 0.160*math.sqrt((A*E)/(m*l+0.33*mb*l))
    return(fnatlong)
    
#define beam mass function
def beammass(A, l, rho):
    mb = A*l*rho
    return(mb)

#define function for maximum bending stress
def bendingstress(F, y, I, l):
    bendstress = (F*l*y)/I
    return(bendstress)

#define function for torsional shear stress
def torsionstresscirc(M, r, t):
    torstresscirc = (M*r)/(math.pi*0.5*(r**4-(r-t)**4))
    return(torstresscirc)

def torsionstressrect(h, w, t, M):
    torstressrect = M/(2*t*(h-t)*(w-t))
    return(torstressrect)

#define function for tensile stress
def tensilestress(A, F):
    tensstress = F/A
    return(tensstress)

#make list of thicknesses to iterate through
thicknesses = []
for j in range(1, 100, 1):
    thicknesses.append(j/10000)

works = []
masses = []

length = 0.2

#loop to iterate through many different combinations
#loop through all materials
for material in materials:
    density = material[0]
    Emod = material[1]
    yieldstress = material[2]
    #loop through thicknesses
    for thickness in thicknesses:
        #print(thickness)
        combo = [material, length, thickness]
        area = arearing(rcirc, thickness)
        mass = beammass(area, length, density)
        areamoment = areamomentring(rcirc, thickness)
        firstareamoment = qring(rcirc, thickness, area)
        latfnat = fnatlateral(Emod*10**9, areamoment, length, mass, msa)
        #print(latfnat)
        longfnat = fnatlongitudinal(Emod*10**9, area, mass, length, msa)
        #print(longfnat)
        bend1 = (bendingstress(fbendlong, rcirc, areamoment, length) + tensilestress(area, fbendlat))/(10**6)
        #print(bend1)
        bend2 = (bendingstress(fbendlat, rcirc, areamoment, length) + tensilestress(area, fbendlat))/(10**6)
        #print(bend2)
        tors = torsionstresscirc(mtors, rcirc, thickness)/(10**6)
        shear1 = shearstress(fbendlat, areamoment, firstareamoment, thickness)/(10**6)
        shear2 = shearstress(fbendlong, areamoment, firstareamoment, thickness)/(10**6)
        shear1max = tors + shear1
        shear2max = tors + shear2
        combo.append(mass)
        masses.append(mass)
        if bend1 < yieldstress and bend2 < yieldstress and shear1max < yieldstress and shear2max < yieldstress and latfnat > latfreq and longfnat > longfreq:
            works.append(combo)
minimum = masses.index(min(masses))
print(works[minimum])

works2 = []
masses2 = []


for material in materials:
    density = material[0]
    Emod = material[1]
    yieldstress = material[2]
    #loop through thicknesses
    for thickness in thicknesses:
        #print(thickness)
        combo = [material, length, thickness]
        area = arearect(widthsa, rcirc*2, thickness)
        mass = beammass(area, length, density)
        areamoment1 = areamomentrect(widthsa, 2*rcirc, thickness)
        areamoment2 = areamomentrect(2*rcirc, widthsa, thickness)
        firstareamoment1 = qrect(2*rcirc, widthsa, thickness, area)
        firstareamoment2 = qrect(widthsa, 2*rcirc, thickness, area)
        latfnat = fnatlateral(Emod*10**9, areamoment2, length, mass, msa)
        #print(latfnat)
        longfnat = fnatlongitudinal(Emod*10**9, area, mass, length, msa)
        #print(longfnat)
        bend1 = (bendingstress(fbendlong, rcirc, areamoment1, length) + tensilestress(area, fbendlat))/(10**6)
        #print(bend1)
        bend2 = (bendingstress(fbendlat, widthsa/2, areamoment2, length) + tensilestress(area, fbendlat))/(10**6)
        #print(bend2)
        tors = torsionstressrect(2*rcirc, widthsa, thickness, mtors)/(10**6)
        shear1 = shearstress(fbendlat, areamoment2, firstareamoment1, thickness)/(10**6)
        shear2 = shearstress(fbendlong, areamoment1, firstareamoment2, thickness)/(10**6)
        shear1max = tors + shear1
        shear2max = tors + shear2
        combo.append(mass)
        masses2.append(mass)
        if bend1 < yieldstress and bend2 < yieldstress and shear1max < yieldstress and shear2max and latfnat > latfreq and longfnat > longfreq:
            works2.append(combo)
minimum = masses2.index(min(masses2))
print(works2[minimum])
print(beammass(arearect(widthsa, rcirc*2, 0.0002), length, 2710))
print(2*beammass(arearect(widthsa, rcirc*2, 0.0002), length, 2710))