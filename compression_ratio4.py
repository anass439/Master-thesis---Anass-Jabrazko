import math
from math import pi

#Data from GX240 and that you can change

ratio_compression = 8.5
Phi = 1
xV_ethanol = 1
xV_gazoline = 0
with_azote=1
rpm = 3600 #3600
#list_power_gx240 = [2000:3800,2100:4000,2200:4250,2300:4500,2400:4700,2500:4850,2600:5000,2700:5100,2800:5300,2900:5400,3000:5500,3100:5600,3200:5650,3300:5700,3400:5800,3500:5900]

n_eff = 0.3 #!!!! Change following simulation
Power_ini = 5900#5900
Power = Power_ini/n_eff # GX240 5.9 kW with an efficiency of 40%
#Power = 5900/0.4 # GX240 5.9 kW with an efficiency of 40%

print("---------Geometry---------")
#Given data for 242 cc
print('Old power: ' + str(Power))
thickness = 0.057

longueur_valve = 0.0328 #m
largeur_valve = 0.0014 #m
volume_valves = 2*longueur_valve*largeur_valve*thickness

Displacement = 2.42*(10**(-4))
Stroke_gx240 = 0.058
Bore = 0.073

largeur_base = (Displacement)/((ratio_compression-1)*Bore*thickness)
largeur_base2 = volume_valves/(Bore*thickness)
largeur_base = largeur_base + largeur_base2

test = ((largeur_base*Bore*thickness)+(Stroke_gx240*Bore*thickness))/(largeur_base*Bore*thickness)
test2=((largeur_base*Bore*thickness)+(Stroke_gx240*Bore*thickness))
print(test)
print(test2)
print("Largeur a apply= ", largeur_base)
Crank_radius = Stroke_gx240/2
print("The crank radius is ", Crank_radius)
print("The connecting rod is " + str((0.096666/0.029)*Crank_radius))
print("Thickness " + str(thickness))
print("position of spark " + str((largeur_base*1000)-1.4))

print("************************")

print("---------Heat transfer coefficient---------")
charac_gaz_velocity = 2*(rpm/60)*Stroke_gx240

print(charac_gaz_velocity)

h = (10.4*0.06*(((charac_gaz_velocity*Bore)/(100*10**(-6)))**(3/4)))/Bore

print("Heat transfer coefficient ", str(h))

print("************************")
print("---------Fluid---------")

Lambda = 1/Phi

MM_ethanol = 46 #C2H6O  (12*2 + 6*1 + 16)
MM_gazoline = 100 #C7H16  (12*7 + 16)

Delta_H_ethanol = 26.9*(10**6)#J/kg
Delta_H_gazoline = 42.9*(10**6) #J/kg

rho_ethanol = 789.3
rho_gazoline = 679.5

if xV_gazoline==1:
    m_dot_ethanol = 0
    m_dot_gazoline = Power/Delta_H_gazoline
    m_dot_fuel = m_dot_gazoline
    print("Le debit massique du fuel est : " + str(m_dot_fuel))

else:
    ratio = (rho_gazoline/rho_ethanol)*(xV_gazoline/xV_ethanol) #adimensionnel
    m_dot_ethanol = Power/(ratio*Delta_H_gazoline + Delta_H_ethanol) # kg/s
    m_dot_gazoline  = ratio*m_dot_ethanol # kg/s
    m_dot_fuel = m_dot_ethanol + m_dot_gazoline #kg/s
    print("Le debit massique du fuel est : " + str(m_dot_fuel))

n_dot_ethanol = m_dot_ethanol/(MM_ethanol*0.001)   #kg/s  /  kg/mol --> mol/s
n_dot_gazoline = m_dot_gazoline/(MM_gazoline*0.001)  #
n_dot_fuel = n_dot_ethanol + n_dot_gazoline

xN_ethanol_fuel = n_dot_ethanol/n_dot_fuel
xN_gazoline_fuel = n_dot_gazoline/n_dot_fuel


y = xN_ethanol_fuel

MM_air = 28.84 * (10**(-3))
MM_fuel = (((1-y)*(16+(12*7)))+(y*((12*2)+6+16)))*(10**(-3))
MM_Oxygen = 16*2
#A_Fst = ((11*Lambda)-(8*Lambda*y))*4.76*(MM_air/MM_fuel)
if with_azote == 1:
    A_Fst = (11-(8*y))*(1+(79/21))*(MM_air/MM_fuel)
    print("A_Fst: " + str(A_Fst))
    A_F = A_Fst*Lambda
    print("A_F: " + str(A_F))
    m_dot_air = A_F * m_dot_fuel
    m_dot_tot = m_dot_fuel + m_dot_air # kg/s
    print("Le debit massique total: " + str(m_dot_tot))
    print("la masse total: " )

    n_dot_air = m_dot_air/MM_air
    n_dot_tot =  n_dot_fuel + n_dot_air

# Result
#m_dot_tot
    xN_ethanol = n_dot_ethanol/n_dot_tot
    print("Fraction molaire ethanol: " + str(xN_ethanol))
    xN_gazoline = n_dot_gazoline/n_dot_tot
    print("Fraction molaire gazoline: " + str(xN_gazoline))
    xN_air = n_dot_air/n_dot_tot
    xN_N2 = 0.79 * xN_air
    print("Fraction molaire Azote: " + str(xN_N2))
    xN_O2 = 0.21 * xN_air
    print("Fraction molaire Oxygene: " + str(xN_O2))

    print("Must be one : " + str(xN_ethanol+xN_gazoline+xN_N2+xN_O2))
else:
    m_dot_oxygen = ((11*Lambda)-(8*Lambda*y))*MM_Oxygen*0.001
    m_dot_tot = m_dot_fuel + m_dot_oxygen# kg/s
    print("Le debit massique total: " + str(m_dot_tot))
    n_dot_oxygen = m_dot_oxygen/MM_Oxygen
    n_dot_tot =  n_dot_fuel + n_dot_oxygen

    xN_ethanol = n_dot_ethanol/n_dot_tot
    print("Fraction molaire ethanol: " + str(xN_ethanol))
    xN_gazoline = n_dot_gazoline/n_dot_tot
    print("Fraction molaire gazoline: " + str(xN_gazoline))
    xN_oxygen = n_dot_oxygen/n_dot_tot
    print("Fraction molaire oxygen: " + str(xN_oxygen))
    print("Must be one : " + str(xN_ethanol+xN_gazoline+xN_oxygen))


print("---------Mixture Materials---------")

print("--REACTANT--")
print("Ethanol coeff:", xN_ethanol_fuel, "and rate:",0.15)
print("Gazoline coeff:", xN_gazoline_fuel, "and rate:",0.25)
print("Oxygen coeff:", (11)-(8*y), "and rate:",1.5)
print("Azote coeff:", ((11)-(8*y))*(79/21), "and rate:", 1.5,"?")
print("")
print("--PRODUCT--")
print("CO2 coeff:", 7-(5*y), "and rate:",0)
print("H2O coeff:", 8-(5*y), "and rate:",0)
print("Azote coeff:", ((11)-(8*y))*(79/21), "and rate:",0)





