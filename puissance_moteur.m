clear all; clc;

%Données

A=1.35*0.8;
C_D = 0.3;
rho = 1.292;
v=[0:1:10];
m=250;
g=9.81;
C_roul = 0.0025;
alpha = 0;%10
delta_t= 5;

v_fin = (30*1000)/3600; %30km/h
v_ini = (20*1000)/3600; %20km/h
v_moy = (v_fin+v_ini)/2; %25km/h

acc = 2;
%---------------------------------------------------------------------------
%---------------------------CALCUL DE PUISSANCE-----------------------------
%---------------------------------------------------------------------------
%P_cst
F_air = 0.5*C_D*A*rho*v.^(2);
F_roul = m*g*C_roul;
F_grav = m*g*sin(deg2rad(alpha));

F_drag= F_air + F_roul + F_grav;
P_cst=F_drag .*v;

%P_acc : tu aurais pu juste utiliser masse*acceleration mdr
%E_kin = (1/2) * (m*(v.^2));
%P_acc = E_kin/delta_t;
P_acc = m*acc*v;

%total
P_watt = P_cst + P_acc;

%total in hp
P_hp = P_watt/750;

%%Plot de La courbe
figure(1)
plot(v,P_hp)
title('Velocity VS Power (de 0 à 30 ou de 10 à 30, il faut la m puissance)')
xlabel('Velocity [m/s]') 
ylabel('Power [hp]')      

%Plot les deux autres puissances
figure(2)

plot(v,P_cst/750)
title('Velocity VS Power in cst speed and acceleration')
xlabel('Velocity [m/s]') 
ylabel('Power [hp]')      
hold on;
plot(v,P_acc/750)
legend('Power at cst speed phase','Power at acceleration phase')

%---------------------------------------------------------------------------
%-------------La stratégie à appliquer pour consommer le moins--------------
%---------------------------------------------------------------------------

energy_liter = 23.6 * 10^(6);%J/l
alpha=0;

temps = ((v_fin-v_ini)/acc);

%E_kin = (1/2) * (m*((v_fin^2)-(v_ini^2)));
%F_air = @(vi) 0.5*C_D*A*rho*vi.^(2);
%F_air_int = integral( F_air , v_ini , v_fin );
%F_air_div_int = integral( F_air_div, v_ini , v_fin );

F_air_int = 0.5*C_D*A*rho*(((v_fin^3)-(v_ini^3))/3);
F_air_div = v_fin-v_ini; 
F_air_moy = F_air_int/F_air_div;

F_roul = m*g*C_roul;

F_drag = F_air_moy + F_roul;

%vitesse = @(vi) vi;
%vitesse_int = int( vi , v_ini , v_fin );
vitesse_int = ((v_fin^2)-(v_ini^2))/2;
vitesse_div = F_air_div;
vitesse_total = vitesse_int/vitesse_div;
distance = (v_ini*temps)+(0.5*acc*(temps^2));

%Consommation lors de la phase d'acceleration
%Il y a bien un moment ou on lache la pedale amenant l'acceleration de a à
%0
Energie_acceleration = (m*acc*(distance))+(F_drag*distance)
consumption = (Energie_acceleration/energy_liter) %l

%--- Partie desacelleration---

F_air_int_d = F_air_int;
F_air_div_d = abs(v_fin-v_ini); 
F_air_moy_d = F_air_int/F_air_div;

F_roul_d = m*g*C_roul;

F_drag_d = F_air_moy_d + F_roul_d; % même valeur que F_drag (logique)
acc_d = -F_drag_d/m;
temps_d = ((v_ini-v_fin)/acc_d);

vecteur_distance = zeros(1,11);
vecteur_vitesse = zeros(1,11);
vecteur_vitesse(1)= v_ini;
vecteur_vitesse_cst = zeros(1,11);
vecteur_vitesse_cst(1)= (v_fin+v_ini)/2;

distance_d = (v_fin * temps_d )+(0.5*acc_d*(temps_d^2));

F_drag_vitesse_cst = 0.5*C_D*A*rho*(6.9)+m*g*C_roul
consommation_vit_constante = F_drag_vitesse_cst*(distance + distance_d)

distance_totale = 0;
consumption_totale = 0;

for element=2:11
    if rem(element,2)==0
        distance_totale=distance_totale + distance;
        vecteur_distance(element)=distance_totale;
        vecteur_vitesse(element)= v_fin;
        vecteur_vitesse_cst(element)= (v_fin+v_ini)/2;
        consumption_totale_technique = consumption_totale+consumption;
    end
    if rem(element,2)~= 0
        distance_totale=distance_totale + distance_d;
        vecteur_distance(element)=distance_totale;
        vecteur_vitesse(element)= v_ini;
        vecteur_vitesse_cst(element)=(v_fin+v_ini)/2;
    end
end
vitesse_cst = (v_fin+v_ini)/2;
tempsss = distance_totale/vitesse_cst;
F_air_new = 0.5*C_D*A*rho*(vitesse_cst^(2));
F_roul_new = m*g*C_roul;
F_drag_new = F_air_new + F_roul_new;

Energy_cst = F_drag_new*distance_totale;
Energy_acc = Energie_acceleration * 5

ratio = ((Energy_acc/Energy_cst)*100)
TEMPS_acceleration = (temps + temps_d)*5
TEMPS_constant = tempsss

%Plot
figure(3)
plot(vecteur_distance,vecteur_vitesse)
title('Position VS Velocity')
xlabel('Position [m]') 
ylabel('Velocity [m/s]')      
hold on;
plot(vecteur_distance,vecteur_vitesse_cst)
hold off;

%---------------------------------------------------------------------------
%-------------------------CALCUL RENDEMENT CYCLE---------------------------
%---------------------------------------------------------------------------
gama = 1.4;
comp_ratio = [1:0.5:50]; %check celui du GX270
thermal_efficiency = 1-(1./(comp_ratio.^(gama-1)));

figure(4)
plot(comp_ratio,thermal_efficiency)
title('compression ratio VS Thermal efficiency')
xlabel('Compression ratio') 
ylabel('Thermal efficiency')

%---------------------------------------------------------------------------
%------------------------------TORQUE--------------------------------------
%---------------------------------------------------------------------------

vitesse_torque = [0:1:10];
mu_s = 0.7;
F_air_torque = 0.5*C_D*A*rho*vitesse_torque.^(2);
F_roul_torque = m*g*C_roul;
F_grav=0;
i=0 
if i==1
    rayon_roue=0.25*ones(1,11);
    F_drag_torque= m*acc; %(avec 1.2 fonctionne mieux);
    Torque_needed =(F_drag_torque .* rayon_roue);
    figure(5)
    plot(vitesse_torque,Torque_needed)
    title('Velocity VS Torque without considering resistive forces')
    xlabel('Velocity [m/s]') 
    ylabel('Torque [Nm]')
else
    rayon_roue=0.25;
    F_drag_torque= m*acc + F_air_torque + F_roul_torque; %(avec 1.2 fonctionne mieux);
    Torque_needed =(F_drag_torque .* rayon_roue);
    figure(5)
    plot(vitesse_torque,Torque_needed)
    title('Velocity VS Torque considering resistive forces')
    xlabel('Velocity [m/s]') 
    ylabel('Torque [Nm]')
end

%---RESULTAT---
ratio
TEMPS_acceleration
TEMPS_constant

