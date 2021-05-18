#Iván Durán Santos
#Ing. en Automatización y Autotrónica
#Materia: Modelado de sistemas
#Profesor: Oscar Ramos Arroyo
clc
clear all
close all
clf

function xp = seg_ord(t,x)
    #Inicio de las matriz de variables de estado
    c = x(1);
    cp = x(2);
    #**********************************************
    R_1 = 730000       ;
    R_2 = 1155000       ;
    C_1 = 0.0000092      ;
    C_2 = 0.000001      ;
    
    wn = sqrt(1/(C_1*C_2*R_1*R_2))      ;

    zeta = (C_2*(R_1+R_2))/((2*wn)*(C_1*C_2*R_1*R_2))       ;
    
    #Matrices
    AA = [0 1; -wn*wn -2*zeta*wn];
    BB = [0; wn*wn];
    
    r(t>=0) = 1; #Entrada (escalón unitario)
    
    #Argumento para la expresión de matrices
    xp = AA*x + BB*r; 
    
endfunction


function M = MT(T,k)
    switch T
        case "TX" #Traslación sobre eje X
            M = [1,0,k;0,1,0;0,0,1];
        case "TY" #Traslación sobre eje Y
            M = [1,0,0;0,1,k;0,0,1];
        case "EX" #Escala en eje de las X
            M = [k,0,0;0,1,0;0,0,1];
        case "EY" #Escala en eje de las Y
            M = [1,0,0;0,k,0;0,0,1];
        otherwise #Matriz de Rotación 
            M = [cos(k),-sin(k),0;sin(k),cos(k),0;0,0,1];
    endswitch
endfunction


#Masuki representación visual nomás 
A = 1.5;
B = 1;
XM = [A,-A,-A,A,A];
YM = [B,B,-B,-B,B];
UM = ones(1,length(XM));
M = [XM;YM;UM];

#Rueduki Nomás es la representación 
D = linspace(0,2*pi,15);
RA = 0.25;
XR = RA*[cos(D),0,cos(2*pi/3),0,cos(4*pi/3)];
YR = RA*[sin(D),0,sin(2*pi/3),0,sin(4*pi/3)];
UR = ones(1,length(XR));
R = [XR;YR;UR];

#Resortuki Igual, solo es representación
XS = [0,1.5,2,3,4,5,6,6.5,8]/8;
YS = [0,0,B,-B,B,-B,B,0,0]/2;
US = ones(1,length(XS));
S = [XS;YS;US];

#Acomodando los elementos antes mencionados
#T = linspace(0,4*pi,50);
#X0 = 1.5;
#X = X0*cos(T);

#T = linspace(0,12,50); #Incremento (vector) de tiempo, pero no tiempo

h = 0.5; ti = 0; tff = 70; #PARA QUE SE APRECIE
#h = 0.001; ti = 0; tff = 30; #MI MUESTRA
T = (ti:h:tff);

X0 = 1;
#X0 = [1,length(T)];

#X = X0*sin(T)

#****************************************************************************
#h = 0.01; ti = 0; tff = 30;
#T = (ti:h:tff);

opciones = odeset('RelTol',1e-3,'InitialStep',h,'MaxStep',h); #Artilujio curioso de Dios
EQUIS_CERO = [0,0]; #Condiciones iniciales 
[T,EQ] = ode45(@seg_ord,T,EQUIS_CERO,opciones); #Método de integracion utilizado por MatLab
X = EQ(:,1)

#****************************************************************************

#transformaciones
P1 = MT("TX",X(1)+A)*M;
P2 = MT("TY",-1)*MT("TX",X(1))*MT("R",-X(1)/RA)*R;
P3 = MT("TY",-1)*MT("TX",X(1)+2*A)*MT("R",-X(1)/RA)*R;
P4 = MT("TX",-A-X0)*MT("EX",X0+A+X(1))*S;

#Primer secuencia
figure(1)
LW = 5;
Q1 = plot(P1(1,:),P1(2,:),'linewidth',LW);
hold on
Q2 = plot(P2(1,:),P2(2,:),'linewidth',LW,'color',[0,0,0]);
Q3 = plot(P3(1,:),P3(2,:),'linewidth',LW,'color',[0,0,0]);
Q4 = plot(P4(1,:),P4(2,:),'linewidth',LW,'color',[0.5,0.5,0.5]);

#Ejes
axis([-A-X0,3*A+X0,-B-RA,B+2*RA])
axis equal
ejes = round(-A-X0:3*A+X0);

#Animazzione
for i=1:length(T);
    P1 = MT("TX",X(i)+A)*M;
    P2 = MT("TY",-1)*MT("TX",X(i))*MT("R",-X(i)/RA)*R;
    P3 = MT("TY",-1)*MT("TX",X(i)+2*A)*MT("R",-X(i)/RA)*R;
    P4 = MT("TX",-A-X0)*MT("EX",X0+A+X(i))*S;
    set(Q1,'Xdata',P1(1,:),'Ydata',P1(2,:))
    set(Q2,'Xdata',P2(1,:),'Ydata',P2(2,:))
    set(Q3,'Xdata',P3(1,:),'Ydata',P3(2,:))
    set(Q4,'Xdata',P4(1,:),'Ydata',P4(2,:))
    set(gca,'xtick',[ejes,X(i)])
    title(num2str(X(i)))
    pause(0.1)
end