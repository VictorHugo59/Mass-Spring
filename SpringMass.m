clear all;
clc;
l0 = 2; %spring length
k=5; %spring constant
[xv,yv]=meshgrid(-4:0.5:4); %create a grid to assign a force value to each point thus converting it to a vector field
X= -k.*(sqrt(xv.^2+yv.^2)-l0).*(xv./(sqrt(xv.^2+yv.^2)));% 
Y= -k.*(sqrt(xv.^2+yv.^2)-l0).*(yv./(sqrt(xv.^2+yv.^2)));% substract l0 (spring length) to the position vector and multiply by its norm (direction)
quiver(xv,yv,X,Y) %vector plot
hold on
m=1;% object mass
tf=12.9; %final time simulation
n=100000; %iterations
h=tf/n; %step
t(1)=0;%initial values time
vx(1)=0;%initial values, x direction velocity
vy(1)=.7;%initial values, y direction velocity
x(1)=2.7;%initial position in x
y(1)=0;%initial position in y
Emtot(1)=0; %we create a initial total mechanic energy in order to calculate it in the loop
for i=1:n
    %solve spring-mass differential equation with euler method
    t(i+1)=t(i)+h;
    vx(i+1)=vx(i)+(-k.*(sqrt(x(i).^2+y(i).^2)-l0).*(x(i)./sqrt(x(i).^2+y(i).^2))/m)*h;
    vy(i+1)=vy(i)+(-k.*(sqrt(x(i).^2+y(i).^2)-l0).*(y(i)./sqrt(x(i).^2+y(i).^2))/m)*h;
    x(i+1)=x(i)+vx(i)*h;
    y(i+1)=y(i)+vy(i)*h;

    %Fuerza disipativa:  /UNCOMMENT this if you want to add a disipative
    %force, it is proportional to its velocoty (friction).
%     f=0.0000491; /dissipative force constant
%     vx(i+1)=vx(i+1)-f.*vx(i+1);
%     vy(i+1)=vy(i+1)-f.*vy(i+1);
    %demostrar Conservación de la energía mecánica./ HERE I demonstrate the
    %conservation of energy law by calculating the velocity's magnitud in
    %each timestep in order to obtain its kinetic energy
    %calculamos magnitud velocidad en cada tiempo del resorte para calcular
    %su energia cinetica
    Vtot=sqrt(vx(i).^2+vy(i).^2);
    %calculamos su magnitud de posición para calcular su energía potencial
    %/calculate position mangitud to obtain potential energy
    magnitudr=sqrt(x(i).^2+y(i).^2);
    %/kinetic energy
    Ec=0.5.*m.*(Vtot.^2);
    %/potential energy
    Ep=0.5*k*(magnitudr-l0).^2;
    %calculamos la energía mecanica en cada tiempo del resorte, si se
    %conserva, siempre debería dar un numero igual
    %/we calculate total energy each timestep, when there is no dissipative
    %force, conservation of energy is achieved, otherwise, the energy will
    %eventually end up being zero.
    Emtot(i)=Ec+Ep;

end 
%/Plot trajectory
plot(x,y)
title("Trayectoria")
hold off
%hacemos t un vector del mismo tamaño que Emtot para graficar
t(n+1) = [];
%/UNCOMMENT this if you want to plot a Energy vs Time graph
%plot(Emtot,t)
% title("Energia Vs Tiempo")
% xlabel("t")
% ylabel("Energía")
%Es posible observar que la energía mecánica, ex aproximadamente igual al
%inicio, en un punto medio del tiempo, y al final. La pequeña variación se
%puede deber aal error del metodo de Euler.
Emtot(1)
Emtot(50000)
Emtot(n)
%/Finally we print the total energy at the start, middle and final time for
%observation purposes.
