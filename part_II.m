%Vellios Georgios-Serafeim AEM:9471 PartB

global m1 m2 L1 Lc1 Lc2 Iz1 Iz2 g m1pro m2pro L1pro Lc1pro Lc2pro Iz1pro Iz2pro gpro eisodosb1_1 eisodosb1_2 timeu eisodosb2_1 eisodosb2_2 timeu2;
m1=6;
m2=4;
L1=0.5;
Lc1=0.2;
Lc2=0.4;
Iz1=0.43;
Iz2=0.05;
g=9.81;

m1pro=6;
m2pro=2;
L1pro=0.5;
Lc1pro=0.35;
Lc2pro=0.1;
Iz1pro=0.05;
Iz2pro=0.02;
gpro=9.81;

% arrays initialization
eisodosb1_1=[];
eisodosb1_2=[];
timeu=[];

eisodosb2_1=[];
eisodosb2_2=[];
timeu2=[];

a;
b1;
b2;

% create diagams for 1st part
function a()

% sole the system equations for 10 s and the given initial values
[t,x] = ode23s( @rhs1, [0 10], [-87 0 167 0 ] );

% Creation of a matrix that will contain q1d, q2d, first and second
% derivatives of q1d and q2d for all values created by solving the input
% system
sze=size(t);
sz=sze(1,1);
A=zeros(sz,6);
for v=1:sz
    Qd=qdesired(t(v,1));
    A(v,:)=Qd;
end

figure(1);
hold on;
title("arm 1");
a1= plot(t, x(:,1));
M1="q1";
xlabel('time');
a2= plot(t,x(:,1)-A(:,1));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(2);
hold on;
title("arm 2");
a1= plot(t, x(:,3));
M1="q2";
xlabel('time');
a2= plot(t,x(:,3)-A(:,2));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(3);
hold on;
title("arm 1");
a1= plot(t, x(:,2));
M1="q1dot";
xlabel('time');
a2= plot(t,x(:,2)-A(:,3));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off

figure(4);
hold on;
title("arm 2");
a1= plot(t, x(:,4));
M1="q2dot";
xlabel('time');
a2= plot(t,x(:,4)-A(:,4));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off
end

% create diagrams for part II
function b1()

global eisodosb1_1 timeu eisodosb1_2;

% solve the system for 10 s
[t,x] = ode23s( @rhs, [0 10], [-87 0 167 0 ] );

% Creation of a matrix that will contain q1d, q2d, first and second
% derivatives of q1d and q2d for all values created by solving the input
% system
sze=size(t);
sz=sze(1,1);
A=zeros(sz,6);
for v=1:sz
    Qd=qdesired(t(v,1));
    A(v,:)=Qd;
end

%Diagrammata. Oi perigrafes e3igoun ti sxediazetai ka8e fora
figure(5);
hold on;
title("arm 1");
a1= plot(t, x(:,1));
M1="q1";
xlabel('time');
a2= plot(t,x(:,1)-A(:,1));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(6);
hold on;
title("arm 2");
a1= plot(t, x(:,3));
M1="q2";
xlabel('time');
a2= plot(t,x(:,3)-A(:,2));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(7);
hold on;
title("arm 1");
a1= plot(t, x(:,2));
M1="q1dot";
xlabel('time');
a2= plot(t,x(:,2)-A(:,3));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off

figure(8);
hold on;
title("arm 2");
a1= plot(t, x(:,4));
M1="q2dot";
xlabel('time');
a2= plot(t,x(:,4)-A(:,4));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off

figure(9);
hold on;
title("input for arm 1");
ylabel("u1");
xlabel('time');
plot(timeu(1,:),eisodosb1_1(1,:));
hold off

figure(10);
hold on;
title("input for arm 2");
ylabel("u2");
xlabel('time');
plot(timeu(1,:),eisodosb1_2(1,:));
hold off
end

function b2()

global eisodosb2_1 timeu2 eisodosb2_2;

eisodosb2_1mikro=[];
timeu2mikro=[];
eisodosb2_2mikro=[];



% solving the system for 10s
[t,x] = ode23s( @rhsB2, [0 10], [-87 0 167 0 ] );

% Creation of a matrix that will contain q1d, q2d, first and second
% derivatives of q1d and q2d for some values created by solving the input
% system
Size=size(timeu2);
for i=1:4000:Size(1,2)
    eisodosb2_1mikro(end+1)=eisodosb2_1(1,i);
    eisodosb2_2mikro(end+1)=eisodosb2_2(1,i);
    timeu2mikro(end+1)=timeu2(1,i);
end

sze=size(t);
sz=sze(1,1);
A=zeros(sz,6);
for v=1:sz
    Qd=qdesired(t(v,1));
    A(v,:)=Qd;
end

figure(11);
hold on;
title("arm 1");
a1= plot(t, x(:,1));
M1="q1";
xlabel('time');
a2= plot(t,x(:,1)-A(:,1));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(12);
hold on;
title("arm 2");
a1= plot(t, x(:,3));
M1="q2";
xlabel('time');
a2= plot(t,x(:,3)-A(:,2));
M2="position error";
legend([a1,a2], [M1, M2]);
hold off

figure(13);
hold on;
title("arm 1");
a1= plot(t, x(:,2));
M1="q1dot";
xlabel('time');
a2= plot(t,x(:,2)-A(:,3));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off

figure(14);
hold on;
title("arm 2");
a1= plot(t, x(:,4));
M1="q2dot";
xlabel('time');
a2= plot(t,x(:,4)-A(:,4));
M2="speed error";
legend([a1,a2], [M1, M2]);
hold off

figure(15);
hold on;
title("Plots for tracking and speed error for arm 1");
xlabel('e1');
ylabel('e1dot');
plot(x(:,1)-A(:,1),x(:,2)-A(:,3));
hold off;

figure(16);
hold on;
title("Plots for tracking and speed error for arm 2");
xlabel('e2');
ylabel('e2dot');
plot(x(:,3)-A(:,2),x(:,4)-A(:,4));
hold off;

figure(17);
hold on;
title("Combined errors");
a1= plot(t,x(:,2)-A(:,3)+10*(x(:,1)-A(:,1)) );
M1="Combined error for arm 1";
xlabel('time');
a2= plot(t,x(:,4)-A(:,4) +10*(x(:,3)-A(:,2)));
M2="Combined error for arm 2";
legend([a1,a2], [M1, M2]);
hold off

figure(18);
hold on;
title("input for arm 1");
ylabel("u1");
xlabel('time');
plot(timeu2mikro(1,:),eisodosb2_1mikro(1,:));
hold off

figure(19);
hold on;
title("input for arm 2");
ylabel("u2");
xlabel('time');
plot(timeu2mikro(1,:),eisodosb2_2mikro(1,:));
hold off
end


function dxA=rhs1(t,x)

Qd=qdesired(t);

q=[x(1); x(3)];
qdot=[x(2); x(4)];
qd=[Qd(1); Qd(2)];
qd_dot=[Qd(3); Qd(4)];
qd_dot2=[Qd(5); Qd(6)];

v=qd_dot2 -20*qdot + 20*qd_dot -100*q + 100*qd;
qdot2=v;

dx_1=qdot(1,1);
dx_2=qdot2(1,1);
dx_3=qdot(2,1);
dx_4=qdot2(2,1);

dxA=[dx_1; dx_2; dx_3; dx_4;];

end


function dxB1=rhs(t,x)
  
global m1 m2 L1 Lc1 Lc2 Iz1 Iz2 g m1pro m2pro L1pro Lc1pro Lc2pro Iz1pro Iz2pro gpro eisodosb1_1 eisodosb1_2 timeu;

H11=m2*(Lc2*Lc2+ 2*L1*Lc2*cos(x(3))+ L1*L1) +Lc1*Lc1*m1 + Iz2 + Iz1;
H12= m2*Lc2*Lc2 +L1*Lc2*m2*cos(x(3)) +Iz2;
H21=H12;
H22=Lc2*Lc2*m2 +Iz2;
H=[H11 H12; H21 H22];

C11=-m2*L1*Lc2*sin(x(3))*x(4);
C12=-m2*L1*Lc2*sin(x(3))*(x(2) + x(4));
C21=m2*L1*Lc2*sin(x(3))*x(2);
C22=0;
C=[C11 C12; C21 C22];

G1=m2*Lc2*g*cos(x(1) + x(3)) +(m2*L1 +m1*Lc1)*g*cos(x(1));
G2=m2*Lc2*g*cos(x(1) +x(3));
G=[G1; G2];

HinvC=H\C;
HinvG=H\G;


H11pro=m2pro*(Lc2pro*Lc2pro+ 2*L1pro*Lc2pro*cos(x(3))+ L1pro*L1pro) +Lc1pro*Lc1pro*m1pro + Iz2pro + Iz1pro;
H12pro= m2pro*Lc2pro*Lc2pro +L1pro*Lc2pro*m2pro*cos(x(3)) +Iz2pro;
H21pro=H12pro;
H22pro=Lc2pro*Lc2pro*m2pro +Iz2pro;
Hpro=[H11pro H12pro; H21pro H22pro];

C11pro=-m2pro*L1pro*Lc2pro*sin(x(3))*x(4);
C12pro=-m2pro*L1pro*Lc2pro*sin(x(3))*(x(2) + x(4));
C21pro=m2pro*L1pro*Lc2pro*sin(x(3))*x(2);
C22pro=0;
Cpro=[C11pro C12pro; C21pro C22pro];

G1pro=m2pro*Lc2pro*gpro*cos(x(1) + x(3)) +(m2pro*L1pro +m1pro*Lc1pro)*gpro*cos(x(1));
G2pro=m2pro*Lc2pro*gpro*cos(x(1) +x(3));
Gpro=[G1pro; G2pro];

Qd=qdesired(t);

q=[x(1); x(3)];
qdot=[x(2); x(4)];
qd=[Qd(1); Qd(2)];
qd_dot=[Qd(3); Qd(4)];
qd_dot2=[Qd(5); Qd(6)];

v=qd_dot2 -20*qdot + 20*qd_dot -100*q + 100*qd;
qdot2= -HinvC*qdot -HinvG +H\(Cpro*qdot +Gpro +Hpro*v);

eisodos=Cpro*qdot +Gpro +Hpro*v;
eisodosb1_1(end+1)=eisodos(1,1);
eisodosb1_2(end+1)=eisodos(2,1);
timeu(end+1)=t;

dx_1=qdot(1,1);
dx_2=qdot2(1,1);
dx_3=qdot(2,1);
dx_4=qdot2(2,1);

dxB1=[dx_1; dx_2; dx_3; dx_4;];
end

function dxB2=rhsB2(t,x)

global m1 m2 L1 Lc1 Lc2 Iz1 Iz2 g m1pro m2pro L1pro Lc1pro Lc2pro Iz1pro Iz2pro gpro eisodosb2_1 eisodosb2_2 timeu2;

H11=m2*(Lc2*Lc2+ 2*L1*Lc2*cos(x(3))+ L1*L1) +Lc1*Lc1*m1 + Iz2 + Iz1;
H12= m2*Lc2*Lc2 +L1*Lc2*m2*cos(x(3)) +Iz2;
H21=H12;
H22=Lc2*Lc2*m2 +Iz2;
H=[H11 H12; H21 H22];

C11=-m2*L1*Lc2*sin(x(3))*x(4);
C12=-m2*L1*Lc2*sin(x(3))*(x(2) + x(4));
C21=m2*L1*Lc2*sin(x(3))*x(2);
C22=0;
C=[C11 C12; C21 C22];

G1=m2*Lc2*g*cos(x(1) + x(3)) +(m2*L1 +m1*Lc1)*g*cos(x(1));
G2=m2*Lc2*g*cos(x(1) +x(3));
G=[G1; G2];

HinvC=H\C;
HinvG=H\G;


H11pro=m2pro*(Lc2pro*Lc2pro+ 2*L1pro*Lc2pro*cos(x(3))+ L1pro*L1pro) +Lc1pro*Lc1pro*m1pro + Iz2pro + Iz1pro;
H12pro= m2pro*Lc2pro*Lc2pro +L1pro*Lc2pro*m2pro*cos(x(3)) +Iz2pro;
H21pro=H12pro;
H22pro=Lc2pro*Lc2pro*m2pro +Iz2pro;
Hpro=[H11pro H12pro; H21pro H22pro];

C11pro=-m2pro*L1pro*Lc2pro*sin(x(3))*x(4);
C12pro=-m2pro*L1pro*Lc2pro*sin(x(3))*(x(2) + x(4));
C21pro=m2pro*L1pro*Lc2pro*sin(x(3))*x(2);
C22pro=0;
Cpro=[C11pro C12pro; C21pro C22pro];

G1pro=m2pro*Lc2pro*gpro*cos(x(1) + x(3)) +(m2pro*L1pro +m1pro*Lc1pro)*gpro*cos(x(1));
G2pro=m2pro*Lc2pro*gpro*cos(x(1) +x(3));
Gpro=[G1pro; G2pro];

Qd=qdesired(t);

q=[x(1); x(3)];
qdot=[x(2); x(4)];
qd=[Qd(1); Qd(2)];
qd_dot=[Qd(3); Qd(4)];
qd_dot2=[Qd(5); Qd(6)];

Lamda=[10 0; 0 10];

qr_dot=qd_dot -Lamda*(q-qd);
qr_dot2=qd_dot2 -Lamda*(qdot-qd_dot);

Ya11=g*L1*cos(x(1));
Ya12=g*m1*cos(x(1));
Ya13=-L1*sin(x(3))*x(4)*qr_dot(1,1) -L1*sin(x(3))*(x(2)+x(4))*qr_dot(2,1) +g*cos(x(1) +x(3));
Ya21=0;
Ya22=0;
Ya23=L1*sin(x(3))*x(2)*qr_dot(1,1) +g*cos(x(1) +x(3));
Ya=[Ya11 Ya12 Ya13; Ya21 Ya22 Ya23];

Yb11=qr_dot2(1,1) +qr_dot2(2,1);
Yb12=qr_dot2(1,1)*L1*L1;
Yb13=m1*qr_dot2(1,1);
Yb14=qr_dot2(1,1);
Yb15=Yb11;
Yb16=2*L1*cos(x(3))*qr_dot2(1,1) +L1*cos(x(3))*qr_dot2(2,1);
Yb21=Yb11;
Yb22=0;
Yb23=0;
Yb24=0;
Yb25=Yb11;
Yb26=L1*cos(x(3))*qr_dot2(1,1);
Yb=[Yb11 Yb12 Yb13 Yb14 Yb15 Yb16; Yb21 Yb22 Yb23 Yb24 Yb25 Yb26];

AbsDtheta1max=[3; 0.25; 2.05];
AbsDtheta2max=[0.993; 3; 0.113; 0.45; 0.13; 2.05];
c=[0.1; 0.1];

s=qdot-qd_dot +Lamda*(q -qd);
Signs=[prosimo(s(1,1)); prosimo(s(2,1))];

a1=(abs(Ya)*AbsDtheta1max +abs(Yb)*AbsDtheta2max +c);
a2=[a1(1,1)*Signs(1,1); a1(2,1)*Signs(2,1)];
u=Hpro*qr_dot2 +Cpro*qr_dot +Gpro-a2;

eisodosb2_1(end+1)=u(1,1);
eisodosb2_2(end+1)=u(2,1);
timeu2(end+1)=t;

qdot2= -HinvC*qdot -HinvG +H\u;

dx_1=qdot(1,1);
dx_2=qdot2(1,1);
dx_3=qdot(2,1);
dx_4=qdot2(2,1);

dxB2=[dx_1; dx_2; dx_3; dx_4;];

end

n
function qd= qdesired(t)
if t<=5
    q1d=-90 +50*(1-cos(0.63*t));
    q1d_dot=31.5*sin(0.63*t);
    q1d_dot2=19.845*cos(0.63*t);
    
    q2d=170 -60*(1-cos(0.63*t));
    q2d_dot=-37.8*sin(0.63*t);
    q2d_dot2= -23.814*cos(0.63*t);
elseif t>5
    q1d=10;
    q1d_dot=0;
    q1d_dot2=0;
    
    q2d=50;
    q2d_dot=0;
    q2d_dot2=0;
end

qd=[q1d; q2d; q1d_dot; q2d_dot; q1d_dot2; q2d_dot2];

end

function gtoux= prosimo(x)
e=0.00001;

if x>=e || x<=-e
    gtoux=x/abs(x);
elseif x<e && x>-e
    gtoux=x/e;
end

end
