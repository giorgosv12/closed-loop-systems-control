%Vellios Georgios Serafeim AEM:9471 PartA
  
% functions take the following inputs:
% initial value for x1
% initial value for x2
% derivative of r
% diagramm ascendind number

%control using step function
 xtime(-2 ,1.5, 0, 1);
 xtime(-2.5, 0.8, 0, 2);
 xtime(1.5, 2, 0, 3);
 xtime(0.2, 1.8, 0, 4);
 xtime(2.5 ,-0.8, 0, 5);
 xtime(2, -2, 0, 6);
 xtime(-0.2, -1.8, 0, 7);
 xtime(-1, -2.5, 0, 8);
 
 % control using 1.2 ramp function
 xtime(-2 ,1.5, 1.2, 15);
 xtime(-2.5, 0.8, 1.2, 16);
 xtime(1.5, 2, 1.2, 17);
 xtime(0.2, 1.8, 1.2, 18);
 xtime(2.5 ,-0.8, 1.2, 19);
 xtime(2, -2, 1.2, 20);
 xtime(-0.2, -1.8, 1.2, 21);
 xtime(-1, -2.5, 1.2, 22);
 
% phase portraits for ramp 1.2
hold on
phase_portrait(-2, 1.5, 1.2, 9);
phase_portrait(-2.5, 0.8, 1.2, 9);
phase_portrait(1.5, 2, 1.2, 9);
phase_portrait(0.2, 1.8, 1.2, 9);
phase_portrait(2.5, -0.8, 1.2, 9);
phase_portrait(2, -2, 1.2, 9);
phase_portrait(-0.2, -1.8, 1.2, 9);
phase_portrait(-1, -2.5, 1.2, 9);
hold off
 
% phase portraits for step function
hold on
phase_portrait(-2, 1.5, 0, 10);
phase_portrait(-2.5, 0.8, 0, 10);
phase_portrait(1.5, 2, 0, 10);
phase_portrait(0.2, 1.8, 0, 10);
phase_portrait(2.5, -0.8, 0, 10);
phase_portrait(2, -2, 0, 10);
phase_portrait(-0.2, -1.8, 0, 10);
phase_portrait(-1, -2.5, 0, 10);
hold off
 
% phase portraits for non linean with step input
hold on
phase_portrait_nonlinear(-2, 1.5, 0, 11);
phase_portrait_nonlinear(-2.5, 0.8, 0, 11);
phase_portrait_nonlinear(1.5, 2, 0, 11);
phase_portrait_nonlinear(0.2, 1.8, 0, 11);
phase_portrait_nonlinear(2.5, -0.8, 0, 11);
phase_portrait_nonlinear(2, -2, 0, 11);
phase_portrait_nonlinear(-0.2, -1.8, 0, 11);
phase_portrait_nonlinear(-1, -2.5, 0, 11);
hold off
 
% phase portraits for non linear with ramp 1.2 input
hold on
phase_portrait_nonlinear(-2, 1.5, 1.2, 12);
phase_portrait_nonlinear(-2.5, 0.8, 1.2, 12);
phase_portrait_nonlinear(1.5, 2, 1.2, 12);
phase_portrait_nonlinear(0.2, 1.8, 1.2, 12);
phase_portrait_nonlinear(2.5, -0.8, 1.2, 12);
phase_portrait_nonlinear(2, -2, 1.2, 12);
phase_portrait_nonlinear(-0.2, -1.8, 1.2, 12);
phase_portrait_nonlinear(-1, -2.5, 1.2, 12);
hold off
 
% phase portraits for non linear with ramp 0.04 input
hold on
phase_portrait_nonlinear(-2, 1.5, 0.04, 13);
phase_portrait_nonlinear(-2.5, 0.8, 0.04, 13);
phase_portrait_nonlinear(1.5, 2, 0.04, 13);
phase_portrait_nonlinear(0.2, 1.8, 0.04, 13);
phase_portrait_nonlinear(2.5, -0.8, 0.04, 13);
phase_portrait_nonlinear(2, -2, 0.04, 13);
phase_portrait_nonlinear(-0.2, -1.8, 0.04, 13);
phase_portrait_nonlinear(-1, -2.5, 0.04, 13);
hold off
 
% phase portraits for non linear with ramp 0.4 input
hold on
phase_portrait_nonlinear(-2, 1.5, 0.4, 14);
phase_portrait_nonlinear(-2.5, 0.8, 0.4, 14);
phase_portrait_nonlinear(1.5, 2, 0.4, 14);
phase_portrait_nonlinear(0.2, 1.8, 0.4, 14);
phase_portrait_nonlinear(2.5, -0.8, 0.4, 14);
phase_portrait_nonlinear(2, -2, 0.4, 14);
phase_portrait_nonlinear(-0.2, -1.8, 0.4, 14);
phase_portrait_nonlinear(-1, -2.5, 0.4, 14);
hold off
  
% creation of extra diagrams present in the Part A report
xtime_nonlinear(1.5, 2, 0, 23);
xtime_nonlinear(0.3, 2, 1.2, 24);
xtime_nonlinear(-2, 1.5, 0, 25);
xtime(1.5, 2, 0.04, 26);
xtime_nonlinear(1.5, 2, 0.04, 27);
xtime_nonlinear(0.3, 2, 0.04, 28);
xtime_nonlinear(1.5, 2, 0.4, 29);
 
 
 
% function that plots the respnse of the functions for 8 seconds
function xtime(initial_x1, initial_x2, r_dot, fig)
 
t=0:0.001:8;
% soliving the system
[t,x] = ode45( @rhs, t, [initial_x1 initial_x2] );

b=num2str(initial_x2);
c=num2str(r_dot);
figure(fig);
hold on
title("arxikes times: (" + a + ", " + b + ") with rdot= " + c);
a1= plot(t, x(:,1));
M1="x1";
xlabel('t');
a2= plot(t, x(:,2));
M2="x2";
legend([a1,a2], [M1, M2]);
hold off
 
    function dxdt=rhs(t,x)
    dxdt_1=x(2);
    dxdt_2= -4*x(1) -x(2) + r_dot;
        
    dxdt=[dxdt_1; dxdt_2];
    end
 
 
end
 
% function for creating the phase portrait
function phase_portrait(initial_x1, initial_x2, r_dot, fig)
 
    function dxdt=rhs(t,x)
    dxdt_1=x(2);
    dxdt_2= -4*x(1) -x(2) + r_dot;
        
    dxdt=[dxdt_1; dxdt_2];
    end
 
[t,x]=ode45(@rhs,[0 8], [initial_x1 initial_x2]);
 
 
x1dot(:,1)=x(:,2);
x2dot(:,1)=-4*x(:,1)-x(:,2) +r_dot;
hold on
figure (fig);
xlabel('x1');
ylabel('x2');
 
quiver(x(:,1),x(:,2),x1dot(:,1),x2dot(:,1));
hold off
end
 
% function for creating the phase portrait for non linear input
function phase_portrait_nonlinear(initial_x1, initial_x2, r_dot, fig)
 
    function dxdt=rhs2(t,x)
    dxdt_1=x(2);
    
    if (x(1)<0.2) && (x(1)>-0.2)
        dxdt_2= -4*0.06*x(1) -x(2) + r_dot;
    else
        dxdt_2= -4*x(1) -x(2) + r_dot;
    end
       
    dxdt=[dxdt_1; dxdt_2];
    end
 
[t,x]=ode45(@rhs2,[0 8], [initial_x1 initial_x2]);
 
 
x1dot(:,1)=x(:,2);
 
if (x(:,1)<0.2) & (x(:,1)>-0.2)
    x2dot(:,1)=-4*0.06*x(:,1)-x(:,2) +r_dot;
else
    x2dot(:,1)=-4*x(:,1)-x(:,2) +r_dot;
end
 
hold on
figure (fig);
xlabel('x1');
ylabel('x2');
 
quiver(x(:,1),x(:,2),x1dot(:,1),x2dot(:,1));
hold off
end
 
% plotting the response time for non linear function
function xtime_nonlinear(initial_x1, initial_x2, r_dot, fig)
 
t=0:0.001:8;
 
[t,x] = ode45( @rhs2, t, [initial_x1 initial_x2] );
a=num2str(initial_x1);
b=num2str(initial_x2);
c=num2str(r_dot);
figure(fig);
hold on
title("(" + a + ", " + b + ") with rdot= " + c);
a1= plot(t, x(:,1));
M1="x1";
xlabel('t');
 
a2= plot(t, x(:,2));
M2="x2";
 
legend([a1,a2], [M1, M2]);
hold off
 
    function dxdt=rhs2(t,x)
    dxdt_1=x(2);
    
    if (x(1)<0.2) && (x(1)>-0.2)
        dxdt_2= -4*0.06*x(1) -x(2) + r_dot;
    else
        dxdt_2= -4*x(1) -x(2) + r_dot;
    end  
    dxdt=[dxdt_1; dxdt_2];
    end
 
end
