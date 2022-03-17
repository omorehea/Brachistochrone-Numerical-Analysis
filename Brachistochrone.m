%{
File: Brachistochrone.m
Author: Owen Morehead
Date: Feb 21, 2022
Purpose: Recover a numerical soluiton to the ordinary differential 
equation with governing boundary conditions
%}

%----- Recovering Numerical Solution -----
%----- Plotting all curves for one endpoint set -----

%use BVP4C
%shift from (x_a,y_a) to (x_ep,y_ep) to avoid singularity on boundary

a1 = 1; b1 = 4; %boundary points y_a and x_b
ep = 0.001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;    %number of grid points for accuracy

xs = ep:1/N:b1; %form evenly spaced integration range

xmesh = linspace(ep,b1,N); 

solinit = bvpinit(xmesh,@(x)guess(x,a1,b1));  %build initial solution for bvp4c
sol = bvp4c(@(x,y)bvpfcn(x,y,a1),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

figure(1); 
subplot(1,2,1)
grid on; hold on;
plot(sol.x,sol.y(1,:), '-','LineWidth',5)
xlabel('x', 'Fontsize', 25, 'Interpreter', 'latex')
y_label = ylabel('y', 'Fontsize', 25, 'Interpreter', 'latex');
title({'Curve of Fastest Descent',sprintf('A = (0,%.0f), B = (%.0f,0)',a1,b1)}, 'Fontsize', 22, 'Interpreter', 'latex')
xlim([0 b1+1]);
ylim([-6 a1+0.1]);
set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle')

%find constant k from y(epsion)
k = ((a1 - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);


%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a1 - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;

%text(.8*b,.55*a,sprintf('A = (0,%.0f)',a), 'Fontsize', 15,'Interpreter', 'latex')
%text(.8*b,.48*a,sprintf('B = (%.0f,0)',b), 'Fontsize', 15,'Interpreter', 'latex')


% ------- Comparing with other defined curves --------
tt_s = [];
syms x


f1 = a1 - (a1/b1)*x;                       %line
f2 = a1*(1-(1-(x/b1-1)^2)^(1/2));         %ellipse
f3 = ((a1-0)/(0-b1)^2)*(x-b1)^2 + 0;       %parabola
%f4 = ((x-3)^2 - 4)/5;  %not monotonically dec parabola, only works for y_a=1, x_b=5




%{
times1(i) = curvetime(f1,a1,b1);
times2(i) = curvetime(f2,a1,b1);
times3(i) = curvetime(f3,a1,b1);
%}


tt_s(1) = othercurve(f1,a1,b1); %calculate total times for each function
tt_s(2) = othercurve(f2,a1,b1);
tt_s(3) = othercurve(f3,a1,b1);
%tt_s(4) = othercurve(f4,a1,b1);

fplot(x,f1,'LineWidth',1.8);
fplot(x,f2,'LineWidth',1.8);
fplot(x,f3,'LineWidth',1.8);
%fplot(x,f4,'LineWidth',1.8);

legend({sprintf('Brachistochrone: %.5f s',totalT),...
    sprintf('Line: %.5fs',tt_s(1)),...
    sprintf('Ellipse: %.5fs',tt_s(2)),...
    sprintf('Parabola: %.5f s',tt_s(3))},'fontsize',18,'Interpreter','latex');
    %sprintf('Non-Monotonic Parabola: %.5f s',tt_s(4))},...
    %'fontsize',18,'Interpreter','latex');

%disp(['Total travel time for line from x = 0 to x = ',num2str(b),' --> T = ',num2str(tt_s(1))])
%disp(['Total travel time for parabola from x = 0 to x = ',num2str(b),' --> T = ',num2str(tt_s(2))])
%disp(['Total travel time for circle from x = 0 to x = ',num2str(b),' --> T = ',num2str(tt_s(3))])
%disp(['Total travel time for other parabola from x = 0 to x = ',num2str(b),' --> T = ',num2str(tt_s(4))])

% ---- Repeat for switched y_a and x_b values

a = b1; b = a1; %boundary points, switch coordinate values from above code
ep = 0.001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;

xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));


sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

subplot(1,2,2)
grid on; hold on;
plot(sol.x,sol.y(1,:), '-','LineWidth',5)
xlabel('x', 'Fontsize', 25, 'Interpreter', 'latex')
y_label = ylabel('y', 'Fontsize', 25, 'Interpreter', 'latex');
%y_label.Position(1) = -80;
title({'Curve of Fastest Descent',sprintf('A = (0,%.0f), B = (%.0f,0)',a,b)}, 'Fontsize', 22, 'Interpreter', 'latex')
xlim([0 b+1]);
ylim([-6 a+0.1]);
set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle')

%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);



%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);
%fx = (1/sqrt(2*9.8))*((1 + yy(2,:).^2)./(a - yy(1,:))).^(1/2);


Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;

%text(.8*b,.55*a,sprintf('A = (0,%.0f)',a), 'Fontsize', 15,'Interpreter', 'latex')
%text(.8*b,.48*a,sprintf('B = (%.0f,0)',b), 'Fontsize', 15,'Interpreter', 'latex')


% ------- Comparing with other defined curves --------
tt_s = [];
syms x


f1 = a - (a/b)*x;                       %line
f2 = a*(1-(1-(x/b-1)^2)^(1/2));         %ellipse
f3 = ((a-0)/(0-b)^2)*(x-b)^2 + 0;       %parabola

%{
times1(i) = curvetime(f1,a,b);
times2(i) = curvetime(f2,a,b);
times3(i) = curvetime(f3,a,b);
%}


tt_s(1) = othercurve(f1,a,b);
tt_s(2) = othercurve(f2,a,b);
tt_s(3) = othercurve(f3,a,b);
%tt_s(4) = othercurve(f4,a,b);


fplot(x,f1,'LineWidth',1.8);
fplot(x,f2,'LineWidth',1.8);
fplot(x,f3,'LineWidth',1.8);
%fplot(x,f4,'LineWidth',1.8);


legend({sprintf('Brachistochrone: %.5f s',totalT),...
    sprintf('Line: %.5fs',tt_s(1)),...
    sprintf('Ellipse: %.5fs',tt_s(2)),...
    sprintf('Parabola: %.5f s',tt_s(3))},...
    'fontsize',18,'Interpreter','latex');




%---------------Functions----------------
%function to calculate total travel time for any provided curve
function curves_t = othercurve(f,a,b)
df = diff(f);
tFunc = sqrt((1+(df)^2)./(2.*9.8.*(a-f)));
curves_t = vpaintegral(tFunc,0,b);
end

% ------- functions for bvp4c ---------

%original 2nd order ODE: y'' = (1+y'^2)/(2y_a - 2y)                
function dydx = bvpfcn(x,y,a)  %2nd order ODE split into two 1st order ODEs
dydx = zeros(2,1);
dydx = [y(2)
       (1+y(2).^2)/(2*(a) - 2*y(1))];   
end
%------------------------

%new boundary condition associated with x shift by distance epsilon
function res = bcfcn(ya,yb,ep,a)  %boundary condition function
res = [ya(2) + (2/(3*ep))*(a - ya(1))
       yb(1)];
end
%------------------------

function g = guess(x,a,b)    %initial guess, linear equation
g = [a-a*x/b 
     -(a/b)];
end



