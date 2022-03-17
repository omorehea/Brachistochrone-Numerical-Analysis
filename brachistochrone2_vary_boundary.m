%{
File: Brachistochrone2.m
Author: Owen Morehead
Date: Feb 21, 2022
Purpose: Plotting Brachistochrone curves for varying 
         starting (0,a) and ending (b,0) points
%}

%use BVP4C
%shift from (x_a,y_a) to (x_ep,y_ep) to avoid singularity on boundary

%-- switch end values to plot more curves --
as = 1:5; bs = 1:5; %boundary points 

all_totalT = [];
ep = 0.001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;
i = 1;
b = 1;

figure(1); subplot(1,2,1);
grid on; hold on;

%run through varying y_a values for plots
for a = as 

xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));

sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

plot(sol.x,sol.y(1,:), '-','LineWidth',3)

%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);

%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;
all_totalT(i) = totalT;
i = i + 1;

end


xlabel('x', 'Fontsize', 24, 'Interpreter', 'latex')
y_label = ylabel('y', 'Fontsize', 24, 'Interpreter', 'latex');
%y_label.Position(1) = -80;
title(sprintf('Varying Starting Points'), 'Fontsize', 22, 'Interpreter', 'latex')
%xlim([0 b]);
%ylim([0 a+0.1]);
%set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle')

%disp(['Total travel time from x = 0 to x = ',num2str(b),' --> T = ',num2str(totalT)])

all_totalT = [];
ep = 0.001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;
i = 1;
a = 1;

subplot(1,2,2);
grid on; hold on;

%vary endpoints x_b for plots
for b = bs 

xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));


sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

plot(sol.x,sol.y(1,:), '-','LineWidth',3)

%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);

%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;
all_totalT(i) = totalT;
i = i + 1;

end


xlabel('x', 'Fontsize', 24, 'Interpreter', 'latex')
y_label = ylabel('y', 'Fontsize', 24, 'Interpreter', 'latex');
%y_label.Position(1) = -80;
title(sprintf('Varying Ending Points'), 'Fontsize', 22, 'Interpreter', 'latex')
%xlim([0 b]);
%ylim([0 a+0.1]);
%set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle')

%------ Funcitons -------

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



