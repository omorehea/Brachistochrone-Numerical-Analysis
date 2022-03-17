%{
File: Brachistochrone4.m
Author: Owen Morehead
Date: Feb 21, 2022
Purpose: Recover a numerical error against the analytical solution.
%}

%------ Comparing Numerical to Analytical Derivation Results -----

%--- analytical results ---
cA = [1.14583,4.81121,13.8379,30.8731,58.5785,...
    99.6192,156.661,232.37,329.413,450.456];

thetaA = [-2.41201,-1.40138,-.968656,-.736425, -.592962,...
    -.495899,-.425978, -.373258,-.332108,-.299105];

times_A = -thetaA.*(cA/(2*(9.8))).^(1/2);

%---- use BVP4C ----
%shift from (x_a,y_a) to (x_ep,y_ep) to avoid singularity on boundary

all_totalT = [];

as = [2:2:10]; %array of y_a values
b = 1;         %x_b boundary value
N = 1000;
i = 1;

eps = [.0001:.001:.01]; %shift numerical integration a distance epsilon away from x = 0

error_num = zeros(10,length(eps)); %empty array for errors
figure(1); grid on; hold on;

for a = as   %loop through values of y_a
i = 1;
for ep = eps %loop through values of epsilon
  
xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));

sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);


%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);

%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to x_b
totalT = T + Q;
%all_totalT(i) = totalT; store total times in array if needed

error_num(a,i) = abs(totalT-times_A(a))/abs(totalT);
i = i+1;

end

plot(eps,error_num(a,:)); %plot error against epsilon

end

%plot(eps,error_num);
xlabel('$\epsilon$', 'Fontsize', 24, 'Interpreter', 'latex')
y_label = ylabel('Error', 'Fontsize', 24, 'Interpreter', 'latex');
title({'Error Between Numerical and Analyitcal Solutions', 'for Fixed Endpoint $B = (1,0)$'}, 'Fontsize', 20, 'Interpreter', 'latex')
legend({sprintf('A = (0,%.0f)',as(1)),...
    sprintf('A = (0,%.0f)',as(2)),...
    sprintf('A = (0,%.0f)',as(3)),...
    sprintf('A = (0,%.0f)',as(4)),...
    sprintf('A = (0,%.0f)',as(5))},...
    'fontsize',18,'Interpreter','latex');

%------- Funcitons --------

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



