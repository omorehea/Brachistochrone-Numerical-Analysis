%{
File: Brachistochrone2.m
Author: Owen Morehead
Date: Feb 21, 2022
Purpose: Record total travell times for all curves for varying 
         x_b (b in code) and y_a (a in code)
%}

%calculate total times for varying endpoints b_s

%use BVP4C
%shift from (x_a,y_a) to (x_ep,y_ep) to avoid singularity on boundary

%all_totalT = [];
%change values to iterate over different endpoint coordinates
as = 4; bs = 4; %boundary points



ep = 0.0001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;
syms x
i = 1;  %iteration holder

%tt_s = zeros(bs,1);
%times1 = sym(zeros(bs,1)); times2 = sym(zeros(bs,1)); 
%times3 = sym(zeros(bs,1)); times4 = sym(zeros(bs,1));
%times1 = zeros(bs,1); times2 = zeros(bs,1); 
%times3 = zeros(bs,1); times4 = zeros(bs,1);

a = 1;
for b = 1:bs 

xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));


sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

%figure(1); grid on; hold on;
%plot(sol.x,sol.y(1,:), '-','LineWidth',3)


%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);

%disp( " ")
%disp(['Time from x = 0 to x = ',num2str(ep),' --> T = ',num2str(T)])

%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;
all_totalTb(i) = totalT;


% ------- Comparing with other defined curves --------


f1 = a - (a/b)*x;                       %line
f2 = a*(1-(1-(x/b-1)^2)^(1/2));          %ellipse
f3 = ((a-0)/(0-b)^2)*(x-b)^2 + 0;       %parabola
%f4 = (-(x-1)/(x+1))^1;      %more...

timesb1(i) = curvetime(f1,a,b);
timesb2(i) = curvetime(f2,a,b);
timesb3(i) = curvetime(f3,a,b);
%times4(i) = curvetime(f4,a,b);

i = i + 1;

end


%calculate total times for varying startingpoints y_a

%use BVP4C
%shift from (x_a,y_a) to (x_ep,y_ep) to avoid singularity on boundary

%all_totalTa = [];

ep = 0.0001;   %shift numerical integration a distance epsilon away from x = 0
N = 10000;
syms x
i = 1;  %iteration holder

%tt_s = zeros(bs,1);
%times1 = sym(zeros(bs,1)); times2 = sym(zeros(bs,1)); 
%times3 = sym(zeros(bs,1)); times4 = sym(zeros(bs,1));
%times1 = zeros(bs,1); times2 = zeros(bs,1); 
%times3 = zeros(bs,1); times4 = zeros(bs,1);

b = 1;
for a = 1:as

xs = ep:1/N:b; %form evenly spaced integration range

xmesh = linspace(ep,b,N);
solinit = bvpinit(xmesh,@(x)guess(x,a,b));


sol = bvp4c(@(x,y)bvpfcn(x,y,a),@(ya,yb)bcfcn(ya,yb,ep,a),solinit);

%figure(1); grid on; hold on;
%plot(sol.x,sol.y(1,:), '-','LineWidth',3)


%find constant k from y(epsion)
k = ((a - sol.y(1,1))/(3*ep/2)^(2/3))^(3);

%find total time from x = 0 to x = ep. Can analytically solve approximated integral
T = (sqrt(2)/sqrt(9.8))*(3/2)^(1/3)*(ep)^(1/3)*k^(1/6);

%disp( " ")
%disp(['Time from x = 0 to x = ',num2str(ep),' --> T = ',num2str(T)])

%numeically integrate from x = ep to x = x_b to find total time
fxx = (1/sqrt(2*9.8))*((1 + sol.y(2,:).^2)./(a - sol.y(1,:))).^(1/2);

Q = trapz(sol.x,fxx);  %integrate fx over range from epsilon to b
totalT = T + Q;
all_totalTa(i) = totalT;




%xlabel('x', 'Fontsize', 24, 'Interpreter', 'latex')
%y_label = ylabel('y', 'Fontsize', 24, 'Interpreter', 'latex');
%y_label.Position(1) = -80;
%title(sprintf('Brachistochrone Curves for Varying Starting Points'), 'Fontsize', 22, 'Interpreter', 'latex')
%xlim([0 b]);
%ylim([0 a+0.1]);
%set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle')

%disp(['Total travel time from x = 0 to x = ',num2str(b),' --> T = ',num2str(totalT)])

%text(.8*b,.55*a,sprintf('A = (0,%.0f)',a), 'Fontsize', 15,'Interpreter', 'latex')
%text(.8*b,.48*a,sprintf('B = (%.0f,0)',b), 'Fontsize', 15,'Interpreter', 'latex')


% ------- Comparing with other defined curves --------


f1 = a - (a/b)*x;                       %line
f2 = a*(1-(1-(x/b-1)^2)^(1/2));          %ellipse
f3 = ((a-0)/(0-b)^2)*(x-b)^2 + 0;       %parabola
%f4 = (-(x-1)/(x+1))^1;      %more...

timesa1(i) = curvetime(f1,a,b);
timesa2(i) = curvetime(f2,a,b);
timesa3(i) = curvetime(f3,a,b);
%times4(i) = curvetime(f4,a,b);

i = i + 1;

end



%Plotting results
%for varying x_b

figure(1);

%bs = 10;

bs_array = 1:bs;

%calculating analytical total times
cB = [1.14583,1.0344,0.668618,0.939986,.680173,.800717,.912421,1.10312,1.298,1.47408];
thetaB = [-2.41501,-3.50837,-8.98682,-8.98682,-15.4483,-15.4506,-15.4495,-15.0865,-14.7086,-14.502];
times_B = -thetaB.*(cB/(2*(9.8))).^(1/2);

subplot(1,2,1);
grid on; hold on;

plot(bs_array,all_totalTb,'LineWidth',4) %numerical brach
%plot(bs_array,times_B,'LineWidth',4)    %analytical brach
plot(bs_array,timesb1,'LineWidth',1.8)   %line
plot(bs_array,timesb2,'LineWidth',1.8)   %ellipse
plot(bs_array,timesb3,'LineWidth',1.8)   %parabola

set(gca,'xtick',[0:1:bs+1])
%set(gca,'ytick',[0:.5:(round(timesb1(bs)) + 2)])


xlabel('$x_b$', 'Fontsize', 28, 'Interpreter', 'latex')
y_label = ylabel('Travel time (s)', 'Fontsize', 24, 'Interpreter', 'latex');
title(sprintf('Varying $x_b$'), 'Fontsize', 22, 'Interpreter', 'latex')
legend("Brachistochrone","Line","Ellipse","Parabola",'fontsize',19,'interpreter','latex');





%Plotting results
%for varying y_a

%as = 15;

as_array = 1:as;

subplot(1,2,2);
grid on; hold on;

plot(as_array,all_totalTa,'LineWidth',4)
plot(as_array,timesa1,'LineWidth',1.8)  %line
plot(as_array,timesa2,'LineWidth',1.8)  %ellipse
plot(as_array,timesa3,'LineWidth',1.8)  %parabola

set(gca,'xtick',[0:1:as+1])
%set(gca,'ytick',[0:.2:(round(timesa1(as)) + 2)])


xlabel('$y_a$', 'Fontsize', 28, 'Interpreter', 'latex')
y_label = ylabel('Travel time (s)', 'Fontsize', 24, 'Interpreter', 'latex');
title(sprintf('Varying $y_a$'), 'Fontsize', 22, 'Interpreter', 'latex')
legend("Brachistochrone","Line","Ellipse","Parabola",'fontsize',19,'interpreter','latex');




%function to calculate total travel time for any provided curve
function curves_t = curvetime(f,a,b)
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


