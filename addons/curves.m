clear
%charger 1
a=0;
%charger 2
b=2;
lamda=1;
initial_phase = pi+0.2;
%intervalamda displayed
x = [a+lamda-0.7:0.001:b-lamda+0.7];
for counter=1:length(x)
    y(counter)=( cos((2*pi*(-a - b + 2*x(counter)))/lamda) + cos((2*pi*(-a - b + 2*x(counter)))/lamda))/((x(counter)-a)*(b-x(counter))) +1/(x(counter)-a)^2 + 1/(b-x(counter))^2;
    
end
hold on;
plot (x,(1/(4*pi)^2)*y,'b')

a=0;
%charger 2
b=2.5;
lamda=1;
initial_phase = pi+0.2;
%intervalamda displayed
x = [a+lamda-0.7:0.001:b-lamda+0.7];
for counter=1:length(x)
    ynew(counter)=( cos((2*pi*(-a - b + 2*x(counter)))/lamda) + cos((2*pi*(-a - b + 2*x(counter)))/lamda))/((x(counter)-a)*(b-x(counter))) +1/(x(counter)-a)^2 + 1/(b-x(counter))^2;
    
end


plot (x,(1/(4*pi)^2)*ynew,'c')
hold off;
ylabel('Power (watts)')
xlabel('Position (m)')
legend('initial positions','after readjastment')