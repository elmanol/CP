%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Lambda Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Initializiations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plane bounds
stopx = 10;
stopy = 10;

%plane discretization step
step = 0.001;

%lambda
lambda = 0.3;

%centers of charger placement areas
x_c_Ti = [ 0.3 5];
y_c_Ti = [ 0 0];

x_c_T = x_c_Ti;
y_c_T = y_c_Ti;


%devices positions
x_c_R=[3 4];
y_c_R=[0 0];

for i = 1:length(x_c_T)
    for j = 1:length(x_c_R)
        distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_power_received = total_power( x_c_T,[1 2],distance,lambda);
y=y_c_T(1);
power1=[];
power2=[];

line = x_c_Ti(1)-lambda/2 : 0.001 : x_c_Ti(1)+lambda/2;

for x=line
    x_c_T(1)=x;
    for i = 1:length(x_c_T)
        for j = 1:length(x_c_R)
            distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
        end
    end
    total_power_received = total_power( x_c_T,[1 2],distance,lambda);
    
    power1=[power1 total_power_received(1,1)];
    power2=[power2 total_power_received(1,2)];
    
end




subplot(3,1,1);
plot(line,power1,"b-");
title('Device 1 power')
xlabel('Charger 1 position')
ylabel('Device 1 Power(Watts)')


subplot(3,1,2);
plot(line,power2,"r-");
title('Device 2 power')
xlabel('Charger 1 position')
ylabel('Device 2 Power(Watts)')

subplot(3,1,3);
plot(line,power1+power2,"g-");
title("Sum of both devices' power")
xlabel('Charger 1 position')
ylabel('Device 2 Power(Watts)')


