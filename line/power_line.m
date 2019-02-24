%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Line Charger Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_iterations=400;
iter_step=100;
more=0;
neigh_size=int16(num_iterations/iter_step);
neigh_gain = zeros(1,neigh_size);

neigh_radius=1;
u=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Initializiations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iterations=1:num_iterations
    if mod(iterations,iter_step)==0 && iterations~=num_iterations
        more=more+1;
    end

%plane bounds
stopx = 10;
stopy = 10;

%plane discretization step
step = 0.05;

%lambda
lambda = 0.3;

%centers of charger placement areas
x_c_T = randi(5,1,10);%[ 0 4 1];
y_c_T = randi(5,1,10);%[ 0 2 5];
x_c_T = [ 2.5 2.5 1 6 4 ];
y_c_T = [ 2.1 4.5 1 5 7];
%devices positions
locDevx=[3.5 3];
locDevy=[0 0];
%or
%number of devices for random placement
nDev=9; %an to allaksw allazw kai to x_c_R
 
%minimum allowed distance from the chargers
minAllowableDistance = 3*lambda/2;

%find random devices positions
[locDev]=locations(nDev, stopx,stopy,minAllowableDistance,x_c_T,y_c_T);
x_c_R = locDev(1,:);
y_c_R = locDev(2,:);
nDev = length(x_c_R);

% [x_c_T,y_c_T, x_c_R,y_c_R,distance] = deployment(5, 10, stopx,lambda)
for i = 1:length(x_c_T)
    for j = 1:nDev
        distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
    end
end
init_total_power_received = sum(total_power( x_c_T,[1:nDev],distance,lambda));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the devices that are in the range of each charger
C=cell(1,length(x_c_T));
for i=1:length(x_c_T)
    for j=1:length(x_c_R)
        if norm([x_c_R(j) y_c_R(j)] - [x_c_T(i) y_c_T(i)])< neigh_radius+more
           C{i}=[C{i} j];
        end
    end 
end


%charger movement step
x_step = 0.01;
num_steps = lambda/x_step;
total_power_received=zeros(length(x_c_T),num_steps+1);%,length(x_c_R));
final_total_power_received = total_power_received;
past_power = ones(length(x_c_T),num_steps+1);
repeats=0;

while(abs(past_power-final_total_power_received)>0.0001)
%for i=1:length(x_c_T)
%   total_power_received=zeros(num_steps+1);%,length(x_c_R));

    i=randi(length(x_c_T),1,1);
    if x_c_T(i)-lambda/2>=0
    for d = 1:nDev
        distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);
    end
        x_start = x_c_T(i)-lambda/2;
    else
        x_start = 0;
    end
    x_c_T(i)=x_start;
    if ~isempty(C{1,i})
        total_power_received(i,1) = sum(total_power( x_c_T,C{1,i},distance,lambda));
    end  
    x=x_start;
    for j=2:num_steps+1     
        x=x+x_step;
        x_c_T(i)=x;
            for d = 1:nDev
                distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);
            end
        if ~isempty(C{1,i})
            total_power_received(i,j) = sum(total_power( x_c_T,C{1,i},distance,lambda));
        end
    end
    [max_pow,max_pos] = max(total_power_received(i,:));
    x_c_T(i)=x_start+x_step*(max_pos-1);
    for d = 1:nDev
        distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);
    end
%end
    past_power = final_total_power_received;
    final_total_power_received = sum(total_power( x_c_T,[1:nDev],distance,lambda));

end
    neigh_gain(more+1) = neigh_gain(more+1) + final_total_power_received;
    
  gain = gain+final_total_power_received - init_total_power_received;

end
plot([1:neigh_size],neigh_gain,"g*");
%gain_mean = gain/iterations;