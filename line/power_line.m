%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Line Charger Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Initializiations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=0;


%plane bounds
stopx = 10;
stopy = 10;

%plane discretization step
step = 0.05;

%lambda
lambda = 0.3;

% %centers of charger placement areas
%     x_c_T = rand(6,1,1)*10;%[ 0 4 1];
%     y_c_T = rand(6,1,1)*10;%[ 0 2 5];
% % x_c_T = [ 2.5 2.5 1 6 4 ];
% % y_c_T = [ 2.1 4.5 1 5 7];
% 
% x_c_Ti = x_c_T;
% y_c_Ti = y_c_T;


%devices positions
% locDevx=[3.5 3];
% locDevy=[0 0];
%or
%number of devices for random placement
nDev=50; %an to allaksw allazw kai to x_c_R

%minimum allowed distance from the chargers
minAllowableDistance = 3*lambda/2;

%charger radius
c_r = 100;

neigh_radius = 100;

%find random devices positions
% [locDev]=locations(nDev, stopx,stopy,minAllowableDistance,x_c_T,y_c_T);
% x_c_R = locDev(1,:);
% y_c_R = locDev(2,:);
% nDev = length(x_c_R);


num_iterations=100;
% iter_step=100;
% % more=0;
% neigh_size=int16(num_iterations/iter_step);
% neigh_gain = zeros(1,neigh_size);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Main  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_power=0;
for iter=1:100

x_c_T = x_c_Tg(iter,:);
y_c_T = y_c_Tg(iter,:);

locDev = [locDevx(iter,:);locDevy(iter,:)];
x_c_R = locDev(1,:);
y_c_R = locDev(2,:);

x_c_Ti = x_c_T;
y_c_Ti = y_c_T;

% for iterations=1:num_iterations

    for i = 1:length(x_c_T)
        for j = 1:nDev
            distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
        end
    end
    init_total_power_received = sum(total_power( x_c_T,1:nDev,distance,lambda));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %find the devices that are in the range of each charger
    C=cell(1,length(x_c_T));
    for i=1:length(x_c_T)
        for j=1:length(x_c_R)
            if norm([x_c_R(j) y_c_R(j)] - [x_c_T(i) y_c_T(i)])< neigh_radius
               C{i}=[C{i} j];
            end
        end 
    end


    %charger movement step
    x_step = 0.01;
    num_steps = lambda/x_step;
    total_power_received=zeros(length(x_c_T),num_steps+1);%,length(x_c_R));
    final_total_power_received = 0;
    past_power = 1;
    repeats=0;

    powerss=[];

    total_power_all_chargers=0;
    points=[];
    for rep=1:100

    %   total_power_received=zeros(num_steps+1);%,length(x_c_R));

        %get a random charger
        i=randi(length(x_c_T),1,1);

        %set initial charger position after checking if it's in the plane
        if x_c_Ti(i)-lambda/2>0

            x_start = x_c_Ti(i)-lambda/2;

        else

            x_start = 0;

        end

        %set the found initial position
        x_c_T(i)=x_start;

        %calculate the new distances, afer the charger replacement
        for d = 1:nDev

            distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);

        end

        %get the power of the devices inside the charger's radius for the first
        %step
        if ~isempty(C{1,i})

            total_power_received(i,1) = sum(total_power( x_c_T,C{1,i},distance,lambda));

        end  
        x=x_start;

        %(j gives us the number of the step) 
        %after making each step and changing the distances, get the power 
        %gained into an array (total_power_received)
        for j=2:num_steps+1 

            %check if the new step gets us out of bounds
            if x+x_step<stopx

                x=x+x_step;
                x_c_T(i)=x;

                for d = 1:nDev
                    distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);
                end

                if ~isempty(C{1,i})
                    %find the power that only the devices inside charger's i 
                    %radius get (from all the chargers)
                    total_power_received(i,j) = sum(total_power( x_c_T,C{1,i},distance,lambda));
                end

            end
            points=[points; x y_c_T(i)];
        end
        %find the step at which the power for the devices inside charger's i
        %radius was maximized
        [max_pow,max_pos] = max(total_power_received(i,:));

        %set charger's position where the power of the devices inside its
        %radius where maximized
        x_c_T(i)=x_start+x_step*(max_pos-1);

        %change the corresponding distances
        for d = 1:nDev
            distance(i,d) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(d) y_c_R(d)]);
        end

    %     past_power = final_total_power_received;
        final_total_power_received = sum(total_power( x_c_T,1:nDev,distance,lambda));

    %     powerss=[powerss;final_total_power_received];
    end


%     neigh_gain(more+1) = neigh_gain(more+1) + final_total_power_received;

    %gain = gain+final_total_power_received - init_total_power_received;


    gain = final_total_power_received - init_total_power_received;
    realtive_gain = gain/init_total_power_received*100
%     fprintf('Realtive gain: : %f .\n', realtive_gain*100);
    
    
    
%     [Pt, Gt, Gr, lamda, k,P_Transfered]=powers( x_c_T,y_c_T,stopx,stopy,step);
%     [X,Y] = meshgrid(0:step:stopx,0:step:stopy);
%     abs(P_Transfered);
%     figure
%     surf(Y,X,abs(P_Transfered));
    final_power(iter)= sum(total_power( x_c_T,1:nDev,distance,lambda));
    x_c_T = x_c_Ti;
    y_c_T = y_c_Ti;
% end
line_power = line_power + total_power( x_c_T,1:nDev,distance,lambda);
end   
hold on
%plot the devices positions
% plot(locDev(1,:), locDev(2,:), 'g*');
% plot(points(:,1), points(:,2), 'r*');



% figure(11)
% plot(x_c_T,y_c_T,'og',  x_c_Ti,y_c_Ti,'ok',  x_c_R,y_c_R,'*r')
% xlabel('x(m)')
% ylabel('y(m)')
% legend('Chargers Final','Chargers Initial','Nodes','Location','northoutside','Orientation','horizontal')
% 

% plot(1:num_iterations,final_power,'-g')         
