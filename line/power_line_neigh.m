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

%centers of charger placement areas
%     x_c_T = randi(5,1,10);%[ 0 4 1];
%     y_c_T = randi(5,1,10);%[ 0 2 5];
% x_c_T = [ 2.5 2.5 1 6 4 ];
% y_c_T = [ 2.1 4.5 1 5 7];
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
c_r = 1;

neigh_radius = c_r;

[x_c_T,y_c_T] = sparse_c(10,stopx,stopy);
[locDev]=locations(nDev, stopx,stopy,minAllowableDistance,x_c_T,y_c_T);

x_c_R = locDev(1,:);
y_c_R = locDev(2,:);

x_c_Ti = x_c_T;
y_c_Ti = y_c_T;

num_iterations=100;
% iter_step=100;
% % more=0;
% neigh_size=int16(num_iterations/iter_step);
% neigh_gain = zeros(1,neigh_size);
neigh_gain=zeros(4,num_iterations);
every_iter_power1=zeros(num_iterations,90);
every_iter_power2=zeros(num_iterations,90);
every_iter_power3=zeros(num_iterations,90);
every_iter_power4=zeros(num_iterations,90);
power_per_rep = zeros(4,90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Main  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iterations=1:num_iterations


%     x_c_T = [ 2.5 2.5 1 6 4 ];
%     y_c_T = [ 2.1 4.5 1 5 7];


    %find random devices positions

    % nDev = length(x_c_R);


    x_c_T = x_c_Ti;
    y_c_T = y_c_Ti;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    powers=[];

    total_power_all_chargers=0;

    for rep=1:90

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

    %     powers=[powers;final_total_power_received];
        power_per_rep(1,rep) =  final_total_power_received;
    end

        gain = final_total_power_received;% - init_total_power_received;
        realtive_gain = gain;%/init_total_power_received;


        neigh_gain(1,iterations)=realtive_gain;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   increased radius +1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_c_T = x_c_Ti;
    y_c_T = y_c_Ti;

    for i = 1:length(x_c_T)
        for j = 1:nDev
            distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
        end
    end
%     init_total_power_received = sum(total_power( x_c_T,1:nDev,distance,lambda));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %find the devices that are in the range of each charger
    C=cell(1,length(x_c_T));
    for i=1:length(x_c_T)
        for j=1:length(x_c_R)
            if norm([x_c_R(j) y_c_R(j)] - [x_c_T(i) y_c_T(i)])< neigh_radius+1
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

    powers=[];

    total_power_all_chargers=0;

    for rep=1:90

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
        power_per_rep(2,rep) = final_total_power_received;
    %     powers=[powers;final_total_power_received];
    end



        gain = final_total_power_received;% - init_total_power_received;
        realtive_gain = gain;%/init_total_power_received;


        neigh_gain(2,iterations)=realtive_gain;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   increased radius +2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_c_T = x_c_Ti;
    y_c_T = y_c_Ti;

    for i = 1:length(x_c_T)
        for j = 1:nDev
            distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
        end
    end
%     init_total_power_received = sum(total_power( x_c_T,1:nDev,distance,lambda));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %find the devices that are in the range of each charger
    C=cell(1,length(x_c_T));
    for i=1:length(x_c_T)
        for j=1:length(x_c_R)
            if norm([x_c_R(j) y_c_R(j)] - [x_c_T(i) y_c_T(i)])< neigh_radius+2
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

    powers=[];

    total_power_all_chargers=0;

    for rep=1:90

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
        power_per_rep(3,rep) = final_total_power_received;
    %     powers=[powers;final_total_power_received];
    end



    gain = final_total_power_received; %- init_total_power_received;
    realtive_gain = gain;%/init_total_power_received;

    neigh_gain(3,iterations)=realtive_gain;







    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   increased radius +3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_c_T = x_c_Ti;
    y_c_T = y_c_Ti;

    for i = 1:length(x_c_T)
        for j = 1:nDev
            distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [x_c_R(j) y_c_R(j)]);
        end
    end
%     init_total_power_received = sum(total_power( x_c_T,1:nDev,distance,lambda));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %find the devices that are in the range of each charger
    C=cell(1,length(x_c_T));
    for i=1:length(x_c_T)
        for j=1:length(x_c_R)
            if norm([x_c_R(j) y_c_R(j)] - [x_c_T(i) y_c_T(i)])< neigh_radius+15
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

    powers=[];

    total_power_all_chargers=0;

    for rep=1:90

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
        power_per_rep(4,rep) = final_total_power_received;
    %     powers=[powers;final_total_power_received];
    end



    gain = final_total_power_received;% - init_total_power_received;
    realtive_gain = gain;%/init_total_power_received;

    neigh_gain(4,iterations)=realtive_gain;
    every_iter_power1(iterations,:) = power_per_rep(1,:);
    every_iter_power2(iterations,:) = power_per_rep(2,:);
    every_iter_power3(iterations,:) = power_per_rep(3,:);
    every_iter_power4(iterations,:) = power_per_rep(4,:);
end


init(1,1:90)=init_total_power_received;

% neigh_gain = neigh_gain/num_iterations;
% gain=[init_total_power_received/100 neigh_gain(1,num_iterations) neigh_gain(2,num_iterations) neigh_gain(3,num_iterations) neigh_gain(4,num_iterations)];
% br=bar(gain);
% br.FaceColor = 'flat';
% br.CData(1,:) = [0 0.4 0.7];
% br.CData(1,:) = [0 0.45 0.74];
% br.CData(2,:) = [0.85 0.33 0.1];
% br.CData(3,:) = [0.93 0.69 0.13];
% br.CData(4,:) = [0.4 0.6 0.7];
% xlabel('Radius length') % x-axis label
% ylabel('Mean Gain') % y-axis label


init(1,1:90)=init_total_power_received;
figure
set(gca,'xticklabel',{'1','2','3','4'})
x=1:90;
plot(...
     x,mean(every_iter_power1,1),'g-',...
     x,mean(every_iter_power2,1),'r-',...
     x,mean(every_iter_power3,1),'k-',...
     x,mean(every_iter_power4,1),'b-','LineWidth',2)
xlabel('Time(Rounds)')
ylabel('Cumulative Power(Watts)')
legend('Radius 1m','Radius 2m','Radius 3m','Unbounded Radius','northoutside','Orientation','horizontal')


init(1,1:90)=init_total_power_received;
figure
set(gca,'xticklabel',{'1','2','3','4'})
x=1:90;
plot(x,init,'c-',...
     x,(every_iter_power1(3,:)),'g-',...
     x,(every_iter_power2(3,:)),'r--',...
     x,(every_iter_power3(3,:)),'k-.',...
     x,(every_iter_power4(3,:)),'b:','LineWidth',2)
xlabel('Time(Rounds)')
ylabel('Cumulative Power(Watts)')
legend('Initial Power','Radius 1m','Radius 2m','Radius 3m','Unbounded Radius','northoutside','Orientation','horizontal')
