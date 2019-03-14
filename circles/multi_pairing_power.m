%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Circle Charger Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Initializiations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plane bounds
stopx = 10;
stopy = 10;

%plane discretization step
step = 0.05;

%lambda
lambda = 0.3;

nPoints = 50;

r = lambda;

minAllowableDistance = r+lambda;
% x_c_T = [ 2.5 2.5];
% y_c_T = [ 4 7];
%centers of charger placement areas
% num_chargers = 5;55
% x_c_T = [ 2.5 2.5 1 6 4 ];
% y_c_T = [ 2.1 4.5 1 5 7];

for i=1:100
    [x_c_Tg(i,:),y_c_Tg(i,:)] = sparse_c(10,stopx,stopy,2*r);
    [locDev]=locations(nPoints, stopx,stopy,minAllowableDistance,x_c_Tg(i,:),y_c_Tg(i,:));
    locDevx(i,:)=locDev(1,:);
    locDevy(i,:)=locDev(2,:);
end
initial_sum_total_power_received=0;

for iter=1:1
    
x_c_T = x_c_Tg(iter,:);
y_c_T = y_c_Tg(iter,:);

locDev = [locDevx(iter,:);locDevy(iter,:)];
%radius of charger placement areas
% r = 2*lambda;

%charger radius
r_c = 2.5;

%number of devices
%nPoints = 20;

%minimum allowed distance from the chargers
%minAllowableDistance = r+lambda;

%set pairing way ('random' or 'closest')
pairing_way='closest';


%grid placement

% x_c_T = [3 6 9 3 6 9];
% y_c_T = [3 3 3 8 8 8];
% 
% [x_c_T,y_c_T] = sparse_c(2,stopx,stopy,2*r);

% x_c_T = x_c_Ti;
% y_c_T = y_c_Ti;
%save chargers initial positions
x_c_Ti = x_c_T;
y_c_Ti = y_c_T;




%find random devices positions
% [locDev]=locations(10, stopx,stopy,minAllowableDistance,x_c_T,y_c_T);
%  locDev = [8.2000    6.4500    8.1000    3.5000    8.7500;
%      7.9500    3.8000    5.3500    9.4000    5.5000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Main  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cell matrix containing the radiuses of the circles for each device 
%for each charger
dCradius = cell(length(locDev),length(x_c_T));

%find devices initial power
% [Pt, Gt, Gr, lamda, k,P_Transfered]=powers( x_c_T,y_c_T,stopx,stopy,step);
% dev_power = zeros(1,length(locDev));
% for i=1:length(locDev)
%     dev_power(i) = P_Transfered(int16(locDev(1,i)/0.05+1), int16(locDev(2,i)/0.05+1));
% end


%find devices in each charger's radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=cell(1,length(x_c_T));
for i=1:length(locDev)
    for j=1:length(x_c_T)
        if norm([locDev(1,i) locDev(2,i)] - [x_c_T(j) y_c_T(j)])< (r+r_c+2)
           C{j}=[C{j} i];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get a random charger x
x=randi([1,length(x_c_T)]);
x=1;
x_range = C{x};
closeToDevArray=[];
allint=[];
%find the centroid of x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(x_range)
    interi=x_range(i);
    coefficients = polyfit([locDev(1,interi) x_c_Ti(x)], [locDev(2,interi) y_c_Ti(x)], 1);
    if  abs(x_c_Ti(x)-locDev(1,interi))>0.001
        a = coefficients (1);
        b = coefficients (2);
        [C_interx, C_intery] = linecirc(a,b,x_c_Ti(x),y_c_Ti(x),r);
    else
        C_intery = [y_c_Ti(x)+r y_c_Ti(x)-r];
        C_interx = [x_c_Ti(x) x_c_Ti(x)];
    end
    
    minc=10^3;
    for j=1:numel(C_interx) 
        clos=norm([C_interx(j) C_intery(j)]-  [x_c_Ti(x) y_c_Ti(x)]);
        if clos<minc
           minc = clos;
           iterator=j;
        end
    end
    closeToDev=[C_interx(iterator),C_intery(iterator)];
    closeToDevArray = [closeToDevArray; closeToDev];
end

polyin1=[];
if size(closeToDevArray,1)>2
    polyin = polyshape(closeToDevArray(:,1),closeToDevArray(:,2));
    [centrx1,centry1] = centroid(polyin);
     x_c_T(x)=centrx1;
     y_c_T(x)=centry1;      
elseif isempty(closeToDevArray)
     x_c_T(x)=x_c_Ti(x);
     y_c_T(x)=y_c_Ti(x);
elseif size(closeToDevArray,1)==2
     x_c_T(x)=(closeToDevArray(1,1)+closeToDevArray(2,1))/2;
     y_c_T(x)=(closeToDevArray(1,2)+closeToDevArray(2,2))/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




closeToDevArray=[];
x_new=-1;
plots=[];
chargers_remained=1:length(x_c_T);
chargers_remained(x)=0;
p=0;
pairs=[];
while sum(chargers_remained==0)~=length(x_c_T)
    
    %gets the next closest charger for pairing 
    %buildin a Minimum spanning tree
    if strcmp(pairing_way,'closest')
        
        mind=10^3;
        for i=1:length(chargers_remained)
            c_iter1=chargers_remained(i);
            distc1=0;
            if c_iter1~=0
                for j=1:length(chargers_remained)
                    c_iter0=chargers_remained(j);
                    if c_iter0==0 
                        distc1 = norm([x_c_T(c_iter1) y_c_T(c_iter1)]- [x_c_T(j) y_c_T(j)]);
                        if distc1 < mind     
                           mind = distc1;
                           x_new = c_iter1;
                           x = j;
                        end
                    end
                end
            end
        end
        
    elseif strcmp(pairing_way,'random')
       
       %randomly find a charger that hasn't been used yet 
       rand_iter = randi(length(chargers_remained),1);
       while chargers_remained(rand_iter)==0
            rand_iter = randi(length(chargers_remained),1);
       end
       
       rand_iter_x = randi(length(chargers_remained),1);
       while chargers_remained(rand_iter_x)==chargers_remained(rand_iter)
            rand_iter_x = randi(length(chargers_remained),1);
       end
       x_new = rand_iter;
       x=rand_iter_x;
    else
        
        print('wrong pairing way selected');
        
    end
    chargers_remained(x_new)=0;
    
    inter = intersect(C{x},C{x_new});

    %sort the chargers according their distance to chargers
    inter_dist = zeros(1,numel(inter));
    for i=1:numel(inter)
         x1 = norm([x_c_T(x) y_c_T(x)]- [locDev(1,inter(i)) locDev(1,inter(i))]);
         x2 = norm([x_c_T(x_new) y_c_T(x_new)]- [locDev(1,inter(i)) locDev(1,inter(i))]);
         x_mean = x1+x2;
         inter_dist(i) = x_mean;        
    end
    [~,pos]=sort(inter_dist,'ascend');
    inter = inter(pos);
    
%     pairs = [pairs; x_c_T(x) y_c_T(x) x_c_T(x_new) y_c_T(x_new)];
    found=[];
    
    inter_pair_num = (numel(inter)*(numel(inter)-1))/2;
    total_errors = zeros(inter_pair_num,4);
    
    %for every device that is in the raius of both chargers 
    allint=[];
    best_int_power=[];
    for i=1:numel(inter)
        
        %get the coordinates of the device
        device = [locDev(1,inter(i)) locDev(2,inter(i))];


        %find the line intersecting both the device and the charger
        %if line is not vertical
        if  abs(x_c_Ti(x_new)-device(1))>0.001
            coefficients = polyfit([device(1) x_c_Ti(x_new)], [device(2) y_c_Ti(x_new)], 1);
            a = coefficients (1);
            b = coefficients (2);
            %find the intersections between the line and the circle
            [C_interx, C_intery] = linecirc(a,b,x_c_Ti(x_new),y_c_Ti(x_new),r);
        else
            %in case the line is vertical
            C_intery = [y_c_Ti(x_new)+r y_c_Ti(x_new)-r];
            C_interx = [x_c_Ti(x_new) x_c_Ti(x_new)];
        end


        %get the closest point on the line (and circle) to the device
        dc1=norm([C_interx(1), C_intery(1)]-device);
        dc2=norm([C_interx(2), C_intery(2)]-device);
        if dc1<=dc2
            closeToDev=[C_interx(1),C_intery(1)];
            %find the line between the two intersections
            [ linex, liney ] = liner(closeToDev, [C_interx(2),C_intery(2)]);
        else
            closeToDev=[C_interx(2),C_intery(2)];
            %find the line between the two intersections
            [ linex, liney ] = liner(closeToDev,[C_interx(1),C_intery(1)]);
        end

    %         %find the line between the two intersections
    %     [ linex, liney ] = liner([C_interx(1),C_intery(1)], [C_interx(2),C_intery(2)]);

        
        %GIA TO POLYGONO
        closeToDevArray = [closeToDevArray; closeToDev];
        
        %find points where d1-d2=k*lambda
        d1=norm( device - closeToDev);
        d2=norm( device - [x_c_T(x) y_c_T(x)]);
        %substract = d1>d2;
        dif = abs(d1-d2);
        moddif= mod(dif,lambda);
                
        %find the first point where d1-d2==k*lambda in the circle
        minr=10^3;
        for k=1:numel(linex)
            dist_ = norm([linex(k) liney(k)]-closeToDev);

            dif = abs(moddif - dist_);
            if  dif < minr %&& norm([linex(k) liney(k)]-[x_c_T(x_new) y_c_T(x_new)])<=r
                minr = dif;
                mulK = [linex(k) liney(k)];        
            end
        end
       
        
        dCradius{inter(i),x_new}=[dCradius{inter(i),x_new} norm(mulK(size(mulK,1),:)- device)];

        %find the rest of the points where d1-d2==k*lambda
        for k=1:numel(linex)
            dist_ = norm([linex(k) liney(k)]-mulK(size(mulK,1),:));
            if dist_> lambda-0.01 && dist_< lambda+0.01 %&& norm([linex(k) liney(k)]-[x_c_Ti(x_new) y_c_Ti(x_new)])<=r
                mulK = [mulK;linex(k) liney(k) ];
                dCradius{inter(i),x_new}=[dCradius{inter(i),x_new} norm(mulK(size(mulK,1),:)- device)];
            end
        end
        
        %if there is only one device in the radius of the two chargers,
        %we must choose the point on the first circle that is the closest 
        %to the device
        if length(inter)==1
            mul=1;
            while mulK(mul,1)<0 || mulK(mul,1)>stopx || mulK(mul,2)<0 || mulK(mul,2)>stopy 
                mul=mul+1;
            end
        	x_c_T(x_new)=mulK(mul,1);
            y_c_T(x_new)=mulK(mul,2);
            continue
        end
        
        
        
            
            %if i>=2, we can create a pair of inter and find their
            %circles' intersections
            if i>=2
                %find the intersection between inter(i) and every previous
                %inter (inter(i-1),inter(i-2)...)
                for j=1:i-1
                    intersections = [];
                    x1  = locDev(1,inter(j));
                    y1  = locDev(2,inter(j));        
                    x2  = locDev(1,inter(i));
                    y2  = locDev(2,inter(i));
                    
                    %get the intersections between inter(i) and inter(j)
                    %into intersections vector
                    for c2 = 1:size(dCradius{inter(i),x_new},2)
                        r2 = dCradius{inter(i),x_new}(c2);
                        for c1 = 1:size(dCradius{inter(j),x_new},2)
                            r1 = dCradius{inter(j),x_new}(c1);            
                            [xout,yout] = circcirc(x1,y1,r1,x2,y2,r2);
                            for s=1:numel(xout)
                                if (norm([xout(s) yout(s)] - [x_c_Ti(x_new) y_c_Ti(x_new)]) <= r && ...
                                	xout(s)>0 && xout(s)<stopx && yout(s)>0 && yout(s)<stopy)
                                    intersections = [intersections; xout(s) yout(s)];
                                    allint = [allint; xout(s) yout(s)];
                                end
                            end
                        end
                    end
                    

                end
       
            end
    end
        allint=allint';
        
        
        for j = 1:size(allint,2)
            x_c_T(x_new) = allint(1,j);
            y_c_T(x_new) = allint(2,j);
            for i = 1:length(x_c_T)
                for z = 1:size(locDev,2)
                    distance_int(i,z) = norm([x_c_T(i) y_c_T(i)] - [locDev(1,z) locDev(2,z)]);
                end
            end 
            int_power = sum(circles_total_power(x_c_T, 1:size(locDev,2), distance_int, lambda));
            best_int_power = [best_int_power int_power];
            
        end

        [~,k]=max(best_int_power);
        if ~isempty(allint)
            x_c_T(x_new)= allint(1,k);
            y_c_T(x_new)= allint(2,k);
        end
       
        
end
































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Initial Power   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(x_c_T)
    for j = 1:size(locDev,2)
        distance(i,j) = norm([x_c_Ti(i) y_c_Ti(i)] - [locDev(1,j) locDev(2,j)]);
    end
end

initial_total_power_received = circles_total_power(x_c_Ti, 1:size(locDev,2), distance, lambda);
initial_sum_total_power_received = sum(initial_total_power_received);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Random Power  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:length(x_c_T)    
%     signx = rand(1)>0.5;
%     signy = rand(1)>0.5;
%     x_c_Tr(i) = x_c_T(i)+r*rand(1)*(-1)^signx;
%     y_c_Tr(i) = y_c_T(i)+r*rand(1)*(-1)^signy;
% end
% 
% for i = 1:length(x_c_T)
%     for j = 1:size(locDev,2)
%         distance(i,j) = norm([x_c_Tr(i) y_c_Tr(i)] - [locDev(1,j) locDev(2,j)]);
%     end
% end
% 
% random_total_power_received = circles_total_power(x_c_Ti, 1:size(locDev,2), distance, lambda);
% random_sum_total_power_received = sum(itial_total_power_received);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    Final Power   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(x_c_T)
    for j = 1:size(locDev,2)
        distance(i,j) = norm([x_c_T(i) y_c_T(i)] - [locDev(1,j) locDev(2,j)]);
    end
end

final_total_power_received = circles_total_power(x_c_T, 1:size(locDev,2), distance, lambda);
final_sum_total_power_received = sum(final_total_power_received);


gain=final_sum_total_power_received-initial_sum_total_power_received;
realtive_gain = gain/initial_sum_total_power_received;
fprintf('Realtive gain: : %f .\n', realtive_gain*100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

power_many_power(1,iter)=final_sum_total_power_received;
initial(1,iter) =  initial_sum_total_power_received;
end %iterations stop


sum_power_many_power = sum(power_many_power);
initial_sum_one= sum(initial);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plots  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%plot the power distribution
% [Pt, Gt, Gr, lamda, k,P_Transfered]=powers( x_c_T,y_c_T,stopx,stopy,step);
% [X,Y] = meshgrid(0:step:stopx,0:step:stopy);
% abs(P_Transfered);
% figure
% surf(Y,X,abs(P_Transfered));


hold on
% if ~isempty(polyin)
%     plot(centrx,centry, 'm*','LineWidth',10);
% end  

% hold on
% plot(allint(:,1), allint(:,2), 'r*');


%%%%%%%%plot the devices positions
plot(x_c_Ti, y_c_Ti, 'g*');
plot(locDev(1,:), locDev(2,:), 'r*');

% % plot the placement areas
for i=1:numel(x_c_Ti)
    th = 0:step:2*pi+step;
    xunit = r * cos(th) + x_c_Ti(i);
    yunit = r * sin(th) + y_c_Ti(i);
    
    %check for points outside the plane
    for j=1:numel(xunit)
        if xunit(j)<0 xunit(j)=0;
        elseif xunit(j)>stopx xunit(j)=stopx;end

        if yunit(j)<0 yunit(j)=0;
        elseif yunit(j)>stopy yunit(j)=stopy;end
    end
    plot(xunit, yunit,"y-",'LineWidth',2);
end
hold on
legend('Chargers','Devices','Charger Movement Areas','Location','northoutside','Orientation','horizontal')
% 
% % % plot the circles for possible positions
% colors=["w-","g-","r-","m-","b-","r-","m-","b-","k-"];
% for k=1:numel(inter)
%     for i=1:numel(dCradius{inter(k),x_new})
%         th = 0:step:2*pi+step;
%         xunit = dCradius{inter(k),x_new}(i) * cos(th) + locDev(1,inter(k));
%         yunit = dCradius{inter(k),x_new}(i) * sin(th) + locDev(2,inter(k));
%         
%         %check for points outside the plane
%         for j=1:numel(xunit)
%             if xunit(j)<0 xunit(j)=0;
%             elseif xunit(j)>stopx xunit(j)=stopx;end
%             
%             if yunit(j)<0 yunit(j)=0;
%             elseif yunit(j)>stopy yunit(j)=stopy;end
%         end
%         plot(xunit, yunit,colors(k),'LineWidth',1.5);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% set(gca,'xticklabel',{'1','2','3','4'})
% x=1:100;
% plot(x,sum_power_many,'g-' ,x, sum_power,'r--')
% xlabel('Time(Rounds)')
% ylabel('Cumulative Power(Watts)')
% legend('Radius 1','Radius 2','Radius 3','Unbounded Radius','northoutside','Orientation','horizontal')
% 
% 
% mean_powers=[sum(sum_power)/100 sum(sum_power_many)/100];
% br=bar(mean_powers);
% br.FaceColor = 'flat';
% br.CData(1,:) = [0 0.4 0.7];
% br.CData(2,:) = [0 0.65 0.14];
% % br.CData(3,:) = [0.85 0.33 0.1];
% % br.CData(4,:) = [0.93 0.69 0.13];
% % br.CData(5,:) = [0.4 0.6 0.7];
% xlabel('Number of intersections') % x-axis label
% ylabel('Power') % y-axis label



% figure
% set(gca,'xticklabel',{'Multi','Centroid','Many'})
% x=1:100;
% plot(x,sum_power_many,'g-' ,x, sum_power,'r--')
% xlabel('Time(Rounds)')
% ylabel('Cumulative Power(Watts)')
% legend('Multi','Centroid','Many','Unbounded Radius','northoutside','Orientation','horizontal')
% 
% 
mean_powers=[initial_sum sum_power_centroids sum_power_multi_r sum_power_multi sum(line_power) sum_power_many_power]/100;
mean_powers=[initial_small centroids_small random_small multi_small sum(line_small) power_small]/100;


br=bar(mean_powers);
br.FaceColor = 'flat';
br.CData(1,:) = [0 0.4 0.7];
br.CData(2,:) = [0 0.65 0.14];
br.CData(3,:) = [0.85 0.33 0.1];
br.CData(4,:) = [0.93 0.69 0.13];
br.CData(5,:) = [0.4 0.6 0.7];
br.CData(4,:) = [0.8 0.5 0.1];
xlabel('Intersections Approach') % x-axis label
ylabel('Power') % y-axis label
legend('Multi','Centroid','Many','Unbounded Radius','northoutside','Orientation','horizontal')
names={'Initial';'Appr. C';'Appr. B';'Appr. A';'Alg. 1';'Appr. D'};
set(gca,'xticklabel',names)
