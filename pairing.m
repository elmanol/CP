%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Charger Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%centers of charger placement areas
x_c_T = [ 2.5 2.5];
y_c_T = [ 2.1 4.5];
x_c_Ti = x_c_T;
y_c_Ti = y_c_T;

%radius of charger placement areas
r = 2;

%charger radius
r_c = 3.5;

%number of devices
nPoints = 4;

%minimum allowed distance from the chargers
minAllowableDistance = r;

%find random devices positions
[locDev]=locations(nPoints, stopx,stopy,minAllowableDistance,x_c_T,y_c_T);
%  locDev = [8.2000    6.4500    8.1000    3.5000    8.7500;
%      7.9500    3.8000    5.3500    9.4000    5.5000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Main  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cell matrix containing the radiuses of the circles for each device 
%for each charger
dC = cell(length(locDev),length(x_c_T));

%find devices initial power
[Pt, Gt, Gr, lamda, k,P_Transfered]=powers( x_c_T,y_c_T,stopx,stopy,step);
dev_power = zeros(1,length(locDev));
for i=1:length(locDev)
    dev_power(i) = P_Transfered(int16(locDev(1,i)/0.05+1), int16(locDev(2,i)/0.05+1));
end

C=cell(1,length(x_c_T));
for i=1:length(locDev)
    for j=1:length(x_c_T)
        if norm([locDev(1,i) locDev(2,i)] - [x_c_T(j) y_c_T(j)])< (r+r_c)
           C{j}=[C{j} i];
        end
    end
end

%find the devices that are in the range of both chargers
inter = intersect(C{1},C{2});

for i=1:numel(inter)
    %dokimastiko
    device1 = [locDev(1,inter(1)) locDev(2,inter(1))];
    
    device = [locDev(1,inter(i)) locDev(2,inter(i))];
    
	d1=norm( device1 - [x_c_T(1) y_c_T(1)]);
    d2=norm( device1 - [x_c_T(2) y_c_T(2)]);
    move=(d1<d2)+1;

    %find the line intersecting both the device and the charger
    coefficients = polyfit([device(1) x_c_Ti(move)], [device(2) y_c_Ti(move)], 1);
    a = coefficients (1);
    b = coefficients (2);
    
    %find the intersections between the line and the circle
    [C_interx, C_intery] = linecirc(a,b,x_c_Ti(move),y_c_Ti(move),r);


    
    %find the line between the two intersections
    [ linex, liney ] = liner([C_interx(1),C_intery(1)], [C_interx(2),C_intery(2)]);

   
    %get the closest point on the line (and circle) to the device
    closeToDev=[C_interx(1),C_intery(1)];
    
    %find points where d1-d2=k*lambda
    d1=norm( device - closeToDev);
    d2=norm( device - [x_c_T(~(move-1)+1) y_c_T(~(move-1)+1)]);
    substract = d1>d2;
    dif = abs(d1-d2);
    moddif= mod(dif,lambda);
    minr=10^3;
    
    %find the first point where d1-d2==k*lambda in the circle
    for k=1:numel(linex)
        dist_ = norm([linex(k) liney(k)]-closeToDev);

        dif = abs(moddif - dist_);
        if  dif < minr && norm([linex(k) liney(k)]-[x_c_T(move) y_c_T(move)])<=r
            minr = dif;
            mulK = [linex(k) liney(k)];        
            difw=dif;
        end
    end
    dC{inter(i),move}=[dC{inter(i),move} norm(mulK(size(mulK,1),:)- device)];
    
    %find the rest of the points where d1-d2==k*lambda
    for k=1:numel(linex)
        dist_ = norm([linex(k) liney(k)]-mulK(size(mulK,1),:));
        if dist_> lambda-0.004 && dist_< lambda+0.004 && norm([linex(k) liney(k)]-[x_c_Ti(move) y_c_Ti(move)])<=r
            mulK = [mulK;linex(k) liney(k) ];
            dC{inter(i),move}=[dC{inter(i),move} norm(mulK(size(mulK,1),:)- device)];
        end
    end
    
%     x_c_T(move)= mulK(2,1);
%     y_c_T(move)= mulK(2,2);

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Plots  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot the power distribution
[Pt, Gt, Gr, lamda, k,P_Transfered]=powers( x_c_T,y_c_T,stopx,stopy,step);
[X,Y] = meshgrid(0:step:stopx,0:step:stopy);
abs(P_Transfered);
figure
surf(Y,X,abs(P_Transfered));
hold on

    
%plot the devices positions

plot(locDev(1,:), locDev(2,:), 'g*');

%plot the placement areas
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

%plot the circles for possible positions
colors=["w-","g-","r-","y-","o-"];
for k=1:numel(inter)
for i=1:numel(dC{inter(k),move})
    th = 0:step:2*pi+step;
    xunit = dC{inter(k),move}(i) * cos(th) + locDev(1,inter(k));
    yunit = dC{inter(k),move}(i) * sin(th) + locDev(2,inter(k));
    
    %check for points outside the plane
    for j=1:numel(xunit)
        if xunit(j)<0 xunit(j)=0;
        elseif xunit(j)>stopx xunit(j)=stopx;end

        if yunit(j)<0 yunit(j)=0;
        elseif yunit(j)>stopy yunit(j)=stopy;end
    end
    plot(xunit, yunit,colors(k),'LineWidth',1.5);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



