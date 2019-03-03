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

%centers of charger placement areas
x_c_T = [ 2.5 2.5];
y_c_T = [ 2.1 4.5];
x_c_Ti = x_c_T;
y_c_Ti = y_c_T;

%radius of charger placement areas
r = 2*lambda;

%charger radius
r_c = 2.5;

%number of devices
nPoints = 12;

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
        if norm([locDev(1,i) locDev(2,i)] - [x_c_T(j) y_c_T(j)])< (r+r_c+2)
           C{j}=[C{j} i];
        end
    end
end









move=1;

%find the devices that are in the range of both chargers
inter = intersect(C{1},C{2});
closeToDevArray=[];

for i=1:numel(inter)
    interi=inter(i);
    coefficients = polyfit([locDev(1,interi) x_c_Ti(2)], [locDev(2,interi) y_c_Ti(2)], 1);
    a = coefficients (1);
    b = coefficients (2);
    [C_interx, C_intery] = linecirc(a,b,x_c_Ti(2),y_c_Ti(2),r);
    minc=10^3;
    for j=1:numel(C_interx) 
        clos=norm([C_interx(j) C_intery(j)]-  [x_c_Ti(2) y_c_Ti(2)]);
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
closeToDevArray=[];



for i=1:numel(inter)
    %dokimastiko
    device1 = [locDev(1,inter(1)) locDev(2,inter(1))];   
    device = [locDev(1,inter(i)) locDev(2,inter(i))];
    
	d1=norm( device1 - [x_c_T(1) y_c_T(1)]);
    d2=norm( device1 - [x_c_T(2) y_c_T(2)]);
    move=(d1>d2)+1;
    move=1;

    %find the line intersecting both the device and the charger
    coefficients = polyfit([device(1) x_c_Ti(move)], [device(2) y_c_Ti(move)], 1);
    a = coefficients (1);
    b = coefficients (2);
    
    %find the intersections between the line and the circle
    [C_interx, C_intery] = linecirc(a,b,x_c_Ti(move),y_c_Ti(move),r);

   
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
    
    
    closeToDevArray = [closeToDevArray; closeToDev];
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
        if  dif < minr %&& norm([linex(k) liney(k)]-[x_c_T(move) y_c_T(move)])<=r
            minr = dif;
            mulK = [linex(k) liney(k)];        
        end
    end
    dC{inter(i),move}=[dC{inter(i),move} norm(mulK(size(mulK,1),:)- device)];
    
    %find the rest of the points where d1-d2==k*lambda
    for k=1:numel(linex)
        dist_ = norm([linex(k) liney(k)]-mulK(size(mulK,1),:));
        if dist_> lambda-0.01 && dist_< lambda+0.01 %&& norm([linex(k) liney(k)]-[x_c_Ti(move) y_c_Ti(move)])<=r
            mulK = [mulK;linex(k) liney(k) ];
            dC{inter(i),move}=[dC{inter(i),move} norm(mulK(size(mulK,1),:)- device)];
        end
    end
        
        if i==2
            intersections = [];
            x1  = locDev(1,inter(1));
            y1  = locDev(2,inter(1));        
            x2  = locDev(1,inter(2));
            y2  = locDev(2,inter(2));

            for c2 = 1:size(dC{inter(2),move},2)
                r2 = dC{inter(2),move}(c2);
                for c1 = 1:size(dC{inter(1),move},2)
                    r1 = dC{inter(1),move}(c1);            
                    [xout,yout] = circcirc(x1,y1,r1,x2,y2,r2);
                    for s=1:numel(xout)
                        if norm([xout(s) yout(s)] - [x_c_Ti(move) y_c_Ti(move)]) <= r
                            intersections = [intersections; xout(s) yout(s)];
                        end
                    end
                end
            end
            errors=zeros(size(intersections,1),2);
        elseif i>=3
            x3p  = locDev(1,inter(i));
            y3p  = locDev(2,inter(i));     
            for c3p = 1:size(dC{inter(i),move},2)
                r3p = dC{inter(i),move}(c3p); 
                for s=1:size(intersections,1)
                   interToDev = norm([x3p y3p]-[intersections(s,1) intersections(s,2)]); 
                   if abs(interToDev-r3p)<lambda/5
                        errors(s,1)= errors(s,1)+1;
                        errors(s,2)= errors(s,2)+ abs(interToDev-r3p);
                   else
                        %error(s,i)=0;
                   end
                end
            end
        end

end

m=max(errors);
minE=10^3;
for er=1:size(errors,1)
    if m(1,1)==errors(er,1) && errors(er,2)<minE
        minE=errors(er,2);
        numE=er;
    end
end   
x_c_T(move)=intersections(numE,1);
y_c_T(move)=intersections(numE,2);


polyin=[];
if size(closeToDevArray,1)>2
    polyin = polyshape(closeToDevArray(:,1),closeToDevArray(:,2));
    [centrx,centry] = centroid(polyin);
end

% x_c_T(move)= 2.65 %centrx%mulK(2,1);
% y_c_T(move)= 2%centry%mulK(2,2);


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
if ~isempty(polyin)
    plot(centrx,centry, 'm*','LineWidth',10);
end  

%plot the devices positions
plot(locDev(1,:), locDev(2,:), 'g.');
plot(locDev(1,inter(1)), locDev(2,inter(1)), 'r.');
plot(locDev(1,inter(2)), locDev(2,inter(2)), 'y.');

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
% colors=["w-","g-","r-","m-","b-","r-","m-","b-","k-"];
% for k=1:numel(inter)
%     for i=1:numel(dC{inter(k),move})
%         th = 0:step:2*pi+step;
%         xunit = dC{inter(k),move}(i) * cos(th) + locDev(1,inter(k));
%         yunit = dC{inter(k),move}(i) * sin(th) + locDev(2,inter(k));
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



