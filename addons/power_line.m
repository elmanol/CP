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
x_c_T = [ 0 4];
y_c_T = [ 0 0];


%devices positions
 locDevx=[3.5 3];
 locDevy=[0 0];
 
[x_c_T,y_c_T, x_c_R,y_c_R,distance] = deployment(5, 10, stopx,lambda)
 total_power_received = total_power( x_c_T,x_c_R,distance,lambda)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[x_c_T(1)];
[~, ~, ~, ~, ~,P_Transfered] = powers( x_c_T, y_c_T,stopx,stopy,step);
y=[P_Transfered(locDevx(1)/0.05+1,0+1)+P_Transfered(locDevx(2)/0.05+1,0+1)];
y1=P_Transfered(locDevx(1)/0.05+1,0+1);
y2=P_Transfered(locDevx(2)/0.05+1,0+1);
for i=0:0.01:2.5
	
    x_c_T(1)=x_c_T(1)+0.01;
    x=[x x_c_T(1)];
    [~, ~, ~, ~, ~,P_Transfered] = powers( x_c_T, y_c_T,stopx,stopy,step);
    y=[y P_Transfered(locDevx(1)/0.05+1,0+1)+P_Transfered(locDevx(2)/0.05+1,0+1)];
    y1=[y1 P_Transfered(locDevx(1)/0.05+1,0+1)];
    y2=[y2 P_Transfered(locDevx(2)/0.05+1,0+1)];
    
end
plot(x,y,"g*-",x,y1,"r-",x,y2,"c*")
