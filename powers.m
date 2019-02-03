function [Pt, Gt, Gr, lamda, k,P_Transfered] = powers( x_coordination_Trans, y_coordination_Trans,stopx,stopy,step)

%search space-moning the unique mote
startx = 0;
%stopx = 5;
stepx = step;

starty = 0;
%stopy = 5;
stepy = step;

%flag for distance from ET to ER
flag=0;

%initialize parameters
Pt = 2;
Gt = 2;
Gr = 1;
lamda = 0.3;
k= 2*pi/lamda;



number_of_Trans = numel(x_coordination_Trans);

%different position of node
countx=0;
 max_value = zeros (int8((stopx-startx)/stepx)+1,1);
 P_Transfered = zeros (int8((stopx-startx)/stepx)+1,int8((stopy-starty)/stepy)+1);
for counterx=startx:stepx:stopx
    countx = countx+1;
    county=0;
    for countery=starty:stepy:stopy
    % for counter=0:0.0001:1
    %     a=2*pi*counter;
    %     r=0.003;
    %     x_coordinations_Receiv = r*cos(a);
    %     y_coordinations_Receiv = r*sin(a);
    county = county+1;
    

        x_coordinations_Receiv = counterx;
        y_coordinations_Receiv = countery;
        
        
        % create costant row vector P
        P = zeros(1,number_of_Trans);
        for i = 1:number_of_Trans
            % find the distances from receiver to every transmitter
            P(1,i) = Pt*Gt*Gr*(lamda/(4*pi))^2;
        end
        
        % create inverse vector R
        R = zeros(number_of_Trans,1);
        distance = zeros(number_of_Trans,1);
        for i = 1:number_of_Trans
            % find the distances from receiver to every transmitter
            distance(i) = sqrt((x_coordinations_Receiv - x_coordination_Trans(i))^2 + (y_coordinations_Receiv - y_coordination_Trans(i))^2);
            R(i) = 1/distance(i);
            if (distance(i)<lamda) flag=1; break; end
        end
        
        %not too close to charger
        if(flag) flag=0; continue; end
        
        
        % create corr matrix
        corr = zeros(number_of_Trans,number_of_Trans);
        for i = 1:number_of_Trans
            for j = 1:number_of_Trans
                
                if(i~=j)
                    corr(i,j)= cos(k*(distance(i)-distance(j)))/distance(i);      
                else 
                    corr(i,j)= R(i); 
                end              
            end
        end
        
        %the total received energy P_T at receiver node
%         P_Transfered(iter)= P*corr*R;
        P_Transfered(countx,county)= P*corr*R;
        
        %countery for
        end
        %counterx for
end
      
% [max_val,max_grammh]=max(abs(P_Transfered));
% [max_val,max_sthlh]=max(max_val);
% P_Transfered(max_grammh(max_sthlh),max_sthlh);



r=0.25;
x=3;
y=3;

maxP=-1;
k=5;
c=1;

%max global radiation
%    [min_val,min_index]=min(abs(P_Transfered));
%     max_value(count) = val;
%max_value = val;
    % abs(P_Transfered(count))
% end


%  
% while c<=k
%     x=round(rand(1)*stopx*20,0)/20
%     y=round(rand(1)*stopy*20,0)/20
%     for i=(x-r):stepx:(x+r)
%       for j=(y-r):stepy:(y+r)
%           if dist([i,j],[x,y])<=r+0.01
%             if P_Transfered(int16(i/0.05+1), int16(j/0.05+1))>maxP
%                 maxP=P_Transfered(int16(i/0.05+1), int16(j/0.05+1));
%                 posMx=i;
%                 posMy=j;
%             end
%           end
%       end
%     end
%     maxP=-1;
%     c=c+1;
%     hold on
%     th = 0:pi/50:2*pi;
%     xunit = r * cos(th) + x;
%     yunit = r * sin(th) + y;
%     h = plot(xunit, yunit,"g");
%     hold on
%     plot(posMx,posMy,"r*")
%     hold on
% end

 
% hold on
%imagesc(startx:stepx:stopx,starty:stepy:stopy,abs(P_Transfered))

%surf(startx:stepx:stopx,starty:stepy:stopy,abs(P_Transfered(1:countx,1:county)));

%contourf(abs(P_Transfered),10);
%count=1;
%for counter=0:0.0001:1
%    a=2*pi*counter;
%    r=0.003;
%    x_coordinations_Receiv(count) = r*cos(a);
%    y_coordinations_Receiv(count) = r*sin(a);
%    count = count+1;
%end

%plot(0:0.0001:1,P_Transfered);

end