function [x_c_T,y_c_T, x_c_R,y_c_R,distance] = deployment(Trans_num, Receiv_num, plane_size,lambda)

x_c_R = 0+(plane_size-0).*rand(Receiv_num,1);
y_c_R = 0+(plane_size-0).*rand(Receiv_num,1);


x_c_R = 0+(plane_size-0).*rand(1,1);
y_c_R = 0+(plane_size-0).*rand(1,1);

ii=2;
while(true)
    if(Receiv_num==1)
        break;
    end
    flag=1;
    x_coordination_Receiv_temp   = 0+(plane_size-0).*rand(1,1);
    y_coordination_Receiv_temp   = 0+(plane_size-0).*rand(1,1);
        for iii=1:length(x_c_R)
         % find the distances from receiver to every deployed receiver
         distanceR2R = sqrt((x_c_R(iii) - x_coordination_Receiv_temp)^2 + (y_c_R(iii) - y_coordination_Receiv_temp)^2);
         if (distanceR2R < lambda)
             flag=0;
             break;
         end
     end
     if (flag==1)
         x_c_R(ii) = x_coordination_Receiv_temp;
         y_c_R(ii) = y_coordination_Receiv_temp;
         ii=ii+1;
         if(ii==Receiv_num+1)
             break;
         end
    end
end







number_of_Receiv = Receiv_num;
i=1;

while(true)
    flag=1;
    x_coordination_Trans_temp   = 0+(plane_size-0).*rand(1,1);
    y_coordination_Trans_temp   = 0+(plane_size-0).*rand(1,1);
    for ii = 1:number_of_Receiv
        % find the distances from receiver to every transmitter
        distance(i,ii) = sqrt((x_c_R(ii) - x_coordination_Trans_temp)^2 + (y_c_R(ii) - y_coordination_Trans_temp)^2);
        if (distance(i,ii)<lambda)
            flag=0;
            break;
        end
    end
    if (flag==1)
        x_c_T(i) = x_coordination_Trans_temp;
        y_c_T(i) = y_coordination_Trans_temp;
        i=i+1;
        if(i==Trans_num+1)
            break;
        end
    end
end



% figure(11)
% plot(x_c_T,y_c_T,'ok',x_c_R,y_c_R,'*r')
% xlabel('x(m)')
% ylabel('y(m)')
% legend('Chargers','Nodes','Location','northoutside','Orientation','horizontal')
% saveas(gcf,'images/deployment.eps','eps');
% saveas(gcf,'images/deployment.png','png');
% %saveas(gcf,'images1/deployment.eps','eps');
% %saveas(gcf,'images1/deployment.png','png');



