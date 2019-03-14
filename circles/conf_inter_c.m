%initial_one sum_power_centroids_r sum_power_multi_one_r sum_power_multi_one sum(line_power_r) sum_power_many_power_one


%     xx = initial_one';% Create Data
%     SEM = std(xx)/sqrt(length(xx));    % Standard Error
%     ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
%     CI_Discrete1(1,:)  = mean(xx) + ts*SEM;    % Confidence Intervals

    xx1 = power_centroids';% Create Data
    SEM = std(xx1)/sqrt(length(xx1));    % Standard Error
    ts = tinv([0.025  0.975],length(xx1)-1);    % T-Score
    CI_Discrete1(:,1) = mean(xx1) + ts*SEM   ;    % Confidence Intervals

    xx2 = power_multi_r';% Create Data
    SEM = std(xx2)/sqrt(length(xx2));    % Standard Error
    ts = tinv([0.025  0.975],length(xx2)-1);    % T-Score
    CI_Discrete1(:,2)  = mean(xx2) + ts*SEM;    % Confidence Intervals
    
    xx2 = power_multi';% Create Data
    SEM = std(xx2)/sqrt(length(xx2));    % Standard Error
    ts = tinv([0.025  0.975],length(xx2)-1);    % T-Score
    CI_Discrete1(:,3)  = mean(xx2) + ts*SEM;    % Confidence Intervals    

    
    xx3 = line_powert';% Create Data
    SEM = std(xx3)/sqrt(length(xx3));    % Standard Error
    ts = tinv([0.025  0.975],length(xx3)-1);    % T-Score
    CI_Discrete1(:,4)  = mean(xx3) + ts*SEM;    % Confidence Intervals

    xx4 = power_many_power';% Create Data
    SEM = std(xx4)/sqrt(length(xx4));    % Standard Error
    ts = tinv([0.025  0.975],length(xx4)-1);    % T-Score
    CI_Discrete1(:,5) = mean(xx4) + ts*SEM;    % Confidence Intervals    


hold on;


colors=['r'; 'b';'g';'m'];
labels = {'Appr. C' 'Appr. B' 'Appr. A' 'Alg. 1' 'Appr. D'};
    
hold on;
h=boxplot(CI_Discrete1,labels,'PlotStyle','traditional','MedianStyle','target', 'colors',colors);


% hLegend = legend(findall(gca,'Tag','Box'), {'Naive','MinDRD','Surf','Graph'},'FontSize', 10);
%        hChildren = findall(get(hLegend,'Children'), 'Type','Line');
       % Set the horizontal lines to the right colors

xlabel('Algorithm','FontSize', 13); % x-axis label
ylabel('Cumulative Power (Watts)','FontSize', 13); % y-axis label

  set(gca,'fontsize',15);








% % boxplot(CI_Discrete1, 'colors','g','MedianStyle','target')
% % set(gca,'XTickLabel',{' '})
% boxplot(CI_Discrete2', 'colors','r','MedianStyle','target')
% %set(gca,'XTickLabel',{' '})
% boxplot(CI_Discrete3','colors','k','MedianStyle','target')
% %set(gca,'XTickLabel',{' '})
% boxplot(CI_Discrete4', 'colors','b','MedianStyle','target')
% %set(gca,'XTickLabel',{' '})
% boxplot(CI_Discrete5', 'colors','y','MedianStyle','target')
% %set(gca,'XTickLabel',{' '})
% 
% set(gca,'XTickLabel',{'CI2','CI3','CI3','CI4'})
% xlabel('Time (Rounds)')
% ylabel('Cumulative Power (Watts)')
% % set(gca,'XTickLabel',{' '})
% 
% % %saveas(gcf,'images/ci.eps','eps');
% % print('images/ci','-depsc','-r0')
% % saveas(gcf,'images/ci.png','png');
% hold off;