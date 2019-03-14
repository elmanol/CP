    
for i=1:90
    xx(:,i) = every_iter_power1(:,i);% Create Data
    SEM = std(xx(:,i))/sqrt(length(xx(:,i)));    % Standard Error
    ts = tinv([0.025  0.975],length(xx(:,i))-1);    % T-Score
    CI_Discrete1(:,i) = mean(xx(:,i)) + ts*SEM;    % Confidence Intervals

    xx(:,i) = every_iter_power2(:,i);% Create Data
    SEM = std(xx(:,i))/sqrt(length(xx(:,i)));    % Standard Error
    ts = tinv([0.025  0.975],length(xx(:,i))-1);    % T-Score
    CI_Discrete2(:,i) = mean(xx(:,i)) + ts*SEM;    % Confidence Intervals

    xx(:,i) = every_iter_power3(:,i);% Create Data
    SEM = std(xx(:,i))/sqrt(length(xx(:,i)));    % Standard Error
    ts = tinv([0.025  0.975],length(xx(:,i))-1);    % T-Score
    CI_Discrete3(:,i) = mean(xx(:,i)) + ts*SEM;    % Confidence Intervals

    
    xx(:,i) = every_iter_power4(:,i);% Create Data
    SEM = std(xx(:,i))/sqrt(length(xx(:,i)));    % Standard Error
    ts = tinv([0.025  0.975],length(xx(:,i))-1);    % T-Score
    CI_Discrete4(:,i) = mean(xx(:,i)) + ts*SEM;    % Confidence Intervals
end




hold on;
boxplot(CI_Discrete1, 'colors','g','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_Discrete2, 'colors','r','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_Discrete3,'colors','k','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_Discrete4, 'colors','b','MedianStyle','target')
set(gca,'XTickLabel',{' '})


%hold off;
xlabel('Time (Rounds)')
ylabel('Cumulative Power (Watts)')
set(gca,'XTickLabel',{' '})

% %saveas(gcf,'images/ci.eps','eps');
% print('images/ci','-depsc','-r0')
% saveas(gcf,'images/ci.png','png');
hold off;