clear
%charger 1
a=0;
%charger 2
b=2;
lamda=1;
initial_phase = pi+0.2;
%intervalamda displayed
x = [a+lamda-0.7:0.001:b-lamda+0.7];
for counter=1:length(x)
    y(counter)=( cos((2*pi*(-a - b + 2*x(counter)))/lamda) + cos((2*pi*(-a - b + 2*x(counter)))/lamda))/((x(counter)-a)*(b-x(counter))) +1/(x(counter)-a)^2 + 1/(b-x(counter))^2;
    
end
hold on;
plot (x,(1/(4*pi)^2)*y,'b')

a=0;
%charger 2
b=2.5;
lamda=1;
initial_phase = pi+0.2;
%intervalamda displayed
x = [a+lamda-0.7:0.001:b-lamda+0.7];
for counter=1:length(x)
    ynew(counter)=( cos((2*pi*(-a - b + 2*x(counter)))/lamda) + cos((2*pi*(-a - b + 2*x(counter)))/lamda))/((x(counter)-a)*(b-x(counter))) +1/(x(counter)-a)^2 + 1/(b-x(counter))^2;
    
end


plot (x,(1/(4*pi)^2)*ynew,'c')
hold off;
ylabel('Power (watts)')
xlabel('Position (m)')
legend('initial positions','after readjastment')
% set(gca,'YTickLabel',{' '})
% set(gca,'Xtick',0+1/4:0.25:2-1/4,'XTickLabel',{'1/4', '2/4', '3/4','1', '5/4', '6/4' , '7/4' })
% %grid on
% %grid minor
% hold on
% 
% %saveas(gcf,'images/toyexample.eps','eps');
% print('figs/toyexample','-depsc','-r0')
% saveas(gcf,'figs/toyexample.png','png');
% 
% plot([2 0],[0 0],'-o',...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)
% 
% plot([5/4],[0],'-s',...    
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)
% 
% legend('vector both chargers ON','initial phase shifted', '(0,0) ON, (2,0) OFF','(0,0) ON, (2,0) OFF','scalar both chargers ON','charger','receiver')
% 
% 
% hold off
% syms initial_phase1
% 
% sinarthsh=( cos(((2*pi*(-a - b + 2*5/4))/lamda)+initial_phase1) + cos(((2*pi*(-a - b + 2*5/4))/lamda)+initial_phase1))/((5/4-a)*(b-5/4)) +1/(5/4-a)^2 + 1/(b-5/4)^2
% 
% prwtiParagwgos = diff(sinarthsh)
% eqn = prwtiParagwgos == 0;
% 
% solx = solve(eqn,initial_phase1)
% initial_phase1 = double(solx)
% 

