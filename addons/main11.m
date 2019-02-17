% oti eixe to main9 sin thn evresi twn opt gwniwn me cholesky

clear
% load a.mat

%deployment
Trans_num  = 7;
Receiv_num = 100;
plane_size = 10;
akriveia = 1;
differentChargerSequences = 100;

%running rounds
time=90;%mexri wra de to xrisimopoiw

%system configuration
Pt=2;
Gt=2;
lambda=0.3;
k=2*pi/lambda;
Gr=1;
A=Gr*(lambda^2/(4*pi));
Zo=119.9169832*pi;
%poses arxikes faseis sto 0
phi_o = zeros(1,Trans_num);
initialPhaseMatrixForDifferentChargerSequence = zeros(differentChargerSequences,Trans_num);
totalPowerForEachRound = zeros(differentChargerSequences,time);
totalPowerForEachRoundRangeDepending = zeros(4,time,differentChargerSequences);
powerReceivedIndividuallyPahseShifted = zeros(differentChargerSequences,Receiv_num);
%edw adi gia 2 thelei Trans_num
dianisma_me_gwnies_dianismatwn = zeros(2,1);%vazw tis gwnies twn dianismatwn se ena node apo tous m chargers
metro_matrix = zeros(Trans_num,Receiv_num);
akoloutheiaChargers = zeros(differentChargerSequences,time);

display('start deployment')

%choose_trasmitter_id = randi([1,Trans_num], 1, time);
transm_initial_state= randi([0,1], 1, Trans_num);
[x_coordination_Trans,y_coordination_Trans,x_coordinations_Receiv,y_coordinations_Receiv,distance] = deployment(Trans_num,...
   Receiv_num, plane_size,lambda);
display('finish deployment')
%##########################################################################
%kalw to kwdika tou proigoumenou paper me 0-1 charger configuration
%##########################################################################
% [operRangePowerOverTime] = final1new(time,Trans_num,x_coordination_Trans,y_coordination_Trans,x_coordinations_Receiv,y_coordinations_Receiv,distance);

%load wrong.mat
% load a.mat

figure(1)
plot(x_coordination_Trans,y_coordination_Trans,'ok',x_coordinations_Receiv,y_coordinations_Receiv,'*r')
xlabel('x(m)')
ylabel('y(m)')
legend('Chargers','Nodes','Location','northoutside','Orientation','horizontal')
saveas(gcf,'images/deployment.eps','eps');
saveas(gcf,'images/deployment.png','png');



%kathe grammi sto distance matrix einai kai diaforetikos transmitter
%kathe stili einai oi diaforetikoi receivers

comm_range_list = [-1 0.4 1];
range=1;

%total power received with nodes' initial phase = 0
[staticPhase,nodesIndividualStaticPower] = powerWithSStaticPower(Trans_num,distance,Receiv_num,k,phi_o,time,Pt,Gt,Zo,A);
tempVarStaticPhaseIndividualPowerSorted = sort(nodesIndividualStaticPower);
staticPhase
sum(nodesIndividualStaticPower)
%mexri wras de to xrisimopoiw
%###################################################################################################################
%[configuration,nodes_in_range,electricField_quadratic] = myquadraticN_distributedAlg_v5 (x_coordination_Trans,...
%    y_coordination_Trans,x_coordinations_Receiv,y_coordinations_Receiv,...
%    choose_trasmitter_id,transm_initial_state, comm_range_list(range) ,time,plane_size,lambda,Pt,Gt,Gr,Zo); 
%###################################################################################################################

%einai mallon logo geitonias
number_of_Receiv=length(x_coordinations_Receiv);
number_of_Trans= length(x_coordination_Trans);

u=zeros(number_of_Trans,number_of_Receiv);
S=zeros(number_of_Trans,number_of_Receiv);


bima = 360/1;
phi_oo = ones(1,bima);
phi_o = ones(1,Trans_num);

uu=zeros(number_of_Trans,number_of_Receiv,bima);
SS=zeros(number_of_Trans,number_of_Receiv);
for i=1:number_of_Trans
    
    for ii=1:number_of_Receiv
        for iii=1:bima
            phi_oo(iii) = iii;%*2*pi/360;
            SS(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);
            uu(i,ii,iii)=sqrt(Zo*SS(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_oo(1,iii))));%));
        end
    end
end
uuu = uu;
testEnergyMatrix = uu;
%choose_trasmitter_id = randi([1,Trans_num], 1, time);
round=1;

% [phi_o] = findInitialPhaseChargerNoSym(bima,uu,Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
% phi_o


%run the simulation for 20 different chargers sequence
for counterSequence = 1:differentChargerSequences
    counterSequence
    choose_trasmitter_id = randi([1,Trans_num], 1, time);
    choose_trasmitter_id_matrix(counterSequence,:) = choose_trasmitter_id;
%%%  choose_trasmitter_id = akoloutheiaChargers(11,:);
    %choose_trasmitter_id = repmat(randperm(Trans_num),1,time/differentChargerSequences);
    phi_o = ones(1,Trans_num);
    for round = 1:time
        counterSequence;
        round;
        %sinarthsh na m epistrepsei arxiki fasi tou charger pou tou lew
        %tic;
%       [phi_o,sinarthsh] = findInitialPhaseCharger(Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
        [phi_o] = findInitialPhaseChargerNoSym(bima,uu,Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
        %toc;
        phi_o;
        for i=1:number_of_Trans
            for ii=1:number_of_Receiv
                S(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                u(i,ii)=sqrt(Zo*S(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_o(1,i))));%*2*pi/bima));
            end
        end

        %a = metro_apo_charger _a
        %b = metro_apo_charger _b
        %theta_a = gwnia dianismatos apo charger a sto node 1
        %theta_b = gwnia dianismatos apo charger b sto node 1

        % for i=1:number_of_Trans
        %     for ii=1:number_of_Receiv
        %         metro_matrix(i,ii) = sqrt( real(u(i,ii))^2 + imag(u(i,ii))^2 );
        %         %gwnia_matrix(i,ii) = atand(imag(u(i,ii))/real(u(i,ii)))  
        %         gwnia_matrix(i,ii) = atan2(imag(u(i,ii)),real(u(i,ii)))*57.2958;
        %     end
        % end

        % a = sqrt( real(u(1,1))^2 + imag(u(1,1))^2 );
        % b = sqrt( real(u(2,1))^2 + imag(u(2,1))^2 );
        % theta_a = atand(imag(u(1,1))/real(u(1,1)));
        % theta_b = atand(imag(u(2,1))/real(u(2,1)));
        % dianisma_me_gwnies_dianismatwn(1) = theta_a ;
        % dianisma_me_gwnies_dianismatwn(2) = theta_b ;
        % aaaaa=max(dianisma_me_gwnies_dianismatwn');
        %disp('fffff')


        for ii=1:number_of_Receiv
        %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
        %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
        %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

            superposition_of_elect_fields(1,ii)=sum(u(:,ii));
            phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
            superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
            superposition_phase(:,ii)=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

            %density of the total transferred power
            ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

            %total power received at point of interest
            total_power_received(1,ii)=ST(1,ii)*A;
        end
        total_power_received_different_phases(1) = sum(total_power_received(1,:));
        totalPowerForEachRound(counterSequence,round) = total_power_received_different_phases(1);
    end
    disp('fffff')
    phi_o2 = phi_o
    kaliteresFaseis(counterSequence,:) = phi_o; 
    total_power_received_different_phases(1)
    powerReceivedIndividuallyPahseShifted(counterSequence,:) = total_power_received(1,:);
    %akoloutheiaChargers(counterSequence,:) = choose_trasmitter_id;
    initialPhaseMatrixForDifferentChargerSequence(counterSequence,:) = phi_o;
    phi_o = zeros(1,Trans_num);
    
end

% na brw tin akoloutheia me to kalo=megalitero power gia na thn balw sto
% plot pou prosegizei to opt
[Mvalue,IIndex] = max(totalPowerForEachRound(:,time))

[a,b]= powerWithSStaticPower(Trans_num,distance,Receiv_num,k,phi_o2,time,Pt,Gt,Zo,A);

figure(2)
plot(1:time,totalPowerForEachRound(1,:),'o',1:time,staticPhase*ones(1,time))



hold on;
for cnt=2:differentChargerSequences
    plot(1:time,totalPowerForEachRound(cnt,:),'o')
end
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
hold off;


%######################################################################################
%###########################   power distance dependance  #############################
%######################################################################################


% figure(2)
% hold on
% yyaxis left
% plot(1:Receiv_num,sort(nodesIndividualStaticPower),'-',1:Receiv_num,sort(powerReceivedIndividuallyPahseShifted(1,:)),'--')
% xlabel('nodes')
% ylabel('power')
% legend('individual nodes received power static initial plase','total power fixed phases','Location','northoutside','Orientation','horizontal')
% hold off




minVector = min(distance);
sortedMinVector = sort(nodesIndividualStaticPower);
for i=1:Receiv_num
    index(i) = find(nodesIndividualStaticPower==sortedMinVector(i));
end
index;


hold on;
yyaxis right
plot(1:Receiv_num,minVector(index(:)))
hold off
    

sumVector = sum(distance);
sortedMinVector = sort(nodesIndividualStaticPower);
for i=1:Receiv_num
    index(i) = find(nodesIndividualStaticPower==sortedMinVector(i));
end
index;



%#######################################################################
%####################  Communicatio Range   ############################
%#######################################################################
for counterSequence = 1:differentChargerSequences
    choose_trasmitter_id_matrix(counterSequence,:);

    rangList=[2,3,4];
    [a,b]=size(rangList);
    rangListChar=['1','1.5','2','2.5'];

    for rangeCounter=1:b
        phi_o = ones(1,Trans_num);
        for round=1:time

            [phi_o]= findInitialPhaseChargerRangeDependingNoSym(bima,rangList(rangeCounter),Trans_num,distance,Receiv_num,k,phi_o,round ...
                ,choose_trasmitter_id_matrix(counterSequence,:),Pt,Gt,Zo,A,akriveia ...
                ,x_coordinations_Receiv,x_coordination_Trans,y_coordinations_Receiv,y_coordination_Trans);


            for i=1:number_of_Trans
                for ii=1:number_of_Receiv
                    S(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                    u(i,ii)=sqrt(Zo*S(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_o(1,i))));%*2*pi/bima));
                end
            end



            for ii=1:number_of_Receiv
                superposition_of_elect_fields(1,ii)=sum(u(:,ii));
                phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
                superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
                superposition_phase(:,ii)=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

                %density of the total transferred power
                ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

                %total power received at point of interest
                total_power_received(1,ii)=ST(1,ii)*A;
            end

            totalPowerForEachRoundRangeDepending(rangeCounter,round,counterSequence)=sum(total_power_received(1,:));



        end
        phi_o;
    %     figure(1)
    %     hold on;
    %     for cnt=1:4
    %         plot(1:time,totalPowerForEachRoundRangeDepending(cnt,:),'o')
    %     end
    %     legend('global','fixed phases','1','1.2','1.6','2','Location','northoutside','Orientation','horizontal')
    %     hold off;
    end
end

for time_counter = 1:time
    xx = totalPowerForEachRoundRangeDepending(1,time_counter,:);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_1(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

    xx = totalPowerForEachRoundRangeDepending(2,time_counter,:);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_2(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

    xx = totalPowerForEachRoundRangeDepending(3,time_counter,:);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_3(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

    xxx = totalPowerForEachRound(:,time_counter);% Create Data
    SEM = std(xxx)/sqrt(length(xxx));    % Standard Error
    ts = tinv([0.025  0.975],length(xxx)-1);    % T-Score
    CI_open(:,time_counter) = mean(xxx) + ts*SEM;    % Confidence Intervals
end


figure(5);
hold on;
boxplot(CI_open, 'colors','g','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_1, 'colors','r','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_2, 'colors','b','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_3,'colors','m','MedianStyle','target')

%hold off;
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
set(gca,'XTickLabel',{' '})
% %saveas(gcf,'images/ci.eps','eps');
% print('images/ci','-depsc','-r0')
% saveas(gcf,'images/ci.png','png');
hold off;


%#######################################################################
%####################  tha vrw ta opr fi apo to christoforos   #########
%#######################################################################
[Y,Z,optFIPower,pinakas_a] = christoforos(Trans_num,distance,Receiv_num,k,phi_o ...
                ,Pt,Gt,Zo,A ...
                ,x_coordinations_Receiv,x_coordination_Trans...
                ,y_coordinations_Receiv,y_coordination_Trans);
% ipologismow twn gwniwn me ton tropo pou eipe o xristoforos apo ta ac
Ac = zeros(1,2*Trans_num);
cholPHI = zeros(1,Trans_num);
[X] = chol(Y);
Ac(1:2*Trans_num) = X(1,:).';
    for counter=1:Trans_num
        if (Ac(counter)<=0)
            if (Ac(1,Trans_num+counter)<=0)
                cholPHI(counter) = 225;
            else
                cholPHI(counter) = 315;
            end
        else
            if (Ac(1,Trans_num+counter)<0)
                cholPHI(counter) = 135;
            else
                cholPHI(counter) = 45;
            end
        end
    end

[OptCholPhase,nodesIndividualOptCholPower] = powerWithSStaticPower(Trans_num,...
distance,Receiv_num,k,cholPHI,time,Pt,Gt,Zo,A);
OptCholPhase


            % ipologismos total opt power me exadlitiki anazitisi
Z;
optFIPower
% % o pinakas_a epistefei me vasi ta cosd
% % o M dinei ola ta pithana apo ta 45 135 225 315
%  M = permn([1:360],Trans_num); 
% for iiii=1:360^Trans_num
%     %optimalFI = pinakas_a(iiii,:);
%     optimalFI = M(iiii,:);
%     for i=1:number_of_Trans
%         for ii=1:number_of_Receiv
%             Sopt(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
%             uopt(i,ii)=sqrt(Zo*Sopt(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(optimalFI(i))));%*2*pi/360));
%         end
%     end
% 
%     for ii=1:number_of_Receiv
%     %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
%     %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
%     %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;
% 
%         superposition_of_elect_fields(1,ii)=sum(uopt(:,ii));
%         phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),...
%             real(superposition_of_elect_fields(1,ii)))*57.2958;
%         superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
%         superposition_phase(:,ii)=atan2(imag(uopt(:,ii)),real(uopt(:,ii)))*57.2958;
% 
%         %density of the total transferred power
%         ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;
% 
%         %total power received at point of interest
%         total_power_received(1,ii)=ST(1,ii)*A;
%     end
%     total_power_received_different_phasesopt(iiii) = sum(total_power_received(1,:));
% 
% end
% [Maximum_total_power_received_different_phasesopt index]=...
%     max(total_power_received_different_phasesopt)     
% %kaliteresFaseis = pinakas_a(index,:)
% kaliteresFaseis = M(index,:)



figure(3)
% plot(1:time,(optFIPower/100)*ones(1,time),'o'...
%     ,1:time,totalPowerForEachRound(differentChargerSequences,:),'-',1:time,staticPhase*ones(1,time))





plot(1:time,optFIPower*ones(1,time),'*k'...
    ,1:time,totalPowerForEachRound(differentChargerSequences,:),'--k',1:time,staticPhase*ones(1,time),'+k'...
    ,1:time,totalPowerForEachRoundRangeDepending(1,:,1),'-k', 1:time,totalPowerForEachRoundRangeDepending(2,:,1),'-.k' ...
    ,1:time,totalPowerForEachRoundRangeDepending(3,:,1),'.k')%,1:time,operRangePowerOverTime,'-')
legend('opt','glob','fix','2','3','4','Location','northoutside','Orientation','horizontal')
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')


%#######################################################################
%####################  4 phases   ############################
%#######################################################################
%opws o apo panw kwdikas gia to power, apla allaxa tis metavlites kai thn
%akriveia gia 4 faseis
akriveia = 90;
for counterSequence = 1:differentChargerSequences
    choose_trasmitter_id_matrix(counterSequence,:);
    
    phi_o4 = ones(1,Trans_num);
    for round = 1:time
        counterSequence;
        round;
        %sinarthsh na m epistrepsei arxiki fasi tou charger pou tou lew
        %tic;
%       [phi_o,sinarthsh] = findInitialPhaseCharger(Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
        [phi_o4] = findInitialPhaseChargerNoSym(bima,uu,Trans_num,distance,Receiv_num,k,phi_o4,round,... 
            choose_trasmitter_id_matrix(counterSequence,:),Pt,Gt,Zo,A,akriveia);
        %toc;
        phi_o4;
        for i=1:number_of_Trans
            for ii=1:number_of_Receiv
                S4(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                u4(i,ii)=sqrt(Zo*S4(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_o4(1,i))));%*2*pi/bima));
            end
        end

        %a = metro_apo_charger _a
        %b = metro_apo_charger _b
        %theta_a = gwnia dianismatos apo charger a sto node 1
        %theta_b = gwnia dianismatos apo charger b sto node 1

        % for i=1:number_of_Trans
        %     for ii=1:number_of_Receiv
        %         metro_matrix(i,ii) = sqrt( real(u(i,ii))^2 + imag(u(i,ii))^2 );
        %         %gwnia_matrix(i,ii) = atand(imag(u(i,ii))/real(u(i,ii)))  
        %         gwnia_matrix(i,ii) = atan2(imag(u(i,ii)),real(u(i,ii)))*57.2958;
        %     end
        % end

        % a = sqrt( real(u(1,1))^2 + imag(u(1,1))^2 );
        % b = sqrt( real(u(2,1))^2 + imag(u(2,1))^2 );
        % theta_a = atand(imag(u(1,1))/real(u(1,1)));
        % theta_b = atand(imag(u(2,1))/real(u(2,1)));
        % dianisma_me_gwnies_dianismatwn(1) = theta_a ;
        % dianisma_me_gwnies_dianismatwn(2) = theta_b ;
        % aaaaa=max(dianisma_me_gwnies_dianismatwn');
        %disp('fffff')


        for ii=1:number_of_Receiv
        %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
        %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
        %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

            superposition_of_elect_fields(1,ii)=sum(u4(:,ii));
            phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
            superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
            superposition_phase(:,ii)=atan2(imag(u4(:,ii)),real(u4(:,ii)))*57.2958;

            %density of the total transferred power
            ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

            %total power received at point of interest
            total_power_received(1,ii)=ST(1,ii)*A;
        end
        total_power_received_different_phases4 = sum(total_power_received(1,:));
        totalPowerForEachRound4(counterSequence,round) = total_power_received_different_phases4;
    end
    disp('fffff')
    phi_o4
    
end

%#######################################################################
%####################  2 phases   ############################
%#######################################################################
%opws o apo panw kwdikas gia to power, apla allaxa tis metavlites kai thn
%akriveia gia 2 faseis
akriveia = 180;
for counterSequence = 1:differentChargerSequences
    counterSequence
    
    phi_o22 = ones(1,Trans_num);
    for round = 1:time
        counterSequence;
        round;
        %sinarthsh na m epistrepsei arxiki fasi tou charger pou tou lew
        %tic;
%       [phi_o,sinarthsh] = findInitialPhaseCharger(Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
        [phi_o22] = findInitialPhaseChargerNoSym(bima,uu,Trans_num,distance,Receiv_num,k,phi_o22,round,...
            choose_trasmitter_id_matrix(counterSequence,:),Pt,Gt,Zo,A,akriveia);
        %toc;
        phi_o22;
        for i=1:number_of_Trans
            for ii=1:number_of_Receiv
                S2(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                u2(i,ii)=sqrt(Zo*S2(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_o22(1,i))));%*2*pi/bima));
            end
        end

        %a = metro_apo_charger _a
        %b = metro_apo_charger _b
        %theta_a = gwnia dianismatos apo charger a sto node 1
        %theta_b = gwnia dianismatos apo charger b sto node 1

        % for i=1:number_of_Trans
        %     for ii=1:number_of_Receiv
        %         metro_matrix(i,ii) = sqrt( real(u(i,ii))^2 + imag(u(i,ii))^2 );
        %         %gwnia_matrix(i,ii) = atand(imag(u(i,ii))/real(u(i,ii)))  
        %         gwnia_matrix(i,ii) = atan2(imag(u(i,ii)),real(u(i,ii)))*57.2958;
        %     end
        % end

        % a = sqrt( real(u(1,1))^2 + imag(u(1,1))^2 );
        % b = sqrt( real(u(2,1))^2 + imag(u(2,1))^2 );
        % theta_a = atand(imag(u(1,1))/real(u(1,1)));
        % theta_b = atand(imag(u(2,1))/real(u(2,1)));
        % dianisma_me_gwnies_dianismatwn(1) = theta_a ;
        % dianisma_me_gwnies_dianismatwn(2) = theta_b ;
        % aaaaa=max(dianisma_me_gwnies_dianismatwn');
        %disp('fffff')


        for ii=1:number_of_Receiv
        %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
        %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
        %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

            superposition_of_elect_fields(1,ii)=sum(u2(:,ii));
            phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
            superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
            superposition_phase(:,ii)=atan2(imag(u2(:,ii)),real(u2(:,ii)))*57.2958;

            %density of the total transferred power
            ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

            %total power received at point of interest
            total_power_received(1,ii)=ST(1,ii)*A;
        end
        total_power_received_different_phases2 = sum(total_power_received(1,:));
        totalPowerForEachRound2(counterSequence,round) = total_power_received_different_phases2;
    end
    disp('fffff')
    phi_o22
    
end

%#######################################################################
%####################  8 phases   ############################
%#######################################################################
%opws o apo panw kwdikas gia to power, apla allaxa tis metavlites kai thn
%akriveia gia 8 faseis
akriveia = 45;
for counterSequence = 1:differentChargerSequences
    counterSequence
    
    phi_o8 = ones(1,Trans_num);
    for round = 1:time
        counterSequence;
        round;
        %sinarthsh na m epistrepsei arxiki fasi tou charger pou tou lew
        %tic;
%       [phi_o,sinarthsh] = findInitialPhaseCharger(Trans_num,distance,Receiv_num,k,phi_o,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
        [phi_o8] = findInitialPhaseChargerNoSym(bima,uu,Trans_num,distance,Receiv_num,k,phi_o8,round,...
            choose_trasmitter_id_matrix(counterSequence,:),Pt,Gt,Zo,A,akriveia);
        %toc;
        phi_o8;
        for i=1:number_of_Trans
            for ii=1:number_of_Receiv
                S8(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                u8(i,ii)=sqrt(Zo*S8(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_o8(1,i))));%*2*pi/bima));
            end
        end

        %a = metro_apo_charger _a
        %b = metro_apo_charger _b
        %theta_a = gwnia dianismatos apo charger a sto node 1
        %theta_b = gwnia dianismatos apo charger b sto node 1

        % for i=1:number_of_Trans
        %     for ii=1:number_of_Receiv
        %         metro_matrix(i,ii) = sqrt( real(u(i,ii))^2 + imag(u(i,ii))^2 );
        %         %gwnia_matrix(i,ii) = atand(imag(u(i,ii))/real(u(i,ii)))  
        %         gwnia_matrix(i,ii) = atan2(imag(u(i,ii)),real(u(i,ii)))*57.2958;
        %     end
        % end

        % a = sqrt( real(u(1,1))^2 + imag(u(1,1))^2 );
        % b = sqrt( real(u(2,1))^2 + imag(u(2,1))^2 );
        % theta_a = atand(imag(u(1,1))/real(u(1,1)));
        % theta_b = atand(imag(u(2,1))/real(u(2,1)));
        % dianisma_me_gwnies_dianismatwn(1) = theta_a ;
        % dianisma_me_gwnies_dianismatwn(2) = theta_b ;
        % aaaaa=max(dianisma_me_gwnies_dianismatwn');
        %disp('fffff')


        for ii=1:number_of_Receiv
        %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
        %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
        %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

            superposition_of_elect_fields(1,ii)=sum(u8(:,ii));
            phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
            superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
            superposition_phase(:,ii)=atan2(imag(u8(:,ii)),real(u8(:,ii)))*57.2958;

            %density of the total transferred power
            ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

            %total power received at point of interest
            total_power_received(1,ii)=ST(1,ii)*A;
        end
        total_power_received_different_phases8 = sum(total_power_received(1,:));
        totalPowerForEachRound8(counterSequence,round) = total_power_received_different_phases8;
    end
    disp('fffff')
    phi_o8
    
end


for time_counter = 1:time
    xx = totalPowerForEachRound2(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_Discrete1(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

    xx = totalPowerForEachRound4(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_Discrete2(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

    xx = totalPowerForEachRound8(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_Discrete3(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

%     xxx = totalPowerForEachRound(:,time_counter);% Create Data
%     SEM = std(xxx)/sqrt(length(xxx));    % Standard Error
%     ts = tinv([0.025  0.975],length(xxx)-1);    % T-Score
%     CI_open(:,time_counter) = mean(xxx) + ts*SEM;    % Confidence Intervals
end


figure(6);
hold on;
boxplot(CI_Discrete1, 'colors','g','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(CI_Discrete2, 'colors','r','MedianStyle','target')
set(gca,'XTickLabel',{' '})
% boxplot(CI_2, 'colors','b','MedianStyle','target')
% set(gca,'XTickLabel',{' '})
boxplot(CI_Discrete3,'colors','m','MedianStyle','target')

%hold off;
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
set(gca,'XTickLabel',{' '})
% %saveas(gcf,'images/ci.eps','eps');
% print('images/ci','-depsc','-r0')
% saveas(gcf,'images/ci.png','png');
hold off;



% figure(3)
% hold on;
% plot(1:time,totalPowerForEachRound2(1,:))
% hold off;

figure(4)
plot(1:time,totalPowerForEachRound2(1,:),'-',1:time,totalPowerForEachRound4(1,:),'.',1:time,totalPowerForEachRound8(1,:),'-.')
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
legend('2','4','8','Location','northoutside','Orientation','horizontal')



%#######################################################################
%####################  brute force me 2 phases   ############################
%#######################################################################
akriveia =180;
phase_2_possible_configurations = permn(1:akriveia:bima,Trans_num);
[rows2,columns2] = size(phase_2_possible_configurations);
totalPowerForEachRound2OPT = zeros(rows2,1);
for posible_conf_counter = 1:rows2
    phi_o2OPT = phase_2_possible_configurations(posible_conf_counter,:);
    for i=1:number_of_Trans
        for ii=1:number_of_Receiv
            S2OPT(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
            u2OPT(i,ii)=sqrt(Zo*S2OPT(i,ii))*exp(-1i*(k*distance(i,ii)+phi_o2OPT(1,i)*2*pi/bima));
        end
    end

    for ii=1:number_of_Receiv
        superposition_of_elect_fields(1,ii)=sum(u2OPT(:,ii));
        phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
        superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
        superposition_phase(:,ii)=atan2(imag(u2OPT(:,ii)),real(u2OPT(:,ii)))*57.2958;

        %density of the total transferred power
        ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

        %total power received at point of interest
        total_power_received(1,ii)=ST(1,ii)*A;
    end
    total_power_received_different_phases2OPT = sum(total_power_received(1,:));
    totalPowerForEachRound2OPT(posible_conf_counter) = total_power_received_different_phases2OPT;
end
disp('fffff')
[OPT_2_phase_power_value OPT_2_phase_power_row] = max(totalPowerForEachRound2OPT);
phi_o2OPT = phase_2_possible_configurations(OPT_2_phase_power_row,:)


%#######################################################################
%####################  brute force me 4 phases   ############################
%#######################################################################
akriveia =90;
phase_4_possible_configurations = permn(1:akriveia:bima,Trans_num);
[rows4,columns4] = size(phase_4_possible_configurations);
totalPowerForEachRound4OPT = zeros(rows4,1);
for posible_conf_counter = 1:rows4
    phi_o4OPT = phase_4_possible_configurations(posible_conf_counter,:);
    for i=1:number_of_Trans
        for ii=1:number_of_Receiv
            S4OPT(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
            u4OPT(i,ii)=sqrt(Zo*S4OPT(i,ii))*exp(-1i*(k*distance(i,ii)+phi_o4OPT(1,i)*2*pi/bima));
        end
    end

    for ii=1:number_of_Receiv
        superposition_of_elect_fields(1,ii)=sum(u4OPT(:,ii));
        phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
        superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
        superposition_phase(:,ii)=atan2(imag(u4OPT(:,ii)),real(u4OPT(:,ii)))*57.2958;

        %density of the total transferred power
        ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

        %total power received at point of interest
        total_power_received(1,ii)=ST(1,ii)*A;
    end
    total_power_received_different_phases4OPT = sum(total_power_received(1,:));
    totalPowerForEachRound4OPT(posible_conf_counter) = total_power_received_different_phases4OPT;
end
disp('fffff')
[OPT_4_phase_power_value OPT_4_phase_power_row] = max(totalPowerForEachRound4OPT);
phi_o4OPT = phase_4_possible_configurations(OPT_4_phase_power_row,:)

%#######################################################################
%####################  brute force me 8 phases   ############################
%#######################################################################
akriveia =45;
phase_8_possible_configurations = permn(1:akriveia:bima,Trans_num);
[rows8,columns8] = size(phase_8_possible_configurations);
totalPowerForEachRound8OPT = zeros(rows8,1);
for posible_conf_counter = 1:rows8
    phi_o8OPT = phase_8_possible_configurations(posible_conf_counter,:);
    for i=1:number_of_Trans
        for ii=1:number_of_Receiv
            S8OPT(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
            u8OPT(i,ii)=sqrt(Zo*S8OPT(i,ii))*exp(-1i*(k*distance(i,ii)+phi_o8OPT(1,i)*2*pi/bima));
        end
    end

    for ii=1:number_of_Receiv
        superposition_of_elect_fields(1,ii)=sum(u8OPT(:,ii));
        phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
        superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
        superposition_phase(:,ii)=atan2(imag(u8OPT(:,ii)),real(u8OPT(:,ii)))*57.2958;

        %density of the total transferred power
        ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

        %total power received at point of interest
        total_power_received(1,ii)=ST(1,ii)*A;
    end
    total_power_received_different_phases8OPT = sum(total_power_received(1,:));
    totalPowerForEachRound8OPT(posible_conf_counter) = total_power_received_different_phases8OPT;
end
disp('fffff')
[OPT_8_phase_power_value OPT_8_phase_power_row] = max(totalPowerForEachRound8OPT);
phi_o8OPT = phase_8_possible_configurations(OPT_8_phase_power_row,:)

figure(4)
hold on;
plot(1:time,OPT_2_phase_power_value*ones(1,time),1:time,OPT_4_phase_power_value*ones(1,time),1:time...
,OPT_8_phase_power_value*ones(1,time))
hold off;




optFIPower


kaliteresFaseis;
Mvalue;
IIndex;
powerForGreedyAlgPlot = zeros(1,Trans_num+1);
% pinakasStavros = zeros(5,Trans_num);
% pinakasMeTisFaseis = zeros(5,Trans_num);
dianismaMeKaliteresFaseis = kaliteresFaseis(IIndex,:);
[GRPhase,nodesIndividualGRPower] = powerWithSStaticPower(Trans_num,...
    distance,Receiv_num,k,dianismaMeKaliteresFaseis,time,Pt,Gt,Zo,A);
powerForGreedyAlgPlot(1,1) = sum(nodesIndividualGRPower);
checkednodes = zeros(1,Trans_num);
for i=1:Trans_num
    pinakasStavros(i,:) = [0 90 180 270 360];
end
for cnt=1:5%oses kai oi times tou stavrou
    pinakasMeTisFaseis(:,cnt) = dianismaMeKaliteresFaseis;
    
end
for i=1:Trans_num
    X = abs(pinakasMeTisFaseis-pinakasStavros);
    for zzz=1:Trans_num
        if(checkednodes(1,zzz) == 1)
            X(zzz,:) = 1000;
        end
    end
    X;
    dianismaMeKaliteresFaseis;
    [Mas Ias] = min(min(X));
    %grammi kombos 
    %sthlh gwnia
    [Ias_row, Ias_col] = find(X==Mas)
    gwnia = Ias_col(1,1);
    if (Ias_col(1,1) == 1)
        gwnia = 5;
    end
    
    dianismaMeKaliteresFaseis(Ias_row(1,1)) = pinakasStavros(1,gwnia);
    checkednodes(1,Ias_row(1,1)) = 1;
    [GRPhase,nodesIndividualGRPower] = powerWithSStaticPower(Trans_num,...
    distance,Receiv_num,k,dianismaMeKaliteresFaseis,time,Pt,Gt,Zo,A);
powerForGreedyAlgPlot(1,i+1) = sum(nodesIndividualGRPower);
    
    for round = 1:time
        [phi_oGreedy] = findInitialPhaseChargerNoSymForGreedy(checkednodes,bima,uu,Trans_num,distance,Receiv_num,k,dianismaMeKaliteresFaseis,round,choose_trasmitter_id,Pt,Gt,Zo,A,akriveia);
            %toc;
            phi_oGreedy;
            for i=1:number_of_Trans
                for ii=1:number_of_Receiv
                    S(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
                    u(i,ii)=sqrt(Zo*S(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(phi_oGreedy(1,i))));%*2*pi/bima));
                end
            end

            %a = metro_apo_charger _a
            %b = metro_apo_charger _b
            %theta_a = gwnia dianismatos apo charger a sto node 1
            %theta_b = gwnia dianismatos apo charger b sto node 1

            % for i=1:number_of_Trans
            %     for ii=1:number_of_Receiv
            %         metro_matrix(i,ii) = sqrt( real(u(i,ii))^2 + imag(u(i,ii))^2 );
            %         %gwnia_matrix(i,ii) = atand(imag(u(i,ii))/real(u(i,ii)))  
            %         gwnia_matrix(i,ii) = atan2(imag(u(i,ii)),real(u(i,ii)))*57.2958;
            %     end
            % end

            % a = sqrt( real(u(1,1))^2 + imag(u(1,1))^2 );
            % b = sqrt( real(u(2,1))^2 + imag(u(2,1))^2 );
            % theta_a = atand(imag(u(1,1))/real(u(1,1)));
            % theta_b = atand(imag(u(2,1))/real(u(2,1)));
            % dianisma_me_gwnies_dianismatwn(1) = theta_a ;
            % dianisma_me_gwnies_dianismatwn(2) = theta_b ;
            % aaaaa=max(dianisma_me_gwnies_dianismatwn');
            %disp('fffff')


            for ii=1:number_of_Receiv
            %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
            %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
            %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

                superposition_of_elect_fields(1,ii)=sum(u(:,ii));
                phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),real(superposition_of_elect_fields(1,ii)))*57.2958;
                superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
                superposition_phase(:,ii)=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;

                %density of the total transferred power
                ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;

                %total power received at point of interest
                total_power_received(1,ii)=ST(1,ii)*A;
            end
            total_power_received_different_phases_greedy(1) = sum(total_power_received(1,:));
            totalPowerForEachRoundGreedy(counterSequence,round) = ...
                total_power_received_different_phases_greedy(1);
            
            for zxzx = 1:number_of_Trans
                if (checkednodes(1,zxzx) == 0)
                    dianismaMeKaliteresFaseis(zxzx) = phi_oGreedy(zxzx);
                end
            end
    end
    
    dianismaMeKaliteresFaseis;
end    
    
dianismaMeKaliteresFaseis;

for zz= number_of_Trans+1:time
    powerForGreedyAlgPlot(zz) = powerForGreedyAlgPlot(number_of_Trans+1);
end


figure(4)
hold on;
plot(1:time,powerForGreedyAlgPlot,'+')
legend('heuristic','Location','northoutside','Orientation','horizontal')
hold off;
powerForGreedyAlgPlot


% %#######################################################################
% %####################  tha vrw ta opt fi apo to christoforos2   #########
% %#######################################################################
% %  einai to idio me to christoforos apla edw dinw san eisodo tis
% %  fixarimenes arxikes faseis
% %pinakas gwniwn stavros
% pinakasStavros = zeros(5,Trans_num);
% for i=1:Trans_num
%     pinakasStavros(:,i) = [0;90;180;270;360];
% end
% fixarismenesArxikesFaseis = zeros(Trans_num,1);
% for trCounter = 1:Trans_num
%     [Z,optFIPower,pinakas_a] = christoforos2(Trans_num,distance,Receiv_num,k,phi_o ...
%                     ,Pt,Gt,Zo,A ...
%                     ,x_coordinations_Receiv,x_coordination_Trans...
%                     ,y_coordinations_Receiv,y_coordination_Trans,fixarismenesArxikesFaseis);
% 
%                 % ipologismos total opt power
% 
% 
%     for iiii=1:2^Trans_num
%         optimalFI = pinakas_a(iiii,:);
%         for i=1:number_of_Trans
%             for ii=1:number_of_Receiv
%                 Sopt(i,ii)=(Gt*Pt)/(4*pi*distance(i,ii)^2);    
%                 uopt(i,ii)=sqrt(Zo*Sopt(i,ii))*exp(-1i*(k*distance(i,ii)+degtorad(optimalFI(i))));%*2*pi/360));
%             end
%         end
% 
%         for ii=1:number_of_Receiv
%         %         superposition_of_elect_fields(ii)=sum(u(:,ii,360));
%         %         superposition_magnitude=abs(superposition_of_elect_fields(ii));
%         %         superposition_phase=atan2(imag(u(:,ii)),real(u(:,ii)))*57.2958;
% 
%             superposition_of_elect_fields(1,ii)=sum(uopt(:,ii));
%             phase_superposition(1,ii) = atan2(imag(superposition_of_elect_fields(1,ii)),...
%                 real(superposition_of_elect_fields(1,ii)))*57.2958;
%             superposition_magnitude(1,ii)=abs(superposition_of_elect_fields(1,ii));
%             superposition_phase(:,ii)=atan2(imag(uopt(:,ii)),real(uopt(:,ii)))*57.2958;
% 
%             %density of the total transferred power
%             ST(1,ii)=(superposition_magnitude(1,ii)^2)/Zo;
% 
%             %total power received at point of interest
%             total_power_received(1,ii)=ST(1,ii)*A;
%         end
%         total_power_received_different_phasesChrisTemp(iiii) = sum(total_power_received(1,:));
%         
%        
% 
%     end
%     [Maximum_total_power_received_different_phasesChrisTemp index]=...
%         max(total_power_received_different_phasesChrisTemp)     
%     kaliteresFaseis = pinakas_a(index,:)
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % gia na dw se kathe vima an gnomai kaliteros alg xristoforos
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     totalPowerAlgCristoforosKatheBima(trCounter)= ...
%         Maximum_total_power_received_different_phasesChrisTemp;
% 
%     %     for jj=1:Trans_num
%     pinakasKaliterwnFasewn = (kaliteresFaseis'*ones(1,5))';
% %     end
%     diaforesGwnies = abs(pinakasStavros-pinakasKaliterwnFasewn);
%     [ii,jj]=find(diaforesGwnies==min(min(abs(pinakasStavros-pinakasKaliterwnFasewn))));
%     fixarismenesArxikesFaseis(jj) = pinakasStavros(ii,1)
%     
%     
% end
% 
% figure(14)
% plot( totalPowerAlgCristoforosKatheBima)
% xlabel('time')
% ylabel('power')
% legend('Chris','Location','northoutside','Orientation','horizontal')
% 
% fixarismenesArxikesFaseis
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % evresi twn 4 diakritwn fasewn (stavros), vriskw elaxisto
% % gwniwn apo 0 90 180 270 kai exw kai agnwsto x
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a = kaliteresFaseis;        
% b = a'*ones(1,4);        
% c = [0 90 180 270];        
% d = c'*ones(1,Trans_num) ;      
% e = abs(b'-d);        
%         
% for i=1:Trans_num
%     arithmosFixarismenouStavrou(i) = min(e(:,i));
%     [iii(i),jjj(i)] = find(e(:,i)==min(e(:,i)));
% end
% 
% 
% % fun = @(x) abs(c(iii(1))+x-kaliteresFaseis(1));
% % for i=12:Trans_num
% %     fun = @(x) abs(mod(c(iii(1))+x-kaliteresFaseis(1),360))...
% %         +abs(mod(c(iii(2))+x-kaliteresFaseis(2),360))...
% %         +abs(mod(c(iii(3))+x-kaliteresFaseis(3),360))...
% %         +abs(mod(c(iii(4))+x-kaliteresFaseis(4),360))...
% %         +abs(mod(c(iii(5))+x-kaliteresFaseis(5),360))...
% %         +abs(mod(c(iii(6))+x-kaliteresFaseis(6),360))...
% %         +abs(mod(c(iii(7))+x-kaliteresFaseis(7),360))...
% %         +abs(mod(c(iii(8))+x-kaliteresFaseis(8),360))...
% %         +abs(mod(c(iii(9))+x-kaliteresFaseis(9),360))...
% %         +abs(mod(c(iii(10))+x-kaliteresFaseis(10),360));
% %      
% % % end
% % 
% % x0 = [0];
% % x = fminsearch(fun,x0)          
% % y = feval(fun,x)        
%         
%         
%         
%         
        
        
        
        
