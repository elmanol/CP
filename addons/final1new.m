for iiiii=1:1
clear


lambda=0.3;
time=90;
Trans_num  = 15;
Receiv_num = 200;
plane_size = 10;
iterations = 20;

Gr=1;
A=Gr*(lambda^2/(4*pi));
Zo=119.9169832*pi;


%gia to k
start=5;
step=5;
stop=45;

%gia to x
startx=30;
stepx=30;
stopx=30;

%choose_trasmitter_id = randi([1,Trans_num], 1, time);
transm_initial_state= randi([0,1], 1, Trans_num);
[x_coordination_Trans,y_coordination_Trans,x_coordinations_Receiv,y_coordinations_Receiv,distance] = deployment(Trans_num, Receiv_num, plane_size,lambda);
comm_range_list = [-1 0.4 0.6 0.8 1 0.3 0.5 0.7 0.8 0.9 1.1 1.2 1.3 1.4 1.5];
m = {'-ok', '-*k', '-+k', '-.k', '-xk', '-^k', '-pk', '--ok'};
nn = {'k-x', 'k--', 'k-', 'k:o', 'k:*'};
zz = {'o', '+', 'x'};
range_list_names = {'open', '0.4', '0.6',  '1'};
%matrix same deployment different pick node sequence iterations x time
%needed for confidence interval 
Pi_matrix_open = zeros(iterations,time);
Pi_matrix_1 = zeros(iterations,time);
Pi_matrix_2 = zeros(iterations,time);
Pi_matrix_3 = zeros(iterations,time);
Pi_matrix_4 = zeros(iterations,time);

%gia tis epipleon embeleies
Pi_matrix_6 = zeros(iterations,time);
Pi_matrix_7 = zeros(iterations,time);
Pi_matrix_8 = zeros(iterations,time);
Pi_matrix_9 = zeros(iterations,time);
Pi_matrix_10 = zeros(iterations,time);
Pi_matrix_11 = zeros(iterations,time);
Pi_matrix_12 = zeros(iterations,time);
Pi_matrix_13 = zeros(iterations,time);
Pi_matrix_14 = zeros(iterations,time);
Pi_matrix_15 = zeros(iterations,time);

%gia ta minimata
message_vector_open=zeros(1,Trans_num);
message_vector2=zeros(1,Trans_num);
message_vector3=zeros(1,Trans_num);
message_vector5=zeros(1,Trans_num);

%gia ta epipleon minimata
message_vector6=zeros(1,Trans_num);
message_vector7=zeros(1,Trans_num);
message_vector8=zeros(1,Trans_num);
message_vector9=zeros(1,Trans_num);
message_vector10=zeros(1,Trans_num);
message_vector11=zeros(1,Trans_num);
message_vector12=zeros(1,Trans_num);
message_vector13=zeros(1,Trans_num);
message_vector14=zeros(1,Trans_num);
message_vector15=zeros(1,Trans_num);

%gia ta configuration analoga to communication range
configuration_open=zeros(1,Trans_num);
configuration_2=zeros(1,Trans_num);
configuration_3=zeros(1,Trans_num);
configuration_5=zeros(1,Trans_num);

for iter=1:iterations
    choose_trasmitter_id = randi([1,Trans_num], 1, time);
    for range=1:length(comm_range_list)
        [configuration,nodes_in_range,electricField_quadratic] = myquadraticN_distributedAlg_v5 (x_coordination_Trans,y_coordination_Trans,x_coordinations_Receiv,y_coordinations_Receiv,choose_trasmitter_id,transm_initial_state, comm_range_list(range) ,time,plane_size,lambda); 
        %open 1st item in comm_range_list
        if (range == 1)
            Pi_matrix_open(iter,1:time) = electricField_quadratic;
            message_vector_open(choose_trasmitter_id)=message_vector_open(choose_trasmitter_id)+nodes_in_range;
            configuration_open=configuration;
        %2nd item in comm_range_list
        elseif (range == 2)
            Pi_matrix_1(iter,1:time) = electricField_quadratic;
            message_vector2(choose_trasmitter_id)=message_vector2(choose_trasmitter_id)+nodes_in_range;
            configuration_2=configuration;
        elseif (range == 3)
            Pi_matrix_2(iter,1:time) = electricField_quadratic;
            message_vector3(choose_trasmitter_id)=message_vector3(choose_trasmitter_id)+nodes_in_range;
            configuration_3=configuration;
        elseif (range == 4)
            Pi_matrix_3(iter,1:time) = electricField_quadratic;
        elseif (range == 5)
            Pi_matrix_4(iter,1:time) = electricField_quadratic;
            message_vector5(choose_trasmitter_id)=message_vector5(choose_trasmitter_id)+nodes_in_range;
            configuration_5=configuration;
        elseif (range == 6)
            Pi_matrix_6(iter,1:time) = electricField_quadratic;
            message_vector6(choose_trasmitter_id)=message_vector6(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 7)
            Pi_matrix_7(iter,1:time) = electricField_quadratic;
            message_vector7(choose_trasmitter_id)=message_vector7(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 8)
            Pi_matrix_8(iter,1:time) = electricField_quadratic;
            message_vector8(choose_trasmitter_id)=message_vector8(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 9)
            Pi_matrix_9(iter,1:time) = electricField_quadratic;
            message_vector9(choose_trasmitter_id)=message_vector9(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 10)
            Pi_matrix_10(iter,1:time) = electricField_quadratic;
            message_vector10(choose_trasmitter_id)=message_vector10(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 11)
            Pi_matrix_11(iter,1:time) = electricField_quadratic;
            message_vector11(choose_trasmitter_id)=message_vector11(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 12)
            Pi_matrix_12(iter,1:time) = electricField_quadratic;
            message_vector12(choose_trasmitter_id)=message_vector12(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 13)
            Pi_matrix_13(iter,1:time) = electricField_quadratic;
            message_vector13(choose_trasmitter_id)=message_vector13(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 14)
            Pi_matrix_14(iter,1:time) = electricField_quadratic;
            message_vector14(choose_trasmitter_id)=message_vector14(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        elseif (range == 15)
            Pi_matrix_15(iter,1:time) = electricField_quadratic;
            message_vector15(choose_trasmitter_id)=message_vector15(choose_trasmitter_id)+nodes_in_range;
            %configuration_5=configuration;
        end
        
        
        %plot(1:time, electricField_quadratic, nn{range});


    end
end


%hold on
communication_overhead_0_4=2*sum(message_vector2)/iterations
communication_overhead_0_6=2*sum(message_vector3)/iterations
communication_overhead_1=2*sum(message_vector5)/iterations

communication_overhead_6=2*sum(message_vector6)/iterations
communication_overhead_7=2*sum(message_vector7)/iterations
communication_overhead_8=2*sum(message_vector8)/iterations
communication_overhead_9=2*sum(message_vector9)/iterations
communication_overhead_10=2*sum(message_vector10)/iterations
communication_overhead_11=2*sum(message_vector11)/iterations
communication_overhead_12=2*sum(message_vector12)/iterations
communication_overhead_13=2*sum(message_vector13)/iterations
communication_overhead_14=2*sum(message_vector14)/iterations
communication_overhead_15=2*sum(message_vector15)/iterations




%find the confident interval
for time_counter = 1:time
    %open 2 x time
    xx = Pi_matrix_open(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_open(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals
    
    %1
    xx = Pi_matrix_1(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_1(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals
    
    %2
    xx = Pi_matrix_2(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_2(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals
    
    %3
    xx = Pi_matrix_3(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_3(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals
    
    %4
    xx = Pi_matrix_4(:,time_counter);% Create Data
    SEM = std(xx)/sqrt(length(xx));    % Standard Error
    ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
    CI_4(:,time_counter) = mean(xx) + ts*SEM;    % Confidence Intervals

end



%evresi mesou orou timwn ilektikou pediou gia ta iterations
mesos_oros_open = zeros(1,time);
mesos_oros_1 = zeros(1,time);
mesos_oros_2 = zeros(1,time);
mesos_oros_3 = zeros(1,time);
mesos_oros_4 = zeros(1,time);
mesos_oros_6 = zeros(1,time);
mesos_oros_7 = zeros(1,time);
mesos_oros_8 = zeros(1,time);
mesos_oros_9 = zeros(1,time);
mesos_oros_10 = zeros(1,time);
mesos_oros_11 = zeros(1,time);
mesos_oros_12 = zeros(1,time);
mesos_oros_13 = zeros(1,time);
mesos_oros_14 = zeros(1,time);
mesos_oros_15 = zeros(1,time);

for time_counter = 1:time
    %size(sum(Pi_matrix_open))
    mesos_oros_open = sum(Pi_matrix_open)./iterations;
    mesos_oros_1 = sum(Pi_matrix_1)./iterations;
    mesos_oros_2 = sum(Pi_matrix_2)./iterations;
    mesos_oros_3 = sum(Pi_matrix_3)./iterations;
    mesos_oros_4 = sum(Pi_matrix_4)./iterations;

    mesos_oros_6 = sum(Pi_matrix_6)./iterations;
    mesos_oros_7 = sum(Pi_matrix_7)./iterations;
    mesos_oros_8 = sum(Pi_matrix_8)./iterations;
    mesos_oros_9 = sum(Pi_matrix_9)./iterations;
    mesos_oros_10 = sum(Pi_matrix_10)./iterations;
    mesos_oros_11 = sum(Pi_matrix_11)./iterations;
    mesos_oros_12 = sum(Pi_matrix_12)./iterations;
    mesos_oros_13 = sum(Pi_matrix_13)./iterations;
    mesos_oros_14 = sum(Pi_matrix_14)./iterations;
    mesos_oros_15 = sum(Pi_matrix_15)./iterations;


end

power=[ A*((mesos_oros_6(1,time)^2))/Zo A*(mesos_oros_1(time).^2)/Zo A*(mesos_oros_7(time).^2)/Zo  A*(mesos_oros_2(time).^2)/Zo  A*(mesos_oros_8(time).^2)/Zo  A*(mesos_oros_9(time).^2)/Zo  A*(mesos_oros_10(time).^2)/Zo  A*(mesos_oros_4(time).^2)/Zo  A*(mesos_oros_11(time).^2)/Zo  A*(mesos_oros_12(time).^2)/Zo A*(mesos_oros_13(time).^2)/Zo A*(mesos_oros_14(time)^2)/Zo A*(mesos_oros_15(time).^2)/Zo A*(mesos_oros_open(time).^2)/Zo];
comm_range=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6];
figure(1);
%hold on;
plot(comm_range,power);
hold off;
xlabel('Communication Range (m)')
ylabel('Cumulative Power (watts)')
%legend(range_list_names,'Location','northwest')
set(gca,'Xtick',0.3:0.1:1.6,'XTickLabel',{'0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9' ,'1' ,'1.1', '1.2', '1.3', '1.4' ,'1.5','open'})
%saveas(gcf,'images/CPperRange.eps','eps');
print('images/CPperRange','-depsc','-r0')
saveas(gcf,'images/CPperRange.png','png');

comm_overhead=[ communication_overhead_6 communication_overhead_0_4 communication_overhead_7 communication_overhead_0_6 communication_overhead_8 communication_overhead_9 communication_overhead_10 communication_overhead_1 communication_overhead_11 communication_overhead_12 communication_overhead_13 communication_overhead_14 communication_overhead_15];
comm_range1=[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 ];
figure(2);
%hold on;
plot(comm_range1,comm_overhead);
hold off;
xlabel('Communication Range (m)')
ylabel('Communication Overhead (messages)')
%legend(range_list_names,'Location','northwest')
%ax.XTickLabel = {'0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9' ,'1' ,'1.1', '1.2', '1.3', '1.4' ,'1.5'};
%saveas(gcf,'images/messages.eps','eps');
print('images/messages','-depsc','-r0')
saveas(gcf,'images/messages.png','png');

figure(3);
power1=[ A*((mesos_oros_6(1,time)^2))/Zo A*(mesos_oros_1(time).^2)/Zo A*(mesos_oros_7(time).^2)/Zo  A*(mesos_oros_2(time).^2)/Zo  A*(mesos_oros_8(time).^2)/Zo  A*(mesos_oros_9(time).^2)/Zo  A*(mesos_oros_10(time).^2)/Zo  A*(mesos_oros_4(time).^2)/Zo  A*(mesos_oros_11(time).^2)/Zo  A*(mesos_oros_12(time).^2)/Zo A*(mesos_oros_13(time).^2)/Zo A*(mesos_oros_14(time)^2)/Zo A*(mesos_oros_15(time).^2)/Zo];
plot(comm_range1,power1./comm_overhead);
xlabel('Communication Range (m)')
ylabel({'Cumulative Power / Communication Overhead'; '(watts/messages)'})
%legend(range_list_names,'Location','northwest')
%ax.XTickLabel = {'0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9' ,'1' ,'1.1', '1.2', '1.3', '1.4' ,'1.5'};
%saveas(gcf,'images/ratio.eps','eps');
print('images/ratio','-depsc','-r0')
saveas(gcf,'images/ratio.png','png');


figure(4);
%hold on;
plot( 1:time, A*(mesos_oros_open.^2)/Zo, ':k', 1:time, A*(mesos_oros_1.^2)/Zo, 'k--', 1:time, A*(mesos_oros_2.^2)/Zo, 'k-', 1:time, A*(mesos_oros_4.^2)/Zo, 'k-.');
hold off;
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
legend(range_list_names,'Location','northwest')
%saveas(gcf,'images/quadratic.eps','eps');
print('images/quadratic','-depsc','-r0')
saveas(gcf,'images/quadratic.png','png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kai ola mazi ta interval gia to distance


figure;
hold on;
boxplot(A*(CI_open.^2)/Zo, 'colors','g','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(A*(CI_1.^2)/Zo, 'colors','r','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(A*(CI_2.^2)/Zo, 'colors','b','MedianStyle','target')
set(gca,'XTickLabel',{' '})
boxplot(A*(CI_4.^2)/Zo,'colors','m','MedianStyle','target')

%hold off;
xlabel('Time (rounds)')
ylabel('Cumulative Power (watts)')
set(gca,'XTickLabel',{' '})
%saveas(gcf,'images/ci.eps','eps');
print('images/ci','-depsc','-r0')
saveas(gcf,'images/ci.png','png');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diaforetikos  arithmos k
GREEDY_mean=zeros(1,((stop-start)/step)+1);
SAMPLE_mean=zeros(1,((stop-start)/step)+1);
MYALG_mean=zeros(1,((stop-start)/step)+1);
OPT_mean=zeros(1,((stop-start)/step)+1);

%pinakas pou tha apothikepsw ton SAMPLE gia ta diaforetika k kai
%diaforetika x     %%%%%%%%   diastaseis: XxK
depend_on_x_SAMPLE=zeros(((stopx-startx)/stepx)+1,((stop-start)/step)+1);

counter44=0;
for k=start:step:stop
    counter44=counter44+1;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %BASELINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %k=20;
    max_min_state=zeros(1,Trans_num);
    %dimiourgia twn sindiasmwn-permutations ektos tis oloi na einai kleistoi
    output=dec2bin(2^Trans_num-1:-1:0)-'0';
    for counter2=1:2^Trans_num%EPANALIPSEIS GIA OLA TA PERMUTATIONS
        transm_initial_state=output(counter2,:);
        if (sum(transm_initial_state)==0) continue; end%an einai ola 0 pigaine sthn epomeni epanalipsi


        %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
        %arxikopoiisi pinakwn
        x_coordination_Trans1=zeros(1,sum(transm_initial_state));
        y_coordination_Trans1=zeros(1,sum(transm_initial_state));
        counter3=1;
        for counter1=1:Trans_num
            if (transm_initial_state(counter1)==1)%an sto initial state einai 1 bale ton sto dianisma sidetagmenwn
                x_coordination_Trans1(counter3)=x_coordination_Trans(counter1);
                y_coordination_Trans1(counter3)=y_coordination_Trans(counter1);
                counter3=counter3+1;
            end
        end

        transm_initial_state;
        %baseline
        individual_power_for_each_receiver=zeros(1,Receiv_num);
        for counter4=1:Receiv_num
            %epistrefei to receiving power tou receiver 
            %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
            x_coordinations_Receiv1=x_coordinations_Receiv(counter4);
            y_coordinations_Receiv1=y_coordinations_Receiv(counter4);
            total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
            individual_power_for_each_receiver(counter4)=(A*total_power_received^2)/Zo;
        end
        clearvars x_coordination_Trans1 y_coordination_Trans1

        %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
        %me ta epomena gia na brw thn max k-ada
        sorted_vector=sort(individual_power_for_each_receiver);
        athroisma_min_k_adas=sum(sorted_vector(1:k));

        if (counter2==1)
            max_min_k_ada=athroisma_min_k_adas;
            max_min_state=transm_initial_state;
        else
            if(max_min_k_ada<athroisma_min_k_adas)
                max_min_k_ada=athroisma_min_k_adas;
                max_min_state=transm_initial_state;
            end
        end


    end
    %paizw me ilectriko pedio
    OPT_max_min_k_ada=max_min_k_ada
    OPT_max_min_state=max_min_state





    GREEDY_max_min_k_ada=zeros(1,iterations);
    SAMPLE_max_min_k_ada=zeros(1,iterations);
    for counter20=1:iterations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %GREEDY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        number_of_chargers_to_check=8;

        %random configuration
        transm_initial_state= randi([0,1], 1, Trans_num);
        for counter9=1:number_of_chargers_to_check
            %random charger
            choose_trasmitter_id = randi([1,Trans_num]);

            for counter5=0:1%2 fores anoigw kleinw ton charger na dw pou einai i max min k-ada

                %tin prwti fora ton exw kleisto kai tin defteri anoikto
                if (counter5==0)
                    transm_initial_state(choose_trasmitter_id)=0;
                    if(sum(transm_initial_state)==0) %an einai oloi kleisti pigaine sthn epomeni epanalipsi
                        transm_initial_state(choose_trasmitter_id)=1;%kai idi thetoume na einai anoixtos
                        max_min_k_ada=0;%an tixei wste stin epomeni epnalipsi na borei na sigrinei me kati  pio katw
                        break;
                    end

                else
                    transm_initial_state(choose_trasmitter_id)=1;
                end

                %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
                x_coordination_Trans1=zeros(1,sum(transm_initial_state));
                y_coordination_Trans1=zeros(1,sum(transm_initial_state));
                counter7=1;
                for counter6=1:Trans_num
                    if (transm_initial_state(counter6)==1)
                        x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                        y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                        counter7=counter7+1;
                    end
                end


                individual_power_for_each_receiver=zeros(1,Receiv_num);
                for counter8=1:Receiv_num
                    %epistrefei to receiving power tou receiver 
                    %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
                    x_coordinations_Receiv1=x_coordinations_Receiv(counter8);
                    y_coordinations_Receiv1=y_coordinations_Receiv(counter8);
                    total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
                    individual_power_for_each_receiver(counter8)=(A*total_power_received^2)/Zo;
                end
                clearvars x_coordination_Trans1 y_coordination_Trans1


                %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
                %me ta epomena gia na brw thn max k-ada
                sorted_vector=sort(individual_power_for_each_receiver);
                athroisma_min_k_adas=sum(sorted_vector(1:k));

                %vriskw thn max min gia otan einai anoixtos i kleistos o enas
                %charger
                if (counter5==0)
                    max_min_k_ada=athroisma_min_k_adas;
                    max_min_state=transm_initial_state;
                else
                    if(max_min_k_ada<athroisma_min_k_adas)
                        max_min_k_ada=athroisma_min_k_adas;
                        max_min_state=transm_initial_state;
                    end
                end


            end

            %evresi tis olikis max min k-adas (krataw to kalitero gia osous exw elegxei)
            if (counter9==1)
                final_max_min_k_ada=max_min_k_ada;
                final_max_min_state=max_min_state;
            else
                if(final_max_min_k_ada<max_min_k_ada)
                    final_max_min_k_ada=max_min_k_ada;
                    final_max_min_state=max_min_state;
                end
            end


        end

        GREEDY_max_min_k_ada(counter20)=final_max_min_k_ada
        GREEDY_max_min_state=final_max_min_state;




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SAMPLE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        counter45=0;
        for x=startx:stepx:stopx
        counter45=counter45+1;
        %x=15;%arithmos k-adwn, configuration
        %dialegw tixaia poies grammes apo to output tha parw
        configuration_matrix= zeros(x,Trans_num);


        %epilegw random x k-ades
        x_random_k_ades_matrix=zeros(x,k);
        for counter11=1:x
            x_random_k_ades_matrix(counter11,:)=randperm(Receiv_num,k);
        end
        x_random_k_ades_matrix;

        %chargers_initial_state= randi([0,1], 1, Trans_num);
        for counter22=1:x

            %dimiourgia tn sidetagmenwn gia receivers sti k ada
            x_coordination_Receiv2=zeros(x,k);
            y_coordination_Receiv2=zeros(x,k);
            for counter26=1:k          
                x_coordination_Receiv2(counter22,counter26)=x_coordinations_Receiv(x_random_k_ades_matrix(counter22,counter26));
                y_coordination_Receiv2(counter22,counter26)=y_coordinations_Receiv(x_random_k_ades_matrix(counter22,counter26));          
            end


            %apo to quadratic to opt configuration
            chargers_initial_state= randi([0,1], 1, Trans_num);
            max_electric_field_quadratic=0;
            opt_configuration_quadratic=zeros(1,time);
            for iter=1:iterations
                choose_trasmitter = randi([1,Trans_num], 1, time);
                [OPT_Transmitters_configuration,electricField_quadratic] = myquadraticN_distributedAlg_v5_return_OPT_configuration (x_coordination_Trans,y_coordination_Trans,x_coordination_Receiv2(counter22,:),y_coordination_Receiv2(counter22,:),choose_trasmitter,chargers_initial_state, -1 ,time,plane_size,lambda); 
                if (max_electric_field_quadratic<electricField_quadratic(time))
                    max_electric_field_quadratic=electricField_quadratic(time);
                    opt_configuration_quadratic=OPT_Transmitters_configuration;
                end
            end
            x_coordination_Receiv2(counter22,:);
            max_electric_field_quadratic;
            opt_configuration_quadratic';

            %pernaw ta configuration apo to quadratic sto matrix
            configuration_matrix(counter22,:)=opt_configuration_quadratic;
            %configuration_matrix(counter10,:)=output(random_conf,:);
            configuration_matrix;

        end

        configuration_matrix;
        %permutation twn chargers
        choose_trasmitter_id_vector=randperm(Trans_num,Trans_num);
        for counter12=1:Trans_num
            %random charger
            choose_trasmitter_id = choose_trasmitter_id_vector(counter12);

            %krataw to athroisma tis energeias apo tis x k ade gia na to sigrinw sto telos
            sumL0=0;
            sumL1=0;

            for counter13=1:x

                %apothikevw gia tin kathe k ada configuration thn energeia apo to 0
                %kai to 1 gia na tin sigrinw kai na to balw sto sumL0 i sumL1
                max_min_state0=0;
                max_min_state1=0;

                %to idio me ton greedy
                for counter5=0:1%2 fores anoigw kleinw ton charger na dw pou einai i max min k-ada

                    %tin prwti fora ton exw kleisto kai tin defteri anoikto
                    if (counter5==0)
                        configuration_matrix(counter13,choose_trasmitter_id)=0;
                        if(sum(configuration_matrix(counter13,:))==0) %an einai oloi kleisti pigaine sthn epomeni epanalipsi
                            configuration_matrix(counter13,choose_trasmitter_id)=1;%kai idi thetoume na einai anoixtos
                            max_min_k_ada=0;%an tixei wste stin epomeni epnalipsi na borei na sigrinei me kati  pio katw
                            break;
                        end

                    else
                        configuration_matrix(counter13,choose_trasmitter_id)=1;
                    end

                    %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
                    x_coordination_Trans1=zeros(1,sum(configuration_matrix(counter13,:)));
                    y_coordination_Trans1=zeros(1,sum(configuration_matrix(counter13,:)));
                    counter7=1;
                    for counter6=1:Trans_num
                        if (configuration_matrix(counter13,counter6)==1)
                            x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                            y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                            counter7=counter7+1;
                        end
                    end


                    individual_power_for_each_receiver=zeros(1,k);
                    for counter38=1:k
                        %epistrefei to receiving power tou receiver 
                        %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
                        x_coordinations_Receiv1=x_coordination_Receiv2(counter22,counter38);
                        y_coordinations_Receiv1=y_coordination_Receiv2(counter22,counter38);
                        total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
                        individual_power_for_each_receiver(counter38)=total_power_received;
                    end
                    clearvars x_coordination_Trans1 y_coordination_Trans1


                    %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
                    %me ta epomena gia na brw thn max k-ada
                    sorted_vector=sort(individual_power_for_each_receiver);
                    athroisma_min_k_adas=sum(sorted_vector(1:k));

                    %vriskw thn max min gia otan einai anoixtos i kleistos o enas
                    %charger
                    if (counter5==0)
                        max_min_k_ada0=athroisma_min_k_adas;
                        max_min_state=configuration_matrix(counter13,:);
                        max_min_state0=max_min_state;
                    else
                        %if(max_min_k_ada<athroisma_min_k_adas)
                            max_min_k_ada1=athroisma_min_k_adas;
                            max_min_state=configuration_matrix(counter13,:);
                            max_min_state1=max_min_state;
                        %end
                    end


                end


                if(max_min_k_ada0>max_min_k_ada1)
                    sumL0=sumL0+max_min_k_ada0-max_min_k_ada1;
                else
                    sumL1=sumL1+max_min_k_ada1-max_min_k_ada0;
                end

                configuration_matrix;

            end

            sumL0;
            sumL1;


            if (sumL0>sumL1)
                configuration_matrix(:,choose_trasmitter_id)=0;
            else
                configuration_matrix(:,choose_trasmitter_id)=1;
            end


            configuration_matrix;

        end

        SAMPLE_max_min_state=configuration_matrix;

        %gia na dw to ilectriko paidio sto configuration pu m bgazei o
        %algorithmos
        %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
        x_coordination_Trans1=zeros(1,sum(configuration_matrix(counter13,:)));
        y_coordination_Trans1=zeros(1,sum(configuration_matrix(counter13,:)));
        counter7=1;
        for counter6=1:Trans_num
            if ( configuration_matrix(1,counter6)==1)
                x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                counter7=counter7+1;
            end
        end

        individual_power_for_each_receiver=zeros(1,Receiv_num);
        for counter8=1:Receiv_num
            %epistrefei to receiving power tou receiver 
            %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
            x_coordinations_Receiv1=x_coordinations_Receiv(counter8);
            y_coordinations_Receiv1=y_coordinations_Receiv(counter8);
            total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
            individual_power_for_each_receiver(counter8)=(A*total_power_received^2)/Zo;
        end
        clearvars x_coordination_Trans1 y_coordination_Trans1


        %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
        %me ta epomena gia na brw thn max k-ada
        sorted_vector=sort(individual_power_for_each_receiver);
        athroisma_min_k_adas=sum(sorted_vector(1:k));
        SAMPLE_max_min_k_ada(counter20)=athroisma_min_k_adas
        depend_on_x_SAMPLE(counter45,counter44)=depend_on_x_SAMPLE(counter45,counter44)+SAMPLE_max_min_k_ada(counter20);
        end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MYALG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find for each mote his opt configuration
    opt_individualConfigurationState=zeros(Receiv_num,Trans_num);
    individualElectricField=zeros(1,Receiv_num);
    for counter40=1:Receiv_num
            %apo to quadratic to opt configuration
            chargers_initial_state= randi([0,1], 1, Trans_num);
            max_electric_field_quadratic=0;
            opt_configuration_quadratic=zeros(1,time);
            for iter=1:iterations
                choose_trasmitter = randi([1,Trans_num], 1, time);
                [OPT_Transmitters_configuration,electricField_quadratic] = myquadraticN_distributedAlg_v5_return_OPT_configuration (x_coordination_Trans,y_coordination_Trans,x_coordinations_Receiv(counter40),y_coordinations_Receiv(counter40),choose_trasmitter,chargers_initial_state, -1 ,time,plane_size,lambda); 
                if (max_electric_field_quadratic<electricField_quadratic(time))
                    max_electric_field_quadratic=electricField_quadratic(time);
                    opt_configuration_quadratic=OPT_Transmitters_configuration;
                end
            end
            x_coordination_Receiv2(counter22,:);
            max_electric_field_quadratic;
            opt_configuration_quadratic';

            opt_individualConfigurationState(counter40,:)=opt_configuration_quadratic;
            individualElectricField(counter40)=max_electric_field_quadratic;
    end


    opt_individualConfigurationState1=opt_individualConfigurationState;
    individualElectricField;
    for counter43=1:iterations
        %gia na xanaparei ta opt kai oxi afta apo thn proigoumeni epanalipsi
        opt_individualConfigurationState=opt_individualConfigurationState1;
        choose_trasmitter_id_vector=randperm(Trans_num,Trans_num);
        for counter41=1:Trans_num
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %otan einai 0 o charger
            choose_trasmitter_id = choose_trasmitter_id_vector(counter41);
            opt_individualConfigurationState(:,choose_trasmitter_id)=0;
            %ta vector pou blepw otan ton kanw 0 i ena kai meta vriskw min k ada
            individual_power_for_each_receiver0=zeros(1,Receiv_num);
            individual_power_for_each_receiver1=zeros(1,Receiv_num);

            for counter42=1:Receiv_num
                %gia na dw to ilectriko paidio sto configuration pu m bgazei o
                %algorithmos
                %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
                x_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(counter42,:)));
                y_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(counter42,:)));
                if(sum(opt_individualConfigurationState(counter42,:))==0)
                    sum(opt_individualConfigurationState(counter42,:))==0
                end
                counter7=1;
                for counter6=1:Trans_num
                    if ( opt_individualConfigurationState(counter42,counter6)==1)
                        x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                        y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                        counter7=counter7+1;
                    end
                end

                if (sum(x_coordination_Trans1)==0)
                    display('ALARM');
                end
                total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv(counter42),y_coordinations_Receiv(counter42));
                individual_power_for_each_receiver0(counter42)=total_power_received;
            end
            clearvars x_coordination_Trans1 y_coordination_Trans1


            %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
            %me ta epomena gia na brw thn max k-ada
            sorted_vector=sort(individual_power_for_each_receiver0);
            athroisma_min_k_adas0=sum(sorted_vector(1:k));
            MYALG_max_min_k_ada0=athroisma_min_k_adas0;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %otan einai 1 o charger
            %choose_trasmitter_id = choose_trasmitter_id_vector(counter41);
            opt_individualConfigurationState(:,choose_trasmitter_id)=1;
            %ta vector pou blepw otan ton kanw 0 i ena kai meta vriskw min k ada

            for counter42=1:Receiv_num
                %gia na dw to ilectriko paidio sto configuration pu m bgazei o
                %algorithmos
                %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
                x_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(counter42,:)));
                y_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(counter42,:)));
                counter7=1;
                for counter6=1:Trans_num
                    if ( opt_individualConfigurationState(counter42,counter6)==1)
                        x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                        y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                        counter7=counter7+1;
                    end
                end


                total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv(counter42),y_coordinations_Receiv(counter42));
                individual_power_for_each_receiver1(counter42)=total_power_received;
                clearvars x_coordination_Trans1 y_coordination_Trans1
            end


            %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
            %me ta epomena gia na brw thn max k-ada
            sorted_vector=sort(individual_power_for_each_receiver1);
            athroisma_min_k_adas1=sum(sorted_vector(1:k));
            MYALG_max_min_k_ada1=athroisma_min_k_adas1;

            %an einai mikrotero kleise ton charger se ola ta configurations
            if (MYALG_max_min_k_ada1<MYALG_max_min_k_ada0)
                opt_individualConfigurationState(:,choose_trasmitter_id)=0;
            end
            opt_individualConfigurationState;

        end


        %gia na dw to ilectriko paidio sto configuration pu m bgazei o
        %my algorithmos
        %dimiourgia tn sidetagmenwn gia osous einai anoiktoi
        x_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(1,:)));
        y_coordination_Trans1=zeros(1,sum(opt_individualConfigurationState(1,:)));
        counter7=1;
        for counter6=1:Trans_num
            if ( opt_individualConfigurationState(1,counter6)==1)
                x_coordination_Trans1(counter7)=x_coordination_Trans(counter6);
                y_coordination_Trans1(counter7)=y_coordination_Trans(counter6);
                counter7=counter7+1;
            end
        end

        individual_power_for_each_receiver=zeros(1,Receiv_num);
        for counter8=1:Receiv_num
            %epistrefei to receiving power tou receiver 
            %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
            x_coordinations_Receiv1=x_coordinations_Receiv(counter8);
            y_coordinations_Receiv1=y_coordinations_Receiv(counter8);
            total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
            individual_power_for_each_receiver(counter8)= A*(total_power_received^2)/Zo;
        end
        clearvars x_coordination_Trans1 y_coordination_Trans1

        %taxinomisi se afxousa kai athroizw ta prwta k kai ta krataw na sigrinw
        %me ta epomena gia na brw thn max k-ada
        sorted_vector=sort(individual_power_for_each_receiver);
        athroisma_min_k_adas=sum(sorted_vector(1:k));
        MYALG_max_min_k_ada(counter43)=athroisma_min_k_adas;
        %MYALG configuration state
        opt_individualConfigurationState(1,:);

    end



    MYALG_max_min_k_ada
    for iiii=1:iterations
       if( MYALG_max_min_k_ada(iiii)>OPT_max_min_k_ada)
           display('ALARM')
       end
    end

    
       counter44
    SEM = std(GREEDY_max_min_k_ada)/sqrt(length(GREEDY_max_min_k_ada));    % Standard Error
    ts = tinv([0.025  0.975],length(GREEDY_max_min_k_ada)-1);    % T-Score
    CI_GREEDY = mean(GREEDY_max_min_k_ada) + ts*SEM;    % Confidence Intervals
    GREEDY_mean(counter44)=mean(GREEDY_max_min_k_ada);

    SEM = std(SAMPLE_max_min_k_ada)/sqrt(length(SAMPLE_max_min_k_ada));    % Standard Error
    ts = tinv([0.025  0.975],length(SAMPLE_max_min_k_ada)-1);    % T-Score
    CI_SAMPLE = mean(SAMPLE_max_min_k_ada) + ts*SEM;    % Confidence Intervals
    SAMPLE_mean(counter44)=mean(SAMPLE_max_min_k_ada);

    SEM = std(MYALG_max_min_k_ada)/sqrt(length(MYALG_max_min_k_ada));    % Standard Error
    ts = tinv([0.025  0.975],length(MYALG_max_min_k_ada)-1);    % T-Score
    CI_MYALG = mean(MYALG_max_min_k_ada) + ts*SEM;    % Confidence Intervals
    MYALG_mean(counter44)=mean(MYALG_max_min_k_ada);

    OPT_mean(counter44)=OPT_max_min_k_ada;

    SAMPLE_max_min_k_ada;

    confidence_intervals=zeros(iterations,3);
    confidence_intervals(:,1)=OPT_max_min_k_ada;
    confidence_intervals(:,2)=GREEDY_max_min_k_ada;
    confidence_intervals(:,3)=SAMPLE_max_min_k_ada;
    confidence_intervals(:,4)=MYALG_max_min_k_ada;

    for counter21=1:4
        xx = confidence_intervals(:,counter21);% Create Data
        SEM = std(xx)/sqrt(length(xx));    % Standard Error
        ts = tinv([0.025  0.975],length(xx)-1);    % T-Score
        CI(:,counter21) = mean(xx) + ts*SEM;    % Confidence Intervals
    end
    
    
    CI_GREEDY_matrix(:,counter44)=CI_GREEDY;
    CI_SAMPLE_matrix(:,counter44)=CI_SAMPLE;
    CI_MYALG_matrix(:,counter44)=CI_MYALG;
    
    
  


end


figure(7)
hold on;
boxplot( CI_GREEDY_matrix, start:step:stop,'colors','r','MedianStyle','target')
set(gca,'XTickLabel',{' '})



%figure   
boxplot( CI_SAMPLE_matrix, start:step:stop,'colors','g','MedianStyle','target')
set(gca,'XTickLabel',{' '})


%figure   
boxplot( CI_MYALG_matrix, start:step:stop,'colors','b','MedianStyle','target')
%set(gca,'Ytick',0:0.1:1.6,'YTickLabel',{'0','0.02', '0.04','0.06', '0.08','0.10', '0.12','0.14', '0.16','0.18','0.2'})
%legend({'OPT ',' GRE ',' SAM ',' FUS '},'Location','northwest')
xlabel('k')
ylabel('Cumulative Power (watts)')
%saveas(gcf,'images/generalCI.eps','eps');
print('images/generalCI','-depsc','-r0')
saveas(gcf,'images/generalCI.png','png');

hold off;


%m = {'-ok', '-*k', '-+k', '-.k', '-xk', '-^k', '-pk', '--ok'};

OPT_mean;
 GREEDY_mean;
 SAMPLE_mean;
 MYALG_mean;
 
figure;
plot(start:step:stop,  OPT_mean, '-ok', start:step:stop, GREEDY_mean,'-xk',start:step:stop,  SAMPLE_mean,'-.+k',start:step:stop, MYALG_mean, '-*k')
%set(gca,'Ytick',0:0.1:1.6,'YTickLabel',{'0','0.02', '0.04','0.06', '0.08','0.10', '0.12','0.14', '0.16','0.18','0.2'})
xlabel('k')
ylabel('Cumulative Power (watt)')
legend({'OPT ',' GRE ',' SAM ',' FUS '},'Location','northwest')
%saveas(gcf,'images/general.eps','eps');
print('images/general','-depsc','-r0')
saveas(gcf,'images/general.png','png');


%SAMPLE depending on x
%depend_on_x_SAMPLE=depend_on_x_SAMPLE/iterations;
%figure
%plot(start:step:stop,depend_on_x_SAMPLE(1,:),start:step:stop,depend_on_x_SAMPLE(2,:),start:step:stop,depend_on_x_SAMPLE(3,:));
%title('SAMPLE depending on x')
%xlabel('k')
%ylabel('cumulative electric field stin min k ada')
%legend({'x=5 ',' x=10 ',' x=15 '},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          POWER BALANCE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comm_range_conf_matrix=zeros(Trans_num,4);
individual_power_for_each_receiver_matrix=zeros(Receiv_num,4);

comm_range_conf_matrix(1:Trans_num,1)=configuration_open;
comm_range_conf_matrix(1:Trans_num,2)=configuration_2;
comm_range_conf_matrix(1:Trans_num,3)=configuration_3;
comm_range_conf_matrix(1:Trans_num,4)=configuration_5;

for counter51=1:4
    counter7=1;
    for counter52=1:Trans_num
        if(comm_range_conf_matrix(counter52,counter51)==1)
            x_coordination_Trans1(counter7)=x_coordination_Trans(counter52);
            y_coordination_Trans1(counter7)=y_coordination_Trans(counter52);
            counter7=counter7+1;
        end
    end
    for counter8=1:Receiv_num
        %epistrefei to receiving power tou receiver 
        %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
        x_coordinations_Receiv1=x_coordinations_Receiv(counter8);
        y_coordinations_Receiv1=y_coordinations_Receiv(counter8);
        total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
        individual_power_for_each_receiver_matrix(counter8,counter51)=A*(total_power_received^2)/Zo;
    end
    clearvars x_coordination_Trans1 y_coordination_Trans1
end
 sum(individual_power_for_each_receiver_matrix(:,1))
sum(individual_power_for_each_receiver_matrix(:,2))
sum(individual_power_for_each_receiver_matrix(:,3))
sum(individual_power_for_each_receiver_matrix(:,4))           


            
figure
ax1 = subplot(4,1,1);
area(sort(individual_power_for_each_receiver_matrix(:,1)))
set(gca,'XTickLabel',{' '})
title('open')
saveas(gcf,'images/PB_range=open.eps','eps');
saveas(gcf,'images/PB_range=open.png','png');
%figure
ax2 = subplot(4,1,2);
area(sort(individual_power_for_each_receiver_matrix(:,2)))
set(gca,'XTickLabel',{' '})
title('0.4')
saveas(gcf,'images/PB_range=0.4.eps','eps');
saveas(gcf,'images/PB_range=0.4.png','png');
%figure
ax3 = subplot(4,1,3);
area(sort(individual_power_for_each_receiver_matrix(:,3)))
set(gca,'XTickLabel',{' '})
title('0.6')
saveas(gcf,'images/PB_range=0.6.eps','eps');
saveas(gcf,'images/PB_range=0.6.png','png');
%figure
ax4 = subplot(4,1,4);
linkaxes([ax1,ax2,ax3,ax4],'xy');
area(sort(individual_power_for_each_receiver_matrix(:,4)))
set(gca,'XTickLabel',{' '})
title('1')
%saveas(gcf,'images/PB_range=1.eps','eps');
print('images/PB_range=1','-depsc','-r0')
saveas(gcf,'images/PB_range=1.png','png');
            
            
			
mean_open=mean(individual_power_for_each_receiver_matrix(:,1)/sum(individual_power_for_each_receiver_matrix(:,1)));
mean_0_4=mean(individual_power_for_each_receiver_matrix(:,2)/sum(individual_power_for_each_receiver_matrix(:,2)));
mean_0_6=mean(individual_power_for_each_receiver_matrix(:,3)/sum(individual_power_for_each_receiver_matrix(:,3)));
mean_1=mean(individual_power_for_each_receiver_matrix(:,4)/sum(individual_power_for_each_receiver_matrix(:,4)));		

total_variation_distance_open=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,1)/sum(individual_power_for_each_receiver_matrix(:,1)))-mean_open))
total_variation_distance_0_4=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,2)/sum(individual_power_for_each_receiver_matrix(:,2)))-mean_0_4))	
total_variation_distance_0_6=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,3)/sum(individual_power_for_each_receiver_matrix(:,3)))-mean_0_6))	
total_variation_distance_1=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,4)/sum(individual_power_for_each_receiver_matrix(:,4)))-mean_1))	

			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   GENERAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


comm_range_conf_matrix(1:Trans_num,1)=OPT_max_min_state;
comm_range_conf_matrix(1:Trans_num,2)=GREEDY_max_min_state;
comm_range_conf_matrix(1:Trans_num,3)=SAMPLE_max_min_state(1,:);
comm_range_conf_matrix(1:Trans_num,4)=opt_individualConfigurationState(1,:);

for counter51=1:4
    counter7=1;
    for counter52=1:Trans_num
        if(comm_range_conf_matrix(counter52,counter51)==1)
            x_coordination_Trans1(counter7)=x_coordination_Trans(counter52);
            y_coordination_Trans1(counter7)=y_coordination_Trans(counter52);
            counter7=counter7+1;
        end
    end
    for counter8=1:Receiv_num
        %epistrefei to receiving power tou receiver 
        %osoi chargers einai kleistoi den tous vazw sta coordination_Trans
        x_coordinations_Receiv1=x_coordinations_Receiv(counter8);
        y_coordinations_Receiv1=y_coordinations_Receiv(counter8);
        total_power_received = baselineAlgorithm (x_coordination_Trans1,y_coordination_Trans1,x_coordinations_Receiv1,y_coordinations_Receiv1);
        individual_power_for_each_receiver_matrix(counter8,counter51)=A*(total_power_received^2)/Zo;
    end
    clearvars x_coordination_Trans1 y_coordination_Trans1
end
            




figure
ax1 = subplot(4,1,1);
area(sort(individual_power_for_each_receiver_matrix(:,1)))
set(gca,'XTickLabel',{' '})
title('OPT')
saveas(gcf,'images/PB_OPT.eps','eps');
saveas(gcf,'images/PB_OPT.png','png');
%figure
ax2 = subplot(4,1,2);
area(sort(individual_power_for_each_receiver_matrix(:,2)))
set(gca,'XTickLabel',{' '})
title('GRE')
saveas(gcf,'images/PB_GRE.eps','eps');
saveas(gcf,'images/PB_GRE.png','png');
%figure
ax3 = subplot(4,1,3);
area(sort(individual_power_for_each_receiver_matrix(:,3)))
set(gca,'XTickLabel',{' '})
title('SAM')
ylabel('           Power (watts)')
saveas(gcf,'images/PB_SAM.eps','eps');
saveas(gcf,'images/PB_SAM.png','png');
%figure
ax4 = subplot(4,1,4);
linkaxes([ax4,ax3,ax2,ax1],'xy');
area(sort(individual_power_for_each_receiver_matrix(:,4)))
set(gca,'XTickLabel',{' '})
title('FUS')
xlabel('Nodes')
%saveas(gcf,'images/PB_FUS.eps','eps');
print('images/PB_FUS','-depsc','-r0')
saveas(gcf,'images/PB_FUS.png','png');

mean_OPT=mean(individual_power_for_each_receiver_matrix(:,1)/sum(individual_power_for_each_receiver_matrix(:,1)));
mean_GRE=mean(individual_power_for_each_receiver_matrix(:,2)/sum(individual_power_for_each_receiver_matrix(:,2)));
mean_SAM=mean(individual_power_for_each_receiver_matrix(:,3)/sum(individual_power_for_each_receiver_matrix(:,3)));
mean_FUS=mean(individual_power_for_each_receiver_matrix(:,4)/sum(individual_power_for_each_receiver_matrix(:,4)));		

total_variation_distance_OPT=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,1)/sum(individual_power_for_each_receiver_matrix(:,1)))-mean_OPT))	
total_variation_distance_GRE=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,2)/sum(individual_power_for_each_receiver_matrix(:,2)))-mean_GRE))	
total_variation_distance_SAM=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,3)/sum(individual_power_for_each_receiver_matrix(:,3)))-mean_SAM))	
total_variation_distance_FUS=(1/2)*sum(abs((individual_power_for_each_receiver_matrix(:,4)/sum(individual_power_for_each_receiver_matrix(:,4)))-mean_FUS))	

%figure
%hold on
%plot(sum(individual_power_for_each_receiver_matrix(1:45,1)),'ok')
%plot(sum(individual_power_for_each_receiver_matrix(1:45,2)),'*r')
%plot(sum(individual_power_for_each_receiver_matrix(1:45,3)),'xg')
%plot(sum(individual_power_for_each_receiver_matrix(1:45,4)),'+b')
%hold off

end





