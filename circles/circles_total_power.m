function total_power_received = circles_total_power( x_c_T,Rec_iter,distance,lambda)
    %initialize parameters
    Pt = 2;
    Gt = 2;
    Gr = 1;
    Zo=119.9169832*pi;
    A=Gr*(lambda^2/(4*pi));
    k= 2*pi/lambda;

    for i=1:length(x_c_T)
        for t=1:length(Rec_iter)
            j=Rec_iter(t);
            S(i,j)=(Gt*Pt)/(4*pi*distance(i,j)^2);    
            u(i,j)=sqrt(Zo*S(i,j))*exp(-1i*(k*distance(i,j)));
        end 
    end

    for t=1:length(Rec_iter)
        j=Rec_iter(t);
        superposition_of_elect_fields(1,j)=sum(u(:,j));        
        superposition_magnitude(1,j)=abs(superposition_of_elect_fields(1,j));

        %density of the total transferred power
        ST(1,j)=(superposition_magnitude(1,j)^2)/Zo;

        %total power received at point of interest
        total_power_received(1,j)=ST(1,j)*A;
    end

    
end