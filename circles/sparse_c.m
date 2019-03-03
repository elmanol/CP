function [allX,allY] = sparse_c(num,stopx,stopy,minAllowableDistance)  


    i=1;
    x=stopx*rand(1,1);
    y=stopy*rand(1,1);
    
    allX(1)=x;
    allY(1)=y;

    minDistance=0;
   
    while(i<num)
        
        x=10*rand(1,1);
        y=10*rand(1,1);
        distances = sqrt((x-allX).^2 + (y - allY).^2);
        minDistance = min(distances);
        
        while(minDistance<minAllowableDistance)
            x=10*rand(1,1);
            y=10*rand(1,1);
            distances = sqrt((x-allX).^2 + (y - allY).^2);
            minDistance = min(distances);
        end
        i=i+1;
        allX(i)=x;
        allY(i)=y;

    end
    
    
end