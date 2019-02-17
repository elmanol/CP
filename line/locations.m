function [locDev]=locations(nPoints, stopx,stopy,minAllowableDistance,x_c_T,y_c_T)

counter = 0;

while(counter<=nPoints)
	% Get a trial point.
    thisX = round(rand(1)*stopx*20,0)/20;
    thisY = round(rand(1)*stopy*20,0)/20;
    
	% See how far is is away from existing keeper points.
	distances = sqrt((thisX-x_c_T).^2 + (thisY - y_c_T).^2);
	minDistance = min(distances);
	if minDistance >= minAllowableDistance
        counter = counter + 1;
		keeperX(counter) = thisX;
		keeperY(counter) = thisY;		
    end
end


%grid on;

locDev=[keeperX;keeperY];
end