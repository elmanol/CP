function [ linex, liney ] = liner(device, center)
    distance=norm(device-center)/0.02;
    p=[ linspace(device(1),center(1),round(distance)); linspace(device(2),center(2),round(distance))];
    linex = p(1,:);
    liney = p(2,:);
end