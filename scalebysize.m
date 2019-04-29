function [ newnet ] = scalebysize( net, sizes )

% Code to scale the weights of a matrix by size of a reef
% Author & copyright: Karlo Hock, University of Queensland. 2019

nrf=length(net);
newnet=zeros(nrf);
for ri=1:nrf
    for rj=1:nrf
        newnet(ri,rj)=net(ri,rj)*sizes(ri,2);
    end
end
end

