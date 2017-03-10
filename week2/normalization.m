function [ xnorm, T ] = normalization( x )
    x(1,1:4)=x(1,1:4)./x(3,1:4);
    x(2,1:4)=x(1,1:4)./x(3,1:4);
    x(3,1:4)=1;
    centroid = mean(x(1:2,1:4)')';
    newx(1,1:4) =  x(1,1:4)-centroid(1);
    newx(2,1:4) = x(2,1:4)-centroid(2);
    dist = sqrt(newx(1,1:4).^2 + newx(2,1:4).^2);
    meandist = mean(dist(:));  
    scale = sqrt(2)/meandist;
    T = [scale   0   -scale*centroid(1)
         0     scale -scale*centroid(2)
         0       0      1      ];
    xnorm = T*x;
end

