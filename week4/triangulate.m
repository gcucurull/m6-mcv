function [ X ] = triangulate( x1, x2, P1, P2, imsize )
%TRIANGULATE Summary of this function goes here
%   Detailed explanation goes here
    x1=homog(x1);
    x2=homog(x2);
    
    x1p1_3t= x1(1)*P1(3,:);
    p1_1t= P1(1,:);
    
    y1p1_3t= x1(2)*P1(3,:);
    p1_2t= P1(2,:);
    
    x2p2_3t=x2(1)*P2(3,:);
    p2_1t=P2(1,:);
    
    y2p2_3t=x2(2)*P2(3,:);
    p2_2t=P2(2,:);
    
    
    A=[ x1p1_3t-p1_1t;
        y1p1_3t-p1_2t;
        x2p2_3t-p2_1t;
        y2p2_3t-p2_2t ];
    
    
    
    [u d v] = svd(A);
    
    X= v(:,size(v,1));
    
    
    
end

