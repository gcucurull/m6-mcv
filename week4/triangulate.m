function [ X ] = triangulate( x1, x2, P1, P2, imsize )
%TRIANGULATE Summary of this function goes here
%   The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

    x1=homog(x1);
    x2=homog(x2);
    
    % build matrix H (slide 8 lecture 7)
    nx = imsize(1);
    ny = imsize(2);
    
    H = [2/nx, 0, -1;
        0, 2/ny, -1;
        0, 0, 1;];
    
    % Transform x, x', P and P' by H
    x1 = euclid(H*x1);
    x2 = euclid(H*x2);
    P1 = H*P1;
    P2 = H*P2;
    
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
    %X = X ./X(end); % Set 4th coordinate to 1
end

