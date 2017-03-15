function [ diff ] = gs_errfunction( P0, Xobs )
%GS_ERRFUNCTION Summary of this function goes here
%   Arguments:
%       - P0: [ Hab(:) ; x(:) ] contains the H matrix (first 9 rows) and the x points
%       - Xobs: [ x(:) ; xp(:) ] contains x points and x'
%   All points are received in euclidean coordinates (2 dimensions)

    % First get the H matrix
    H = reshape(P0(1:9), [3,3]);

    % get x1 and x2 contained in Xobs
    n_points = size(Xobs,1) / 2;
    x1 = Xobs(1:n_points); % x
    x1 = reshape(x1, [2,size(x1,1)/2]); % from nx1 to 2x(n/2)
    x2 = Xobs(n_points+1:end); % x'
    x2 = reshape(x2, [2,size(x2,1)/2]); % from nx1 to 2x(n/2)
    
    % transfer error
    % from euclidean coord to homogeneous
%     x1 = [x1 ; ones(1,size(x1,2))];
%     x2 = [x2 ; ones(1,size(x2,2))];
%     diff = l2_dist(euclid(x1)-euclid(inv(H)*x2))+l2_dist(euclid(x2)-euclid(H*x1));
    
    % reprojection error
    xhat = P0(9+1:end);
    xhat = reshape(xhat, [2,size(xhat,1)/2]); % from nx1 to 2x(n/2)
    xhat = [xhat ; ones(1,size(xhat,2))]; % from euclidean to homogeneous
    xhatp = H*xhat;
    diff = l2_dist(x1-euclid(xhat))+l2_dist(x2-euclid(xhatp));

end

function distance = l2_dist(vec)
    distance = (sum(vec.^2));
end
