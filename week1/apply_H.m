function [ out_im ] = apply_H( im, T)
% Receives an image 'im' and a matrix transofrmation 'T'
    
%%% 1. %%% compute output image size
%     corners = [1,1,size(im,1),size(im,1);
%                 1, size(im,2), 1, size(im,2);
%                 1, 1, 1, 1];
            
    corners = [1,   size(im,2),              1,             size(im,2);
                1,  1,     size(im,1),                      size(im,1);
                1, 1, 1, 1];
            
    t_corners = T*corners;
    % normalize
    t_corners(:,1) = t_corners(:,1) / t_corners(3,1);
    t_corners(:,2) = t_corners(:,2) / t_corners(3,2);
    t_corners(:,3) = t_corners(:,3) / t_corners(3,3);
    t_corners(:,4) = t_corners(:,4) / t_corners(3,4);
    t_corners = round(t_corners);
    
    % compute distances between transformer corner points
    x_dist = abs(min(t_corners(1,:))-max(t_corners(1,:))); % first row are x coordinates
    y_dist = abs(min(t_corners(2,:))-max(t_corners(2,:))); % second row are j coordinates
    
    % compute displacement
    x_disp = min(t_corners(1,:));
    y_disp = min(t_corners(2,:));
            
    out_im=zeros(y_dist, x_dist, 3);
    T = inv(T);
    
    from = im;
    to = out_im;
        
%%% 2. %%% Compute homography
    for i=y_disp:size(out_im,1)
        for j=x_disp:size(out_im,2)
            v1 = [j;i;1]; % point in target image
            v2 = T*v1;
            v3 =v2/v2(3,1);
            v3 = round(v3); % point in source image
            if (v3(1,1)<=0 || v3(2,1)<=0 || v3(2,1) > size(im,1) || v3(1,1) > size(im,2))
               continue;
            else
                % copy from source to target
                to(i-y_disp+1,j-x_disp+1,:) = from(v3(2,1),v3(1,1),:);
            end
        end
    end
    
    out_im = to;

end

