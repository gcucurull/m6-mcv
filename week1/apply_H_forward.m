function [ out_im ] = apply_H_forward( im, T)
% Receives an image 'im' and a matrix transofrmation 'T'
    
%%% 1. %%% compute output image size
    %T = inv(T);
    corners = [1,1,size(im,1),size(im,1);
                1, size(im,2), 1, size(im,2);
                1, 1, 1, 1];
            
    t_corners = T*corners;
    % normalize
    t_corners(:,1) = t_corners(:,1) / t_corners(3,1);
    t_corners(:,2) = t_corners(:,2) / t_corners(3,2);
    t_corners(:,3) = t_corners(:,3) / t_corners(3,3);
    t_corners(:,4) = t_corners(:,4) / t_corners(3,4);
    t_corners = round(t_corners);
    
    % compute distances between transformer corner points
    i_dist = abs(min(t_corners(1,:))-max(t_corners(1,:))); % first row are i coordinates
    j_dist = abs(min(t_corners(2,:))-max(t_corners(2,:))); % second row are j coordinates
    
    % compute displacement
    i_disp = min(t_corners(1,:));
    j_disp = min(t_corners(2,:));
            
    out_im=zeros(i_dist, j_dist, 3);
    %T = inv(T);
        
%%% 2. %%% Compute homography
    for i=1:size(im,1)
        for j=1:size(im,2)
            v1 = [i;j;1]; % point in current image
            v2 = T*v1;
            v3 =v2/v2(3,1);
            v3 = round(v3); % point in target image
            if (v3(1,1)-i_disp+1<=0 || v3(2,1)-j_disp+1<=0 || v3(1,1)-i_disp+1 > size(out_im,1) || v3(2,1)-j_disp+1 > size(out_im,2))
               continue;
            else
                % copy from source to target
                %out_im(i-i_disp+1,j-j_disp+1,:) = im(v3(1,1),v3(2,1),:);
                out_im(v3(1,1)-i_disp+1,v3(2,1)-j_disp+1,:) = im(i,j,:);
            end
        end
    end
    
    %out_im = to;

end

