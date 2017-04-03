function [ disparity ] = stereo_computation( left_im, right_im, min_disp, max_disp, ws, cost )
%STEREO_COMPUTATION Summary of this function goes here
% The input parameters are 5:
%   - left image
%   - right image
%   - minimum disparity
%   - maximum disparity
%   - window size (e.g. a value of 3 indicates a 3x3 window)
%   - matching cost (the user may able to choose between SSD and NCC costs)

    disparity = zeros(size(left_im));    
    [h,w] = size(left_im);
    
    % pad is the distance at each side of the window
    pad = floor(ws/2);
    
    left_im = padarray(left_im, [pad pad]);
    right_im = padarray(right_im, [pad pad]);
    for row=1+pad:h+pad
        for col=1+pad:w+pad
            %[row-pad, row+pad, col-pad, col+pad]
            left_patch = left_im(row-pad:row+pad, col-pad:col+pad);
            %size(left_patch)
            min_col = max(1+pad, col-max_disp);
            max_col = min(w+pad, col+max_disp);
            
            best_ssd = Inf; 
            for i = min_col:max_col
                right_patch = right_im(row-pad:row+pad, i-pad:i+pad);
                ssd = sum(sum( (left_patch-right_patch).^2 ));
                if ssd < best_ssd
                    best_ssd = ssd;
                    best_index = i;
                end
            end
            disparity(row-pad, col-pad) = abs(col-best_index);
        end        
    end
end

