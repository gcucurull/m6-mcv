% Vanishing Points detection algorithm as
% featured in Finding Vanishing Points via Point Alignments in Image 
% Primal and Dual Domains Jose Lezama, Rafael Grompone von Gioi, 
% Gregory Randall and Jean-Michel Morel. In CVPR 2014.
% 
% This version uses the algorithm by Figueiredo and Jain, Unsupervised 
% learning of finite mixture models, to quickly obtain cluster candidates.
%
%
% Copyright (c) 2013-2014 Jose Lezama <jlezama@gmail.com>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
clear
close all

img_in =  'test.jpg'; % input image
folder_out = '.'; % output folder

manhattan = 1;
acceleration = 0;

focal_ratio = 1;

params.PRINT = 1;
params.PLOT = 1;

horizon = detect_vps(img_in, folder_out, manhattan, acceleration, focal_ratio, params);
