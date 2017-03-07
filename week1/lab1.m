%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
s = 0.5;
theta = 30;
t1 = 5;
t2 = 10;
H = [s*cosd(theta), -s*sind(theta), t1 ; s*sind(theta), s*cosd(theta), t2; 0, 0, 1];

I2 = apply_H(I, H);
size(I2);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
clear
I=imread('Data/0005_s.png'); % we have to be in the proper folder
A = [0 1; 1 1];
t1 = 2;
t2 = 3;
H = [A(1,1), A(1,2), t1 ; A(2,1), A(2,2), t2; 0, 0, 1]

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
T = [1 0 t1;0 1 t2; 0 0 1];
[U,D,V] = svd(A);
Rtheta = U*V';
Rphi = V';

Rtheta = [Rtheta(1,1), Rtheta(1,2), 0 ; Rtheta(2,1), Rtheta(2,2), 0; 0,0,1];
Rphi = [Rphi(1,1), Rphi(1,2), 0 ; Rphi(2,1), Rphi(2,2), 0; 0,0,1];
Rphit = Rphi';
Rphi_t = [Rphit(1,1), Rphit(1,2), 0 ; Rphit(2,1), Rphit(2,2), 0; 0,0,1];
D = [D(1,1), D(1,2), 0 ; D(2,1), D(2,2), 0 ; 0,0,1];
H_decomp = T*(Rtheta*Rphit*D*Rphi);
a=2;

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
diff = round(H_decomp)-H;
diff = sum(diff(:));
if diff == 0
    disp('matrix are equal')
else
    disp('not equal')
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I2_decomp = apply_H(I, H_decomp);
figure; imshow(I); figure; imshow(uint8(I2_decomp));

diff = I2-I2_decomp;
if diff == 0
    disp('images are equal')
else
    disp('images are not equal')
end



%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

clear
I=imread('Data/0005_s.png'); % we have to be in the proper folder

corners = [1,1,size(I,1),size(I,1);
           1, size(I,2), 1, size(I,2);
           1, 1, 1, 1];
            
target = [1,  120,    size(I,1)-100,size(I,1)-50;
          1,    1300,   1, size(I,2)-200;
          1, 1, 1, 1];
      
% compute homography matrix
points = {};
for i=1:4
    points{i} = [corners(1,i), corners(2,i), 1, 0, 0, 0, -corners(1,i)*target(1,i), -corners(2,i)*target(1,i);
         0,  0, 0, corners(1,i), corners(2,i), 1, -corners(1,i)*target(2,i), -corners(2,i)*target(2,i)];
end

A = [points{1}; points{2}; points{3}; points{4}];
b = [target(1,1), target(2,1), target(1,2), target(2,2), target(1,3), target(2,3), target(1,4), target(2,4)]';

h = A\b;

H = [h(1), h(2), h(3);
    h(4), h(5), h(6);
    h(7), h(8), 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
a = cross(l1, l2);
a = a/a(3,1); % point at infinity obtained from first pair of lines...
b = cross(l3, l4);
b = b/ b(3,1); % point at infinity obtained from second pair of lines...
l = cross(a, b); % required imaged line at infinity...
H = [1 0 0; 0 1 0; l(1, 1)/l(3,1) l(2, 1)/l(3,1) 1];


I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
p1trans=H*p1;
p1trans = p1trans / p1trans(3);
p2trans=H*p2;
p2trans = p2trans / p2trans(3);
p3trans=H*p3;
p3trans = p3trans / p3trans(3);
p4trans=H*p4;
p4trans = p4trans / p4trans(3);
p5trans=H*p5;
p5trans = p5trans / p5trans(3);
p6trans=H*p6;
p6trans = p6trans / p6trans(3);
p7trans=H*p7;
p7trans = p7trans / p7trans(3);
p8trans=H*p8;
p8trans = p8trans / p8trans(3);
lr1 = cross(p1trans,p2trans);
lr2 = cross(p3trans,p4trans);
lr3 = cross(p5trans, p6trans);
lr4 = cross(p7trans, p8trans);

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

%% 5. OPTIONAL: Affine Rectification of the left facade of image 0000

%% 6. OPTIONAL: Metric Rectification of the left facade of image 0000

%% 7. OPTIONAL: Affine Rectification of the left facade of image 0001

%% 8. OPTIONAL: Metric Rectification of the left facade of image 0001


