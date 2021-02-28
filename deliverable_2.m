%Load image , and convert it to gray-scale
I = imread('im2.jpg');
im = rgb2gray(I);
imshow(im,'InitialMagnification', 'fit');

%resize
N=10;
[rows, columns, numColorChannels] = size(im);
numOutputRows = round(rows/N);
numOutputColumns = round(columns/N);
im = imresize(im, [numOutputRows, numOutputColumns]);
figure

%Harris function
c = myDetectHarrisFeatures(im);

 function corners = myDetectHarrisFeatures(I) 

% create X and Y Sobel filters
h = [-1 0 1; -1 0 1; -1 0 1];
v = h';

% using imfilter to get our gradient in each direction
filtered_x = filter2(h, I);
filtered_y = filter2(v, I);

% store the values in our output variables, for clarity
Ix = filtered_x;
Iy = filtered_y;
sigma=1;
f = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
Ix2 = filter2(f,Ix.^2);
Iy2 = filter2(f,Iy.^2);
Ixy = filter2(f,Ix.*Iy);

%set empirical constant between 0.04-0.06
k = 0.04;

num_rows = size(I,1);
num_cols = size(I,2);

% create a matrix to hold the Harris values
R = zeros(num_rows, num_cols);
Rmax=0;

for i=1:num_rows
    for j=1:num_cols
        M=([Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)]);
        R(i,j)=det(M)-k*(trace(M))^2;
        if(R(i,j)>Rmax)
            Rmax=R(i,j);
        end
    end
end

for i = 2:num_rows-1
    for j = 2:num_cols-1
        %threshold=0.005
        if(R(i,j)>0.005*Rmax &&R(i,j)> R(i-1,j-1) && R(i,j) > R(i-1,j) && R(i,j) > R(i-1,j+1) && R(i,j) > R(i,j-1) && R(i,j) > R(i,j+1) && R(i,j) > R(i+1,j-1) && R(i,j) > R(i+1,j) &&(( R(i,j) > R(i+1,j+1))|| R(i,j)<0))
        result(i,j) = 1;
        end
    end
end
[rows, cols] = find(result==1);

corners=[rows, cols];
imshow(I,'InitialMagnification', 'fit') 
hold on
plot(cols,rows,'rs','MarkerSize',5);

 end