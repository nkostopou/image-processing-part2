%Load image , and convert it to gray-scale
x = imread('im2.jpg');
x = rgb2gray(x);
x = imresize(x, 0.2);

%Gaussian filter
sigma=3.75;
x_s=imgaussfilt(x, sigma);

%edge detector
img=edge(x_s,'sobel'); % exei kalitera apotelesmata sta Hough Lines se antithesi me tous upoloipous edge detectrors
figure
imshow(img,'InitialMagnification', 'fit')

%Hough Function
[h, l, r]= myHoughTransform(img, 1, pi/180, 15);


function [H, L, res] = myHoughTransform(img_binary, Drho, Dtheta, n)

I = img_binary;
[rows, cols] = size(I);

%max distance
rho_maximum = (sqrt(rows^2 + cols^2)) ;
rhoScale = -rho_maximum:Drho:rho_maximum;

%theta to degrees 
thetaDeg = Dtheta * 180/pi;
theta_maximum = 90;
thetaScale = -theta_maximum:thetaDeg:theta_maximum - 1;

H = zeros(length(rhoScale), length(thetaScale));

for row = 1:rows
    for col = 1:cols
        if I(row, col) > 0
            x = col - 1;
            y = row - 1;
            for theta_ = 1 : length(thetaScale)
                tempR = x*cos(thetaScale(theta_)*pi/180) + y*sin(thetaScale(theta_) * pi/180);
                tempR = round((tempR + rho_maximum)/Drho)+1;
                H(tempR,theta_) = H(tempR,theta_) + 1;

            end
        end
    end
end
H
figure
imshow(H,[],'XData',thetaScale,'YData', rhoScale, 'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
title('Hough Matrix');
axis on, axis normal, hold on;
colormap(gca,hot)

%%%PEAKS%%%%

rhos = zeros(n, 1);
thetas = zeros(n, 1); 
H_nms = H; 

% boundary problem
H_padded= padarray(H_nms, [1, 1], 'replicate'); 

% size of hough accumulator
[rows, cols] = size(H_nms); 

% to account for padding
for i = 2:rows-1 
    for j = 2:cols-1
        if any(find((H_padded(i-1:i+1, j-1:j+1) > H_padded(i,j)))) > 0 
              H_nms(i-1,j-1) = 0;
        end
    end
end
for i = 1:(n)
    maxIdx = max(H_nms(:)); 
    [rhoMaxIdx, thetaMaxIdx] = find(H_nms==maxIdx);
    rhos(i) = rhoMaxIdx(1);
    thetas(i) = thetaMaxIdx(1);
    H_nms(rhoMaxIdx(1), thetaMaxIdx(1)) = 0;
end


L=[rhos ,thetas];
L(L(:,1)==0, : )=[];
plot(thetaScale(L(:,2)),rhoScale(L(:,1)),'s','color','blue');


%%%LINES%%%

figure()
x = imread('im2.jpg');
x = imresize(x, 0.2);
imshow(x,'InitialMagnification','fit');
hold on;
size(L,1);
d = [];
for i = 1:size(L,1)
    rhoTemp = rhoScale(L(i,1));
    theTemp = thetaScale(L(i,2));
    if theTemp == 0
        x1 = rhoTemp;
        x2 = rhoTemp;
        y1 = 1;
        y2 = size(I,1);
    else
        x1 = 1;
        x2 = size(I,2);
        y1 = (rhoTemp - x1*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
        y2 = (rhoTemp - x2*cos(theTemp*pi/180)) / sin(theTemp*pi/180);
    end
    d(i) = sqrt((y2-y1)^2+ (x2-x1)^2);
    plot([x1,x2],[y1,y2],'r','LineWidth',2);
    title('Image with hough lines');
end
    res=(size(x,1)*size(x,2))-ceil(sum(d));
    res
    
    
end
