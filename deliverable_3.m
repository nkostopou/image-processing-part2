I = imread('im2.jpg');
N = 10;
[rows, columns, numColorChannels] = size(I);
numberOfRows = round(rows/N);
numberOfColumns = round(columns/N);
image = imresize(I, [numberOfRows numberOfColumns]);

an = 213*(pi/180);
myrotateImage = myImgRotation(image, an);

function rotImg = myImgRotation(Im, angle)
[m,n,p]=size(Im);
mm = sqrt(m^2+n^2);
nn = sqrt(m^2+n^2);
for t=1:mm
   for s=1:nn
      i = uint16((t-mm/2)*cos(angle)+(s-nn/2)*sin(angle)+m/2);
      j = uint16(-(t-mm/2)*sin(angle)+(s-nn/2)*cos(angle)+n/2);
      if i>0 && j>0 && i<=m && j<=n           
         Irotate(t,s,:)=Im(i,j,:);
      end
   end
end
rotImg = Irotate;
figure;
imshow(Irotate,'InitialMagnification','fit');

end