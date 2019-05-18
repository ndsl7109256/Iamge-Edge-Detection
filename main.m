GIm=imread('123.png');

GIm=rgb2gray(GIm);

numofpixels=size(GIm,1)*size(GIm,2);


figure,imshow(GIm);

title('Original Image');

HIm=uint8(zeros(size(GIm,1),size(GIm,2)));

freq=zeros(256,1);

probf=zeros(256,1);

probc=zeros(256,1);

cum=zeros(256,1);

output=zeros(256,1);


%freq counts the occurrence of each pixel value.

%The probability of each occurrence is calculated by probf.


for i=1:size(GIm,1)

    for j=1:size(GIm,2)

        value=GIm(i,j);

        freq(value+1)=freq(value+1)+1;

        probf(value+1)=freq(value+1)/numofpixels;

    end

end


sum=0;

no_bins=255;


%The cumulative distribution probability is calculated. 

for i=1:size(probf)

   sum=sum+freq(i);

   cum(i)=sum;

   probc(i)=cum(i)/numofpixels;

   output(i)=round(probc(i)*no_bins);

end

for i=1:size(GIm,1)

    for j=1:size(GIm,2)

            HIm(i,j)=output(GIm(i,j)+1);

    end

end

figure,imshow(HIm);

title('Histogram equalization');
 
BW1 = edge(GIm,'Sobel')'
figure,imshow(BW1);
title('Sobel');

BW2 = edge(GIm,'Canny');
figure,imshow(BW2);
title('Canny');

BW3 = edge(GIm,'Prewitt');
figure,imshow(BW3);
title('Prewitt');

BW4 = edge(GIm,'Roberts');
figure,imshow(BW4);
title('Roberts');

BW5 = edge(GIm,'log');
figure,imshow(BW5);
title('log');

BW6 = edge(GIm,'zerocross');
figure,imshow(BW6);
title('zerocross');

BW7 = edge(GIm,'approxcanny');
figure,imshow(BW7);
title('approxcanny');
%The result is shown in the form of a table

