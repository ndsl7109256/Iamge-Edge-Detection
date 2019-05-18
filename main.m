clear;

sss=[1 2 3];
sum(sss)

GIm=imread('123.png');

GIm=rgb2gray(GIm);

%I = imnoise(GIm, 'salt & pepper', 0.01);
% median filtering
I = GIm;
[r c] = size(I);
Rep = zeros(r + 2, c + 2);
for x = 2 : r + 1
    for y = 2 : c + 1
        Rep(x,y) = I(x - 1, y - 1);
    end
end
B = zeros(r, c);
for x = 1 : r
    for y = 1 : c
        for i = 1 : 3
            for j = 1 : 3
                q = x - 1;     w = y -1;
                array((i - 1) * 3 + j) = Rep(i + q, j + w);
            end
        end
        B(x, y) = median(array(:));
    end
end
figure, imshow(I);
title('noise');
figure, imshow(uint8(B));
title('solved');
 GIm = I;

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


% A=imread('plane.jpg');
% B=rgb2gray(A);
% C=double(B);
i = 0;
j = 0;
k = 0;

X= double(HIm); 
height = size(X, 1); 
width = size(X, 2); 
channel = size(X, 3);
lenaOutput = X;
Gx = [1 2 1; 0 0 0; -1 -2 -1];
Gy = Gx';
for i = 2 : height-1
   for j = 2 : width-1        
           tempLena = X(i - 1 : i + 1, j - 1 : j + 1);
           aa = Gx.*tempLena;
           xx = 0;
           v= aa(:);
           for k =1:9
           xx = xx + v(k);          
           end
           bb= Gy.*tempLena;
           yy = 0;
           v= bb(:);
           for k =1:9
           yy = yy + v(k);          
           end
           pixValue =sqrt(xx.^2+ yy.^2);
          % pixValue =(x-y);
           lenaOutput(i, j) = pixValue;
   end
end
lenaOutput = uint8(lenaOutput); figure; imshow(abs(lenaOutput),[]); title(' Sobel Edge Detection');

