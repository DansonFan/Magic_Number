% -- Dingxin Fan 2023

clc, clear
I = imread("path");
I1=rgb2gray(I);
p=FastPeakFind(I,55,fspecial("disk",5),3,1); %Adjust this for diff image size
image(I);
hold on
x = p(1:2:end);
y = p(2:2:end);
plot(x,y,'ro')
