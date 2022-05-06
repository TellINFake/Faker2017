I1 = imread('fig_1.jpg');
I2= imread('fig_2.jpg');
I3 = imread('fig_3.jpg');
I4 = imread('fig_4.jpg');

figure();
subplot('Position',[0.02 0.65 0.3 0.3]);
subplot(2,2,1); imshow(I1);
hold on;

subplot('Position',[0.35 0.65 0.3 0.3]);
subplot(2,2,2); imshow(I2);
hold on;

subplot('Position',[0.02 0.3 0.3 0.3]);
subplot(2,2,3); imshow(I3);
hold on;

subplot('Position',[0.35 0.3 0.3 0.3]);
subplot(2,2,4); imshow(I4);
hold on;