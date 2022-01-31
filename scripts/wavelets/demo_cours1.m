close all
clear all

fid = fopen('foreman.qcif','r');
[compY1,compU1,compV1]=yuv_readimage(fid,'qcif');
[compY2,compU2,compV2]=yuv_readimage(fid,'qcif');

figure(1)
subplot(121)
imshow(uint8(compY1))
subplot(122)
imshow(uint8(compY2))

pause

figure(2)
subplot(121)
imshow(uint8(compY1))
subplot(122)
imagesc(compY2-compY1), axis square

pause

figure(3)
compY1 = double(compY1);
It = (compY1 -mean(compY1 (:)))/256;
Ito = wt2d(It,'haar',1);
subplot(121), imagesc(It),colormap(gray), axis square
subplot(122), imagesc(Ito), axis square

pause

figure(3)
subplot(121), imagesc(It),colormap(gray), axis square
subplot(122), imagesc(Ito(1:end/2,1:end/2)), axis square

pause
[y,fs] = wavread('tada.wav');
code='haar'; stages=1;
out = wt(y',code,stages);
t=(0:length(y)-1)/fs;
figure(4)
subplot(211)
plot(t,y);
subplot(212)
plot(t,out);

pause
sound(y,fs)

pause
sound(out(1:end/2),fs/2)