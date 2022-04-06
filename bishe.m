clc;clear;close all;
%%
% QPSK函数调制，成功加入卷积码函数，交织函数,成功加入RGB三通道，恢复彩色
% https://zhuanlan.zhihu.com/p/139381223     交织
%% 参数
M=4;
k = log2(M);
L=7;                %卷积码约束长度
tblen=6*L;          %Viterbi译码器回溯深度

%% 彩色图-灰度图
img_jpg = imread('lena.jpg');
figure(1)
imshow(img_jpg)
title('原图')
                            %img_bmp = rgb2gray(img_jpg);%转灰度图
                            %figure(2)
                            %imshow(img_bmp)
                            %title('灰度图')
R = img_jpg(:,:,1);
G = img_jpg(:,:,2);
B = img_jpg(:,:,3);

figure(2)
imshow(R)
title('R通道');
%imwrite(R,'LenaR.jpg');
figure(3)
imshow(G)
title('G通道');
%imwrite(G,'LenaG.jpg');
figure(4)
imshow(B)
title('B通道');

[a,b] = size(R);
%% 灰度图-比特
img_bmp_double_R=im2double(R);%double类型
img_bit_R=de2bi(img_bmp_double_R(:).*255);%比特
img_bit_list_R=img_bit_R(:).';

img_bmp_double_G=im2double(G);%double类型
img_bit_G=de2bi(img_bmp_double_G(:).*255);%比特
img_bit_list_G=img_bit_G(:).';

img_bmp_double_B=im2double(B);%double类型
img_bit_B=de2bi(img_bmp_double_B(:).*255);%比特
img_bit_list_B=img_bit_B(:).';
%% 卷积码编码
trellis = poly2trellis(7,[171 133]);
code_data_R=convenc(img_bit_list_R,trellis);
code_data_G=convenc(img_bit_list_G,trellis);
code_data_B=convenc(img_bit_list_B,trellis);
%% 交织
date_scramble_R = matintrlv(code_data_R,k,length(code_data_R)/k);
date_scramble_G = matintrlv(code_data_G,k,length(code_data_G)/k);
date_scramble_B = matintrlv(code_data_B,k,length(code_data_B)/k);

date_scramble_jz = matintrlv(img_bit_list_R,k,length(img_bit_list_R)/k);%不经过卷积码编码，直接交织
%% QPSK调制
A_R = reshape(date_scramble_R,2,length(date_scramble_R)/2)';
img_bit_tx_R = bi2de(A_R)';%2进制转10进制
img_bit_tx_PSK_R = pskmod(img_bit_tx_R,M,pi/M);
img_bit_tx_ofdm_R = ifft(img_bit_tx_PSK_R);

A_wcl = reshape(img_bit_list_R,2,length(img_bit_list_R)/2)';%不经过卷积、交织处理
img_bit_tx_wcl = bi2de(A_wcl)';%2进制转10进制
img_bit_tx_PSK_wcl = pskmod(img_bit_tx_wcl,M,pi/M);
img_bit_tx_ofdm_wcl = ifft(img_bit_tx_PSK_wcl);

A_jj = reshape(code_data_R,2,length(code_data_R)/2)';%只经过卷积编码处理
img_bit_tx_jj = bi2de(A_jj)';%2进制转10进制
img_bit_tx_PSK_jj = pskmod(img_bit_tx_jj,M,pi/M);
img_bit_tx_ofdm_jj = ifft(img_bit_tx_PSK_jj);

A_jz = reshape(date_scramble_jz,2,length(date_scramble_jz)/2)';%只经过交织处理
img_bit_tx_jz = bi2de(A_jz)';%2进制转10进制
img_bit_tx_PSK_jz = pskmod(img_bit_tx_jz,M,pi/M);
img_bit_tx_ofdm_jz = ifft(img_bit_tx_PSK_jz);








A_G = reshape(date_scramble_G,2,length(date_scramble_G)/2)';
img_bit_tx_G = bi2de(A_G)';%2进制转10进制
img_bit_tx_PSK_G = pskmod(img_bit_tx_G,M,pi/M);
img_bit_tx_ofdm_G = ifft(img_bit_tx_PSK_G);

A_B = reshape(date_scramble_B,2,length(date_scramble_B)/2)';
img_bit_tx_B = bi2de(A_B)';%2进制转10进制
img_bit_tx_PSK_B = pskmod(img_bit_tx_B,M,pi/M);
img_bit_tx_ofdm_B = ifft(img_bit_tx_PSK_B);
%% 
EbNo = 0:3;
SNR = EbNo + k;
SNR1=10.^(SNR/10);

img_bit_deqpsk_R = zeros(1,length(img_bit_list_R));
img_bit_deqpsk_G = zeros(1,length(img_bit_list_G));
img_bit_deqpsk_B = zeros(1,length(img_bit_list_B));
img_bit_channel_R_wcl = zeros(1,length(img_bit_list_R));
for i = 1:length(EbNo)
%% 通过AWGN信道，并fft。
    img_bit_channel_R = awgn(img_bit_tx_ofdm_R,SNR(i),'measured');
    img_bit_rx_PSK_R = fft(img_bit_channel_R);
    
    wcl = awgn(img_bit_tx_ofdm_wcl,SNR(i),'measured');%不经过任何处理通过信道
    img_bit_rx_PSK_wcl = fft(wcl);

    jj = awgn(img_bit_tx_ofdm_jj,SNR(i),'measured');%只经过卷积码编码处理通过信道
    img_bit_rx_PSK_jj = fft(jj);
    
    jz = awgn(img_bit_tx_ofdm_jz,SNR(i),'measured');%只经过交织编码处理通过信道
    img_bit_rx_PSK_jz = fft(jz);   
    
    
    
    
    
    img_bit_channel_G = awgn(img_bit_tx_ofdm_G,SNR(i),'measured');
    img_bit_rx_PSK_G = fft(img_bit_channel_G);
    
    img_bit_channel_B = awgn(img_bit_tx_ofdm_B,SNR(i),'measured');
    img_bit_rx_PSK_B = fft(img_bit_channel_B);
%% QPSK解调

demodulation_data_R = pskdemod(img_bit_rx_PSK_R,M,pi/M);   
De_data1_R = reshape(demodulation_data_R,[],1);
De_data2_R = de2bi(De_data1_R);
De_Bit_R = reshape(De_data2_R',1,[]);

demodulation_data_wcl = pskdemod(img_bit_rx_PSK_wcl,M,pi/M);%未处理   
De_data1_wcl = reshape(demodulation_data_wcl,[],1);
De_data2_wcl = de2bi(De_data1_wcl);
De_Bit_wcl = reshape(De_data2_wcl',1,[]);

demodulation_data_jj = pskdemod(img_bit_rx_PSK_jj,M,pi/M);%只经过卷积编码处理   
De_data1_jj = reshape(demodulation_data_jj,[],1);
De_data2_jj = de2bi(De_data1_jj);
De_Bit_jj = reshape(De_data2_jj',1,[]);

demodulation_data_jz = pskdemod(img_bit_rx_PSK_jz,M,pi/M);%只经过交织编码处理   
De_data1_jz = reshape(demodulation_data_jz,[],1);
De_data2_jz = de2bi(De_data1_jz);
De_Bit_jz = reshape(De_data2_jz',1,[]);









demodulation_data_G = pskdemod(img_bit_rx_PSK_G,M,pi/M);   
De_data1_G = reshape(demodulation_data_G,[],1);
De_data2_G = de2bi(De_data1_G);
De_Bit_G = reshape(De_data2_G',1,[]);

demodulation_data_B = pskdemod(img_bit_rx_PSK_B,M,pi/M);   
De_data1_B = reshape(demodulation_data_B,[],1);
De_data2_B = de2bi(De_data1_B);
De_Bit_B = reshape(De_data2_B',1,[]);
%% 解交织
rx_date_jiejiaozhi_R = matdeintrlv(De_Bit_R,k,length(De_Bit_R)/k);
rx_date_jiejiaozhi_G = matdeintrlv(De_Bit_G,k,length(De_Bit_G)/k);
rx_date_jiejiaozhi_B = matdeintrlv(De_Bit_B,k,length(De_Bit_B)/k);

rx_date_jiejiaozhi_jz = matdeintrlv(De_Bit_jz,k,length(De_Bit_jz)/k);%只进行交织、解交织
%% 信道译码（维特比译码）
  trellis = poly2trellis(7,[171 133]);
  rx_c_de_R = vitdec(rx_date_jiejiaozhi_R,trellis,tblen,'trunc','hard');   %硬判决
  rx_c_de_G = vitdec(rx_date_jiejiaozhi_G,trellis,tblen,'trunc','hard');   %硬判决
  rx_c_de_B = vitdec(rx_date_jiejiaozhi_B,trellis,tblen,'trunc','hard');   %硬判决
  
  rx_c_de_jj = vitdec(De_Bit_jj,trellis,tblen,'trunc','hard');   %只进行卷积码编码、译码
%% 通过信道后解调的比特流还原灰度图

    img_bit_rx_R=reshape(rx_c_de_R,[],8);%8列
    img_dec_R=reshape(bi2de(img_bit_rx_R),a,b);
    img_bmp_recover_R=uint8(img_dec_R);

   
    img_bit_rx_G=reshape(rx_c_de_G,[],8);%8列
    img_dec_G=reshape(bi2de(img_bit_rx_G),a,b);
    img_bmp_recover_G=uint8(img_dec_G);

   
    
    img_bit_rx_B=reshape(rx_c_de_B,[],8);%8列
    img_dec_B=reshape(bi2de(img_bit_rx_B),a,b);
    img_bmp_recover_B=uint8(img_dec_B);
    
    rgb_image = cat(3,img_bmp_recover_R,img_bmp_recover_G,img_bmp_recover_B);
    figure(5)
    subplot(1,length(EbNo)/1,i);imshow(rgb_image)
    snr = num2str(SNR(i));
    tit3=strcat('snr=',snr,'dB');
    title(tit3);

%% 计算误码率
    [errorBit,BER(i)] = biterr(img_bit_list_R, rx_c_de_R,log2(M));  %计算BER
    [errorSym,SER(i)] = symerr(img_bit_list_R, rx_c_de_R);          %计算SER
    
    [errorBit1,BER1(i)] = biterr(img_bit_list_R,De_Bit_wcl,log2(M));  %计算未经过处理BER
    [errorSym1,SER1(i)] = symerr(img_bit_list_R,De_Bit_wcl);          %计算未经过处理SER
   
    [errorBit2,BER2(i)] = biterr(img_bit_list_R,rx_c_de_jj,log2(M));  %计算只经过卷积码编码处理BER
    [errorSym2,SER2(i)] = symerr(img_bit_list_R,rx_c_de_jj);          %计算只经过卷积码编码处理SER
    
    [errorBit3,BER3(i)] = biterr(img_bit_list_R,rx_date_jiejiaozhi_jz,log2(M));  %计算只经过交织编码处理BER
    [errorSym3,SER3(i)] = symerr(img_bit_list_R,rx_date_jiejiaozhi_jz);          %计算只经过交织编码处理SER   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[errorBit_G,BER_G(i)] = biterr(img_bit_list_G, rx_c_de_G,log2(M));  %计算BER
    %[errorSym_G,SER_G(i)] = symerr(img_bit_list_G, rx_c_de_G);          %计算SER
    %
    %[errorBit_B,BER_B(i)] = biterr(img_bit_list_G, rx_c_de_G,log2(M));  %计算BER
    %[errorSym_B,SER_B(i)] = symerr(img_bit_list_G, rx_c_de_G);          %计算SER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
ser1=2*qfunc(sqrt(2*SNR1)*sin(pi/M));%理论误码率
%ser1 = serawgn(EbNo,'psk',M,'nondiff');%没有这个函数
%ber1=1/log2(M)*ser1;                 %理论误比特率
ber1 = berawgn(EbNo,'psk',M,'nondiff');%理论误比特率

rgb_image = cat(3,img_bmp_recover_R,img_bmp_recover_G,img_bmp_recover_B);
figure(6)
imshow(rgb_image)
title('恢复后的彩色图')

scatterplot(img_bit_rx_PSK_R);

figure(8)
semilogy(EbNo,BER1,"-bo", EbNo,SER1,"-b*",EbNo,ser1,"-ko",EbNo,ber1,"-k*")
legend("未经过处理通过AWGN信道BER","未经过处理通过AWGN信道SER","理论误码率","理论误比特率")
title("QPSK在AWGN信道下的性能")
xlabel("信噪比（dB）")
ylabel("误符号率和误码率")

figure(9)
%semilogy( EbNo,SER,"-r*",  EbNo,SER2,"-g*", EbNo,SER3,"-m*",EbNo,ser1,"-ko")
%legend("卷积编码+交织AWGN信道SER","卷积编码AWGN信道SER","交织编码AWGN信道SER","理论误码率")
semilogy( EbNo,BER,"-r*",  EbNo,BER2,"-g*", EbNo,BER3,"-m*",EbNo,ber1,"-ko")
legend("卷积编码+交织AWGN信道BER","卷积编码AWGN信道BER","交织编码AWGN信道BER","理论误码率")
title("QPSK在AWGN信道下的性能")
xlabel("信噪比（dB）")
ylabel("误码率")



