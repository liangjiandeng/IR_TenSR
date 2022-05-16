% IR-TenSR for Hyperspectral image and multispectral image fusion
% Copyright(c) 2022 Ting Xu
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for Hyperspectral image super-
% resolution from a pair of low-resolution hyperspectral image and a high-
% resolution multispectral image.
% 
% if you use this code, Please cite the following paper:
%  Ting Xu, Ting-Zhu Huang, Liang-Jian Deng, and Naoto Yokoya, An Iterative Regularization Method based on Tensor
%  Subspace Representation for Hyperspectral Image Super-Resolution (IR-TenSR),IEEE TGRS, 2022

clear
clc

addpath(genpath('functions'));
addpath(genpath('data'));
addpath(genpath('tSVD'));

data={'indian'};
aa=load('Noi_Indian.mat');
aa=aa.Noi_Indian;
aa=aa(1:128,1:128,:);
aaa(:,:,1:98)=aa(:,:,5:102);
aaa(:,:,99:133)=aa(:,:,113:147);
aaa(:,:,134:182)=aa(:,:,166:214);
cc=double(aaa);
cc=max(cc,0);
x1= max(max(max(cc)));
dd=cc/x1;
S=double(dd);
[M, N, L]=size(S);      
sz=[M N];

spec0 = load( 'AVIRIS_spec.ascii');
spec = spec0(:,1); % Center wavelength of AVIRIS
wavelength = [450 520; 520 600; 630 690; 760 900; 1550 1750; 2080 2350]; % ETM/Landsat
m_band = size(wavelength,1); % Number of MSI bands
F = zeros(m_band,182);

for b=1:6
    b_i = find(spec>wavelength(b,1),1);
    b_e = find(spec<wavelength(b,2),1,'last');
    F(b,b_i:b_e) = 1/(b_e+1-b_i);
end
F=F(:,1:L);

kernel_type     =    {'uniform_blur'};
s0=2;sf=4;
para = Parameters_setting( sf, kernel_type,data, sz,s0);
S_bar = hyperConvert2D(S);
%% simulated LR-HSI with noise
SNRm=25;
hyper= para.H(S_bar);
HSI_clean=hyperConvert3D(hyper,M/sf,N/sf);
sigmah =sqrt(sum(hyper(:).^2)/(10^(SNRm/10))/numel(hyper));
rng(10,'twister')
hyper= hyper+sigmah*randn(size(hyper));
HSI=hyperConvert3D(hyper,M/sf,N/sf);

%% simulated HR-MSI with noise 
Y=F*S_bar;
SNRm=25; 
sigmah1 =sqrt(sum(Y (:).^2)/(10^(SNRm/10))/numel(Y ));
rng(10,'twister')
MSI = hyperConvert3D((Y+sigmah1*randn(size(Y))), M, N);

MSI_clean = hyperConvert3D(Y, M, N);
MSI0=MSI;
HSI0=HSI;
MSI_clean0=MSI_clean;
HSI_clean0=HSI_clean;
%% Proposed IR_TenSR
MSI_h={};HSI_h={};X={};t4=[]; MSI=MSI0;HSI=HSI0;MSI_clean=MSI_clean0;HSI_clean=HSI_clean0;                                                       
for i=1:para.IR   
     MSI_h{i}=MSI_clean;  
     HSI_h{i}=HSI_clean;
     t0=clock;
     X{i} = IR_TenSR(HSI,MSI,F,para);
     t4(i)=etime(clock,t0); 
     S_bar=hyperConvert2D(X{i});
     hyper= para.H(S_bar);
     HSI1=hyperConvert3D(hyper,M/sf,N/sf);
     MSI1 = hyperConvert3D(F*S_bar, M, N);
     HSI=HSI_h{i}-HSI1;
     MSI=MSI_h{i}-MSI1;
     MSI_clean=MSI;
     HSI_clean=HSI; 
 end
    ZZ=0;
 for ii=1:para.IR
      ZZ=ZZ+X{ii};
 end 
 [psnr,rmse, ergas, sam, ssim] = quality_assessment(double(im2uint8(S)), double(im2uint8(ZZ)), 0, 1.0/sf);
 











 