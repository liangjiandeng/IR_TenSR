% IRTenSR for Hyperspectral image and multispectral image fusion
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
%  Ting Xu, Ting-Zhu Huang, Liang-Jian Deng, and Naoto Yokoya,  An Iterative Regularization Method based on Tensor
%  Subspace Representation for Hyperspectral Image Super-Resolution (IR-TenSR),IEEE TGRS, 2022

clear
clc

addpath(genpath('functions'));
addpath(genpath('data'));
addpath(genpath('tSVD'));

data={'pavia'};
S=imread('.\data\original_rosis.tif');
S=S(1:256,1:256,11:end);
S=double(S);
S=S/max(S(:));
[M, N, L]=size(S);

sz=[M N];
F=load('.\data\R.mat');
F=F.R;
F=F(:,11:end);
for band = 1:size(F,1)
     div = sum(F(band,:));
     for i = 1:size(F,2)
          F(band,i) = F(band,i)/div;
     end
end
kernel_type     =    {'uniform_blur'};
s0=2;sf=16;
para = Parameters_setting( sf, kernel_type,data, sz,s0);
S_bar = hyperConvert2D(S);
%% simulated LR-HSI with noise
SNRm=30;
hyper= para.H(S_bar);
HSI_clean=hyperConvert3D(hyper,M/sf,N/sf);
sigmah =sqrt(sum(hyper(:).^2)/(10^(SNRm/10))/numel(hyper));
rng(10,'twister')
hyper= hyper+sigmah*randn(size(hyper));
HSI=hyperConvert3D(hyper,M/sf,N/sf);

%% simulated HR-MSI with noise 
Y=F*S_bar;
SNRm=30; 
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
 











 