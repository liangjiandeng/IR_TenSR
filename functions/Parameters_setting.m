function   para   =  Parameters_setting( sf, kernel_type,data,sz,s0 )
if strcmp(kernel_type, 'uniform_blur')
    psf        =    ones(sf)/(sf^2);
elseif strcmp(kernel_type, 'Gaussian_blur')
    psf        =    fspecial('gaussian',2*sf+1,2);
end
 para.fft_B      =    psf2otf(psf,sz);
 para.B     =    ifft2(para.fft_B) ;
 para.fft_BT     =    conj(para.fft_B);
 para.H          =    @(z)H_z(z, para.fft_B, sf, sz,s0 );
 para.HT         =    @(y)HT_y(y, para.fft_BT, sf, sz,s0);

if strcmp(data, 'pavia')
para.beta=5*1e-7;para.mu=1e-4;para.gama=1e-1;para.lamada=5*1e-3;para.k_subspace=3;
para.K=2; para.patsize=21;para.Int_Iter=7;para.Out_Iter=[17];para.sf=sf;para.overlap = floor(para.patsize/2);para.IR=6;  
elseif strcmp(data, 'indian')
para.beta=1e-5;para.mu=5*1e-7;para.gama=1e-5;para.lamada=1e-7;para.k_subspace=8;
para.K=32; para.patsize=9;para.Int_Iter=2;para.Out_Iter=[100];para.sf=sf;para.overlap = floor(para.patsize/2);para.IR=7;
end



   
  




    