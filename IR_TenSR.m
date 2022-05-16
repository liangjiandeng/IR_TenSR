function X=IR_TenSR(HSI,MSI,R,para)
L=size(R,2);
FBm=para.fft_B;
FBmC  = conj(FBm);
FBs  = repmat(FBm,[1 1 L]);
FBCs=repmat(FBmC,[1 1 L]);   
[nr, nc,~]=size(MSI);
n_dr=nr/para.sf;
n_dc=nc/para.sf;

%intialize dictionary B and coefficient C by t-svd
HSI_up1=imresize(HSI, para.sf,'bicubic');
HSI_up1_twist=permute(HSI_up1,[1,3,2]);
[U,~,~]=tsvd(HSI_up1_twist); 
B= U(:,1:para.k_subspace,:);%intialize B
C3= tprod(tran(B),HSI_up1_twist);%intialize C0 PAM
       
HSI_UP=zeros(nr,nc,L);
HSI_UP(1:para.sf:end,1:para.sf:end,:)=HSI;
HHH=ifft2((fft2(HSI_UP).*FBCs));
HSI_BS=hyperConvert2D(HHH);
MSI3=Unfold(MSI,size(MSI),3);
S3_fixed=HSI_BS+R'*MSI3;

S1=R'*R+(para.mu+para.beta)*eye(L);
[P,Sigma]=eig(S1);
Sigma=reshape(diag(Sigma),[1 1 L]);
InvLbd=1./repmat(Sigma,[ para.sf*n_dr  para.sf*n_dc 1]);
B2Sum=PPlus(abs(FBs).^2./( para.sf^2),n_dr,n_dc);
InvDI=1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Sigma,[n_dr n_dc 1]));
A=ipermute(tprod(B,C3),[1,3,2]);
X_new=HSI_up1_twist;

for iter=1:para.Out_Iter
    X_old=X_new;
    C_last=C3;
%% update A
     T_BC=ipermute(tprod(B,C3),[1,3,2]);
     T_BC2=Unfold(T_BC,size(T_BC),3);
     A2= T_BC2;
     A11=Unfold(A,size(A),3);
     S3_change=para.mu*A2+para.beta*A11;
     S3=S3_fixed+S3_change;
     S30=fft2(reshape((P\S3)',[nr nc L])).*InvLbd;
     temp  = PPlus_s(S30/(para.sf^2).*FBs,n_dr,n_dc);
     invQUF = S30-repmat(temp.*InvDI,[para.sf para.sf 1]).*FBCs;
     VXF = P*reshape(invQUF,[nr*nc L])';
     A3= reshape(real(ifft2(reshape(VXF',[nr nc L]))),[nr*nc L])';
     A=hyperConvert3D(A3,nr,nc);
%% Update B 
    [U2,~,V2] = tsvd(tprod(C3,tran(permute(A,[1,3,2])))+(para.beta/para.mu)*tran(B));
    B= tprod(V2,tran(U2));
%% Update C by ADMM
     C0=permute((para.mu*tprod(tran(B),permute(A,[1,3,2]))+para.beta*C3)/(para.mu+para.beta),[3,2,1]);
     O=zeros(size(C0));
     C=C0;
    for i=1:para.Int_Iter
        N0=C-O/para.gama;
        [aa, Z_block, predenoised_blocks,  bparams]=NLMATCHING(N0,para);
        for mn=1:max(aa)
            gg=find(aa==mn);
            XES=predenoised_blocks(gg,:,:,:);
            [a, b, c, d ]=size(XES);
            XES = reshape(XES,[a b c*d]);
            [N11] = prox_TNN(XES,para.lamada/(para.gama*(para.mu+para.beta)));
            N11=reshape(N11,[a b c d]); 
            Z_block(:,:,:,gg)=permute(N11,[3 4 2 1]);
         end
         N1= JointBlocks2(Z_block, bparams);  
         C=(C0+para.gama*N1+O)/(1+para.gama);   
         O=O+para.gama*(N1-C);
   end
  C3=ipermute(C,[3,2,1]); 
%% Update reconstructed X
     X_new = tprod(B,C3); 
     Rel_Cha = norm(C3(:)-C_last(:))/norm(C_last(:));
     if Rel_Cha < 0.01      
        break;
     end  
end
X=ipermute(X_new,[1,3,2]);
end
