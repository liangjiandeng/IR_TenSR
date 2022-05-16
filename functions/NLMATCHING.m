function [aa, Z_block, predenoised_blocks,  bparams]=NLMATCHING(X,para)
bparams.block_sz = [para.patsize, para.patsize];
bparams.overlap_sz=[para.overlap, para.overlap];
[nr, nc,~]=size(X);
bparams.sz=[nr nc];
sz=[nr nc];
step=para.patsize-para.overlap;
sz1=[1:step:sz(1)- bparams.block_sz(1)+1];
 sz1=[sz1 sz(1)- bparams.block_sz(1)+1];
sz2=[1:step:sz(2)- bparams.block_sz(2)+1];
sz2=[sz2 sz(2)- bparams.block_sz(2)+1];
bparams.block_num(1)=length(sz1);
bparams.block_num(2)=length(sz2);

predenoised_blocks = ExtractBlocks1(X, bparams);
Z_block=zeros(bparams.block_sz(1), bparams.block_sz(2),para.k_subspace, bparams.block_num(1)* bparams.block_num(2));
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
if para.K==1
    aa=ones(nr*nc,1);
else
  [aa]=fkmeans(Y2,para.K);
end
predenoised_blocks=permute(predenoised_blocks,[4 3 1 2]);
end