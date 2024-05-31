clear all; close all; clc;

load ESAR.mat;
I = I(1:50,1:50,:);
fsize = '3';
for ii=1:3
    for jj=ii:3
        covariance(ii,jj,:,:) = imfilter(I(:,:,ii).*conj(I(:,:,jj)), fspecial('average',str2num(fsize)),'replicate');
        covariance(jj,ii,:,:) = conj(covariance(ii,jj,:,:));
    end
end

I_filtered = Stc_McSAR(covariance,3,21,1);

I_in = Pauli_C(I);
I_out = Pauli_C(I_filtered_1);

figure();
ax(1)=subplot(131); imshow(I_in);
ax(2)=subplot(132); imshow(I_out);
linkaxes(ax,'xy');