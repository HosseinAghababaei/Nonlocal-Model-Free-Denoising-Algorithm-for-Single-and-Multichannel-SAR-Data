function [FC, alpha] = Stc_McSAR(Cint, Ptch_size, win_size, Look)
% Cint, initial covariance matrix
% Ptch_size is the size of local window, 
% win_size is the size of search window
% look represents the number of look 

scale = norm(mean(mean(Cint,4),3));
Cint = Cint./scale;

%__________________________________________________________________________
%% Parameter Setting 
[~,N,Lazi,Lrng] = size(Cint);

Edge = linspace(0,100,1500);
[fdW1, edgeW1,fdW2, edgeW2] = weight_fcn(Edge,Look,N);
[~, indx_cmaxf1] = max(fdW1) ;
[~, indx_cmaxf2] = max(fdW2) ;

[mask(:,:,1), mask(:,:,2), mask(:,:,3),mask(:,:,4)] = mask_window(win_size);
Lmsk = length(find(mask(:,:,4)==1)) ;
Lrpc = 12;

%__________________________________________________________________________
%% Estimation of Tr and Inv_cov

Tr_C = zeros(Lazi,Lrng);
Inv_C = zeros(N,N,Lazi,Lrng);
% In = eye(N);
for ii=1:Lazi
    for jj=1:Lrng
        Tr_C(ii,jj) = abs(trace(Cint(:,:,ii,jj)))/N;
        Inv_C(:,:,ii,jj) = pinv(Cint(:,:,ii,jj),1e-6);
%         Inv_C(:,:,ii,jj) = In/Cint(:,:,ii,jj);
    end
end

%__________________________________________________________________________
%% Start

dispN = round(Lazi*Lrng/10); 
[ii,jj] = meshgrid(1:Lazi,1:Lrng);
ii = ii(:);   
jj = jj(:);

%__________________________________________________________________________
%% Filtering Step

Cv = zeros(N,N,Lazi*Lrng);
alpha = zeros(Lazi*Lrng,1);

parfor kk=1:Lazi*Lrng

    %______________________________________________________________________
    % Patch generation around central pixel (p)
    xx_1 = max(1,ii(kk)-Ptch_size);    xx_2 = min(Lazi,ii(kk)+Ptch_size);
    yy_1 = max(1,jj(kk)-Ptch_size);    yy_2 = min(Lrng,jj(kk)+Ptch_size);
    
    %______________________________________________________________________
    % Window generation around the pixel p 
    xp_1 = max(1,ii(kk)-win_size);     xp_2 = min(Lazi,ii(kk)+win_size);
    yp_1 = max(1,jj(kk)-win_size);     yp_2 = min(Lrng,jj(kk)+win_size);
    
    
    
    Xp = Cint(:,:,xp_1:xp_2,yp_1:yp_2);
    Xp_inv = Inv_C(:,:,xp_1:xp_2,yp_1:yp_2);
    tr_p = Tr_C(xp_1:xp_2,yp_1:yp_2);
    
    wx = size(Xp,3);         wy = size(Xp,4);
    df1 = (2*win_size+1)-wx;    df2 = (2*win_size+1)-wy;
    
    cnt_p = Cint(:,:,ii(kk),jj(kk));
    Xp = cat(4,Xp,repmat(cnt_p,1,1,size(Xp,3),df2));
    Xp = cat(3,Xp,repmat(cnt_p,1,1,df1,size(Xp,4)));
    
    cnt_pinv = Inv_C(:,:,ii(kk),jj(kk));
    Xp_inv = cat(4,Xp_inv,repmat(cnt_pinv,1,1,size(Xp_inv,3),df2));
    Xp_inv = cat(3,Xp_inv,repmat(cnt_pinv,1,1,df1,size(Xp_inv,4)));
    
    cnt_trp = Tr_C(ii(kk),jj(kk));  
    tr_p = [tr_p repmat(repmat(cnt_trp,wx,1),1,df2)];
    tr_p = [tr_p;repmat(repmat(cnt_trp,1,(2*win_size+1)),df1,1)];
    cnt_trp = repmat(cnt_trp,Lrpc,1);  
    
    %______________________________________________________________________
    % Similarity of pixel p and pixel q in the patch
    Edgex = linspace(min(tr_p(:)),max(tr_p(:)),1500);
    
    Hbin = histc([tr_p(mask(:,:,1));cnt_trp],Edgex);
    Fx11 = cumsum(Hbin); Fx11 = Fx11./Fx11(end);
    Hbin = histc([tr_p(mask(:,:,2));cnt_trp],Edgex);
    Fx12 = cumsum(Hbin); Fx12 = Fx12./Fx12(end);
    Hbin = histc([tr_p(mask(:,:,3));cnt_trp],Edgex);
    Fx13 = cumsum(Hbin); Fx13 = Fx13./Fx13(end);
    Hbin = histc([tr_p(mask(:,:,4));cnt_trp],Edgex);
    Fx14 = cumsum(Hbin); Fx14 = Fx14./Fx14(end);
    
    Lx = xx_2-xx_1+1;
    Ly = yy_2-yy_1+1;
    Lpq = Lx*Ly;
    D_1 = zeros(Lx,Ly);
    D_2 = D_1*0;
    [aa,bb] = meshgrid(xx_1:xx_2,yy_1:yy_2); 
    aa = aa(:);
    bb = bb(:);
    
    for pp =1:Lpq
        
        % Window generation around the pixel q 
        xq_1 = max(1,aa(pp)-win_size);     xq_2 = min(Lazi,aa(pp)+win_size);
        yq_1 = max(1,bb(pp)-win_size);     yq_2 = min(Lrng,bb(pp)+win_size);
        
        Xq = Cint(:,:,xq_1:xq_2,yq_1:yq_2);
        Xq_inv = Inv_C(:,:,xq_1:xq_2,yq_1:yq_2);
        tr_q = Tr_C(xq_1:xq_2,yq_1:yq_2);
        
        wx = size(Xq,3);         wy = size(Xq,4);
        df1 = (2*win_size+1)-wx;    df2 = (2*win_size+1)-wy;
        
        cnt_q = Cint(:,:,aa(pp),bb(pp));
        Xq = cat(4,Xq,repmat(cnt_q,1,1,size(Xq,3),df2));
        Xq = cat(3,Xq,repmat(cnt_q,1,1,df1,size(Xq,4)));
        
        cnt_qinv = Inv_C(:,:,aa(pp),bb(pp));
        Xq_inv = cat(4,Xq_inv,repmat(cnt_qinv,1,1,size(Xq_inv,3),df2));
        Xq_inv = cat(3,Xq_inv,repmat(cnt_qinv,1,1,df1,size(Xq_inv,4)));
        
        cnt_trq = Tr_C(aa(pp),bb(pp));
        tr_q = [tr_q repmat(repmat(cnt_trq,wx,1),1,df2)];
        tr_q = [tr_q;repmat(repmat(cnt_trq,1,(2*win_size+1)),df1,1)];
        cnt_trq = repmat(cnt_trq,Lrpc,1);
        
       
        Hbin = histc([tr_q(mask(:,:,1));cnt_trq],Edgex);
        Fx21 = cumsum(Hbin); Fx21 = Fx21./Fx21(end);
        Hbin = histc([tr_q(mask(:,:,2));cnt_trq],Edgex);
        Fx22 = cumsum(Hbin); Fx22 = Fx22./Fx22(end);
        Hbin = histc([tr_q(mask(:,:,3));cnt_trq],Edgex);
        Fx23 = cumsum(Hbin); Fx23 = Fx23./Fx23(end);
        Hbin = histc([tr_q(mask(:,:,4));cnt_trq],Edgex);
        Fx24 = cumsum(Hbin); Fx24 = Fx24./Fx24(end);
        
        D_1(aa(pp)-xx_1+1,bb(pp)-yy_1+1) = 0.25*(max(abs(Fx11-Fx21))+max(abs(Fx12-Fx22))+...
                                                 max(abs(Fx13-Fx23))+max(abs(Fx14-Fx24)) );
        
        
       
        tr_r1 = zeros(Lmsk,4);
        tr_r2 = zeros(Lmsk,4);
        for mm =1:4
            Xpp = Xp(:,:,mask(:,:,mm));
            Xpp_inv = Xp_inv(:,:,mask(:,:,mm));
           
            Xqq = Xq(:,:,mask(:,:,mm));
            Xqq_inv = Xq_inv(:,:,mask(:,:,mm));
           
            dif = Lmsk-size(Xqq,3);
            Xpp = cat(3,Xpp,repmat(cnt_p,1,1,dif));
            Xpp_inv = cat(3,Xpp_inv,repmat(cnt_pinv,1,1,dif));
            Xqq = cat(3,Xqq,repmat(cnt_q,1,1,dif));
            Xqq_inv = cat(3,Xqq_inv,repmat(cnt_qinv,1,1,dif));
            for zz=1:Lmsk
                tr_r1(zz,mm) = abs(trace( Xpp(:,:,zz)*Xqq_inv(:,:,zz) ))/N;
                tr_r2(zz,mm) = abs(trace( Xqq(:,:,zz)*Xpp_inv(:,:,zz) ))/N;
            end
        end
        tr_r1 = [tr_r1; repmat(abs(trace(cnt_p*cnt_qinv))/N,Lrpc,4)];
        tr_r2 = [tr_r2; repmat(abs(trace(cnt_q*cnt_pinv))/N,Lrpc,4)];
      
        Edge = linspace(min([tr_r1(:);tr_r2(:)]),max([tr_r1(:);tr_r2(:)]),1500);
        
        Hbin = histc(tr_r1(:,1),Edge);
        Fr11 = cumsum(Hbin); Fr11 = Fr11./Fr11(end);
        Hbin = histc(tr_r2(:,1),Edge);
        Fr21 = cumsum(Hbin); Fr21 = Fr21./Fr21(end);
        
        Hbin = histc(tr_r1(:,2),Edge);
        Fr12 = cumsum(Hbin); Fr12 = Fr12./Fr12(end);
        Hbin = histc(tr_r2(:,2),Edge);
        Fr22 = cumsum(Hbin); Fr22 = Fr22./Fr22(end);
        
        Hbin = histc(tr_r1(:,3),Edge);
        Fr13 = cumsum(Hbin); Fr13 = Fr13./Fr13(end);
        Hbin = histc(tr_r2(:,3),Edge);
        Fr23 = cumsum(Hbin); Fr23 = Fr23./Fr23(end);
        
        Hbin = histc(tr_r1(:,4),Edge);
        Fr14 = cumsum(Hbin); Fr14 = Fr14./Fr14(end);
        Hbin = histc(tr_r2(:,4),Edge);
        Fr24 = cumsum(Hbin); Fr24 = Fr24./Fr24(end);
        
        D_2(aa(pp)-xx_1+1,bb(pp)-yy_1+1) = 0.25*(max(abs(Fr11-Fr21))+max(abs(Fr12-Fr22))+...
                                                 max(abs(Fr13-Fr23))+max(abs(Fr14-Fr24)) );
    end
    D_1 = D_1.*(sqrt(Lmsk+Lrpc)+0.12+(0.11/sqrt(Lmsk+Lrpc)));
    D_2 = D_2.*(sqrt(Lmsk+Lrpc)+0.12+(0.11/sqrt(Lmsk+Lrpc)));
      
    %______________________________________________________________________
    D_1(ii(kk)-xx_1+1,jj(kk)-yy_1+1) = edgeW1(indx_cmaxf1);
    D_2(ii(kk)-xx_1+1,jj(kk)-yy_1+1) = edgeW2(indx_cmaxf2);
    D_1(isnan(D_1)==1) = edgeW1(end);
    D_2(isnan(D_2)==1) = edgeW2(end);   
%     D_1(find(D_1<edgeW1(1)))=edgeW1(1);
%     D_2(find(D_2<edgeW2(1)))=edgeW2(1);
    
    
    wM1 = D_1*0;
    wM2 = D_1*0;
    for ll=1:Lpq
        wM1(ll)=fdW1(sum(edgeW1<= abs(D_1(ll)) ));
        wM2(ll)=fdW2(sum(edgeW2<= abs(D_2(ll)) ));
    end
    wM = abs((wM1.*wM2).^(1/2));
    wM  = wM./sum(wM(:));
    
    Clcl = Cint(:,:,xx_1:xx_2,yy_1:yy_2);
    Clcl = Clcl(:,:,:);
    wM = wM(:);
    R = 0; varx = 0;
    for zz=1:Lpq
        R = R + Clcl(:,:,zz).*wM(zz);
        varx = varx + (diag(Clcl(:,:,zz)).^2).*wM(zz);
    end
    Cv(:,:,kk) = R;
    
    %______________________________________________________________________
        
    varx = max(0,varx - diag(R).^2) ;
    alpha(kk) = max(max(0,(varx-(1/Look)*(diag(R).^2)))./(varx+eps));

    %______________________________________________________________________
    % disp 
%         parfor_progress;
    if rem(kk, dispN)==0
       kk    
    end
    %______________________________________________________________________
end


FC = Cint.*0;
for kk = 1:Lazi*Lrng
    FC(:,:,ii(kk),jj(kk)) = Cv(:,:,kk);
end
FC = FC.*scale;
alpha = reshape(alpha,Lrng,Lazi)';



end