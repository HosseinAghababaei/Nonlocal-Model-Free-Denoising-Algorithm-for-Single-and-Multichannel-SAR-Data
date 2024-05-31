function [fdW1, edgeW1,fdW2, edgeW2] = weight_fcn(Edge,L,N)


gma = min(1,L/N)^(1/3);
K = 1e+4;
Nr = 15^2;
dM_r1  =zeros(K,1);
dM_r2  =zeros(K,1);

parfor k=1:K
    
    tr_1 = zeros(Nr,1);
    tr_2 = zeros(Nr,1);
    tr_3 = zeros(Nr,1);
    tr_4 = zeros(Nr,1);
    for nn=1:Nr
        x1 = (randn(L,N) + 1i*randn(L,N))*eye(N)./sqrt(2);
        X1 = (x1'*x1)/L;
        X1 = gma*(X1-diag(diag(X1)))+diag(diag(X1));
        
        x2 = (randn(L,N) + 1i*randn(L,N))*eye(N)./sqrt(2);
        X2 = (x2'*x2)/L;        
        X2 = gma*(X2-diag(diag(X2)))+diag(diag(X2));
               
        tr_1(nn) = abs(trace(X1));
        tr_2(nn) = abs(trace(X2));
        tr_3(nn) = abs(trace(X1^(-1/2)*X2*X1^(-1/2)));
        tr_4(nn) = abs(trace(X2^(-1/2)*X1*X2^(-1/2)));
    end
    tr_1 = tr_1./N;
    tr_2 = tr_2./N;
    tr_3 = tr_3./N;
    tr_4 = tr_4./N;
    
    Hbin = histc(tr_1,Edge);
    Fx1 = cumsum(Hbin); Fx1 = Fx1./Fx1(end);
    HbinIv = histc(tr_2,Edge);
    Fx2 = cumsum(HbinIv); Fx2 = Fx2./Fx2(end);
    
%     Hbin12 = histc([tr_1 tr_2],Edge);
%     Fx12 = cumsum(Hbin12); Fx12 = Fx12./Fx12(end);   
%     X = ((Fx1-Fx2).^2)./(Fx12.*(1-Fx12));
%     X(isnan(X))=[];
%     X(isinf(X))=[];
%     dM_r(k) = (sqrt(Nr)+0.12+(0.11/sqrt(Nr)))*sqrt( mean( X ));
    
    dM_r1(k)=(sqrt(Nr)+0.12+(0.11/sqrt(Nr)))*max(abs(  Fx1-Fx2  ));
%    k

    Hbin = histc(tr_3,Edge);
    Fx1 = cumsum(Hbin); Fx1 = Fx1./Fx1(end);
    HbinIv = histc(tr_4,Edge);
    Fx2 = cumsum(HbinIv); Fx2 = Fx2./Fx2(end);
    dM_r2(k)=(sqrt(Nr)+0.12+(0.11/sqrt(Nr)))*max(abs(  Fx1-Fx2  ));  
end



[fdW1, edgeW1]=ksdensity(dM_r1(:));
fdW1=interp(fdW1,10);
edgeW1=interp(edgeW1,10);


[fdW2, edgeW2]=ksdensity(dM_r2(:));
fdW2=interp(fdW2,10);
edgeW2=interp(edgeW2,10);


% figure, plot(edgeW1,fdW1,'b','LineWidth',2)
% hold on 
% plot(edgeW2,fdW2,'r','LineWidth',2)
% xlim([0 4])


fdW1= [0 abs(fdW1)];
edgeW1 = [0 edgeW1];
fdW2= [0 abs(fdW2)];
edgeW2 = [0 edgeW2];


% FdW1 = cumsum(fdW1); FdW1 = FdW1./FdW1(end);
% FdW2 = cumsum(fdW2); FdW2 = FdW2./FdW2(end);

end