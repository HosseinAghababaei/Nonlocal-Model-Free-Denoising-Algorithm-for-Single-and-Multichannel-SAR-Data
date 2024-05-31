function Pauli = Pauli_C(C)

Pauli(:,:,1) = 10*log( 0.5*abs( squeeze(C(1,1,:,:))+squeeze(C(3,3,:,:))-squeeze(C(1,3,:,:))-squeeze(C(3,1,:,:)) ).^2 );
Pauli(:,:,2) = 10*log( 2*abs( squeeze(C(2,2,:,:)) ).^2 );
Pauli(:,:,3) = 10*log( 0.5*abs( squeeze(C(1,1,:,:))+squeeze(C(3,3,:,:))+squeeze(C(1,3,:,:))+squeeze(C(3,1,:,:)) ).^2 );
    

% if nargin >1
%     
%     for kk = 1:3
%        xx = Pauli(:,:,kk);
%        xx = (255./(max(max(xx))-min(min(xx))))*(xx-min(min(xx)));
%        Pauli(:,:,kk) = histeq( uint8(xx) );
%     end
%     clear xx kk
%     
% end
%         

for kk=1:3
    xx = Pauli(:,:,kk);
    xx = (255./(max(max(xx))-min(min(xx))))*(xx-min(min(xx)));
    Pauli(:,:,kk) = histeq( uint8(xx) );
end
clear xx kk 
Pauli = uint8(Pauli);


end