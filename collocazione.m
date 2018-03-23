function [vortex, C_point] = collocazione(wing)
% calcola i punti di controllo (collocazione) e la griglia di votrici
% sfalsata di 1/4 di corda di ogni pannello
 
   vortex(:,:,[1,3]) = wing(:,:,[1,3])+.25*...
                     [wing(2:end,:,[1,3])-wing(1:end-1,:,[1,3]) ; ...
                      wing(end,:,[1,3])-wing(end-1,:,[1,3])];
                  
   vortex(:,:,2) = wing(:,:,2);

   C_point(:,:,2) = wing(1:end-1,1:end-1,2)+.5*...
                         (wing(2:end,2:end,2)-wing(1:end-1,1:end-1,2));
                     
                     
   point_3_4 = wing(1:end-1,:,:) + 3/4*(wing(2:end,:,:)-wing(1:end-1,:,:)); 
     
   dx = point_3_4(:,2:end,1)-point_3_4(:,1:end-1,1);
   dy = point_3_4(:,2:end,2)-point_3_4(:,1:end-1,2);
   dz = point_3_4(:,2:end,3)-point_3_4(:,1:end-1,3);
   alfax = atan2(dx,dy);
   alfaz = atan2(dz,dy);
   
   C_point(:,:,1) = point_3_4(:,1:end-1,1) + (C_point(:,:,2)-...
                                    wing(1:end-1,1:end-1,2)).*tan(alfax);
   C_point(:,:,3) = point_3_4(:,1:end-1,3) + (C_point(:,:,2)-...
                                    wing(1:end-1,1:end-1,2)).*tan(alfaz);


%    C_point(:,:,1) = wing(1:end-1,1:end-1,1)+3/4*...
%                      (wing(2:end,2:end,1)-wing(1:end-1,1:end-1,1));
%                  
%    C_point(:,:,3) = wing(1:end-1,1:end-1,3)+3/4*...
%                         (wing(2:end,2:end,3)-wing(1:end-1,1:end-1,3));
%                  
end