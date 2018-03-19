function [vortex, C_point] = collocazione(wing)
% calcola i punti di controllo (collocazione) e la griglia di votrici
% sfalsata di 1/4 di corda di ogni pannello
 
   vortex(:,:,[1,3]) = wing(:,:,[1,3])+.25*...
                     [wing(2:end,:,[1,3])-wing(1:end-1,:,[1,3]) ; ...
                      wing(end,:,[1,3])-wing(end-1,:,[1,3])];
   vortex(:,:,2) = wing(:,:,2);
                     
   C_point(:,:,[1,3]) = wing(1:end-1,1:end-1,[1,3])+.75*...
                     (wing(2:end,2:end,[1,3])-wing(1:end-1,1:end-1,[1,3]));
   C_point(:,:,2) = wing(1:end-1,1:end-1,2)+.5*...
                         (wing(2:end,2:end,2)-wing(1:end-1,1:end-1,2));
                     
end