function V = induced_semiinf_speed(point_start,point,dir,gamma)
% calcola la velocità indotta da un vortice semiinfinito che parte dal
% punto point_start e viaggia allieato al vettore dir

    r = point-point_start;
    
    % if sqrt(sum(r1.^2,3)) == 0 || sqrt(sum(r2.^2,3)) == 0
    % V = zeros(1,1,3);
    % else
    [riga,colonna] = find(sqrt(sum(r.^2,3))==0);
    dir = repmat(dir,2,1,1);
    num =cross(dir,r,3);
    den = sqrt(sum(r.^2,3)).*(sqrt(sum(r.^2,3))-dot(dir,r,3));
    % correzione per tenere conto del core viscoso
    
    V = gamma/(4*pi).*num./(den);
    
    if ~isempty(riga)
        V(riga,colonna,:) = zeros(size(riga,1),size(colonna,1),3);
    end

end