function VV = span_induced(vortex,gamma,point)
    
   % assegno il valore di gamma corrispondente al segmento
    if numel(gamma ==1) % coefficienti d'influenza
       gamma_c = ones(size(vortex,1)-1,size(vortex,2)-1);
    else
        gamma_c = gamma;
    end

    
    line_e = reshape(vortex(1:end-1,2:end,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)-1),1,3);
    line_s = reshape(vortex(1:end-1,1:end-1,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)-1),1,3);
    gamma_c = reshape(gamma_c,numel(gamma_c),1);
    line = [line_e,line_s];
    
    r1 = point-line(:,2,:);
    r2 = point-line(:,1,:);
    
    %cerco h == 0 per imporre poi velocità nulla
    [riga1,colonna1] = find(sqrt(sum(r1.^2,3))==0);
    [riga2,colonna2] = find(sqrt(sum(r2.^2,3))==0);
    [riga3,colonna3] = find(sqrt(sum(cross(r1,r2,3).^2,3))==0);

    num =(sqrt(sum(r1.^2,3))+sqrt(sum(r2.^2,3))).*cross(r1,r2,3);
    den = sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3)).*...
          (sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3))+dot(r1,r2,3));
    % correzione per tenere conto del core viscoso
    k = 1e-3; %parametro modificabile
    corr1 = k^2*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)).^2)./ ...
                sqrt(sum(cross(r1,r2,3).^2,3)).^2;
    corr2 = 2*k*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)))./ ...
                sqrt(sum(cross(r1,r2,3).^2,3));
    corr = 1 + corr1 + corr2;
    
    V = gamma_c/(4*pi).*num./(den.*corr);
    
    if ~isempty(riga1)
        V(riga1,colonna1,:) = zeros(size(riga1,1),size(colonna1,1),3);
    end
    if ~isempty(riga2)
        V(riga2,colonna2,:) = zeros(size(riga2,1),size(colonna2,1),3);
    end   
    if ~isempty(riga3)
        V(riga3,colonna3,:) = zeros(size(riga3,1),size(colonna3,1),3);
    end 
    
    if numel(gamma)==1
       VV = reshape(V,size(vortex,1)-1,size(vortex,2)-1,3);
    else
        VV = sum(V,1);
    end
end
