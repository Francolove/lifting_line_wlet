function VV_bu = bu_wake_induced(vortex,gamma,point,dir)
    
    % assegno il valore di gamma corrispondente al segmento
    if numel(gamma ==1) % coefficienti d'influenza
       gamma_c_bu = ones(size(vortex,1)-1,size(vortex,2));
    else
        gamma_c_bu = gamma; gamma_c_bu(:,2:end) = gamma(:,1:end-1)-gamma(:,2:end);
        gamma_c_bu(:,end+1) = gamma(:,end);
         for i = 2:size(gamma,1)
             gamma_c_bu(i,:) = gamma_c_bu(i,:)+gamma_c_bu(i-1,:);
         end
    end
    
    gamma_c_bu = gamma_c_bu(end,:);
    
    start_point = reshape(vortex(end,:,:),size(vortex,2),1,3);
    gamma_c_bu = reshape(gamma_c_bu,numel(gamma_c_bu),1);
    
    r = point-start_point;
    
    [riga,colonna] = find(sqrt(sum(r.^2,3))==0);
    dir = repmat(dir,size(vortex,2),1,1);
    num =cross(dir,r,3);
    den = sqrt(sum(r.^2,3)).*(sqrt(sum(r.^2,3))-dot(dir,r,3));
    % correzione per tenere conto del core viscoso
    
    V_bu = gamma_c_bu/(4*pi).*num./(den);
    
    if ~isempty(riga)
        V_bu(riga,colonna,:) = zeros(size(riga,1),size(colonna,1),3);
    end

    
    if numel(gamma)==1
       V_bu = reshape(V_bu,1,size(vortex,2),3);
       VV_bu = ricostruzione_ind_vortice(V_bu);
       VV_bu = repmat(VV_bu,size(vortex,1)-1,1,1);
    else
        VV_bu = sum(V_bu,1);
    end
end


function VV = ricostruzione_ind_vortice(V)
        VV(1,1,:) = V(1,1,:) - V(1,2,:);
        VV(1,2:size(V,2)-1,:) = V(1,3:end,:) + V(1,2:end-1,:);
end