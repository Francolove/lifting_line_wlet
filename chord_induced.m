function VV = chord_induced(vortex,gamma,point,dir)
    
    % assegno il valore di gamma corrispondente al segmento
    if numel(gamma ==1) % coefficienti d'influenza
       gamma_c = ones(size(vortex,1)-1,size(vortex,2));
    else
        gamma_c = gamma; gamma_c(:,2:end) = gamma(:,1:end-1)-gamma(:,2:end);
        gamma_c(:,end+1) = gamma(:,end);
         for i = 2:size(gamma,1)
             gamma_c(i,:) = gamma_c(i,:)+gamma_c(i-1,:);
         end
    end
    
    gamma_c_bu = gamma_c(end,:);
    
    line_e = reshape(vortex(1:end-1,:,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)),1,3);
    line_s = reshape(vortex(2:end,:,:),(size(vortex,1)-1)*...
                                        (size(vortex,2)),1,3);
       
                                    
    start_point = reshape(vortex(end,:,:),size(vortex,2),1,3);
    
    
    gamma_c_bu = reshape(gamma_c_bu,numel(gamma_c_bu),1);                                
    gamma_c = reshape(gamma_c,numel(gamma_c),1);
    
    line = [line_e,line_s];
    
    r1 = point-line(:,2,:);
    r2 = point-line(:,1,:);
    
    r = point-start_point;
    

    %cerco h == 0 per imporre poi velocità nulla
    [riga1,colonna1] = find(sqrt(sum(r1.^2,3))==0);
    [riga2,colonna2] = find(sqrt(sum(r2.^2,3))==0);
    [riga3,colonna3] = find(sqrt(sum(cross(r1,r2,3).^2,3))==0);
    [riga,colonna] = find(sqrt(sum(r.^2,3))==0);

    num =(sqrt(sum(r1.^2,3))+sqrt(sum(r2.^2,3))).*cross(r1,r2,3);
    den = sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3)).*...
          (sqrt(sum(r1.^2,3)).*sqrt(sum(r2.^2,3))+dot(r1,r2,3));
    % correzione per tenere conto del core viscoso
    k = 0;1e-3; %parametro modificabile
    corr1 = k^2*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)).^2)./ ...
                sqrt(sum(cross(r1,r2,3).^2,3)).^2;
    corr2 = 2*k*(sqrt(sum(r1.^2,3))-sqrt(sum(r2.^2,3)))./ ...
                sqrt(sum(cross(r1,r2,3).^2,3));
    corr = 1 + corr1 + corr2;
    
    V = gamma_c/(4*pi).*num./(den.*corr);
       
    dir = repmat(dir,size(vortex,2),1,1);
    
    num_bu =cross(dir,r,3);
    den_bu = sqrt(sum(r.^2,3)).*(sqrt(sum(r.^2,3))-dot(dir,r,3));
       
    V_bu = gamma_c_bu/(4*pi).*num_bu./(den_bu);
    
    if ~isempty(riga)
        V_bu(riga,colonna,:) = zeros(size(riga,1),size(colonna,1),3);
    end
    
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
       V = reshape(V,size(vortex,1)-1,size(vortex,2),3);
       VV = ricostruzione_ind_vortice(V);
       
       V_bu = reshape(V_bu,1,size(vortex,2),3);
       VV_bu = ricostruzione_ind_vortice_bu(V_bu);
       VV_bu = repmat(VV_bu,size(vortex,1)-1,1,1);
       VV = VV+VV_bu;
       
    else
        VV = sum(V,1) + sum(V_bu,1);
    end
    
end


function VV = ricostruzione_ind_vortice(V)
    for i = 1:size(V,1)
        VV(i,1,:) = +sum(V(i:end,1,:),1) - sum(V(i:end,2,:),1);
        VV(i,2:size(V,2)-1,:) = -sum(V(i:end,3:end,:),1) + ...
                                                 sum(V(i:end,2:end-1,:),1);
%         VV(i,size(V,2)-1) = sum(V(i:end,end) + sum(V(i:end,end-1)
    end
end



function VV = ricostruzione_ind_vortice_bu(V)
        VV(1,1,:) = +V(1,1,:) - V(1,2,:);
        VV(1,2:size(V,2)-1,:) = -V(1,3:end,:) + V(1,2:end-1,:);
end
