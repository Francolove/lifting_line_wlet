function V = induced_speed(u_vortex, point, dir, gamma)
% calcola la velocit� indotta da un vortice a ferro di cavallo (che pu� non
% essere per forza "rettilineo" i.e.(segue la curvatura della corda)

    % tratto di aperura
    chord = [u_vortex(1,2,:) u_vortex(1,1,:)]; %[fine partenza]
    
    % tratto lungo la corda
    line_left = [u_vortex(1:end-1,1,:) u_vortex(2:end,1,:)];
    
    line_right = [u_vortex(2:end,2,:) u_vortex(1:end-1,2,:)];
    
    %assemblo i tratti
    line = [chord; line_left; line_right];
    
    gamma_u = gamma*[1;ones(size(line_left,1),1);ones(size(line_right,1),1)];
    
    Vl = induced_line_speed(line,gamma_u,point);
    
    % tratto_semiinfinito
    p_inf = [u_vortex(end,1,:); u_vortex(end,2,:)];
    
    gamma_inf = gamma*[1;-1]; % uno dei due deve essere opposto al primo
    
    V_inf = induced_semiinf_speed(p_inf,point ,dir,gamma_inf);
    
    V = sum([Vl;V_inf],1);
end