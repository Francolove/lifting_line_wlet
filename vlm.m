%% ala
clear all

naca = 0012;
wing_root = 1;
wing_tip = 1;
wing_span = 5;
n_chord = 50;
n_span = 50;
wing_sweep = deg2rad(0);
wing_twist = deg2rad(0);


% winglet
w_let_root = wing_tip;
w_let_tip = 1;
height = 2;
Radius = .5;
n_height = 20;
cant = deg2rad(15);
sweep = deg2rad(0);
toe_out = deg2rad(5);
up = 1; % 1 --> winglet verso l'alto
         % -1 --> winglet verso il basso


%% code
Wing = build_wing(wing_root,wing_tip,wing_span,n_chord,n_span,naca,...
                                                    wing_sweep,wing_twist);

% W_let = build_winglet(w_let_root,w_let_tip,height,Radius,n_chord,...
%                     n_height,cant,sweep,toe_out,naca,Wing,wing_twist,up);
%  Wing = assemble_wing(Wing,W_let);
 Wing2 = Wing; Wing2(:,:,2) = - Wing(:,:,2);
 Wing = assemble_wing(Wing2(:,end:-1:1,:),Wing);

[vortex,p_controllo] = collocazione(Wing);
[N,T] = versori(Wing);



% %% plot ala
%
% figure(100)
% surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3))
% xlabel('chord')
% ylabel('span')
% axis equal
% hold on
% plot3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),'r*')
% quiver3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),...
%                                           N(:,:,1),N(:,:,2),N(:,:,3),'k')
% quiver3(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),...
%                                           T(:,:,1),T(:,:,2),T(:,:,3),'y')
%% risoluzione sistema lineare
% tic
% [At,bt] = influence(vortex,p_controllo,V_inf(1,1,:),N,1);
% toc
alpha = deg2rad(5);linspace(deg2rad(0),deg2rad(12),12);
for i = 1:numel(alpha)
V_inf = zeros(size(Wing,1)-1,size(Wing,2)-1,3);
V_inf(:,:,1) = cos(alpha(i));
V_inf(:,:,3) = sin(alpha(i));7
dummy = NaN;
A = induced(vortex,dummy,p_controllo,V_inf(1,1,:),N);
B = -dot(V_inf,N,3)';
b = B(:);
%%
gamma_ll = A\b;
gamma_ll = reshape(gamma_ll,size(p_controllo,1),size(p_controllo,2));
% surf(Wing(:,:,1),Wing(:,:,2),Wing(:,:,3),gamma_ll)
% % surf(gamma_ll')
% axis equal
% xlabel('chord')
% ylabel('span')
% colorbar
%
middle_point_span = (vortex(1:end-1,1:end-1,:)+vortex(1:end-1,2:end,:))/2;
dspan = vortex(1:end-1,2:end,:)-vortex(1:end-1,1:end-1,:);
middle_point_chord = (vortex(1:end-1,:,:)+vortex(2:end,:,:))/2;
dchord = vortex(1:end-1,:,:)-vortex(2:end,:,:);

tic
mid_span_speed = induced(vortex,gamma_ll,middle_point_span,V_inf(1,1,:),N);
mid_chord_speed = induced(vortex,gamma_ll,middle_point_chord,V_inf(1,1,:),N);
toc

quiver3(middle_point_span(:,:,1),middle_point_span(:,:,2),middle_point_span(:,:,3),V_inf(1,1,1)+mid_span_speed(:,:,1),V_inf(1,1,2)+mid_span_speed(:,:,2),V_inf(1,1,3)+mid_span_speed(:,:,3),'r')
hold on
quiver3(middle_point_chord(:,:,1),middle_point_chord(:,:,2),middle_point_chord(:,:,3),V_inf(1,1,1)+mid_chord_speed(:,:,1),V_inf(1,1,2)+mid_chord_speed(:,:,2),V_inf(1,1,3)+mid_chord_speed(:,:,3),'b')
figure
quiver3(middle_point_chord(:,:,1),middle_point_chord(:,:,2),middle_point_chord(:,:,3),dchord(:,:,1),dchord(:,:,2),dchord(:,:,3),'b')
hold on
quiver3(middle_point_span(:,:,1),middle_point_span(:,:,2),middle_point_span(:,:,3),dspan(:,:,1),dspan(:,:,2),dspan(:,:,3),'r')

%tic
%[V,coeffxcd] = induced(vortex,gamma_ll,p_controllo,V_inf(1,1,:),N);
%toc
%%
%cp = 1-(V(:,:,1)./sqrt(sum(V_inf.^2,3))).^2;

% ogni sezione produce una portanza pari a -rho*V_inf*Gamma
% va adimensionalizzata con (Scusa Franco) la pressione dinamica all'inf
% -rho * V_inf * Gamma * dy
% ------------------------- = 2 * Gamma * dy / V_inf
%  .5 * rho * V_inf.^2

cl = gamma_ll*2.*(Wing(1:end-1,2:end,2)-Wing(1:end-1,1:end-1,2))/(wing_root*(Wing(1,end,2)-Wing(1,1,2)));


dl = Wing(1:end-1,2:end,:)-Wing(1:end-1,1:end-1,:);

cfs = gamma_ll.*cross(V_inf(1,1,:)+mid_span_speed,dspan,3)/(wing_root*(Wing(1,end,2)-Wing(1,1,2)));

gamma_c = gamma_ll; gamma_c(:,2:end) = -gamma_ll(:,1:end-1)+gamma_ll(:,2:end);
gamma_c(:,end+1) = gamma_ll(:,end); % numel(segmenti_riga) = numel(vortici_riga)+1

% per ogni colonna il segmento j-esimo avra
% intensità pari alla somma della sua gamma e della gamma del
% segmento che lo precede
for j = 2:size(gamma_ll,1)
    gamma_c(i,:) = gamma_c(j,:)+gamma_c(j-1,:);
end

cfc = gamma_c.*cross(V_inf(1,1,:)+mid_chord_speed,-dchord,3)/(wing_root*(Wing(1,end,2)-Wing(1,1,2)));



%cd = coeffxcd.*(Wing(1:end-1,2:end,2)-Wing(1:end-1,1:end-1,2))/(wing_root*(Wing(1,end,2)-Wing(1,1,2)));

%cp2 = 1-(sqrt(sum(V.^2,3))./sqrt(sum(V_inf(1,1,:).^2,3))).^2;
CL(i) = sum(sum(cfc(:,:,3),1),2)+sum(sum(cfs(:,:,3),1),2);
CD(i) = sum(sum(cd));

end
% figure(400)
% subplot(2,1,1)
% surf(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),cl)
% axis equal
% shading interp
% colorbar
% subplot(2,1,2)
% surf(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),cp2)
% axis equal
% shading interp
% colorbar
% linkaxes
%
% % andamento cp in corda alle varie sezioni in apertura
%
% figure(500)
% for i = 1:size(gamma_ll,2)
%     subplot(1,2,1)
%     plot3(p_controllo(:,i,1),p_controllo(:,i,2),cl(:,i))
%     hold on
%     subplot(1,2,2)
%     hold on
%     plot3(p_controllo(:,i,1),p_controllo(:,i,2),cp2(:,i))
% end
