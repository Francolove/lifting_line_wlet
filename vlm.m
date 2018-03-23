%% ala
clear all
close all

naca = 2312;
wing_root = 1;
wing_tip = 1;
wing_span = 10;
n_chord = 50;
n_span = 30;
wing_sweep = deg2rad(0);
wing_twist = deg2rad(0);

% winglet
w_let_root = wing_tip;
w_let_tip = .75;
height = 2;
Radius = .5;
n_height = 12;
cant = deg2rad(0);
sweep = deg2rad(20);
toe_out = deg2rad(-15);
up = 1; % 1 --> winglet verso l'alto
         % -1 --> winglet verso il basso
         

%% code
Wing = build_wing(wing_root,wing_tip,wing_span,n_chord,n_span,naca,...
                                                    wing_sweep,wing_twist);
% 
% W_let = build_winglet(w_let_root,w_let_tip,height,Radius,n_chord,...
%                     n_height,cant,sweep,toe_out,naca,Wing,wing_twist,up);
%  Wing = assemble_wing(Wing,W_let); 
%  Wing2 = Wing; Wing2(:,:,2) = - Wing(:,:,2);
%  Wing = assemble_wing(Wing2(:,end:-1:1,:),Wing);                                                

[vortex,p_controllo] = collocazione(Wing);
[N,T] = versori(Wing);

V_inf = zeros(size(Wing,1)-1,size(Wing,2)-1,3);
V_inf(:,:,1) = 10;
V_inf(:,:,3) = 1;
         
%% plot ala

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
tic
A = induced(vortex,1,p_controllo,V_inf(1,1,:),N);
B = dot(V_inf,N,3)';
b = B(:);
toc
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
tic
V = induced(vortex,gamma_ll,Wing,V_inf(1,1,:),N);
toc
%%
%cp = 1-(V(:,:,1)./sqrt(sum(V_inf.^2,3))).^2;

% ogni sezione produce una portanza pari a -rho*V_inf*Gamma
% va adimensionalizzata con (Scusa Franco) la pressione dinamica all'inf
% -rho * V_inf * Gamma * dy
% ------------------------- = 2 * Gamma * dy / V_inf
%  .5 * rho * V_inf.^2

cp = -gamma_ll*2.*(sqrt(sum((Wing(1:end-1,2:end,:)-Wing(1:end-1,1:end-1,:)).^2,3)))./sqrt(sum(V_inf.^2,3));
figure(400)
surf(p_controllo(:,:,1),p_controllo(:,:,2),p_controllo(:,:,3),cp)
axis equal
shading interp
colorbar

% andamento cp in corda alle varie sezioni in apertura

figure(500)
hold on
for i = 1:size(gamma_ll,2)
    plot3(p_controllo(:,i,1),p_controllo(:,i,2),cp(:,i))
end