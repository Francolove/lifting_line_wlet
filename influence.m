function [A,b,V] = influence(vortex,point,dir,N,gamma)
% costrusice la matrice dei coefficienti d'induzione della griglia di
% vortici vortex nei punti point
    A = [];
    b = [];
    V = zeros(size(point));
    if numel(gamma) ==1
        gamma = repmat(gamma,size(vortex,1)-1,size(vortex,2)-1);
    end
    parfor i = 1:size(point,1)
        for j = 1:size(point,2)
            dummy =[];
            for k = 1:size(vortex,1)-1
                for L= 1:size(vortex,2)-1
                    dummy(k,L,:) = induced_speed(vortex(k:end,L:L+1,:),...
                                             point(i,j,:), dir,gamma(k,L));
                end
            end
            V = V+dummy;
            dummyA = dot(dummy,repmat(N(i,j,:),size(dummy,1),...
                                                       size(dummy,2),1),3);
            dummyA = reshape(dummyA,1,numel(dummyA));
            A = [A;dummyA];
            b = [b;-dot(N(i,j,:),dir)];
        end
    end  
end
