function V = induced(vortex,gamma,points,dir,N)

    if numel(gamma)==1
        V = zeros(size(points,1)*size(points,2));
    else
        V = zeros([size(gamma),3]);
    end
    pp = size(points,2);
    for j = 1:size(points,2)
        for i = 1:size(points,1)
            vs = span_induced(vortex,gamma,points(i,j,:));
            vc = chord_induced(vortex,gamma,points(i,j,:),dir);
%             vbu = bu_wake_induced(vortex,gamma,points(i,j,:),dir);
            if numel(gamma)==1
                dummy = vs+vc;
                dummy = reshape(dummy,1,size(dummy,1)*size(dummy,2),3);
                V((i-1)*pp+j,:) = dot(dummy,repmat(N(i,j,:)...
                                                   ,1,size(dummy,2),1),3);
            else
                V(i,j,:) = vs+vc;
            end
            
        end
    end





end