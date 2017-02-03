function A_out=resize(A_in,factor,mode)

% Call: A_out=resize(A_in,factor);
% * A_in is the input array
% * A_out is the output array
% * factor is the factor by which the sides of the input array cells are refined
%   or coarsened: factor > 1 = The array is refined; factor < 1 = The array is coarsened
%
% KSC (2009)

[y,x]=size(A_in);

if nargin<3
    mode=1;
end

try
    % Refine:
    if factor > 1
        for i=1:y
            for j=1:x
                A_out((i-1)*factor+1:i*factor,((j-1)*factor+1:j*factor))=A_in(i,j);
            end
        end
    end

    % Coarsen:
    if factor < 1
        [y,x]=size(ones(round(y*factor),round(x*factor)));
        bin=round(1/factor);
        for i=1:y
            for j=1:x
                if mode==1
                    A_out(i,j)=mean(mean(A_in((i-1)*bin+1:i*bin,((j-1)*bin+1:j*bin))));
                else
                    A_out(i,j)=sum(sum(A_in((i-1)*bin+1:i*bin,((j-1)*bin+1:j*bin))));
                end 
            end
        end
    end
    
    if factor==1
        A_out=A_in;
    end

catch
    disp('Something went wrong: An odd scaling factor was possibly chosen')
end
