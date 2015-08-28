% setup_wavelet_matrix: setup matrix of wavelet to be used for convolution
%
% Call: 
%   W=setup_wavelet_matrix(data,t,nt_out);
%   input:
%     data [1,ns]: wavelet defined at ns time samples
%     t [1,ns]: time of wavelet defined at ns times sample
%               The center of the wavelet should be at t=0;
%     nt_out: size of the output wavelet matrix
%
%   outout:
%     W: [nt_out,nt_out]
%
function W=setup_wavelet_matrix(data,t,nt_out);

it0=find(t==0);
if isempty(it0)
    it0=find(abs(t)==min(abs(t)));
    it0=it0(1);
end

W=zeros(nt_out,nt_out);
nd=length(data);

%%
for i=1:nt_out;
    i_arr=i+[1:1:nd]-it0;
    i_use=find( (i_arr>0) & (i_arr<=nt_out) );
    W(i,i_arr(i_use))=data(i_use);
end
    
    














