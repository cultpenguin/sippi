function d_ireal=sippi_nscore(d_real,o_nscore)

[nd,nr]=size(d_real);
if isstruct(o_nscore)
    for id=1:nd
        d_ireal(id,:)=nscore(d_real(id,:),o_nscore);
    end
else
    for id=1:nd
        d_ireal(id,:)=nscore(d_real(id,:),o_nscore{id});
    end
end
