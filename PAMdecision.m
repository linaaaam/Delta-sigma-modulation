function [ out ] = PAMdecision( in, constellation)

N1 = length(in);
out = zeros(1,N1);

Symbol = constellation;
N2 = length(Symbol);

for ii = 1:N1
    Dis = zeros(1,N2);
    for jj = 1:N2
        Dis(jj) = abs(in(ii)-Symbol(jj)).^2;
        [~, min_ind]=min(Dis);
        out(ii) = Symbol(min_ind);
    end
end

end

