function [ output ] = RLS_Train( input, TrainSeq, delay, lambda, Up2 )
% Training sequence based RLS equalization
taps = 2*delay + 1;
iter = 1;
TrainLen = length(TrainSeq);
DataLen = length(input);
output = zeros(1, DataLen);
w = zeros(1, taps);
% w(delay + 1) = 1;
error = [];
input = [zeros(1,delay) input zeros(1,delay)];
Delta = 0.1*eye(taps,taps);
%% training
for kk =1 : iter
    
for ii = 1 : TrainLen
    x = input( (ii-1)*Up2+1: (ii-1)*Up2+taps );
    err = TrainSeq(ii) - conj(w) * x.';
    
    G = Delta*x.' / (lambda + conj(x) *Delta*x.');
    Delta = 1/lambda*(Delta - G*conj(x)*Delta);
    w = w + G.'.*conj(err);
%     w = w + step * err .* x;
    error = [error err];
end

end

%     figure;plot(abs(error));


%% equalization
for jj = 1 : DataLen
    in = input(jj : jj + taps - 1);
    output(jj) = conj(w) * in.';
end

end

