function [ Out ] = DSM1_NC( In, AmpRatio, Thres, Order )
%% 1st-order delta-sigma modulator

U = In/AmpRatio;
Thre = Thres;
V =  ones(1,length(U));
X1 = 0; X2 = 0;
X = 0;
err = 0;
switch Order
    case 1
        
    for ii = 2:length(U)
         e = U(ii) - V(ii-1);
         X = X + e;
    
         V(ii) = 2*(X>0)-1;
    end

    case 2
        
    for ii = 2:length(U)
         e = U(ii) - 1*V(ii-1);
         X1 = X1 + e;
         
         X = X1 - err;
%          X = X*0 + eee;
        if X > Thre
             a = 3;
        elseif (X <= Thre)&&(X > 0)
             a = 2;
        elseif (X <= 0)&&(X > -Thre)
        a = 1;
       else a = 0;
        end
          V(ii) = a*2-3;
          
         err = (V(ii) - X);
     end

end

Out = V;
% scatterplot(Out);
end

