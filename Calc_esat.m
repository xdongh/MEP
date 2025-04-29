function es = Calc_esat(T)
% This fucntion is calculate saturated vapor pressure with
% Clasusius-Clapeyron relations assuming lambda = constant
A = 2.53*10^11; %[Pa]
B = 5.42*10^3;  %[k]

es = A .* exp(-B./T);

% [m,n] = size(T);
% 
% es = NaN(m,n);
% 
% for i = 1 : m
%     for j = 1 : n
%         es(i,j) = A * exp(-B/T(i,j));
%     end
% end

end

