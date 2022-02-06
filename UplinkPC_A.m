function pc = UplinkPC_A(t)
% e=1,for interference limited case
T = 10.^(t./10);
pc = exp( -sqrt(T).*(pi/2 - atan(T.^(-1/2))));
end