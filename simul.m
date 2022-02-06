function pu = simul(t,e,a)
% Simulation of cellular UL system
% input  
%   t: the SINR threshold in dB 
%   e: the FPC factor varing from 0 to 1
%   a: path loss index a>4

Nt = length(t);
N = 3000;    % number of simulations
Thrsh = 10.^(t./10);  % convert to linear values
p = 1;  % Mobile user transmit power
lambda = 1e-5;  % density of BSs
S = 2000;  % Side of simulation square windows (m)


% output
covered = zeros(N, Nt); %initialize covered BS number

for i =1:N
    % generate number of BSs in simulation window
    Nbs = poissrnd(lambda*S^2);
    G = zeros(1,Nbs);  %initialize the channel gain vector;
    U = zeros(1,Nbs);  %initialize the dictance vector: between UE and typical BS;
    R = zeros(1,Nbs);  %initialize the dictance vector: between UE and serving BS;
    % Build a null association table
    A = cell(Nbs,2);
%     for j = 1:Nbs
%         A{j,1} = [];
%         A{j,2} = [];
%     end
    % Set typical BS position as origin
    A{1,1} = [0,0];
    for j = 2:Nbs
        % generate random location for BS j
        x = rand*(2*S) - S;
        y = rand*(2*S) - S;
        A{j,1} = [x,y];
    end
    % Fill col 2 with UEs
    while sum(sum(cellfun(@isempty, A)))>0
        % generate random location for UE
        u = rand*(2*S) - S;
        v = rand*(2*S) - S;
        % calculate the Nbs distances between the UE and every BS j
        A1 = cell2mat(A(:,1));
        % n is the BS's index which  nearest to the UE
        [n,~] = dsearchn(A1,[u v]);
        if isempty(A{n,2})
            A{n,2}=[u,v];
        end
    end
    
    %calculate SIR at typical BS , update covered
    % generate channel fading gain G (Rayleigh fading)
    G(1) = exprnd(1);
    % calculate the distances between typical BS and tagged UE
    R(1) = pdist2(A{1,1},A{1,2});
    NUM = p*G(1)*R(1)^(-a*(1-e));% numerator of SIR
    DENOM = 0;   % denominator of SIR
    
    for j=2:Nbs
        % generate channel fading gain G (Rayleigh fading)
        G(j) = exprnd(1);
        % calculate the distances between UE j and serving BS 
        R(j) = pdist2(A{j,1},A{j,2});
         % calculate the distances between UE j and typical BS
        U(j) = pdist2(A{1,1},A{j,2});
        DENOM = DENOM + p*G(j)*R(j)^(a*e)*U(j)^(-a);
    end
    
    SIR = NUM/DENOM;
    covered(i,:) = SIR>Thrsh;
end
pu = sum(covered)/N;
end
    
    
