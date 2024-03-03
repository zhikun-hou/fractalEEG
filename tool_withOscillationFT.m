

% [command]
% O(A,f,phi)

% see ft_fractalSim

function [FT] = tool_withOscillationFT(FT,command)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);

    % =======================================

    Params = split(command(3:end-1),',');

    A       = str2double(Params(1));
    f       = str2double(Params(2));

    phi     = Params(3);
    if(any(strcmp(phi,{'random','Random'})))
        Phi = rand(1,N_trial) * pi;
    else
        Phi = ones(1,N_trial) * str2double(phi);
    end

    % =======================================

    Names = string(FT.label{:});
    Names = [Names,Names+'-WithOscillation'];

    FT.label = cellstr(Names);

    % =======================================

    T = FT.time{1};

    for i=1:N_trial
        X = FT.trial{i};
        for j=1:N_channel
            Y = X(j,:);
            
            FT.trial{i}(j+N_channel,:) = Y + A*sin(2*pi*f*T+Phi(i));
        end
    end

end
