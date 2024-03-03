

% [command]
% O(A,f,phi)

% see ft_fractalSim

function [EEG] = tool_withOscillationEEG(EEG,command)

    N_trial   = EEG.trials;
    N_channel = EEG.nbchan;

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

    EEG.nbchan = EEG.nbchan * 2;

    % =======================================

    T = EEG.times;

    for i=1:N_trial
        for j=1:N_channel
            EEG.data(j*2,:,i) = EEG.data(j,:,i) + A*sin(2*pi*f*T+Phi(i));
        end
    end

end
