% [command]
% N(mu,sigma)

% see ft_fractalSim

function [EEG] = tool_withGaussianEEG(EEG,command)

    N_trial   = EEG.trials;
    N_channel = EEG.nbchan;

    % =======================================

    Params = split(command(3:end-1),',');

    mu     = str2double(Params(1));
    sigma  = str2double(Params(2));
    
    % =======================================

    EEG.nbchan = EEG.nbchan * 2;

    % =======================================

    for i=1:N_trial
        
        Noise = randn(1,1025)*sigma+mu;

        for j=1:N_channel
            
            EEG.data(j*2,:,i) = EEG.data(j,:,i) + Noise;
        end
    end

end
