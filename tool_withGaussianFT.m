% [command]
% N(mu,sigma)

% see ft_fractalSim

function [FT] = tool_withGaussianFT(FT,command)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);

    % =======================================

    Params = split(command(3:end-1),',');

    mu     = str2double(Params{1});
    sigma  = str2double(Params{2});
    
    % =======================================

    Names = string(FT.label{:});
    Names = [Names,Names+'-WithGaussian'];

    FT.label = cellstr(Names);

    % =======================================

    for i=1:N_trial
        X = FT.trial{i};
        
        Noise = randn(1,1025)*sigma+mu;

        for j=1:N_channel
            Y = X(j,:);
            
            FT.trial{i}(j+N_channel,:) = Y + Noise;
        end
    end

end
