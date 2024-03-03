
% ========================================================================= 
% [Aim]
    % simulate fractal time series by 1-D fractional Brownian motion
    % and wrap these data as the format of EEGLab

% [Input] cfg
% a struct for arguments
% cfg.N
    % each series generated will be treated as a trial in EEG 
% cfg.Hurst
    % Hurst exponent, should be a scalar or a vector which length==cfg.N
    % and range from 0~1

% [Output] FT

% =========================================================================

function [EEG] = eeg_fractalSim(cfg)

    EEG = [];
    EEG.srate  = 1024;
    EEG.npnts  = 1025;
    EEG.times  = 0:1/1024:1;
    
    % ======================================

    if(not(isfield(cfg,'N')))
        error("[eeglab_fractalSim] 'N' can't be empty.");
    end
    if(mod(cfg.N,1)~=0)
        error("[eeglab_fractalSim] cfg.N should be an integar.");
    end

    
    if(not(isfield(cfg,'Hurst')))
        error("[eeglab_fractalSim] 'Hurst' can't be empty.");
    end
    if(min(cfg.Hurst)<0||max(cfg.Hurst)>1)
        error("[eeglab_fractalSim] Hurst exponent shoud range from 0 to 1.");
    end
    N_hurst = numel(cfg.Hurst);
    if(N_hurst==1)
        cfg.Hurst = ones(1,cfg.N) * cfg.Hurst;
    elseif(N_hurst~=cfg.N)
        error("[eeglab_fractalSim] When Hurst is a vector, it length should be same as trial number N.");
    end


    if(not(isfield(cfg,'with')))
        cfg.with = ["none"];
    elseif(not(isstring(cfg.with)))
        error("[eeglab_fractalSim] cfg.with should be string or string vector.");
    end

    % ======================================

    EEG.trials = cfg.N;
    EEG.nbchan = 1;

    for i=1:cfg.N
        [X,T] = tool_fbm(cfg.Hurst(i),1024,1); % See [1], just rename the function fbm1d
        EEG.data(1,:,i) = X;
    end

    % ======================================
    
    for k=1:numel(cfg.with)
        command = char(cfg.with(k));

        if(any(strcmp(command,{'none','None','NONE'})))
            % Do Nothing
        elseif(startsWith(command,'N(')) % white noise like N(0.1,1)
            EEG = tool_withGaussianEEG(EEG,command);
        elseif(startsWith(command,'O(')) % oscillation like O(A,f,phi)
            EEG = tool_withOscillationEEG(EEG,command);
        else
            ft_error("[ft_fractalSim] Unknown cfg.with command.");
        end
    end

end

% =========================================================================

% [Author]
% ZhiKun Hou,
% State Key Laboratory of Cognitive Neuroscience and Learning,
% Beijing Normal University (BNU)

% [Contact]
% www.github.com/zhikun-hou     # suggest
% zhikun.hou@mail.bnu.edu.cn
% bnu.zhikun.hou@gmail.com

% [Reference]
% [1] Zdravko Botev (2024). Fractional Brownian motion generator (https://www.mathworks.com/matlabcentral/fileexchange/38935-fractional-brownian-motion-generator), MATLAB Central File Exchange.

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly
