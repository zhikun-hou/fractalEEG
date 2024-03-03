
% ========================================================================= 
% [Aim]
    % simulate fractal time series by 1-D fractional Brownian motion
    % and wrap these data as the format of FieldTrip

% [Input] cfg
% a struct for arguments
% cfg.N
    % each series generated will be treated as a trial in FT 
% cfg.Hurst
    % Hurst exponent, should be a scalar or a vector which length==cfg.N
    % and range from 0~1
% cfg.with
    % string or string array
    % N(mu,sigma) for white noise
    % O(A,f,phi)  for oscillations
    % see FT.label in generated data for more detail.

% [Output] FT

% =========================================================================

function [FT] = ft_fractalSim(cfg)

    FT = [];
    FT.fsample  = 1024;
    
    % ======================================

    ft_checkopt(cfg,'N','doublescalar');
    ft_checkopt(cfg,'Hurst','doublevector');

    if(mod(cfg.N,1)~=0)
        ft_error("[ft_fractalSim] cfg.N should be an integar.");
    end


    N_hurst = numel(cfg.Hurst);
    if(min(cfg.Hurst)<0||max(cfg.Hurst)>1)
        ft_error("[ft_fractalSim] Hurst exponent shoud range from 0 to 1.");
    end
    if(N_hurst==1)
        cfg.Hurst = ones(1,cfg.N) * cfg.Hurst;
    elseif(N_hurst~=cfg.N)
        ft_error("[ft_fractalSim] When Hurst is a vector, it length should be same as trial number N.");
    end

    % ======================================

    FT.label = {'Raw'};

    for i=1:cfg.N
        [X,T]            = tool_fbm(cfg.Hurst(i),1024,1); % See [1], just rename the function fbm1d
        FT.time{i}       = T;
        FT.trial{i}(1,:) = X;
    end

    % ======================================
    
    cfg.with = ft_getopt(cfg,'with',["none"]);
    if(not(isstring(cfg.with)))
        ft_error("[ft_fractalSim] cfg.with should be string or string vector.");
    end

    % ======================================

    for k=1:numel(cfg.with)
        command = char(cfg.with(k));

        if(any(strcmp(command,{'none','None','NONE'})))
            % Do Nothing
        elseif(startsWith(command,'N(')) % white noise like N(0.1,1)
            FT = tool_withGaussianFT(FT,command);
        elseif(startsWith(command,'O(')) % oscillation like O(A,f,phi)
            FT = tool_withOscillationFT(FT,command);
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
