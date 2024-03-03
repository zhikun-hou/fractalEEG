
% ========================================================================= 
% [Aim]
    % simulate fractal time series by 1-D fractional Brownian motion
    % and wrap these data as the format of FieldTrip

% [Input] cfg
% a struct for arguments
% cfg.N
    % each series generated will be treated as a channel in FT 
% cfg.Hurst
    % Hurst exponent, should be a scalar or a vector which length==cfg.N
    % and range from 0~1

% [Output] FT

% =========================================================================

function [FT] = ft_fractalSim(cfg)

    FT = [];
    FT.fsample  = 1024;
    FT.label    = {'1'};
    FT.time{1}  = 0:1/1024:1;
    
    % ======================================

    ft_checkopt(cfg,'N','doublescalar');
    ft_checkopt(cfg,'Hurst','doublevector');


    if(mod(cfg.N,1)~=0)
        ft_error("[ft_fractalSim] cfg.N should be an integar.");
    end
    FT.trial{1} = zeros(cfg.N,1025);


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

    for i=1:cfg.N
        [X,T] = fbm1d(cfg.Hurst(i),1024,1);
        FT.trial{1}(i,:) = X;
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
