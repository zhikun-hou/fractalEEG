
% ========================================================================= 
% [Aim]
    % analyze fractal dimension of [TimeSeries] by [Adaptive Fractal Analysis] method

% [Input] X
% time series data
    % [Must] be a vector, can be (1,N_time) or (N_time,1)

% [Input] cfg
% a struct for arguments
% cfg.visualize
    % whether to show the figure or not
    % default: false
% cfg.epsilon
    % a vector to determine which winsize to used for fit local trend
    % [Must] be odd integar to construct overlapped slide window

% [Suggest]
    % Data in actual world can't be perfect fractal,
        % commonly it only obey 1/f distribution in specific scales.
    % So use struct(epsilon=3:2:numel(X),'visualize',true) to check,
        % which part of scales can make log(E)-log(F) curve seems like a line.
    % Then you can use a more sparse setting for your large dataset.
        % for example, struct(epsilon=2*[5:50]+1)

% [Output] D
    % a scalar

% =========================================================================

function [D] = tfd_AFA(cfg,X)
    
    N_time    = numel(X);
    N_epsilon = numel(cfg.epsilon);
    
    if(N_epsilon<2)
        error("[tfd_AFA] cfg.epsilon should be a vector");
    elseif(any(mod(cfg.epsilon,2)~=1))
        error("[tfd_AFA] cfg.epsilon shoud be odd integar");
    end

    % ======================================

    % keep X is always a row vector
    % since time axis in eeglab and fieldtrip are always the 2nd dimension
    if(iscolumn(X))
        X = X';
    end

    % E: short for Epsilon, list of winsize to fit local trend
    E  = cfg.epsilon;
    
    % F: Fluctuation calculated at correspond winsize 
    F = zeros(N_epsilon,1);
    
    % =======================================

    for k=1:N_epsilon
        N_winsize = E(k); % total length of the window
        N_flanker = (N_winsize-1)/2; % half length of the window
        N_window = ceil((N_time-2*N_flanker) / N_flanker); % number of windows with tail dropped
        
        % weights to smooth curve
        W_right = [(0:N_flanker)./N_flanker];
        W_left = 1-W_right;

        % new length of smoothed curve with tail dropped
        N_smooth = (N_window+1)*N_flanker+1; % flanker-1-flanker，其中(1-flanker)与(flanker-1)重叠
        
        % ==========================

        % Fit local trend by overlapped windows

        Trend = zeros(N_window,N_winsize);
        for i=1:N_window
            Idx = 1+N_flanker*(i-1) : 1+N_flanker*(i+1);
            Coeff = polyfit(Idx,X(Idx),1);
            Trend(i,:) = polyval(Coeff,Idx);
        end
        
        % ==========================

        % Weighted sum the overlapped part to build a smooth curve 

        X_smooth = zeros(1,N_smooth);

        % First window and last window can't be smoothed by weighted sum
        
        X_smooth(1:N_flanker) = Trend(1, 1:N_flanker);
        X_smooth(N_smooth-N_flanker+1:end) = Trend(N_window, (N_flanker+1)+1:end); % skip middle point

        for i=1:N_window-1 % with middle point of a window
            Idx = N_flanker*i+1 : N_flanker*(i+1)+1;
            X_smooth(Idx) = W_left .* Trend(i,N_flanker+1:end) + W_right .* Trend(i+1,1:N_flanker+1);
        end

        % =========================

        % Calculate the fluctuation of data from their trend by standard deviation
        
        X_raw = X(1:N_smooth); % drop tail
        F(k) = mean( (X_raw-X_smooth).^2 )^(1/2);

    end

    % =====================================

    % F ~ E^H

    Coeff = polyfit(log(E),log(F),1);
    Hurst = Coeff(1); % slope of the line
    D = 2 - Hurst;

    if(cfg.visualize)
        loglog(E,F,'o');
        hold on;
        loglog(E,exp(polyval(Coeff,log(E))));
        xlabel("log(E)");
        ylabel("log(F)");
        title("Dim="+D+" Hurst="+Hurst+" Method=AFA");
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
% [1] Fractal Analysis of Time-Series Data Sets: Methods and Challenges
% [2] A tutorial introduction to adaptive fractal analysis
% [3] Facilitating Joint Chaos and Fractal Analysis of Biosignals through Nonlinear Adaptive Filtering

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly
