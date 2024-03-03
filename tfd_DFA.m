
% ========================================================================= 
% [Aim] 
    % analyze fractal dimension of [TimeSeries] by [Detrend Fluctuation Analysis] method

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
    % So use struct(epsilon=1:10:numel(X),'visualize',true) to check,
        % which part of scales can make log(E)-log(F) curve seems like a line.
    % Then you can use a more sparse setting for your large dataset.
        % for example, struct(epsilon=2.^[3:6])

% [Output] D
    % a scalar

% =========================================================================

function [D,H] = tfd_DFA(cfg,X)

    N_time  = numel(X);
    N_epsilon = numel(cfg.epsilon);

    if(N_epsilon<2)
        error("[tfd_DFA] cfg.epsilon should be a vector");
    elseif(min(cfg.epsilon)<3)
        error("[tfd_DFA] local trend should be fitted by at least 3 points");
    elseif(max(cfg.epsilon)>N_time/10)
        error("[tfd_DFA] local trend shouldn't be fitted by too many points, or it would be global trend.");
    end

    % =====================================

    % Commonly EEG signals are random walk process,
        % so we just perform

    X = X-mean(X);
    
    % =====================================

    E = cfg.epsilon;
    F = zeros(1,N_epsilon);

    for k=1:N_epsilon
        N_lag = E(k);
        
        % Since tail points not enough for a window are dropped,
            % we need balance results by drop head points.

        N_remain = mod(N_time,N_lag);
        if(N_remain==0)
            X_clip = buffer(X,N_lag); % (N_lag,N_window)
        else
            Clips_droptail = buffer(X(1:end-N_remain),N_lag);
            Clips_drophead = buffer(X(1+N_remain:end),N_lag);
            X_clip = [Clips_droptail,Clips_drophead]; % (N_lag,2*N_window)
        end

        % =====================================

        % Detrend raw data for each window
    
        N_window = size(X_clip,2); % number of windows
        Idx = [1:N_lag]';
        
        for i=1:N_window
            Coeff = polyfit(Idx,X_clip(:,i),1);
            X_trend = polyval(Coeff,Idx); % local trend
            X_clip(:,i) = X_clip(:,i) - X_trend; % detrend
        end
        
        % =====================================

        % Calculate Fluctuations
    
        F(k) = mean(  mean(X_clip.^2,1)  ,2)^(1/2);
    end

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
        title("Dim="+D+" Hurst="+Hurst+" Method=FA");
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
% [1] Multifractal detrended fluctuation analysis of nonstationary time series
% [2] Fractal Dimension of Self-Affine Signals: Four Methods of Estimation
% [3] Introduction to Multifractal Detrended Fluctuation Analysis in Matlab

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly
