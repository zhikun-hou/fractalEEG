
% ========================================================================= 
% [Aim]
    % analyze fractal dimension of [TimeSeries] by [Fluctuation Analysis] method

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

function [D,H] = tfd_FA(cfg,X)

    N_epsilon = numel(cfg.epsilon);

    if(N_epsilon<2)
        error("[tfd_FA] cfg.epsilon should be a vector");
    elseif(max(cfg.epsilon)>100)
        error("[tfd_FA] cfg.epsilon are not suggest to be very large since it would be influenced by local trend. 100 is just a random border to raise this notice.");
    end

    % =====================================

    % In Jianbo Gao's paper, he preprocess data with cumsum function,
        % since X is an increment process in his paper,
        % Y=cumsum(X-mean(X)) would be a random walk process.
    % However, EEG data are commonly a random walk process,
        % so we didn't perform cumsum here.

    X = X-mean(X);
   
    % =====================================

    E  = cfg.epsilon;
    F = zeros(1,N_epsilon);

    % In Jianbo Gao's paper, he use offset=0 (see tfd_Higuchi),
        % by modify his method using higuchi's idea,
        % we use different m to downsample.

    for k=1:N_epsilon
        N_lag = E(k);

        F_m = zeros(1,N_lag);

        for m=1:N_lag % for each offset            
            Y = downsample(X,N_lag,m-1);

            % get standard deviation of fluctuations in specific lag
            F_m(m) = mean(diff(Y).^2) ^(1/2);

        end

        % take average from different offset

        F(k) = mean(F_m);
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
% [1] Detection of Low Observable Targets Within Sea Clutter by Structure Function Based Multifractal Analysis

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly

