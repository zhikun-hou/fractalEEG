
% ========================================================================= 
% [Aim]
    % analyze fractal dimension of [TimeSeries] by [Dubuc's] method

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

function [D] = tfd_Dubuc(cfg,X)
    
    N_time    = numel(X);
    N_epsilon = numel(cfg.epsilon);

    if(N_epsilon<2)
        error("[tfd_Dubuc] cfg.epsilon should be a vector");
    elseif(min(cfg.epsilon)<1 || max(cfg.epsilon)>N_time/2)
        error("[tfd_Dubuc] cfg.epsilon too small or too large.");
    end

    % =====================================

    E = cfg.epsilon;
    V = zeros(1,N_epsilon);

    Upper = zeros(1,N_time);
    Lower = zeros(1,N_time);
    Idx_total = 1:N_time;

    for k=1:N_epsilon
        N_flanker = E(k); % half length of neighbor window

        for i=1:N_time

            % get neighbors and their lower/upper bound

            Idx_neighbor = intersect(Idx_total,i-N_flanker:i+N_flanker);
            X_neighbor   = X(Idx_neighbor);

            Upper(i) = max(X_neighbor);
            Lower(i) = min(X_neighbor);
        end

        % Accumulate length in unit time and measure by unit length

        V_t  = Upper-Lower;
        V(k) = sum(V_t) / N_flanker^2;
    end

    % V ~ E^(-D)

    Coeff = polyfit(log(E),log(V),1);
    D = -Coeff(1); % slope of the line
    Hurst = 2 - D;

    if(cfg.visualize)
        loglog(E,V,'o');
        hold on;
        loglog(E,exp(polyval(Coeff,log(E))));
        xlabel("log(E)");
        ylabel("log(V)");
        title("Dim="+D+" Hurst="+Hurst+" Method=Dubuc");
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
% [2] Evaluating the fractal dimension of profiles

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly
