
% ========================================================================= 
% [Aim]
    % analyze fractal dimension of [TimeSeries] by [Higuchi] method
    % sometimes FD calculated by higuchi method are named as HFD

% [Input] X
% time series data
    % [Must] be a vector, can be (1,N_time) or (N_time,1)

% [Input] cfg
% a struct for arguments
% cfg.visualize
    % whether to show the figure or not
    % default: false
% cfg.epsilon
    % a vector to determine downsample ratio (skip N_winsize points for each N_winsize+1 points)
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

function [D] = tfd_Higuchi(cfg,X)
    
    N_time  = numel(X);
    N_epsilon = numel(cfg.epsilon);

    if(N_epsilon<2)
        error("[tfd_Higuchi] cfg.epsilon should be a vector");
    elseif(any(mod(cfg.epsilon,1)~=0))
        error("[tfd_Higuchi] cfg.epsilon should be integar");
    elseif(max(cfg.epsilon)>N_time/2)
        error("[tfd_Higuchi] cfg.epsilon > N_time/2");
    end

    % ======================================

    % E: short for Epsilon, list of downsample ratio
    E = cfg.epsilon;

    % L: Length of curves with different downsample ratio
    L = zeros(1,N_epsilon);

    % [Notice] "curve length" in higuchi method is a bad name
    % time series is a 1-D wave point actually
    % so t-axis is omitted when calculate "curve length"
    
    % ======================================

    for k=1:N_epsilon
        N_winsize = E(k);

        % Since we keep 1 point for each N_winsize points,
            % we can get N_winsize curves for difference [offset]
            % from 0 to N_winsize-1

        % ======================================
        
        % Now we compute curve length for each curve

        L_m = zeros(1,N_winsize);

        for m=1:N_winsize % for each offset

            Y = downsample(X,N_winsize,m-1);

            % For different epsilon we drop different length of tail,
                % so we have to use unit_time_normalizer to keep them comparable.
            % Since we have N_time points, total length of time is N_time-1.
            % Since we have floor(N_time-m)/N_winsize windows, we can get
                % actual length by *N_winsize.

            unit_time_normalizer =  (N_time-1) / ( floor((N_time-m)/N_winsize) * N_winsize);

            L_m(m) = sum(abs(diff(Y))) * unit_time_normalizer / N_winsize;

            % For last part of the formula, we use "/N_winsize",
                % it serves as a measure of unit length.

        end
        
        % ======================================

        % Now we get average L of such N_winsize curves

        L(k) = mean(L_m);
    end

    % ======================================

    % L ~ E^(-D)

    Coeff = polyfit(log(E),log(L),1);
    D = -Coeff(1); % negative slope
    Hurst = 2-D;


    if(cfg.visualize)
        loglog(E,L,'o');
        hold on;
        loglog(E,exp(polyval(Coeff,log(E))));
        xlabel("log(E)");
        ylabel("log(L)");
        title("Dim="+D+" Hurst="+Hurst+" Method=Higuchi");
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

% [Ref]
% [1] Approach to an irregular time series on the basis of the fractal theory
% [2] Fractal Dimension Feature as a Signature of Severity in Disorders of Consciousness: An EEG Study
% [3] Comparison of Entropy and Complexity Measures for the Assessment of Depth of Sedation
% [4] Fractal Dimension of Self-Affine Signals: Four Methods of Estimation

% [Environment]
% Matlab    R2022a
% EEGLAB    2022.1
% FieldTrip 20230118

% [Alert] Check the difference between your version and mine if code runs incorrectly
