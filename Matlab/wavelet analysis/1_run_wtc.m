clc;clear;warning off;
% time y x1 x2 x3 x4
%  -  -  -  -  -  -
%  -  -  -  -  -  -
global mcount
mcount = 1000; %1 the count of the MonteCarlo
path_file = '\dataset'; %2 the path of the dataset (relative path)
path_out = '\statistic-wtc'; %3 the path that outputs the result
% section one.prepar
path = pwd;
path_out1 = '\awc_pwc';
path_out2 = '\pahse_period';
mkdir(strcat(path, path_out));
mkdir(strcat(path, path_out, path_out1));
mkdir(strcat(path, path_out, path_out2))
% section two.batch run wavelet coherence
file = dir(strcat(path, path_file, '\*.csv'));
for i = 1:size(file, 1)
    filename = file(i).name;
    ds = readtable(strcat(path, path_file, '\',filename));
    statistic = core_wtc(ds);
    detail_info = detail_wtc(ds);
    writetable(statistic, strcat(path, path_out, path_out1, '\', filename), 'WriteRowNames', true) % awc and pwc
    writetable(detail_info, strcat(path, path_out, path_out2, '\', filename), 'WriteRowNames', true) % awc and pwc
end


function container=core_wtc(data_arr)
global mcount
name_y = data_arr(:, 3:end).Properties.VariableNames;
name_x = ["AWC", "PASC"];
container = table('Size', [numel(name_y), 2],...
    'VariableTypes',["double", "double"],...
    'VariableNames', name_x,...
    'RowNames', name_y);
x_ = data_arr(:, 2).Variables;
for j = 3:size(data_arr, 2)
    y_ = data_arr(:, j).Variables;
    name = data_arr(:, j).Properties.VariableNames;
    [rsq,~,~,~,wtcsig,~] = wtc(x_, y_, 'MonteCarloCount', mcount); % 3 the count of MonteCarlo
    awc = mean(rsq, 'all');
    pwc = sum(wtcsig > 1, 'all') / numel(wtcsig);
    container(name, 'AWC') = num2cell(awc);
    container(name, 'PASC') = num2cell(pwc);
end
end

function container2table=detail_wtc(data_arr)
global mcount
x_ = data_arr(:, 2).Variables;
for j = 3:size(data_arr, 2)
    date = data_arr.time;
    y_ = data_arr(:, j).Variables;
    % main -------------
    [Wxy,~,~,~,~] = xwt(x_, y_);
    [~,period,~,~,wtcsig, t] = wtc(x_, y_, 'MonteCarloCount', mcount); % "MonteCarloCount" denotes the count of the monte carlo
    [coherence, ptile]  =globalcoher(x_, y_, 'MonteCarloCount', mcount);
    loc = find(diff(sign(diff(coherence)))==-2)+1; % find the local maximum
    c1 = coherence(loc);
    p1 = ptile(loc);
    p1 = p1';
    loc = loc((c1 - p1)>0); % the local maximum at 0.05 level
%     coherence95 = coherence(loc);
    resonance = period(loc);
    
    pp = 1;
    iterator = [];
    for k=resonance
        wtcsig_one_scale = wtcsig((abs(period-k) == min(abs(period-k))), :);
        
        Wxy_one_scale = Wxy((abs(period-k) == min(abs(period-k))), :);
%         Rsq_one_scale = Rsq((abs(period-k) == min(abs(period-k))), :);
        t_significant = t(wtcsig_one_scale >= 1);  % "t_significant" denotes the significant time period at specific scale.
%         awc = Rsq_one_scale(wtcsig_one_scale >= 1); % average of wavelet coherence at specific scale
        t_significant = unique(t_significant);
        m=1; n=1;
        container = [];
        
        for l=1:length(t_significant)
            if l==length(t_significant)
                container(m, n) = t_significant(l);
            elseif t_significant(l + 1) - t_significant(l) == 1
                container(m, n) = t_significant(l);
                n = n + 1;
            else
                container(m ,n) = t_significant(l);
                n = 1;
                m = m + 1;
            end
        end
        
        if isempty(container) || size(container, 2) <= 1
            continue
        else
            container = container(container(:, 2) ~= 0, :);
        end
        
        for l=1:size(container, 1)
            event = container(l, :);
            event = event(event ~= 0);
            
            wxy = Wxy_one_scale(ismember(t_significant, event'));
            wxy = angle(wxy);
            meantheta = anglemean(wxy);
            degree = rad2deg(meantheta);
            lag = k * meantheta / (2 * pi); % convet the phase angle into time domain
            %         iterator(pp, 1) = string(round(k/12, 1)); if sample resolution is monthly, suggest divided by 12 that convert into yearly.
            iterator{pp, 1} = round(k, 1);
            iterator{pp, 2} = strcat(string(date(1)), '_', string(date(end)));
            iterator{pp, 3}= round(degree, 1);
            iterator{pp, 4} = round(lag, 1);
            pp = pp + 1;
        end
        
        
    end
    try
        container2table = table('VariableNames', {'Coherence period (t_significant)', 't_significants', 'Phase difference', 'Leads/lags'},...
            'Size', size(iterator),...
            'VariableTypes', {'double', 'string', 'double', 'double'}) ;
    catch
        a = 1
    end
    container2table(:, :) = iterator;
end

end

