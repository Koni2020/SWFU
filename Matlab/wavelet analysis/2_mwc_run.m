clc;clear;warning off;
% time y x1 x2 x3 x4
%  -  -  -  -  -  -
%  -  -  -  -  -  -
global mcount
mcount = 1; %1 the count of the MonteCarlo
path_file = '\dataset'; %2 the path of the dataset (relative path)
path_out = '\statistic-pwc'; %3 the path that outputs the result
% section one.prepar
path = pwd;
path_out1 = '\awc_pasc';
path_awc_pasc = '\statistic-wtc';
mkdir(strcat(path, path_out));
mkdir(strcat(path, path_out, path_out1));
assert(exist(strcat(path, path_awc_pasc), 'file') == 7, 'Please run (run_wtc.m) first.');
file = dir(strcat(path, path_file, '\*.csv'));
for i = 1:size(file, 1)
    filename = file(i).name;
    ds = readtable(strcat(path, path_file, '\',filename));
    pasc = readtable(strcat(path, path_awc_pasc, '\awc_pasc\',filename), 'ReadRowNames', true);
    pwc_statistic = core_pwc(ds, pasc);
    writetable(pwc_statistic, strcat(path, path_out, path_out1, '\', filename), 'WriteRowNames', true) % awc and pwc
end


function vector=z_standardize(vector)
vector = (vector - mean(vector)) ./ std(vector);
end

function container=core_pwc(data_arr, wtc_pasc)
global mcount
candidate_index = wtc_pasc.Properties.RowNames;
candidate_index = string(candidate_index);
first_index = candidate_index(wtc_pasc(:, 2).Variables == max(wtc_pasc(:, 2).Variables));
candidate_index(candidate_index == first_index) = [];
x = data_arr(:, 2).Variables;
y = data_arr(:, first_index).Variables;
x = z_standardize(x);
y = z_standardize(y);
n = 1;
combine = first_index;
database = [x, y];
while ~isempty(candidate_index)
    
    for k = 1:length(candidate_index)
        disp(strcat("Searching: ",combine, '-', candidate_index(k)))
        data = data_arr(:, candidate_index(k)).Variables;
        try % two
            [Rsq,~,~,~,pwcsig] = pwc([database, data], 'MonteCarloCount', mcount);
        catch
            break
        end
        pwcsig_per=length(pwcsig(pwcsig>=1))/numel(pwcsig);
        middle{k, 1} = candidate_index(k);
        middle{k, 2} = pwcsig_per;
        middle{k, 3} = mean(Rsq(pwcsig > 1));
        middle{k, 4} = [database, data];
        
    end
    if ~exist("middle")
        break
    end
    sig = middle{:, 2};
    pos = find(sig == max(sig));
    apwc = middle{pos, 3};
    
    database = middle{pos, 4};
    combine = strcat(combine, '-', middle{pos,1});
    disp(combine)
    container_best{n, 1} = combine;
    container_best{n, 2} = apwc;
    container_best{n, 3} = max(sig);
    disp(strcat('The best: ', combine))
    n = n + 1;
    candidate_index(candidate_index == middle{pos, 1}) = [];
    clear("middle")
    clear("arc")
    clear("databased")
    
end

container = table('VariableNames', {'Combination','AWC', 'PASC'},...
    'Size', size(container_best),...
    'VariableTypes', {'string','double', 'double'});
container(:, :) = container_best;
end
