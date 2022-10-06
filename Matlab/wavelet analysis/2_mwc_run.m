clc;clear;warning off;
% time y x1 x2 x3 x4
%  -  -  -  -  -  -
%  -  -  -  -  -  -
global mcount
mcount = 1; %1 the count of the MonteCarlo
path = pwd;
path_file = '\dataset'; %2 the path of the dataset (relative path)
path_out = '\statistic-mwc'; %3 the path that outputs the result
path_awc_pasc = '\statistic-wtc';
assert(exist(strcat(path, path_awc_pasc), 'file') == 7, 'Please run (run_wtc.m) first.');
% section one.prepar
path_out1 = '\awc_pwc';
mkdir(strcat(path, path_out));
mkdir(strcat(path, path_out, path_out1));
% section two.batch run wavelet coherence
file = dir(strcat(path, path_file, '\*.csv'));

for i = 1:size(file, 1)
    filename = file(i).name;
    ds = readtable(strcat(path, path_file, '\',filename));
    pasc = readtable(strcat(path, path_awc_pasc, '\awc_pasc\',filename), 'ReadRowNames', true);
    mwc_statistic = core_mwc(ds, pasc);
end

function vector=z_standardize(vector)
vector = (vector - mean(vector)) ./ std(vector);
end

function container=core_mwc(data_arr, wtc_pasc)
global mcount
candidate_index = wtc_pasc.Properties.RowNames;
candidate_index = string(candidate_index);
first_index = candidate_index(wtc_pasc(:, 2).Variables == max(wtc_pasc(:, 2).Variables));
candidate_index(candidate_index == first_index) = [];
x = data_arr(:, 2).Variables;
y = data_arr(:, first_index).Variables;
x = z_standardize(x);
y = z_standardize(y);
database = [x, y];
combine = first_index;
n = 1;
while ~isempty(candidate_index)
    
    for k = 1:length(candidate_index)
        disp(strcat("Searching: ",combine, '+', candidate_index(k)))
        data = data_arr(:, candidate_index(k)).Variables;
        try % two
            [Rsq,~,~,~,mwcsig] = mwc([database, data], 'MonteCarloCount', mcount);
        catch
            break
        end
        mwcsig_per=length(mwcsig(mwcsig>=1))/numel(mwcsig);
        middle{k, 1} = candidate_index(k);
        middle{k, 2} = mwcsig_per;
        middle{k, 3} = mean(Rsq(mwcsig > 1));
        middle{k, 4} = [database, data];
        
    end
    if ~exist("middle")
        break
    end
    sig = middle{:, 2};
    pos = find(sig == max(sig));
    apwc = middle{pos, 3};
    
    database = middle{pos, 4};
    combine = strcat(combine, '+', middle{pos,1});
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
