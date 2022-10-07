from pathlib import Path

import numpy as np
import xarray as xr
from datetime import timedelta

# params
max_nan_consecutive = 31  # 最大连续缺测
max_nan_rate = 0.05  # 最大缺测占比
time_period = [1961, 2020]  # 提取时间段

path = Path.cwd()
filename_output = 'observation_interpolation' + '_' + str(max_nan_consecutive) + '_' + str(
    int(max_nan_rate * 100)) + '_' + str(time_period[0]) + '-' + str(time_period[-1])
path_out_root = path.joinpath(filename_output)
path_out_root.mkdir(exist_ok=True)
ds = xr.open_dataset('china_observation.nc')
n = ds.time.size  # the number of sample
for i in ds.variables:
    if i in ['id', 'time']:
        continue
    else:
        path_out = path_out_root.joinpath(i)
        path_out.mkdir(exist_ok=True)
        data_arr = ds[i]
        years = data_arr['time.year']
        data_arr = data_arr.sel(time=((years >= time_period[0]) & (years <= time_period[1])))
        data_arr = data_arr.interpolate_na(dim='time', max_gap=timedelta(days=int(np.floor(n * max_nan_rate))), limit=max_nan_consecutive)
        data_arr = data_arr.dropna(dim='id')
        if data_arr.size == 0:
            print(f'{i} Too many nan, continue.')
            continue
        else:
            for j in data_arr.id:
                data_arr_ = data_arr.sel(id=j)
                data_arr_ = data_arr_.to_series()
                data_arr_.to_csv(path_out.joinpath(str(j.values) + '.csv'))
                print(f'{i}, {j.values}, successful.')

