from pathlib2 import Path
import fiona
import rasterio
import rasterio.mask
import numpy as np
def clip_raster(src_filename, dst_filename, mask_filename):
    shapefile = fiona.open(mask_filename, 'r')
    shapes = [feature['geometry'] for feature in shapefile]
    raster = rasterio.open(src_filename)
    raster_extract, transform_extract = rasterio.mask.mask(raster, shapes, crop=True, nodata=254)
    out_meata = raster.meta
    out_meata.update({'driver': 'GTiff',
                      'height': raster_extract.shape[1],
                      'width': raster_extract.shape[2],
                      'transform':transform_extract,
                     'nodata': 254})
    dest = rasterio.open(dst_filename, 'w', **out_meata)
    dest.write(raster_extract)

if __name__ == '__main__':

    '''
    批量按掩膜裁剪
    mask文件夹下放掩膜：shapefile
    raster文件夹下放栅格：tif
    
    '''
    pathCwd = Path.cwd()
    pathMask = pathCwd.joinpath('mask')
    pathRaster = pathCwd.joinpath('raster')
    pathOut = pathCwd.joinpath('result')
    pathOut.mkdir(exist_ok=True)
    for mask in pathMask.glob('*.shp'):
        pathOut1 = pathOut.joinpath(mask.stem)
        pathOut1.mkdir(exist_ok=True)
        for raster in pathRaster.glob('*.tif'):
            clip_raster(raster.__str__(), pathOut1.joinpath(raster.name).__str__(), mask.__str__())
            print('裁剪', raster.stem, ' 成功')
