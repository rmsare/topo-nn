import os
import numpy as np
import fiona, rasterio, shapely
import matplotlib

import rasterio.features

from shapely.geometry import shape
from spatail_utils import intersects

import matplotlib.pyplot as plt

from copy import copy

from progressbar import ETA, Bar, Percentage, ProgressBar



def make_chip(feature, raster):
    """
    Extracts an image chip of area around feature from input raster.
    """

    if intersects(feature, raster):
        geometry = shape(feature['geometry'])
        ul = raster.index(*geometry.bounds[0:2])
        lr = raster.index(*geometry.bounds[2:4])
        
        window = ((lr[0], ul[0] + 1), (ul[1], lr[1] + 1))

        try:
            data = raster.read(1, window=window)
        except ValueError:
            raise ValueError("Feature bounding box {} not contained in raster bounds".format(window))

        coords = map(str, ul + lr)
        name = '_'.join(coords)

        return data, name
    else:
        raise ValueError("Feature and raster do not intersect")

def calculate_slope(z, dx, dy):
    """
    Calculates topographic slope in x and y directions.
    """

    pass

def calculate_curvature(z, dx, dy):
    """
    Calculates topographic curvature in x and y directions.
    """

    pass

def hillshade(z, az=315, elev=45, dx=0.5, dy=0.5):
    """
    Produces a grayscale hillshade image of input elevation grid.
    """

    ls = matplotlib.colors.LightSource(azdeg=az, altdeg=elev)
    hs = ls.hillshade(z, vert_exag=1, dx=dx, dy=dy)      
    return hs

def download_imagery(feature):
    """
    Downloads image chips around a feature using available satellite imagery  
    """
    
    pass

def divide_raster(filename, size=1024, directory='divided/'):

    with rasterio.open(filename) as dataset:
        ny, nx = dataset.shape
        for i in range(ny / size):
            for j in range(nx / size):
                r0 = i * size
                r1 = (i + 1) * size
                c0 = j * size
                c1 = (j + 1) * size

                # XXX: changed to affine 
                transform = list(copy(dataset.transform))
                shifted_transform = list(copy(dataset.transform))
                shifted_transform[0] += c0
                shifted_transform[3] -= r0

                window = ((r0, r1), (c0, c1))
                z = dataset.read(1, window=window)
                z[z < 0] = np.nan # XXX: quick hack to deal with nodata mis-specification
                fn = '{}_{}.tif'.format(r0, c0)

                with rasterio.open(directory + fn, 
                        'w', 
                        driver='GTiff',
                        transform=shifted_transform,
                        width=size,
                        height=size,
                        count=1,
                        dtype=rasterio.float32) as out:
                    out.write(z.astype(np.float32), 1)


def create_training_images(z, name):
    """
    Saves training images using provided image chips 
    """

    hs = hillshade(z)

    plt.figure()
    plt.imshow(hs, alpha=1, cmap='gray')
    plt.axis('off')
    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    plt.margins(0,0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())

    plt.savefig('data/hs/' + name + '.jpg', dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    
    plt.imshow(z, alpha=0.5, cmap='terrain', vmin=600, vmax=800)
    plt.savefig('data/color/' + name + '.jpg', dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)

    #slope_x, slope_y = calculate_slope(z)
    #plt.figure()
    #plt.imshow(slope_x, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/slope/' + name + '_x.jpg', dpi=300, bbox_inches='tight')

    #plt.figure()
    #plt.imshow(slope_y, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/slope/' + name + '_y.jpg', dpi=300, bbox_inches='tight')

    #curv = calculate_curvature(z)
    #plt.figure()
    #plt.imshow(curv, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/curv/' + name + '.jpg', dpi=300, bbox_inches='tight')

    #for image in images:
        #plt.figure()
        #plt.imshow(image['data'], alpha=1, origin='lower')
        #plt.axis('off')
        #plt.savefig('data/sat/' + name +  '_' + image['id'] + '.jpg', dpi=300, bbox_inches='tight')

    plt.close('all')

def save_labelled_chips(path, out_path, label_values=None):
    """
    Save labelled pixels corresponding to labels of image chips
    """
    
    tif_path = 'data/divided/'
    subdirectories = os.listdir(path)

    if not label_values:
        label_values = np.arange(1, len(subdirectories) + 1)

    for d, val in zip(subdirectories, label_values):
        img_path = os.path.join(path, d)
        files = os.listdir(img_path)
        for f in files:
            f = f.replace('jpg', 'tif')
            infile = os.path.join(tif_path, f)
            with rasterio.open(infile, 'r') as r:
                profile = r.profile
                profile.update(dtype=np.uint8, nodata=0)
                arr = val * np.ones(r.shape)
            outfile = os.path.join(out_path, f)
            with rasterio.open(outfile, mode='w', **profile) as out:
                out.write(arr.astype(np.uint8), 1)

def main(data_file, chips_path):

    print('splitting raster...')
    divide_raster(data_file, size=256, directory=chips_path)

    print('saving images...')
    files = os.listdir(chips_path)

    for f in files:
        with rasterio.open(chips_path + f, 'r') as dataset:
            fn = f.replace('.tif', '')
            z = dataset.read(1)
            if np.sum(np.isnan(z)) == 0:
                create_training_images(z, fn)


if __name__ == "__main__":

    main('data/output_be.tif', 'data/divided/')
