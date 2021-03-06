import os
import numpy as np
import fiona, rasterio, shapely
import matplotlib

import rasterio.features

from shapely.geometry import shape

from spatial_utils import intersect

import matplotlib.pyplot as plt

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

def create_training_images(chips, name):
    """
    Saves training images using provided image chips 
    """

    z = chips[0]
    images = chips[1::]

    hs = hillshade(z)
    plt.figure()
    plt.imshow(hs, alpha=1, cmap='gray', origin='lower')
    plt.axis('off')
    plt.savefig('data/hs/' + name + '.png', dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)

    #slope_x, slope_y = calculate_slope(z)
    #plt.figure()
    #plt.imshow(slope_x, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/slope/' + name + '_x.png', dpi=300, bbox_inches='tight')

    #plt.figure()
    #plt.imshow(slope_y, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/slope/' + name + '_y.png', dpi=300, bbox_inches='tight')

    #curv = calculate_curvature(z)
    #plt.figure()
    #plt.imshow(curv, alpha=1, cmap='RdBu', origin='lower')
    #plt.axis('off')
    #plt.savefig('data/curv/' + name + '.png', dpi=300, bbox_inches='tight')

    #for image in images:
        #plt.figure()
        #plt.imshow(image['data'], alpha=1, origin='lower')
        #plt.axis('off')
        #plt.savefig('data/sat/' + name +  '_' + image['id'] + '.png', dpi=300, bbox_inches='tight')

    plt.close('all')

def create_chips_from_features(features_path, data_path):
    
    features = fiona.open(features_path, 'r')

    files = os.listdir(data_path)
    files = [f for f in files if 'tif' in f and 'xml' not in f]
    images = []

    pbar = ProgressBar(widgets=[Percentage(), ' ', Bar(), ' ', ETA()], maxval=len(files))
    pbar.start()
    
    for i, f in enumerate(files):
        with rasterio.open(data_path + f) as dataset:
            for feature in features:
                try:
                    z, name = make_chip(feature, dataset)
                    #images = download_imagery(feature)
                    chips = [z] + images
                    create_training_images(chips, name=name)
        pbar.update(i+1)

