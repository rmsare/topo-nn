"""
Convenience functions for working raster and vector data
"""

import fiona, rasterio, shapely

import rasterio.features

from shapely.geometry import shape


def intersects(feature, dataset):
    """
    Returns True if vector feature and raster bounding boxes intersect.
    """

    feature_left, feature_bottom, feature_right, feature_top = rasterio.features.bounds(feature)
    raster_left, raster_bottom, raster_right, raster_top = dataset.bounds

    intersect_bottom = feature_bottom <= raster_top and feature_bottom >= raster_bottom
    intersect_top = feature_top >= raster_bottom and feature_top <= raster_top
    intersect_left = feature_left <= raster_right and feature_left >= raster_left
    intersect_right = feature_right >= raster_left and feature_right <= raster_right

    intersect_x = intersect_left or intersect_right
    intersect_y = intersect_bottom or intersect_top

    if intersect_x and intersect_y:
        return True

    return False
