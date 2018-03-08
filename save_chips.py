import sys

from imagery_utils import create_chips_from_features


if __name__ == "__main__":
    
    data_path = '/media/rmsare/GALLIUMOS/data/ot_data/tif/'
    faults_path = '/media/rmsare/GALLIUMOS/data/qfaults/qfaults_norcal.shp'

    create_chips_from_features(faults_path, data_path)
