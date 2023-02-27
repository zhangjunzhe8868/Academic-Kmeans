# algorithm of scwm (k-means)

import gc
import numpy as np
from osgeo import gdal

num_class = 3
MaxIteration = 2
minChangeThreshold = 0.1

# def read_img(filename):
#    dataset=gdal.Open(filename)
#
#    im_width = dataset.RasterXSize
#    im_height = dataset.RasterYSize
#
#    im_geotrans = dataset.GetGeoTransform()
#    im_proj = dataset.GetProjection()
#    im_data = dataset.ReadAsArray(0,0,im_width,im_height)
#
#    del dataset
#    return im_proj,im_geotrans,im_data

img = gdal.Open(r"C:\Program Files\Exelis\ENVI53\data\qb_boulder_msi")

im_width = img.RasterXSize
im_height = img.RasterYSize
im_band = img.RasterCount
im_geotrans = img.GetGeoTransform()
im_proj = img.GetProjection()
im_data = img.ReadAsArray(0, 0, im_width, im_height)
band_data = np.reshape(im_data, (im_band, int(im_data.size / im_band)))

del img
gc.collect()

data_temp1 = np.loadtxt(r"D:\cv_dis.txt", delimiter=',')
data_temp2 = np.loadtxt(r"D:\cv_sam.txt", delimiter=',')

print(data_temp1.shape)
print(data_temp2.shape)

cluster_center = np.zeros((im_band, num_class, MaxIteration+1))
mean_index = np.zeros((im_band))
stddev_index = np.zeros((im_band))
dis_center = np.zeros((num_class))
plot_cluster = np.zeros((im_height, im_width))

for i in range(im_band):
    b = np.arange(len(band_data[i, :]))
    index = b[np.where(band_data[i, :] > 0)]
    mean_index[i] = np.mean(band_data[i, index])
    stddev_index[i] = np.std(band_data[i, index])

interation_last = 0.0
for interation in range(MaxIteration):
    if interation == 0:
        for j in range(num_class):
            if num_class != 0:
                cluster_center[:, j, interation] = (mean_index +
                                                    stddev_index *
                                                    (2.0 * j / (num_class - 1)
                                                     - 1.0))
            else:
                cluster_center[:, j, interation] = 0
    print(cluster_center)
    for n in range(im_height):
        for m in range(im_width):
            T_SA = np.sum(im_data[:, n, m] ** 2)
            for i in range(num_class):
                dis = np.sum((im_data[:, n, m] -
                              cluster_center[:, i, interation]) ** 2)
                T_AB = np.sum(im_data[:, n, m] *
                              cluster_center[:, i, interation])
                T_SB = np.sum(cluster_center[:, i, interation] ** 2)
                dis_center[i] = (0.5) * np.sqrt(dis) + (0.5) * \
                                (1.0 - T_AB / (np.sqrt(T_SA * T_SB)))
            b = np.arange(len(dis_center))
            index = b[np.where(dis_center == min(dis_center))]
            plot_cluster[n, m] = index[0] + 1
    new = np.reshape(plot_cluster, (plot_cluster.size))
    temp = interation + 1
    if temp <= MaxIteration:
        for numclass in range(num_class):
            b = np.arange(len(new))
            index = b[np.where(new == numclass + 1)]
            for band in range(im_band):
                cluster_center[band, numclass,
                               interation+1] = np.mean(band_data[band, index])

    change = 0.0
    change_sum = 0.0

    for i in range(num_class):
        if np.sum(cluster_center[:, i, interation]) != 0:
            change = np.absolute(np.sum(cluster_center[:, i, interation+1] -
                                        cluster_center[:, i, interation]) /
                                 np.sum(cluster_center[:, i, interation]))
        else:
            change = 0.0
        change_sum = change_sum + change
        print(change_sum)

    if change_sum < minChangeThreshold:
        interation_last = interation
        interation = MaxIteration + 1
    else:
        interation_last = MaxIteration

del cluster_center, dis_center, band_data, im_data
gc.collect()

# dis=np.reshape(data_temp1,(1,data_temp1.size))
# sam=np.reshape(data_temp2,(1,data_temp2.size))
# for i in range(num_class):
#     temp=np.reshape(plot_cluster,(plot_cluster.shape[0]*plot_cluster.shape[1]))
#     b=np.arange(len(temp))
#     position=b[np.where(temp == i+1)]
#     print('dis_weight',np.mean(dis[position]))
#     print('sam_weight',np.mean(sam[position]))

# outfile=r"d:\test.tif"
# out=ga.SaveArray(plot_cluster,outfile,format = "GTiff")
# out=None


def write_img(filename, im_proj, im_geotrans, im_data):

    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    if len(im_data.shape) == 3:
        im_band, im_height, im_width = im_data.shape
    else:
        im_band, (im_height, im_width) = 1, im_data.shape

    driver = gdal.GetDriverByName("ENVI")
    dataset = driver.Create(filename, im_width, im_height, im_band, datatype)

    dataset.SetGeoTransform(im_geotrans)
    dataset.SetProjection(im_proj)

    if im_band == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)
    else:
        for i in range(im_band):
            dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset


write_img(r"d:\test.dat", im_proj, im_geotrans, plot_cluster)


print('********')
print(interation_last+1)
print('********')
print(change_sum)
