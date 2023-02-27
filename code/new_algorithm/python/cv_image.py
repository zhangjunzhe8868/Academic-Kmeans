# compute the weight of two kernals

import gc
import numpy as np
from osgeo import gdal

img = gdal.Open(r"C:\Program Files\Exelis\ENVI53\data\qb_boulder_msi")

im_width = img.RasterXSize
im_height = img.RasterYSize
im_band = img.RasterCount
im_geotrans = img.GetGeoTransform()
im_proj = img.GetProjection()
im_data = img.ReadAsArray(0, 0, im_width, im_height) / 10000

print(im_data.dtype.name)

del img
gc.collect()

# ROI_data=np.loadtxt(r"C:\Users\zhang\Dropbox\roi.txt")
# ROI1=ROI_data[0:625,7:]/10000
# ROI2=ROI_data[625:1250,7:]/10000
# ROI3=ROI_data[1250:1875,7:]/10000
# ROI4=ROI_data[1875:2500,7:]/10000
# ROI5=ROI_data[2500:3125,7:]/10000
# ROI6=ROI_data[3125:3750,7:]/10000

ROI1 = im_data[:, 0:10, 0:10]
ROI2 = im_data[:, 25:35, 25:35]
ROI3 = im_data[:, 100:110, 100:110]
ROI4 = im_data[:, 175:185, 175:185]
ROI5 = im_data[:, 250:260, 250:260]
ROI6 = im_data[:, 325:335, 325:335]

temp = [ROI1, ROI2, ROI3, ROI4, ROI5, ROI6]
ROI = np.asarray(temp)
print(ROI.shape)

# ROI=np.zeros((6,im_band,100))
# ROI[0,:,:]=np.reshape(ROI1,(im_band,100))
# ROI[1,:,:]=np.reshape(ROI2,(im_band,100))
# ROI[2,:,:]=np.reshape(ROI3,(im_band,100))
# ROI[3,:,:]=np.reshape(ROI4,(im_band,100))
# ROI[4,:,:]=np.reshape(ROI5,(im_band,100))
# ROI[5,:,:]=np.reshape(ROI6,(im_band,100))
num_roi = ROI.shape[0]

del ROI1, ROI2, ROI3, ROI4, ROI5, ROI6
gc.collect()

data_roi = np.zeros((num_roi, im_band))
for j in range(num_roi):
    for i in range(im_band):
        data_roi[j, i] = np.mean(ROI[j, i, :])

sam_center = np.zeros((num_roi))
dis_center = np.zeros((num_roi))
mean_index_sam = np.zeros((im_height, im_width))
stddev_index_sam = np.zeros((im_height, im_width))
mean_index_dis = np.zeros((im_height, im_width))
stddev_index_dis = np.zeros((im_height, im_width))
cv_dis = np.zeros((im_height, im_width))
cv_sam = np.zeros((im_height, im_width))

for a in range(im_height):
    for b in range(im_width):
        T_SB = np.sum(im_data[:, a, b] ** 2)
        for j in range(num_roi):
            T_SA = np.sum(im_data[:, a, b] * data_roi[j, :])
            T_SC = np.sum(data_roi[j, :] ** 2)
            sam_center[j] = 1.0 - (T_SA / (np.sqrt(T_SB * T_SC)))
            distance = np.sum((im_data[:, a, b] - data_roi[j, :]) ** 2)
            dis_center[j] = np.sqrt(distance)

        mean_index_sam[a, b] = np.mean(sam_center)
        stddev_index_sam[a, b] = np.std(sam_center)
        mean_index_dis[a, b] = np.mean(dis_center)
        stddev_index_dis[a, b] = np.std(dis_center)
        cv_dis[a, b] = stddev_index_dis[a, b] / mean_index_dis[a, b]
        cv_sam[a, b] = stddev_index_sam[a, b] / mean_index_sam[a, b]

print(np.mean(cv_dis), np.mean(cv_sam))

np.savetxt('cv_dis.txt', cv_dis, delimiter=',', fmt='%10.5f')
np.savetxt('cv_sam.txt', cv_sam, delimiter=',', fmt='%10.5f')

# mean_index_sam=np.zeros((im_width,im_height))
# stddev_index_sam=np.zeros((im_width,im_height))
# mean_index_dis=np.zeros((im_width,im_height))
# stddev_index_dis=np.zeros((im_width,im_height))
# cv_dis= np.zeros((im_width,im_height))
# cv_sam= np.zeros((im_width,im_height))
#
# for a in range(im_width):
#     for b in range(im_height):
#         mean_index_sam[a,b]=np.mean(sam_center[a,b,:])
#         stddev_index_sam[a,b]=np.std(sam_center[a,b,:])
#         mean_index_dis[a,b]=np.mean(dis_center[a,b,:])
#         stddev_index_dis[a,b]=np.std(dis_center[a,b,:])
#
#         cv_dis[a,b]= stddev_index_dis[a,b]/mean_index_dis[a,b]
#         cv_sam[a,b]= stddev_index_sam[a,b]/mean_index_sam[a,b]
#
# np.savetxt('cv_dis.txt',cv_dis,delimiter=',',fmt='%10.5f')
# np.savetxt('cv_sam.txt',cv_sam,delimiter=',',fmt='%10.5f')
