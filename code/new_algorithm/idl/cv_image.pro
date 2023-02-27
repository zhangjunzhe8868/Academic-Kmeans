PRO cv_image, infile, roi, cv_dis, cv_sam
 
 ENVI_OPEN_FILE, infile, r_fid=fid, /INVISIBLE, /NO_INTERACTIVE_QUERY, /NO_REALIZE
 ENVI_FILE_QUERY, fid, dims=dims, nl=nl, ns=ns, nb=nb, data_ignore_value=nan
 data = make_array(dims[2]+1,dims[4]+1,nb,/float)
 for i=0, nb-1 do begin
   data[*,*,i] =  ENVI_GET_DATA(fid=fid, dims=dims, pos=i)
 endfor

 ENVI_RESTORE_ROIS, roi
 roi_ids = ENVI_GET_ROI_IDS(fid=fid, roi_colors=roi_colors)
 num_roi = n_elements(roi_ids)
 ENVI_GET_ROI_INFORMATION, roi_ids

 data_roi = make_array(num_roi,nb,/float)
 for j=0, num_roi-1 do begin
   for i=0, nb-1 do begin
     a = ENVI_GET_ROI_DATA(roi_ids[j], fid=fid, pos=i)
     temp = mean(a)
     data_roi[j,i] = temp 
   endfor
 endfor
      
  ;分配变量
  sam_center_org = dblarr(ns,nl,num_roi)
  sam_center = dblarr(ns,nl,num_roi)
  dis_center_org = dblarr(ns,nl,num_roi)
  dis_center = dblarr(ns,nl,num_roi)
  mean_index_sam = dblarr(ns,nl)
  stddev_index_sam = dblarr(ns,nl)
  mean_index_dis = dblarr(ns,nl)
  stddev_index_dis = dblarr(ns,nl) 
  cv_dis = FLTARR(ns,nl)
  cv_sam = FLTARR(ns,nl)
    
  ;计算每个波段欧氏距离和光谱角余弦的的平均值和标准差
  for a=0, ns-1 do begin
   for b=0, nl-1 do begin
    t_sb = total(data[a,b,*]^2.0)
    for j=0, num_roi-1 do begin    
         ;计算每个像元与聚类中心的光谱余弦值sam
         t_sa = total(data[a,b,*]*data_roi[j,*])
         t_sc = total(data_roi[j,*]^2.0)
         sam_center_org[a,b,j] = t_sa/(sqrt(t_sb*t_sc))   
         sam_center[a,b,j] = 1.0-sam_center_org[a,b,j]

         ;计算每个像元与聚类中心的光谱距离dis
         distance = total((data[a,b,*]-data_roi[j,*])^2.0)
         dis_center[a,b,j] = sqrt(distance)
    endfor
   endfor
  endfor

  for a=0, ns-1 do begin
    for b=0, nl-1 do begin
         mean_index_sam[a,b] = mean(sam_center[a,b,*])
         stddev_index_sam[a,b] = stddev(sam_center[a,b,*])
         mean_index_dis[a,b] = mean(dis_center[a,b,*])
         stddev_index_dis[a,b] = stddev(dis_center[a,b,*])
    
         cv_dis[a,b] = stddev_index_dis[a,b]/mean_index_dis[a,b]
         cv_sam[a,b] = stddev_index_sam[a,b]/mean_index_sam[a,b]
         
    endfor
  endfor
  
  ENVI_FILE_MNG, id=fid, /remove
  print,'done'
end