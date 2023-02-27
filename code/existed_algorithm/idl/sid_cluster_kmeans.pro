;确定初始参数,聚类的类别总数num_class，迭代次数MaxIteration，迭代阈值minChangeThreshold 
PRO sid_cluster_kmeans, num_class, MaxIteration, minChangeThreshold 

  ;选择需要聚类的图像
  infile = 'd:\input'
  if (infile eq '') then return
  ENVI_OPEN_FILE, infile, r_fid=fid1
  ENVI_FILE_QUERY, fid1, dims=dims, nb=nb

  ;选择聚类后聚类图像的存储路径
  outfile = 'd:\sid'
  
  ;为聚类中心等变量分配空间
  cluster_center = DBLARR(nb,num_class,MaxIteration+1)
  mean_index = DBLARR(nb)
  stddev_index = DBLARR(nb)
  data = make_array(dims[2]+1,dims[4]+1,nb,/float)
  dis_center = dblarr(num_class)
  plot_cluster = intarr(dims[2]+1,dims[4]+1)
  
  ;计算每个波段的平均值和标准差
  for i=0, nb-1 do begin
    data[*,*,i] = ENVI_GET_DATA(fid=fid1, dims=dims, pos=i)
    data_band = ENVI_GET_DATA(fid=fid1, dims=dims, pos=i)
    index = where(data_band gt 0.0,count)
    mean_index[i] = mean(data_band[index])
    stddev_index[i] = stddev(data_band[index])
  endfor
  count = 0l
  nsum = 0.0d
  interation_last = 0l
  
  ;***********************************************
  
  ;循环迭代聚类
  for interation=0, MaxIteration-1 do begin
    ;初始类别中心的选择
    if interation eq 0 then begin
        for j=0, num_class-1 do begin
          if num_class ne 0 then begin
            cluster_center[*,j,interation] = mean_index+stddev_index*(2.0*j/(num_class-1)-1.0)
          endif else begin
            cluster_center[*,j,interation] = 0
          endelse
        endfor
    endif
    
     ;图像的迭代聚类
     for n=0, dims[2] do begin
      for j=0, dims[4] do begin
        t_sa = total(data[n,j,*])
        p = data[n,j,*]/t_sa
        for i=0, num_class-1 do begin
            ;计算每个像元与聚类中心的光谱信息离散度sid
            t_sb = total(cluster_center[*,i,interation])           
            q = cluster_center[*,i,interation]/t_sb
            d_a = p*alog10(abs(p/q))
            d_b = q*alog10(abs(q/p))
         dis_center[i] = total(d_a)+total(d_b)
        endfor
        ;sid值越小，距离越近
        index = where(dis_center eq min(dis_center),count1)
        plot_cluster[n,j] = index[0]+1
       endfor
      endfor
    
    ;更新聚类中心
    temp = interation+1
    if temp le MaxIteration then begin
     for numclass=0, num_class-1 do begin
        ;计算每个类别的象元个数
        index = where(plot_cluster eq numclass+1,ncount)
         for band=0, nb-1 do begin
          data1 = ENVI_GET_DATA(fid=fid1, dims=dims, pos=band)
          ;计算每个类别每个波段所有象元的平均值
          cluster_center[band,numclass,interation+1] = mean(data1[index])
         endfor
     endfor 
    endif
       
    ;判断迭代是否结束（以新集群中心和旧集群中心的变化率判断）
    change = 0.0D
    change_sum = 0.0D
    total_c = 0.0D
    
    for i=0, num_class-1 do begin
      if total(cluster_center[*,i,interation]) ne 0.0 then begin
        change = abs(total(cluster_center[*,i,interation+1]-cluster_center[*,i,interation])/total(cluster_center[*,i,interation]))
      endif else begin
        change = 0.0
      endelse
      change_sum = change_sum+change
    endfor
    
    if change_sum LT minChangeThreshold then begin
      interation_last = interation
      interation = MaxIteration+1
    endif else begin
      interation_last = MaxIteration
    endelse
  endfor
  
;******************************************************
  
  ;输出最后的迭代次数和变化的阈值
  print, interation_last
  print, '*********'
  print, change_sum

  mapinfo = ENVI_GET_MAP_INFO(fid=fid1)
  ENVI_WRITE_ENVI_FILE, plot_cluster, out_name=outfile, map_info=mapinfo
  
;  ;写聚类图像
;  OPENW, lun, outfile, /append, /get_lun
;  WRITEU, lun, plot_cluster
;  CLOSE, lun
;  free_lun, lun
;   
;  ;写聚类图像的头文件，波段数为1，数据类型为整型
;  zhu_write_envi_header, infile, outfile, bands=1, data_type=2
  
  print, 'done'
  
end