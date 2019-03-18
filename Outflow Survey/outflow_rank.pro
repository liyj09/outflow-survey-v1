pro outflow_rank

cd,'/home/lee/W3/'

;if ~keyword_set(region) then region='BFS52'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC1'
;if ~keyword_set(Lvrange) then Lvrange=[2,9]
;distance=2
;if ~keyword_set(region) then region='GGMC2'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
if ~keyword_set(region) then region='region_C_III'
if ~keyword_set(Lvrange) then Lvrange=[-40,-20]
distance=2.0
;if ~keyword_set(region) then region='GGMC4'
;if ~keyword_set(Lvrange) then Lvrange=[-3,5]
;distance=2
;if ~keyword_set(region) then region='lynds'
;if ~keyword_set(Lvrange) then Lvrange=[-3,3]
;distance=0.4
;if ~keyword_set(region) then region='west'
;if ~keyword_set(Lvrange) then Lvrange=[-1,4]
;distance=0.6
;if ~keyword_set(region) then region='swallow'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='horn'
;if ~keyword_set(Lvrange) then Lvrange=[12,18]
;distance=3.8
;if ~keyword_set(region) then region='remote'
;if ~keyword_set(Lvrange) then Lvrange=[18,28]
;distance=8.5

cd,'/home/lee/W3/'+region+'/candidates/'

;readcol,'blue_out.cat',bsn,bpeak1,bpeak2,bw0,bw1,bnote,format='I,F,F,F,F,A',stringskip='#'
;readcol,'red_out.cat',rsn,rpeak1,rpeak2,rw0,rw1,rnote,format='I,F,F,F,F,A',stringskip='#'
readcol,'outflowcat.cat',num,glon,glat,bluev0,bluev1,redv0,redv1,tag,bnum,rnum,format='I,F,F,F,F,F,F,A,I,I',stringskip='#'

readcol,'./score/lobe_blue.cat',bsn,bpeak1,bpeak2,bw0,bw1,brank1,brank2,brank3,format='I,F,F,F,F,I,I,I',stringskip='#'
readcol,'./score/lobe_red.cat',  rsn,rpeak1,rpeak2,rw0,rw1,rrank1,rrank2,rrank3,format='I,F,F,F,F,I,I,I',stringskip='#'

openw,para,/get_lun,'lobe_info.cat'
printf,para,'#osn','tag','lsn','Glon','Glat','w0','w1','spec','pv','con','l_rank','bi_rank',$
  format='(3a4,a9,a7,2a8,3a6,2a8)'
;;printf,para,'#####################################################################'
brankA=0
brankB=0
brankC=0
brankD=0
rrankA=0
rrankB=0
rrankC=0
rrankD=0
for i=0, n_elements(num)-1 do begin
  j=where(bsn eq bnum[i])
  if j ne -1 then begin
    bmark=brank1[j]+brank2[j]+brank3[j]
    case bmark of
      9: b_rank='A+'
      8: b_rank='A-'
      7: b_rank='B+'
      6: b_rank='B-'
      5: b_rank='C+'
      4: b_rank='C-'
      3: b_rank='D'
    endcase
    if strmid(b_rank,0,1) eq 'A' then brankA++
    if strmid(b_rank,0,1) eq 'B' then brankB++
    if strmid(b_rank,0,1) eq 'C' then brankC++
    if strmid(b_rank,0,1) eq 'D' then brankD++
  endif
  k=where(rsn eq rnum[i])
  if k ne -1 then begin
    rmark=rrank1[k]+rrank2[k]+rrank3[k]
    case rmark of
      9: r_rank='A+'
      8: r_rank='A-'
      7: r_rank='B+'
      6: r_rank='B-'
      5: r_rank='C+'
      4: r_rank='C-'
      3: r_rank='D'
    endcase
    if strmid(r_rank,0,1) eq 'A' then rrankA++
    if strmid(r_rank,0,1) eq 'B' then rrankB++
    if strmid(r_rank,0,1) eq 'C' then rrankC++
    if strmid(r_rank,0,1) eq 'D' then rrankD++
  endif
  if tag[i] eq 'D' then begin
    mark=round(mean([bmark,rmark]))
    case mark of
      9: rank='A+'
      8: rank='A-'
      7: rank='B+'
      6: rank='B-'
      5: rank='C+'
      4: rank='C-'
      3: rank='D'
    endcase
  endif else if tag[i] eq 'B' then rank=b_rank else if tag[i] eq 'R' then rank=r_rank
  
  if j ne -1 then printf,para,num[i],'B',bsn[j],bpeak1[j],bpeak2[j],bw0[j],bw1[j],brank1[j],brank2[j],brank3[j],b_rank,rank,$
    format='(i4,a4,i4,f9.3,f7.3,2f8.1,3i6,2a8)' $
  else printf,para,num[i],'B','-','-','-','-','-','-','-','-','-',rank,format='(i4,2a4,a9,a7,2a8,3a6,2a8)'
  if k ne -1 then printf,para,'    ','R',rsn[k],rpeak1[k],rpeak2[k],rw0[k],rw1[k],rrank1[k],rrank2[k],rrank3[k],r_rank,'',$
    format='(2a4,i4,f9.3,f7.3,2f8.1,3i6,2a8)' $
  else printf,para,'    ','R','-','-','-','-','-','-','-','-','-',format='(3a4,a9,a7,2a8,3a6,2a8)'
  
endfor
sum=float(brankA+rrankA+brankB+rrankB+brankC+rrankC+brankD+rrankD)
perA=num2str(((brankA+rrankA)/sum)*100,format='(I2)')+'%'
perB=num2str(((brankB+rrankB)/sum)*100,format='(I2)')+'%'
perC=num2str(((brankC+rrankC)/sum)*100,format='(I2)')+'%'
perD=num2str(((brankD+rrankD)/sum)*100,format='(I2)')+'%'
printf,para,'#All:',n_elements(num),'blue:',n_elements(bsn),' red:',n_elements(rsn),'A:'+perA,brankA+rrankA,'B:'+perB,brankB+rrankB,'C:'+perC,brankC+rrankC,'D:'+perD,brankD+rrankD,$
  FORMAT='(7(A8,I3))'

free_lun,para

end