pro download_wise

;if ~keyword_set(region) then region='BFS52'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC1'
;if ~keyword_set(Lvrange) then Lvrange=[2,9]
;distance=2
;if ~keyword_set(region) then region='GGMC2'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
;if ~keyword_set(region) then region='GGMC3'
;if ~keyword_set(Lvrange) then Lvrange=[5,11]
;distance=2
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
if ~keyword_set(region) then region='region_E_I'
if ~keyword_set(Lvrange) then Lvrange=[-20,-13]
distance=0.6

cd,'/media/alpha/W3/'+region+'/candidates'
L_lim=0.6
def_rmsU=0.25
def_rmsL=0.3
outsz=(8./distance)>5.5
outsz=outsz<7

readcol,'blue_out.cat',bsn,bpeak1,bpeak2,bw0,bw1,format='I,F,F,F,F',stringskip='#'
readcol,'red_out.cat',rsn,rpeak1,rpeak2,rw0,rw1,format='I,F,F,F,F',stringskip='#'
readcol,'outflowcat.cat',num,glon,glat,bluev0,bluev1,redv0,redv1,tag,bnum,rnum,format='I,F,F,F,F,F,F,A,I,I',stringskip='#'

if ~file_test('./wisedata/') then spawn,'mkdir ./wisedata'

for i=0, n_elements(num)-1 do begin
 ; if i ne 40 then continue
  source=num2str(bnum[i])+'-'+num2str(rnum[i])
  wiseBname='wise_4.6_'+source+'.fits'
  wiseGname='wise_12_'+source+'.fits'
  wiseRname='wise_22_'+source+'.fits'
  wiseRGBfits=[wiseBname,wiseGname,wiseRname]
  wisesurvey=['WISE 4.6','WISE 12','WISE 22']
  foreach wisename, wiseRGBfits, index do begin
    if ~file_test('./wisedata/'+wisename) then begin
      ;print,' Downloading '+wisename+' from SURVEY '+wisesurvey[index]
      shellcmd='skvbatch_wget file='+wisename+" position='"+string([glon[i],glat[i]],format='(f7.3,",",f6.3)')+"' Survey='"+wisesurvey[index]+$
        "' Coordinates='Galactic' Projection='Car' Pixels=300 Size="+string([outsz/30.0,outsz/30.0],format='(f7.5,",",f7.5)')
 ;     print,shellcmd
      spawn,shellcmd
      print,num2str(num[i])+': ',num2str((i+(index+1)/3.)*100/n_elements(num),format='(f5.1)')+'%'+' downloaded.'
      ;fileinfo=file_info(wisename)
      ;while fileinfo.size lt 100000 do begin
      ;  spawn, 'rm ./'+wisename
      ;  spawn,shellcmd
      ;  fileinfo=file_info(wisename)
      ;endwhile
      spawn,'mv wise*.fits ./wisedata/'
    endif
  endforeach
endfor

end