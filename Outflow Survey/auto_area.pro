function walk,arr,pos,footprint,count,limit
sz=size(arr)
cc=count
x=pos[0]
y=pos[1]
all=[[x,y+1],[x-1,y],[x,y-1],[x+1,y]]
block=$
(all[0,0] eq -1 || all[0,0] eq sz[1] || all[1,0] eq -1 || all[1,0] eq sz[2] || (where(all[*,0] eq footprint))[0] ne -1 || arr[all[0,0]>0<(sz[1]-1),all[1,0]>0<(sz[2]-1)] lt limit) && $
(all[0,1] eq -1 || all[0,1] eq sz[1] || all[1,1] eq -1 || all[1,1] eq sz[2] || (where(all[*,1] eq footprint))[0] ne -1 || arr[all[0,1]>0<(sz[1]-1),all[1,1]>0<(sz[2]-1)] lt limit) && $
(all[0,2] eq -1 || all[0,2] eq sz[1] || all[1,2] eq -1 || all[1,2] eq sz[2] || (where(all[*,2] eq footprint))[0] ne -1 || arr[all[0,2]>0<(sz[1]-1),all[1,2]>0<(sz[2]-1)] lt limit) && $
(all[0,3] eq -1 || all[0,3] eq sz[1] || all[1,3] eq -1 || all[1,3] eq sz[2] || (where(all[*,3] eq footprint))[0] ne -1 || arr[all[0,3]>0<(sz[1]-1),all[1,3]>0<(sz[2]-1)] lt limit)
if block then begin
repeat begin
xx=(footprint[cc-1])[0]
yy=(footprint[cc-1])[1]
allf=[[xx,yy+1],[xx-1,yy],[xx,yy-1],[xx+1,yy]]
cc--
if cc eq 0 then begin
return,-1
endif
blockf=$
(allf[0,0] eq -1 || allf[0,0] eq sz[1] || allf[1,0] eq -1 || allf[1,0] eq sz[2] || (where(allf[*,0] eq footprint))[0] ne -1 || arr[allf[0,0]>0<(sz[1]-1),allf[1,0]>0<(sz[2]-1)] lt limit) && $
(allf[0,1] eq -1 || allf[0,1] eq sz[1] || allf[1,1] eq -1 || allf[1,1] eq sz[2] || (where(allf[*,1] eq footprint))[0] ne -1 || arr[allf[0,1]>0<(sz[1]-1),allf[1,1]>0<(sz[2]-1)] lt limit) && $
(allf[0,2] eq -1 || allf[0,2] eq sz[1] || allf[1,2] eq -1 || allf[1,2] eq sz[2] || (where(allf[*,2] eq footprint))[0] ne -1 || arr[allf[0,2]>0<(sz[1]-1),allf[1,2]>0<(sz[2]-1)] lt limit) && $
(allf[0,3] eq -1 || allf[0,3] eq sz[1] || allf[1,3] eq -1 || allf[1,3] eq sz[2] || (where(allf[*,3] eq footprint))[0] ne -1 || arr[allf[0,3]>0<(sz[1]-1),allf[1,3]>0<(sz[2]-1)] lt limit)
endrep until blockf eq 0
for i=0,3 do begin
noutf=$  
(allf[0,i] eq -1 || allf[0,i] eq sz[1] || allf[1,i] eq -1 || allf[1,i] eq sz[2] || $
(where(allf[*,i] eq footprint))[0] ne -1 || arr[allf[0,i]>0<(sz[1]-1),allf[1,i]>0<(sz[2]-1)] lt limit)
if noutf eq 0 then return,allf[*,i]
endfor
endif else begin
for i=0,3 do begin
nout=$  
(all[0,i] eq -1 || all[0,i] eq sz[1] || all[1,i] eq -1 || all[1,i] eq sz[2] || $
(where(all[*,i] eq footprint))[0] ne -1 || arr[all[0,i]>0<(sz[1]-1),all[1,i]>0<(sz[2]-1)] lt limit)
if nout eq 0 then return,all[*,i]
endfor
endelse
end

pro auto_area,arr,pos,limit=limit,mask=mask
sz=size(arr)
mask=make_array(sz[1],sz[2],value=0)
footprint=list()
if ~keyword_set(limit) then limit=arr[pos[0],pos[1]]/2.0
temp=[pos[0],pos[1]]
count=0
while min(temp) ne -1 do begin
  footprint.add,temp
  mask[temp[0],temp[1]]=1
  temp=ccwalk(arr,temp,footprint,count,limit)
  count++
endwhile
return
end
