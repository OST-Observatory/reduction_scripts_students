;
; December 2012
; ------------
; stack add an image series to one image and makes a dark and flatfield correction.
; stack use the program add_images written by Nadine Giese and was modified by Marcel Pietschmann
; for use within the astrophysics lab course at Potsdam University.
; The source code may be modified, reused, and distributed as long as it 
; retains a reference to the original author(s).
;
; May 2014: bug in the cross corelation eliminated which caused too small
;           shifts between consecutive images (rh)
;
; For problems and questions: nadine.giese@gmx.net or mpietsch@astro.physik.uni-potsdam.de
;

pro stack

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; configuration - please only edit here :-)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Images should be sorted into directories, folders with images should
; not contain further subfolders
; * Images of the object
; * Darkframes for the object images and flatfield
; * flatfields


; darks

; if no darks are used for the flatfields, the paths can be empty strings
; or whitespaces, e.g. darksFlat = ' '

; path to darkframes for flats
darksFlat = 'dark-flats/'
; path to darkframes for pictures
darksImg = 'darks/60s'


; flats

; path to flats
flats = 'flats'


; images

; path to visual images
img = 'images'



; correlation stuff - please only edit this section if problems occure
;                     the following values represent the assumed maximum 
;                     displacement between two images in pixel
;                   - default values reduced to 60 pixel to avoid false 
;                     idetification of artefacts from the cross crorelation 
;                     which lead to wrong shifts between the images (rh)
;                   - for very long image series it might is more appropriated 
;                     to split the time series rather than increasing the 
;                     assumed maximum shift, since the latter often lead to
;                     the idenfication of cross-corelation artefacts
xMax = 60
yMax = 60


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; main program
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fileExt=['.fit','.FIT']
; create output directory if not present
if(file_test('output/') eq 0) then spawn,'mkdir output'

; find out if flatdarks are used
checkflatdark = max(strlen(file_search(darksflat+'*'+fileExt)))
print, checkflatdark
if checkflatdark lt 1 then message,/info,'WARNING: no flat darks given'

; find out if darks are used
checkdark = max(strlen(file_search(darksImg+'*'+fileExt)))
if checkdark lt 1 then message,/info,'WARNING: no darks given'

; find out if flats are used
checkflats = max(strlen(file_search(flats+'*'+fileExt)))
if checkflats lt 1 then message,/info,'WARNING: no flats given'

; find out if images are used
checkimg = max(strlen(file_search(img+'*'+fileExt)))
if checkimg lt 1 then begin
  message,/info,'WARNING: no images given -> EXIT'
  stop
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read darkframes and calculate mean
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; masterdark
if checkdark ne 0 then begin
    fileList = File_Search(darksImg,'*')  

    ; size of images
    imgSize = size(readfits(fileList[0],hnull,/SILENT))
    nx = imgSize[1]
    ny = imgSize[2]
    headerSize = N_elements(hnull)
    
    masterdark = lonarr(nx,ny)
    for i=0,len(fileList)-1 do begin
        masterdark = masterdark + long(readfits(fileList[i],/SILENT))
    endfor
    masterdark = masterdark/len(fileList)
    writefits, 'output/masterdark.fit', masterdark
endif

; darks for flatfield
if checkflatdark ne 0 then begin
    fileList = File_Search(darksFlat,'*')  

    ; size of images
    imgSize = size(readfits(fileList[0],hnull,/SILENT))
    nx = imgSize[1]
    ny = imgSize[2]
    headerSize = N_elements(hnull)
    
    ; flatdark mean 
    flatdark = lonarr(nx,ny)
    for i=0,len(fileList)-1 do begin
        flatdark = flatdark + long(readfits(fileList[i],/SILENT))
    endfor
    flatdark = flatdark/len(fileList)
    writefits, 'output/flatdark.fit', flatdark
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; process flatfields 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if checkflats ne 0 then begin
    fileList = File_Search(flats,'*')
    if checkflatdark lt 1 then begin
        ; size of images
        imgSize = size(readfits(fileList[0],hnull,/SILENT))
        nx = imgSize[1]
        ny = imgSize[2]
        headerSize = N_elements(hnull)
    endif

    flat = lonarr(nx,ny)

    for i=0,len(fileList)-1 do begin
        flat = flat + long(readfits(fileList[i],/SILENT))
    end
    flat = flat/len(fileList)
    if checkflats ne 0 then flat = (flat-flatdark) > 0
    flat = flat/mean(flat)
    writefits, 'output/flat.fit', flat
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; process image files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; read images and estimate shifts via cross correlation
fileList = File_Search(img,'*')

if checkflatdark lt 1 or checkdark lt 1 then begin
  ; size of images
  imgSize = size(readfits(fileList[0],hnull,/SILENT))
  nx = imgSize[1]
  ny = imgSize[2]
  headerSize = N_elements(hnull)
endif

addImage = lonarr(nx,ny,len(fileList))
visual = lonarr(nx,ny,len(fileList))
addImagecorr = intarr(2,len(fileList))
header = strarr(headerSize,len(fileList))
print, STRING(13B)
print, 'Displacement'
print, '   Image', '       x', '       y','   Path'
print, '------------------------', '   --------'
for i=0,len(fileList)-1 do begin
  addImage[*,*,i] = long(readfits(fileList[i],h,/SILENT))
  if max(addImage[*,*,i]) eq 65535 then begin
     message,/info,'WARNING: Image '+strtrim(i,2)+' is overexposed!'
  endif

  if (checkflats lt 1) and (checkdark lt 1) then begin
;     visual[*,*,i] = long(readfits(fileList[i],h,/SILENT)) > 0.
    visual[*,*,i] = addImage[*,*,i]
  endif else if checkflats lt 1 then begin
;     visual[*,*,i] = long(readfits(fileList[i],h,/SILENT)-masterdark) > 0.
    visual[*,*,i] = addImage[*,*,i] - masterdark > 0.
  endif else if checkdark lt 1 then begin
;     visual[*,*,i] = long(readfits(fileList[i],h,/SILENT))/flat > 0.
    visual[*,*,i] = addImage[*,*,i] / flat > 0.
  endif else begin
;     visual[*,*,i] = long(readfits(fileList[i],h,/SILENT)-masterdark)/flat > 0.
    visual[*,*,i] = (addImage[*,*,i] - masterdark) / flat > 0.
  endelse
  addImagecorr[*,i] = correl_images(visual[*,*,0],visual[*,*,i],xMax,yMax)
  header[*,i] = h 
  if i eq 0 then begin
    print,i, '       0', '       0','   ',fileList[i]
  endif else  print, i, addImagecorr[0,i], addImagecorr[1,i],'   ',fileList[i]
end

; maximum and minimum shifts
addImagecorr[0,0] = 0
addImagecorr[1,0] = 0

minShiftx = min(addImagecorr[0,*])
maxShiftx = max(addImagecorr[0,*])

minShifty = min(addImagecorr[1,*])
maxShifty = max(addImagecorr[1,*])

deltaShiftx = maxShiftx - minShiftx
deltaShifty = maxShifty - minShifty

; define new arrays with cut size
addImagecut = lonarr(nx-deltaShiftx,ny-deltaShifty,len(fileList))
addImagePic = lonarr(nx-deltaShiftx,ny-deltaShifty)

; cut pictures and add them
for i=0,len(fileList)-1 do begin
    xs = 0 + addImagecorr[0,i] - minShiftx
    xe = nx- 1 -(deltaShiftx-(addImagecorr[0,i]-minShiftx))
    ys = 0 + addImagecorr[1,i] - minShifty
    ye = ny- 1 -(deltaShifty-(addImagecorr[1,i]-minShifty))
    addImagecut[*,*,i] = addImage[ xs:xe , ys:ye , i]
    addImagePic = addImagePic + addImagecut[*,*,i]
end

; free memory
addImage = 0
addImagecut = 0
addImagecorr = 0

; free other stuff
fileList = 0

writefits, 'output/image.fit', addImagePic

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; function definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;; cross correlation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function correl_images, imagea, imageb, xMax, yMax

; Idea and further information:
; http://en.wikipedia.org/wiki/Phase_correlation

lx = (size(imagea))[1]
ly = (size(imagea))[2]

imagea_cut = float(imagea);
imageb_cut = float(imageb);

imafft = fft(imagea_cut)
imbfft = fft(imageb_cut)
imbfftcc = conj(complex(imbfft))
fftcc = imafft*imbfftcc
fftcc = fftcc/abs(fftcc)
cc = fft(fftcc,-1)
cc[0,0] = 0.

; commented out since it seems to be unnecessary and it definitely is the
; reason why some shifts might not be recognized (rh)
; cc = abs(cc)

for i = xMax,lx-xMax-1 do begin
  for j = 0,ly-1 do begin
    cc[i,j] = 0
  end
end
for i = 0,lx-1 do begin
  for j = yMax,ly-yMax-1 do begin
    cc[i,j] = 0
  end
end


maxCorr = max(cc,location)
ind = array_indices(cc,location)

if ind[0] gt lx/2. then ind[0] = (ind[0]-1)-lx+1 else ind[0] = ind[0] - 1
if ind[1] gt ly/2. then ind[1] = (ind[1]-1)-ly+1 else ind[1] = ind[1] - 1

return, ind

; free memory
imagea_cut = 0
imageb_cut = 0
imbfftcc = 0
fftcc = 0
cc = 0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;; list length ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function len, list

listsize = size(list)
n = listsize[1]
return, n

end

