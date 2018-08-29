;pro transit, file, teles, inidate, dura

;;input parameters
;file='wasp_planet.list'; file with the inforamtion about planet,
;position and ephemeris. The file with planetary data has
;;name ra dec TO(days) Period(days) width(minuts)


;;teles the telescope used =0 for WHT 
                         ;; =1 for VLT
;;                         = 2 for Liverpool
;;;;;;;;;;;;;;;;;;;;;;;;;;=3 for SALT
;;;;;;;;;;;;;;;;;;;;;;;;;;=4 for NTT
;;;;;;;;;;;;;;;;;;;;;;;;;;=5 for greek
;;;;;;;;;;;;;;;;;;;;;;;;;;=6 for hawaii ;;;faulkes north
;;;;;;;;;;;;;;;;;;;;;;;;;;=7 for sliding spring observatory australia



;inidate=2454869.27229d0; start date to calculate the transits
;dura=100 period of time where the transits will be calculated


; can also run with the parameters set
file='wasp_planet.list'
;file='try'
teles=4

inidate=2455652.86252d0;april2011
dura=200

;;teles=4;;2



;;GET_JULDATE Get the current Julian date as a double precision scalar
;;JDCNV Convert from calendar date to Julian date.
;;JULDATE Convert from calendar date to reduced Julian date.
;CALDAT Return the calendar date and time given Julian date. (ITTVIS Library)
;JULDAY Calculate the Julian Day Number for a given month, day, and year. (ITTVIS Library)
;;;JD for the 16 january =2454848.30139



;;;;Start of the program.
;It calulates any type of transit but could bee adapated to get only
;full transits

 
ncols=10
np=100
d=dblarr(ncols-1,np)
name=strarr(np)
readtable,file , ncols,np, name, d, n

planet=name[0:n-1]

ra=reform(((d[2,0:n-1]/60.+d[1,0:n-1])/60.+d[0,0:n-1])*15.)
dec=reform(((d[5,0:n-1]/60.+d[4,0:n-1])/60.+d[3,0:n-1]))
to=reform(d[6,0:n-1]); in to of ephemris days
period=reform(d[7,0:n-1]);; period in days
wid=reform(d[8,0:n-1]/(24.*60.));; transform minuts to days

;;29°15.67S 70°43.88\u2032W
;;-29.261167  70.731333

;;observatorio
if teles eq 0 then begin
    geocord=[28.760722,-17.8816,2332]
endif else begin
    if teles eq 1 then begin
        geocord=[-24.6259,-70.4032,2648]
    endif else begin
        if teles eq 2 then begin  
            geocord=[28.76254,-17.879192,2344.]
        endif else begin
            if teles eq 3 then begin  
                geocord=[-32.3783, 20.8117,1770.]
            endif else begin
                if teles eq 4 then  begin 
                    geocord=[ -29.261167,  -70.731333  ,2347.]
                endif else begin
                   if teles eq 5 then  begin 
                      geocord=[ 37.98, 22.22  ,2340.]
                   endif else begin
                      if teles eq 6 then  begin 
                         geocord=[19.8250, -155.46,4215]
                      endif else begin
                         if teles eq 7 then  geocord=[-31.273333,149.06444 ,1165]
                      endelse
                    endelse
                endelse
            endelse
        endelse
    endelse
 endelse


 




lat=geocord[0]
lon=geocord[1]
altobs=geocord[2]



;For each planet get how many transit within the time


ncycles=fix((inidate-to)/period)

ntran=fix(dura/period)

;;estimate the number of transits
ntotal=total(ntran)
namelist=strarr(ntotal)
table=dblarr(7,3,ntotal)
maxalt=fltarr(ntotal)
minalt=fltarr(ntotal)
minhou=fltarr(ntotal)
maxhou=fltarr(ntotal)
midhou=fltarr(ntotal)
l=0
for j=0,n-1 do begin ;; for each object
    tntran=indgen(ntran[j])+1
    time=dblarr(3,Ntran[j])
    time[1,*]=to[j]+reform(period[j]*(ncycles[j]+tntran))
;print,format='(F30, F30, F30)',time[1,*]

    time[0,*]=time[1,*]-wid[j]/2.
    time[2,*]=time[1,*]+wid[j]/2.

    FOR i=0L,ntran[j]-1 DO BEGIN ; number of transits
      
        For k=0,2 do begin ;;ingress, midle, egress
            eq2hor, ra[j], dec[j],reform(time[k,i]), alt, az, ha, LAT=lat , LON=lon ,ALTITUDE=altobs, /double
            table[1,k,l]=alt
            table[2,k,l]=ha
 ;print, ra, dec, reform(time[k,i]), alt
        ENDFOR

        maxalt[l]=max(table[1,*,l])
        minalt[l]=min(table[1,*,l])
;;save the results to a array that has to include
;name time0 time1 time2 altitute0 altitude1 altitude2 hourangle0
;hourangle2

        namelist[l]=planet[j]
        table[0,0:2,l]=time[0:2,i]

        CALDAT,table[0,1,l] , Month, Day, Year, Hour, Minute, Second
        table[3,0:2,l]=[ Month, Day, Year]
        table[4,0:2,l]=[Hour, Minute, Second];;;mid transit time

        CALDAT,table[0,0,l] , Month, Day, Year, Hour, Minute, Second
        table[5,0:2,l]=[Hour, Minute, Second];;begining
        CALDAT,table[0,2,l] , Month, Day, Year, Hour, Minute, Second
        table[6,0:2,l]=[Hour, Minute, Second];;end

        minhou[l]=min(table[4:6,0,l]+table[4:6,1,l]/60.); miminum of the 3 times
        maxhou[l]=max(table[4:6,0,l]+table[4:6,1,l]/60.);; maximum of the 3 times
        midhou[l]=table[4,0,l]+table[4,1,l]/60.; mid transit time in hours
;print,format='(F30, F30, F30)',  


        l=l+1
;print, j, l, ntran[j]
    endfor


endfor



;;;the hour angle is the time that pass since the star crossed the
;;;local meridian, so since it reach the higesh position.
;the distance of a star to the local meredian 
;is somehow a distance to the solar noon which mean
;;;it can be used as a measure of the distance to the sun
;;; from the results it appears that from ha 0:180 its nightime and
;;;from 180 to 360 is the other side of the sun nighttime but i cannot
;;;be sure
;;only works for ojects above the ground
;; this is not interely accurate

 ;; and (min(hourang[*,i]) le 180) 




;;and( ( minhou le 6) or( maxhou ge 20.5)),number) for summer( ( minhou le 6.5) or( maxhou ge 19.5)),


;order=sort(table[0,1,*])
;junk=order
;for i=0, ntotal-1 do begin

;index=where( (maxalt gt 20) and (minha le 180) ,number)
if (teles eq 3) then begin 
    index=where( (maxalt gt 30) and ( ( minhou le 4) or( maxhou ge 17.)),number) 
endif else begin
 if teles eq 5 then begin
        index=where( (maxalt gt 15) and ( ( minhou lt 1.5) or ( maxhou ge 19.5)),number);;; careful with this is almost only full transits!
    endif else begin
    if teles eq 6 then begin
        index=where( (maxalt gt 25) and ( ( minhou gt 4) and ( maxhou le 16.5)),number);;; careful with this is almost only full transits!
    endif else begin
    if teles eq 7 then begin
        index=where( (maxalt gt 25) and ( ( minhou gt 9) and ( maxhou le 20.)),number);;; careful with this is almost only full transits!  
    endif else begin
        if (teles eq 4) then  begin
       ;;     index=where( (maxalt gt 25) and ( ( minhou gt 0.) and  ( maxhou le 9.0)),number) 
              index=where( (maxalt gt 25) and( ( minhou le 10.3) or ( maxhou ge 23.0)),number)
        endif else  index=where( (maxalt gt 25) and( ( minhou le 6.5) or( maxhou ge 19.5)),number)
    endelse
endelse
 endelse
endelse







;index=indgen(l)
;number=l

order=sort(table[0,1,index])
junk=index[order]


;; for i=0, number-1 do begin


;; print,format='(A10,I6,"/",I2,"/",I4,I5,":",I2,I6,":",I2,I4,":",I2,F10.2,F10.2,F10.2,F10.2,F10.2 )', namelist[junk[i]],table[3,1,junk[i]], table[3,0,junk[i]], table[3,2,junk[i]],table[4,0,junk[i]], table[4,1,junk[i]],table[5,0,junk[i]], table[5,1,junk[i]],table[6,0,junk[i]], table[6,1,junk[i]], table[1,0,junk[i]], table[1,1,junk[i]], table[1,2,junk[i]], table[2,0,junk[i]], table[2,2,junk[i]]

;; endfor


openw,unit,'transits_apr11ntt.txt',/get_lun
;openw,unit,'transits_ntt.txt',/get_lun
printf,unit,'#   Object      DATE        T0    Observing Period  Alt_beg    Alt_To   Alt_end ' 
for i=0, number-1 do begin
    printf, unit,format='(A10,I6.2,"/",I2.2,"/",I4,I5.2,":",I2.2,I6.2,":",I2.2,I4.2,":",I2.2,F10.2,F10.2,F10.2)', namelist[junk[i]],table[3,1,junk[i]], table[3,0,junk[i]], table[3,2,junk[i]],table[4,0,junk[i]], table[4,1,junk[i]],table[5,0,junk[i]], table[5,1,junk[i]],table[6,0,junk[i]], table[6,1,junk[i]], table[1,0,junk[i]], table[1,1,junk[i]], table[1,2,junk[i]]
endfor
free_lun, unit



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;full transits only work for the NTT!!!!!!!!!!!!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop

;;;;full transit
;;;not working
;;;;:????????????????
;;;index=where( (minalt gt 20) and( ( midhou le 10. ) and ( midhou ge 23.5)),number)
index=where( (minalt gt 20) and (maxhou le 10.3) and (minhou le 10.0),number)


order=sort(table[0,1,index])
junk=index[order]
openw,unit,'transitsfull_apr11ntt.txt',/get_lun
printf,unit,'#   Object      DATE        T0    Observing Period  Alt_beg    Alt_To   Alt_end ' 
for i=0, number-1 do begin
    printf, unit,format='(A10,I6.2,"/",I2.2,"/",I4,I5.2,":",I2.2,I6.2,":",I2.2,I4.2,":",I2.2,F10.2,F10.2,F10.2)', namelist[junk[i]],table[3,1,junk[i]], table[3,0,junk[i]], table[3,2,junk[i]],table[4,0,junk[i]], table[4,1,junk[i]],table[5,0,junk[i]], table[5,1,junk[i]],table[6,0,junk[i]], table[6,1,junk[i]], table[1,0,junk[i]], table[1,1,junk[i]], table[1,2,junk[i]]
endfor
free_lun, unit

stop

;;; also write the hour angle

openw,unit,'transitshour.txt',/get_lun
printf,unit,'#   Object      DATE        T0    Observing Period  Alt_beg    Alt_To   Alt_end      Hour Angle' 
for i=0, number-1 do begin
    printf, unit,format='(A10,I6.2,"/",I2.2,"/",I4,I5.2,":",I2.2,I6.2,":",I2.2,I4.2,":",I2.2,F10.2,F10.2,F10.2,F10.2,F10.2 )', namelist[junk[i]],table[3,1,junk[i]], table[3,0,junk[i]], table[3,2,junk[i]],table[4,0,junk[i]], table[4,1,junk[i]],table[5,0,junk[i]], table[5,1,junk[i]],table[6,0,junk[i]], table[6,1,junk[i]], table[1,0,junk[i]], table[1,1,junk[i]], table[1,2,junk[i]], table[2,0,junk[i]], table[2,2,junk[i]]
endfor
free_lun, unit




;; times in JD

 ;DAYCNV, time[1,0], YR, MN, DAY, HR
;print,  YR, MN, DAY, HR


;;; its not working for south objcts also....not exacty the epemeris
;;; of the sun

end
