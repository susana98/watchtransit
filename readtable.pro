;;; reads the tables with the bvrimk colors for brown swarfs, ijhk tables for dwarfs and giants

;;;# bvrimk.txt has the calibration calibration of MK spectral types
;;;# taken allens astrophysics quatities

;;	ijhk_dwarfs.txt','ijhk_giants.txt were taken from 
;;;	JHKLM photometry - Standard systems, passbands, and intrinsic colors
;;Bessell, M. S.; Brett, J. M. 1988

;;;this is just like any open.pro but  the first colum is a string with the name of the spectral type  and is separatly writen to name and the number table is writen to d

PRO readtable, file,ncols,nrows, name, d, nd


;------- Default directory with data ---
	afile=file
;------- Opening the file ---
        openr,unit,afile ,/get_lun ;
  	print,'  Reading file: ',afile
;------- Skcping the comments ---
	cskpcom='#'
	point_lun,unit,0
	while(cskpcom eq '#') do begin
		tmpskpcom=fstat(unit)
		readf,unit,format='(a1)',cskpcom
	endwhile
	point_lun,unit,tmpskpcom.cur_ptr
;-------
	d=DBLARR(ncols-1,nrows)
        name=strarr(nrows)
	nd=0L
        xr=''
	while (not EOF(unit)) do begin             
            readf,unit, xr
            data=strsplit(xr,/extract)
            name[nd]=data[0]
            d[*,nd]=double(data[1:ncols-1])
           nd=nd+1 
        ENDWHILE


        free_lun, unit
	print,nd,' points!'
    END
