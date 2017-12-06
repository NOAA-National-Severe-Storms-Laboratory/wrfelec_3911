      program PLOT_INTERP_LIGHTNING_DATA

!      Developer: Alexandre Fierro (Alex.Fierro@noaa.gov) - December 2014.
!
!      Code to interpolate lat/lon point flash data onto a given
!      domain of a Mercator/equidistant lat/lon grid for
!      WRF [map_proj='mercator' in WPS]. This code also features an option to inflate the lightning
!      densities onto a coarser footprint on the same WRF grid.
!
!      This software must utilize the same environment variables for NETCDF than
!      the version of the WRF build. I.e., ${NETCDF} must be the path of
!      netcdf build on the local machine/computer cluster.
!
!      compile code as follows: 
!
!      NETCDF 3.6.X: ifort -CB -traceback -132 LIGHT_INTERP_WRF_GRID_2.f -L${NETCDF}/lib -I${NETCDF}/include -lnetcdf -o LIGHT_INTERP_WRF_GRID_2.exe
!      NETCDF 4: ifort -CB -traceback -132 LIGHT_INTERP_WRF_GRID_2.f -L${NETCDF}/lib -I${NETCDF}/include -lnetcdf -lnetcdff -o LIGHT_INTERP_WRF_GRID_2.exe
!
!     Important: The first record of the raw lightning data must equate the starting time of the simulation (on the grid: e.g., d01, d02) 
!                minus one lightning accumulation interval. 
!     (e.g., if model is started at 12Z on d01 using lightning density rates of 10 min (600 s), then the raw lightning data records must start at 1150Z).
!
!     the format of the lightning data is assumed to be of the following
!     format: day,hour,min,sec, lat, lon
!     with time = 0,0,0,0 being the starting time of the simulation domain where the DA is taking place.
!
!     The code below is set up for  the parent domain d01.
!     To perform this operation on inner nests (e.g., d02, d03 etc), the lat/lon
!     data on the WRF grid (XLAT(:,:) and XLONG(:,:) arrays 
!     are obtained by reading the first wrfout_d0[2-3] output file instead of wrfinput_d01. 
!     These should be available in most cases as
!     control run without lightning assimilation often is performed first.
!
!     For simplicity, it is assumed that wrfinput_d01 or wrfout_d0[2-3]*
!     and the file containing the raw lightning 'point flash' data
!     are located in the same (current) directory 
!
!     This code won't work for an inner moving nest as lat/lon of
!     the observations must be domain relative. Code can be extended by
!     including a loop on wrfout_d0[2-3] files to update xlon and xlat as the domain
!     moves before interpolating the point flash data.  
!
!     Reminder; Important: This code only supports Mercator/Equidistant projections.

      use netcdf

      implicit none

      integer its,i,j,k,t,memberid,it,npts,npts2
      integer a,n,b
      INTEGER tmp,len,system  

!     %%%%%%%% START input parameters %%%%%%%      

      real, parameter :: start_grid=6.*3600. ! start time (seconds) of simulation on grid (e.g., d01, d02) relative to time of first record in data
      real, parameter :: timeint=1800. ! time interval in sec for accumulation (must equate light_nudge_int namelist parameter)
      real, parameter :: stime=0. ! start time in sec of lightning assimilation (must equate light_nudge_times namelist parameter). Relative of model time on domain
      real, parameter :: etime=7200. ! end time in sec of lightning assimilation (must equate light_nudge_timee namelist parameter). Relative of model time on domain
      
! Example: if domain of interest is spawned at 06Z and the lightning DA is applied between 07-08UTC using a 30 min DA interval then: 
!      start_grid=6.*3600
!      timeint=1800.
!      stime=3600.
!      etime=7200.
!      the raw lightning data must start at the *earliest* with:
!      day,hour,min,sec, lat,         lon
!       0,  05, 30, 00, lat value, lon value
!      accumulation intervals are constructed as follows:
!              530Z   6Z   630Z   7Z   730Z   8Z   
!               [.....[.....[.....[.....[.....[
!                   \      \     \     \     \
! intervals(i)=       1     2     3     4     5       

      character path*150 /''/ ! location of raw lightning data file and wrfinput_d01 for d01 or wrfout_d02 file for inner nest
      character lightdata_fname*150 /'WWLLN_clean.txt'/ ! name of raw lightning data file (fortran formatted ASCII). Data assumed to have the following format for each record: day,timeh,timem,times,lat,lon 
      real, parameter :: buff = 20.  ! Halo zone (in km) near boundary where the gridded lightning densities are forced to zero to avoid potential CFL problems  

      logical :: Ichoose_9km=.false. ! flag to turn on inflation of the interpolated data onto a coarser footprint on the same grid (here assumed to be 9-km) 
      real, parameter :: dxkminterp=9.0 ! optional:  target interpolation resolution of lightning footprints on the WRF grid -assume here as 9-km=Anticipated pseudo-GLM resolution
      character grid_id*3    /'d01'/ ! suffix of lightning files assumed to be light.d0'grid_id'.out.'timeint'.txt' - which are assumed in the WRF code-do not change except for grid_id = here d01. If first inner nest is chosen then change d01 to d02 etc.
!     %%%%%%%% END input parameters %%%%%%%      
   
      integer :: maxtime ! number of time intervals in assimilation window = number of gridded, time-parsed lightning files
      integer :: starttime, endtime
      
      integer :: nx, ny, dxkm,nor ! domain size and horizontal resolution - read in wrfinput_d01 or wrfout_d02 file for inner nest
      real, allocatable ::  xlon(:,:),xlat(:,:) ! lat/lon coordinates of WRF grid read in wrfinput_d01 or wrfout_d02 file for inner nest 
           
      real,allocatable ::  intervals(:),lat(:),lon(:),time(:)
      real :: timeh,timem,times,day


      real, allocatable ::  grid_light(:,:,:) ! gridded, time-parsed lightning data array

      character command*160
      character memberid2*10
      character fname*150,fname2*150

! NETCDF stuff
      integer ncid,variddat
      integer cmode, status
      integer nxID,nyID
      
! GRADS stuff      
      integer irec10,irec20,irec30,irec40
      real, allocatable :: gga0(:,:,:,:,:), flash(:,:,:), light9km(:,:,:,:)
      integer, parameter :: nfld=1

       maxtime=   1 + NINT((etime-stime)/timeint)
       starttime= 1 + NINT((stime)/timeint) 
       endtime=   1 + NINT((etime)/timeint)
   
! READ
       fname=trim(path)//trim(lightdata_fname) ! raw data file - use 'sort -u old_file > lightdata_fname' to remove duplicate lines if applicable.
     
!       print*,fname
        command='wc -l < '//trim(lightdata_fname)//'> /tmp/foo.txt'

        len=system(command)
        open(unit=44,file='/tmp/foo.txt',status='unknown',form='formatted')
! '
        read(44,*) len
        print*, 'Number of point flash records=',len
        tmp=system('rm -rf /tmp/foo.txt')
        close(44)
      
       allocate(intervals(len))
       allocate(lat(len))
       allocate(lon(len))
       allocate(time(len))
       intervals(:)=0 ; time(:)=0. ; lat(:)=0. ; lon(:)=0. 

       open(unit=44,file=fname,status='unknown',form='formatted')
 
       do i=1,len
       read(44,*) day,timeh,timem,times,lat(i),lon(i)
       time(i)=INT(day*86400+timeh*3600+timem*60+times)
       enddo

      close(44)
      
!     Get grid specs from wrfinput_d01 or wrfout_d02,3 etc file (no striding):

      fname=trim(path)//'wrfinput_d01' ! for D01 
!      fname=trim(path)//'wrfout_d02_2012-08-28_12:00:00'  ! for D02 etc
!      print*,'filename is: ',fname
      cmode = NF90_NOWRITE
      status = nf90_open(fname,cmode,ncid)
      if (status /= nf90_noerr) call handle_err(status)

! get dimension IDs
      status = nf90_inq_dimid(ncid, "west_east", nxID)
      status = nf90_inq_dimid(ncid, "south_north", nyID)
!      status = nf90_inq_dimid(ncid, "bottom_top",  nzID)
! get dimension values
      status = nf90_inquire_dimension(ncid, nxID, len=nx)
      status = nf90_inquire_dimension(ncid, nyID, len=ny)
     
!     add +1 because of staggering      
      print*, 'nx, ny dims of interp domain with staggering: ',nx+1,ny+1
      
      allocate( xlat(nx,ny) )
      allocate( xlon(nx,ny) )

!     add +1 because of staggering      
      allocate(grid_light(maxtime,nx+1,ny+1))  
      grid_light(:,:,:)=0


      status = NF90_GET_ATT(ncid, NF90_GLOBAL, 'DX',dxkm)
      if (status /= nf90_noerr) call handle_err(status)

      dxkm=1.e-3*dxkm ! assume dx=dy

! get xlat  ID
      status = NF90_INQ_VARID(ncid,'XLAT',variddat)
      if (status /= nf90_noerr) call handle_err(status)
! get xlat array
      status = NF90_GET_VAR(ncid,variddat,xlat)
      if (status /= nf90_noerr) call handle_err(status)
! get xlon  ID
      status = NF90_INQ_VARID(ncid,'XLONG',variddat)
      if (status /= nf90_noerr) call handle_err(status)
! get xlon array
      status = NF90_GET_VAR(ncid,variddat,xlon)
      if (status /= nf90_noerr) call handle_err(status)
      
! open empty light.grid_id.out files and assign a unique unit number to each before interpolation

       DO memberid=1,maxtime

!              if (memberid.lt.10) then
!                write (memberid2,'(i1.1)') memberid
!              else if (memberid.lt.100) then
!                write (memberid2,'(i2.2)') memberid
!              else  if (memberid.lt.1000) then
!                write (memberid2,'(i3.3)') memberid
!              else
                write (memberid2,'(i5.5)') memberid
!              endif

!       open empty light files for interp

       fname=trim(path)//'light.'//trim(grid_id)//'.out.'//trim(memberid2)//'.txt'
       print*,fname
       open(unit=memberid,file=fname,status='unknown',form='formatted')
! '
       ENDDO

       nor =INT(buff/dxkm) ! buffer zone around boundary where grid_light=0


       DO i=1,len 

!         fits in time (in seconds) into a time interval - must start at 1
 
           intervals(i)=( int(time(i)/timeint)+1 ) - ( int((start_grid-timeint)/timeint))

!          write(0,*) i,intervals(i),starttime,endtime

!!!!!!!!!!!  Interpolation/projection of point flash data onto WRF lat/lon Equidistant grid !!!!!!!!!!!!!!!

           IF(intervals(i).ge.starttime.and.intervals(i).le.endtime) THEN ! are we in the DA interval?


!     Brute force interpolation (search) (slow).
!           outer : do a=nor,nx-nor ! assume gridded lightning to be zero with buff km from  boundaries
!           do b=nor,ny-nor ! assume gridded lightning to be zero within buff km from  boundaries
!           if (intervals(i)-(starttime).eq.1) print*,xlat(1,b),xlon(a,1)
!            if(  (lat(i) .ge. xlat(1,b)) .and. (lat(i) .lt. xlat(1,b+1))  ) then
!            if(  (360+lon(i) .ge. 360+xlon(a,1)) .and. (360+lon(i) .lt. 360+xlon(a+1,1))  ) then ! LON are negative !
!             grid_light(intervals(i)-(starttime-1),a,b)=grid_light(intervals(i)-(starttime-1),a,b)+1. 
!             exit outer
!            endif
!            endif
!           enddo
!           enddo outer
!
!        if (MAXVAL(grid_light(intervals(i)-(starttime),:,:)).gt.0) print*,'MAXVAL',MAXVAL(grid_light(intervals(i)-(starttime),:,:))


!     Interpolation using Fortran's intrinsic function (fast)


           if ( (lat(i).ge.xlat(1,1)) .and. (lat(i).le.xlat(1,ny)) ) then
           if ( (lon(i).ge.xlon(1,1)) .and. (lon(i).le.xlon(nx,1)) ) then
             
             b=minloc((abs(lat(i)-xlat(1,:))),1)
             a=minloc((abs(lon(i)-xlon(:,1))),1)


! Implies lightning densities 2pts near domain boundary to zero to avoid promoting W there:
             IF ((a.gt.nor+1) .and. (b.gt.nor+1)) THEN
             IF ((a.lt.nx-nor) .and. (b.lt.ny-nor)) THEN

!            accounts for obs points in [b:b+1] to stay there, except:

             if (b.lt.ny) then
             if (lat(i).gt.0.5*(xlat(1,b)+xlat(1,b+1))) b=b-1
             elseif (b.eq.ny) then
             if (lat(i).gt.0.5*(xlat(1,b-1)+xlat(1,b))) b=b-1              
             endif
        
             if (a.lt.nx) then
             if (lon(i).gt.0.5*(xlon(a,1)+xlon(a+1,1))) a=a-1
             elseif (a.eq.nx) then
             if (lon(i).gt.0.5*(xlon(a-1,1)+xlon(a,1))) a=a-1
             endif

!       AVOID MINS that are large --> points outside domain will be
!       assigned border indices
!
             grid_light(intervals(i)-(starttime-1),a,b)=grid_light(intervals(i)-(starttime-1),a,b)+1 
!'
             ENDIF
             ENDIF
!'
            endif
            endif


           ENDIF
   

       ENDDO
      
      DEALLOCATE(lat,lon,time,intervals)
 
!     This loop inflates the interp data onto a grid of horizontal grid spacing given by dxkminterp - assumed to be pseudo-GLM=9 km:

      npts=ceiling(dxkminterp/dxkm)

      npts2=ceiling(npts*0.5)

      allocate(light9km(maxtime,nx+1,ny+1,1))
      light9km(:,:,:,:)=0.

      nor=INT(buff/dxkminterp)

      IF(Ichoose_9km) THEN
 
      DO t=1,maxtime
 
      do j=1, ny+1 
      do i=1, nx+1 

        if (i.gt.nor+npts2.and.i.lt.(nx+1)-npts2-nor) then
         if (j.gt.nor+npts2.and.j.lt.(ny+1)-npts2-nor) then 
            if (grid_light(t,i,j).ge.1) then

              do a=1,npts ! goes from i-8 to i+9 aka 8+9+1 = 18 points
                do b=1,npts
              
                 light9km(t,(i-npts2)+a,(j-npts2)+b,1)=light9km(t,(i-npts2)+a,(j-npts2)+b,1)+grid_light(t,i,j)+1

                enddo
              enddo

            endif
          endif
        endif 

      enddo
      enddo
      
      ENDDO

      ENDIF

!      write parsed lightning data for assimilation       
       DO i=1,maxtime 
       do a=1,nx+1
       do b=1,ny+1
        IF(Ichoose_9km) THEN
        write (i,*) light9km(i,a,b,1) 
        ELSE
        write (i,*) grid_light(i,a,b)
        ENDIF
       enddo
       enddo
       write (0,*) maxval(grid_light(i,:,:)) 
       ENDDO
       

       deallocate(grid_light,light9km)


!      close all light files

       DO memberid=1,maxtime
       close(memberid)
       ENDDO
  
!      Grads plotting stuff        
       open(unit=33,file='LIGHT.ctl',form='formatted')

       write(33,*)'DSET LIGHT.dat'
       write(33,*)'TITLE FLASH'
       write(33,*)'UNDEF -9.99E33'
       write(33,*) 'XDEF', nx, 'LEVELS'
       DO a=1,nx
       write(33,*) xlon(a,1)  
       ENDDO
       write(33,*) 'YDEF', ny, 'LEVELS'
       DO b=1,ny
       write(33,*) xlat(1,b)  
       ENDDO
       write(33,*) 'ZDEF    1   LEVELS  1000'
       write(33,*)   'TDEF',maxtime,'LINEAR 00:00Z01JAN2015 10mn'
       write(33,*)   'VARS 1'
       write(33,*)'Z1 1   1001    FLASH'
       write(33,*)'ENDVARS'

       close(33)


      deallocate(xlat,xlon)

      irec10=1
      irec20=0
      irec30=0
      irec40=0

       DO it=1,maxtime
!              if (it.lt.10) then
!                write (memberid2,'(i1.1)') it
!              else if (it.lt.100) then
!                write (memberid2,'(i2.2)') it
!              else  if (it.lt.1000) then
!                write (memberid2,'(i3.3)') it
!              else
                write (memberid2,'(i5.5)') it
!              endif

       fname=trim(path)//'light.'//trim(grid_id)//'.out.'//trim(memberid2)//'.txt'

       open(unit=150,file=fname,status='unknown',form='formatted')

       print*, 'PLOTTING', fname 

       allocate(flash(nx+1,ny+1,1))
       allocate(gga0(1,1,1,ny,nx)) ! to match xlon and xlat dims - not on staggered grid !

       do i=1,nx+1
       do j=1,ny+1
       read(150,*) flash(i,j,1)
       enddo
       enddo
       print*,'MAX',maxval(flash(:,:,1)),it
       close(150)

       print*,'read data ok, time=',it 
     
      do k=1,1
      do j=1,ny
      do i=1,nx
      gga0(1,1,k,j,i)=flash(i,j,k)
      enddo
      enddo
      enddo

      call grads_write(gga0,nx,ny,1,nfld,irec10,irec20,irec30,irec40)

      deallocate (flash,gga0)

      print*, 'GRADS DONE. it=',it

      ENDDO

       stop
       END

! ################################################
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
! ################################################
      subroutine handle_err(errcode)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
    
      if(errcode /= nf90_noerr) then
        print *, 'Error: ', trim(nf90_strerror(errcode))
       stop "Stopped"
      endif
      end subroutine handle_err

!   ==========SUBROUTINE for GRADS===========

      subroutine grads_write(a,nx,ny,nz,nfld,irec1,irec2,irec3,irec4)
      integer nx,ny,nfld,nz
      parameter (ntimes=1) 
      real*4 a(ntimes,nfld,nz,ny,nx)

      open(29,file='LIGHT.dat',form='binary',access='direct',recl=4*nx*ny)
      do 10 lt=1,ntimes
       do la=1,nfld
      do i=1,nz
      itot=irec1+irec2+irec3+irec4
      write(29,rec=itot)((a(lt,la,i,j,k),k=1,nx),j=1,ny)
       irec4=irec4+1
      enddo
       irec4=irec4-1
      irec3=irec3+1
      enddo
      irec3=irec3-1
      irec2=irec2+1
 10   continue
!     irec2=irec2-1
      return
      end subroutine grads_write
