program Track_density

!
! This code is a copy of /home/chen/codes/cir_disttr/Track_density_all_seasons.f90

! Compile with:
! module load compiler/gcc-7.3 development/netcdf
! gfortran $(nf-config --fflags) $(nf-config --flibs) XXX.f90 -o XXX.Abs
!
! This program calculated the number of tracks passing a grid cell (with repeated entried of the same track
! being counted as one) (Neu et al., 2013)

use netcdf
implicit none

character(len=110) fileout
integer ix, jx
real, allocatable, dimension (:) :: vlon, vlat
real, allocatable, dimension (:,:) :: trackDen1, trackDen2, trackDen3, trackDen4
real, allocatable, dimension (:,:) :: lat, lon
integer, allocatable, dimension (:,:) :: yy,mm,date
integer, allocatable, dimension (:) :: pmax
logical, allocatable, dimension (:,:) :: Trackpass1, Trackpass2, Trackpass3, Trackpass4
real R,dist,timediff, mindis, total_p_InDom, Tlat, Tlon, Dlat, Dlon
integer i, j, p, t, num, point
integer nummax, px, topnum, num_previous, tp, seg_point
integer line, linemax, lineheader, vix, vjx, vi, vj, II, JJ, sgd

!parameter (linemax = 257119, lineheader = 1, nummax = 24602, px=556)
parameter (linemax = 1833497, lineheader = 1, nummax = 24604, px=556)

!nummax:total number of tracked storms; px:max points per storm
!topnum : Not sure what it is, but it doesn't seem to be used in the code

! Output data
  R=2.5e5 !in meters
  fileout='/pampa/cloutier/Track_density_1979_2020_all_seasons_250km.nc'

! Define a virtual domain covers the entire domain of the filtered ETCs

! Lon: 0 to 359.5 (vix=359.5/0.5+1)
! Lat: -20 to 90  (vjx=110/0.5+1)

  vix=720
  vjx=221
allocate (vlat(vjx))
allocate (vlon(vix))
allocate (trackDen1(vix, vjx))
allocate (trackDen2(vix, vjx))
allocate (trackDen3(vix, vjx))
allocate (trackDen4(vix, vjx))
allocate (Trackpass1(vix, vjx))
allocate (Trackpass2(vix, vjx))
allocate (Trackpass3(vix, vjx))
allocate (Trackpass4(vix, vjx))

! The virtual domain covers all storm tracks:
! with 0.5 degree interval
do i= 1, vix
do j= 1, vjx
   trackDen1(i,j) = 0.
   trackDen2(i,j) = 0.
   trackDen3(i,j) = 0.
   trackDen4(i,j) = 0.
   vlon(i) = 0.+float(i-1)*0.5
   vlat(j) = -20.+float(j-1)*0.5
enddo
enddo
!---

! Input data
allocate (yy(nummax,px))
allocate (mm(nummax,px))
allocate (date(nummax,px))
allocate (lat(nummax,px))
allocate (lon(nummax,px))
allocate (pmax(nummax))

pmax = 1

! input data
open (10, FILE='/pampa/cloutier/ERA5_all_stormtracks.txt', STATUS='old')

print*,'Start reading the .txt file...'

do line = 1, linemax

  if (line .le. lineheader) then
    print*,'Start reading the .txt file...'
    read(10,*)    !header
  else
    read(10,100)  num, point, yy(num,point), mm(num,point),date(num,point), lat(num, point), lon(num, point)
    pmax(num)=point
  endif
enddo

!print *,'num =', num, 'point =', point, 'yy =', yy(num,point), 'mm =', mm(num,point)

!print *, 'date =', date(num,point), 'lat =', lat(num,point), 'lon =', lon(num,point)

close(10)
!nummax=num

print*,'Finish reading the .txt file...'

!=============================================================

! Now calculating the track density (the number of tracks, swept by a given radius R,
! per grid point in all time)

do num = 1 , nummax
! 
   if (yy(num,point).ne.0) then
      print*, yy(num,point)
   endif 
  ! print *, 'yy = ', yy(num,point)
   do vi = 1, vix
   do vj = 1, vjx
      Trackpass1(vi,vj) = .False.
      Trackpass2(vi,vj) = .False.
      Trackpass3(vi,vj) = .False.
      Trackpass4(vi,vj) = .False.
   enddo
   enddo
   !
   IF (pmax(num).gt.1) THEN
   do point = 1, pmax(num)

      Tlat = lat(num,point)
      Tlon = lon(num,point)
      ! Find the cloest (II) and (JJ) for the storm center (Tlat, Tlon):
      ! INT(A): If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
      ! then INT(A) is the integer whose magnitude is the largest integer that does 
      ! not exceed the magnitude of A and whose sign is the same as the sign of A.
         JJ = int((Tlat+20)/0.5)+1 
         II = int(Tlon/0.5)+1
         if ( mod(Tlat,0.5).gt.0.25 ) JJ = JJ+1
         if ( mod(Tlon,0.5).gt.0.25 ) II = II+1  !find the "CLOEST" II and JJ not the "smallest" integer
      !!   print*,'point:',point,' II,JJ=',II,JJ

      ! sgd: searching distance in virtual grid points (must be larger than R)
        sgd= int(R/1.1e5/0.5)+20  ! 20 is the buffer as 1.1e5 is just an approx. 
                                  ! to convert km to lat/lon grid length
         do i = max(1,II-sgd), min(II+sgd,vix) !Searching i, j that is within the distance of R to II,JJ
         do j = max(1,JJ-sgd), min(JJ+sgd,vjx)

            call calc_dist(vlat(j),vlon(i),Tlat,Tlon,dist)
            if (dist.le.R) then
                   if (mm(num,1).eq.6 .or. mm(num,1).eq.7 .or. mm(num,1).eq.8 ) then
                   Trackpass1(i,j) = .True.
                   elseif (mm(num,1).eq.9 .or. mm(num,1).eq.10 .or. mm(num,1).eq.11 ) then
                   Trackpass2(i,j) = .True.
                   elseif (mm(num,1).eq.12 .or. mm(num,1).eq.1 .or. mm(num,1).eq.2 ) then
                   Trackpass3(i,j) = .True.
                   elseif (mm(num,1).eq.3 .or. mm(num,1).eq.4 .or. mm(num,1).eq.5 ) then
                   Trackpass4(i,j) = .True.
                   endif
            endif
         enddo
         enddo

      ! Track path (linearly interpolated)
      !
      if (point.gt.1) then

        Dlat = lat(num,point)-lat(num,point-1)
        Dlon = lon(num,point)-lon(num,point-1)
        if (abs(Dlat).ge.abs(Dlon)) then
           seg_point = int(abs(Dlat)/0.2) ! 0.2 is adjustable
        else
           seg_point = int(abs(Dlon)/0.2)
        endif
      !
      !  print*,'Tlat(point-1)=',lat(num,point-1),'Tlat(point)=',lat(num,point),'Dlat=',Dlat
      !  print*,'Tlon(point-1)=',lon(num,point-1),'Tlon(point)=',lon(num,point),'Dlon=',Dlon

        if (seg_point.ge.1) then
          do tp = 1, seg_point
          !
             Tlon = lon(num,point-1)+float(tp)*Dlon/float(seg_point)
             Tlat = lat(num,point-1)+(float(tp)/float(seg_point))*Dlat
             JJ = int((Tlat+20)/0.5)+1
             II = int(Tlon/0.5)+1
             if ( mod(Tlat,0.5).gt.0.25 ) JJ = JJ+1
             if ( mod(Tlon,0.5).gt.0.25 ) II = II+1

             sgd= int(R/1.1e5/0.5)+20
             do i = max(1,II-sgd), min(II+sgd,vix)
             do j = max(1,JJ-sgd), min(JJ+sgd,vjx)
               call calc_dist(vlat(j),vlon(i),Tlat,Tlon,dist)
               if (dist.le.R) then
                   if (mm(num,1).eq.6 .or. mm(num,1).eq.7 .or. mm(num,1).eq.8 ) then
                   Trackpass1(i,j) = .True.
                   elseif (mm(num,1).eq.9 .or. mm(num,1).eq.10 .or. mm(num,1).eq.11 ) then
                   Trackpass2(i,j) = .True.
                   elseif (mm(num,1).eq.12 .or. mm(num,1).eq.1 .or. mm(num,1).eq.2 ) then
                   Trackpass3(i,j) = .True.
                   elseif (mm(num,1).eq.3 .or. mm(num,1).eq.4 .or. mm(num,1).eq.5 ) then
                   Trackpass4(i,j) = .True.
                   endif
               endif
             enddo
             enddo
          !
          enddo ! tp loop
        endif ! seg_point.ge.1
        !--------------------
      endif  
   enddo !for point loop
   !print*, 'yy = ', yy(num,point)
   ENDIF
   !
   do vi = 1, vix
   do vj = 1, vjx
      if (Trackpass1(vi,vj)) then
          trackDen1(vi,vj) = trackDen1(vi,vj)+1.
      elseif (Trackpass2(vi,vj)) then
          trackDen2(vi,vj) = trackDen2(vi,vj)+1.
      elseif (Trackpass3(vi,vj)) then
          trackDen3(vi,vj) = trackDen3(vi,vj)+1.
      elseif (Trackpass4(vi,vj)) then
          trackDen4(vi,vj) = trackDen4(vi,vj)+1.
      endif
   enddo
   enddo
   !
   !print *, 'yy = ', yy(num,point)
enddo !for num

!=============================================================

100 format (i6, 7x, i4, i4, i2, i4, 5x, f6.2, 4x, f7.2)
!50 format (i6, 7x, i4)
!100 format (3x, i5, 2x, i3, 11x, i4, i2, i4, 3x, f5.2, 2x, f6.2)

! When we read the data, it's like a cursor moving through the numbers. 

! For example, we have : 
! storm  lifetime   datetime  latitude  longitude
! 24602       124 2020123123     77.25     294.25
! The cursor only reads the data once it's to the right of it. For example, 
! with i5, the cursor will go through 2, 4, 6, 0 and then will be at 2. 
! Since it's on the 2, it will not read it, we need to move it another time 
! (with i6). 
! Now, we skip 7 characters with i7. The cursor will now be on 1. 
! We read the next 4 characters for lifetime. After i4, the cursor is on the 2 
! of datetime, and so on. 
! To understand the above, highlight the numbers with the cursor
! For the float, we need to give the space for the whole number 
! (so f6 if we have 4 digits + the numeral point) and the precision
! (so 2 if there are 2 numbers after the numeral point). 
!=============================================================

call writegrid(fileout, vlon, vlat, trackDen1, trackDen2, trackDen3, trackDen4, vix, vjx)
contains
!stop
!=====================================================================
subroutine writegrid(fileout, var1, var2, var3, var4, var5, var6, ix,jx)
  use netcdf
  implicit none
  character(len=110) fileout
  real, allocatable, dimension(:) :: var1, var2
  real, allocatable, dimension(:,:) :: var3, var4, var5, var6
  integer, dimension(2) :: dimids
  integer ix,jx
  integer ncid5, varid1,varid2,varid3,varid4,varid5, varid6
  integer i_dimid, j_dimid, t_dimid
  integer i_varid, j_varid, t_varid
!
!Creat the netCDF file
!  
  call check(nf90_create(fileout, NF90_CLOBBER, ncid5))
!
!Define the dimensions
!
  call check(nf90_def_dim(ncid5,'longitude', ix, i_dimid))
  call check(nf90_def_dim(ncid5,'latitude', jx, j_dimid))
  call check(nf90_def_var(ncid5,'longitude', NF90_FLOAT, i_dimid, i_varid))
  call check(nf90_def_var(ncid5,'latitude', NF90_FLOAT, j_dimid, j_varid))
  dimids = (/i_dimid, j_dimid/)
  call check(nf90_def_var(ncid5,'trackDen_JJA', NF90_FLOAT, dimids, varid3))
  call check(nf90_def_var(ncid5,'trackDen_SON', NF90_FLOAT, dimids, varid4))
  call check(nf90_def_var(ncid5,'trackDen_DJF', NF90_FLOAT, dimids, varid5))
  call check(nf90_def_var(ncid5,'trackDen_MAM', NF90_FLOAT, dimids, varid6))
  call check(nf90_enddef(ncid5))
  call check(nf90_put_var(ncid5,  i_varid, var1))
  call check(nf90_put_var(ncid5,  j_varid, var2))
  call check(nf90_put_var(ncid5,  varid3, var3))
  call check(nf90_put_var(ncid5,  varid4, var4))
  call check(nf90_put_var(ncid5,  varid5, var5))
  call check(nf90_put_var(ncid5,  varid6, var6))
  call check(nf90_close(ncid5))

end subroutine writegrid
!----------------------------------------------------------
subroutine calc_dist (lat1, lon1, lat2, lon2, dist12)
!
! Author
! Katja Winger (Oct 2008)
!
! Description
! Subroutine to calculate the distance in meters between
! two points (lat1,lon1) and (lat2,lon2) on the sphere
!

  implicit none

  real  lat1, lon1, lat2, lon2, dist12
  real  rlat1, rlon1, rlat2, rlon2
  real  x1(3), x2(3)
  real  d, pi, a

! Define some constants

! PI = 3.14159...
  pi = 2.*ASIN(1.)

! Mean radius of the Earth in meters
  a  = 6371.22E3

! Convert to radians
  rlon1=lon1
  rlon2=lon2
  if ( rlon1 .lt. 0 ) rlon1 = 360. + rlon1
  if ( rlon2 .lt. 0 ) rlon2 = 360. + rlon2

  rlat1=lat1*pi/180.
  rlon1=rlon1*pi/180.
  rlat2=lat2*pi/180.
  rlon2=rlon2*pi/180.

! Locate points in Cartesian space

  x1(1)=COS(rlat1)*COS(rlon1)
  x1(2)=COS(rlat1)*SIN(rlon1)
  x1(3)=SIN(rlat1)

  x2(1)=COS(rlat2)*COS(rlon2)
  x2(2)=COS(rlat2)*SIN(rlon2)
  x2(3)=SIN(rlat2)

! Find shortest distance in Cartesian space (divided by Earth's radius)

  d=SQRT((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)

! Find distance following Earth's surface

  dist12 = 2.*a*ASIN(d/2.)

  return

end subroutine calc_dist
!----------------------------------------------------------
subroutine check(status)
implicit none
integer, intent ( in) :: status

 if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
 end if

end subroutine check

end program Track_density
