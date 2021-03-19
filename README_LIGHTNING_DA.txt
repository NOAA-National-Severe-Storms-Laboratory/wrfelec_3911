The following namelist options (under &physics) are used for the lightning DA Qv nudging:

# Lightning Qv Nudging
rconfig   integer  nudge_lightning        namelist,physics      max_domains   0       rh       "flag for lightning Qv nudging"  ""      ""
rconfig   integer  nudge_light_times      namelist,physics      max_domains   0       rh       "start time in sec (domain relative) for lightning Qv nudging"  ""      ""
rconfig   integer  nudge_light_timee      namelist,physics      max_domains   7200    rh       "end time in sec (domain relative) for lightning Qv nudging"  ""      ""
rconfig   integer  nudge_light_int        namelist,physics      max_domains   600     rh       "time interval in sec of input lightning data files"  ""      ""
rconfig   character  path_to_files        namelist,physics      1             "~/WRFV3/"    rh   "path on local machine of input lightning data files"  ""      ""

Example:

&physics
 mp_physics                          = 26,     26,   
 nssl_cccn                           = 0.5e+08, 0.5e+08,
 nudge_lightning                     = 0, 1,
 nudge_light_times                   = 7200,  0,
 nudge_light_timee                   = 14400, 10800,
 nudge_light_int                     = 600,   3600,
 path_to_files                       = "/home/hpc/afierro/WRFV36/dyn_em/run1/"
 nssl_ipelec                         = 0, 0,

For the pre-processing of the lightning obs, follow the instructions commented out at the top of the stand-alone Fortran code LIGHT_INTERP_WRF_GRID_2.f

The WRF-ELEC release should include phys/module_mp_nudge_light.F with a description of SUBROUTINE nudge_light(nx,ny,nz,qg,qv,gridlight,temp,press,dz,dx,dy,dt), which is called in module_microphysics_driver.F as follows:

! LIGHTNING NUDGING -

      IF (nudge_lightning.eq.1) THEN
     
      if (maxval(gridlight(its:ite,jts:jte)).gt.0) then
        write(0,*) 'TILE MAXVAL GRIDLIGHT',maxval(gridlight(its:ite,jts:jte))
      endif


     CALL nudge_light(ite-its+1,jte-jts+1,kte-kts+1,qg_curr(its:ite,kts:kte,jts:jte),   &
                   & qv_curr(its:ite,kts:kte,jts:jte),gridlight(its:ite,jts:jte),  &
                   & t8w(its:ite,kts:kte,jts:jte),p8w(its:ite,kts:kte,jts:jte),  &
                   & dz8w(its:ite,kts:kte,jts:jte),dx,dy,dt)

      ENDIF
! END LIGHTNING NUDGING


There are key tuning DA params in phys/module_mp_nudge_light.F that you can modify once the system works [see the 2012 paper for details]:

         real ::  RH= 0.95 ! minimum RH for nudging: larger RH -> larger grid volume activated for nudging
         real ::  QG_TRESH=0.003 ! threshold nudging on graupel mass to further constraint nudging volume.
         real ::  t_bot=273.15 !top of nudging layer
         real ::  t_up= 253.15 ! bottom of nudging layer
         real ::  MAX_RH= 1.02
         real ::  BB= 0.25 ! controls the slope of TANH curve
         real ::  CC= 0.25 ! controls spacing of QG curves in eq 1 o
