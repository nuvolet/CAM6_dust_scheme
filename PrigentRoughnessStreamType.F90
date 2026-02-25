
module PrigentRoughnessStreamType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in the Prigent et al. (1997) roughness length streams file
  ! Created by Danny M. Leung 30 Nov 2023
  !
  ! !USES
  use shr_strdata_mod , only : shr_strdata_type, shr_strdata_create  ! dmleung
  use shr_kind_mod   , only: r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod        , only: mpicom, masterproc
  use clm_varctl     , only: iulog
  use abortutils     , only: endrun
  use decompMod      , only: bounds_type
  !use ch4varcon      , only: finundation_mtd

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: prigentroughnessstream_type
     real(r8), pointer, private :: prigent_rghn  (:)         ! Prigent et al. (1997) roughness length (m)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: CalcDragPartition ! Calculate drag partitioning based on input streams
      procedure, public :: UseStreams      ! If streams will be used

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type prigentroughnessstream_type


  ! ! PRIVATE DATA:

  type, private :: streamcontrol_type
     character(len=CL)  :: stream_fldFileName_prigentroughness   ! data Filename
     character(len=CL)  :: stream_meshfile_prigentroughness      ! mesh Filename
     character(len=CL)  :: prigentroughnessmapalgo               ! map algo
     character(len=CL)  :: stream_fldFileName_ch4finundated   ! Filename
     character(len=CL)  :: ch4finundatedmapalgo               ! map algo
     character(len=CL)  :: fldList                            ! List of fields to read
  contains
     procedure, private :: ReadNML     ! Read in namelist
  end type streamcontrol_type

  type(streamcontrol_type), private :: control        ! Stream control data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !    
   ! Initialize the ch4 finundated stream object
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar, get_curr_date
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_mpi_mod      , only : shr_mpi_bcast
   use ndepStreamMod    , only : clm_domain_mct
   use domainMod        , only : ldomain
   use decompMod        , only : bounds_type, gsmap_lnd_gdc2glo
   use mct_mod          , only : mct_ggrid, mct_avect_indexra
   use shr_strdata_mod  , only : shr_strdata_type, shr_strdata_create
   use shr_strdata_mod  , only : shr_strdata_print, shr_strdata_advance
   use spmdMod          , only : comp_id, iam
   !use ch4varcon        , only : finundation_mtd_h2osfc
   !use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   !
   ! arguments
   implicit none
   class(prigentroughnessstream_type) :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: ig, g, n            ! Indices
   type(mct_ggrid)    :: dom_clm          ! domain information 
   type(shr_strdata_type) :: sdat_rghn         ! input data stream
   !integer            :: index_ZWT0       = 0 ! Index of ZWT0 field
   !integer            :: index_F0         = 0 ! Index of F0 field
   !integer            :: index_P3         = 0 ! Index of P3 field
   !integer            :: index_FWS_TWS_A  = 0 ! Index of FWS_TWS_A field
   !integer            :: index_FWS_TWS_B  = 0 ! Index of FWS_TWS_B field
   integer            :: index_Z0a         = 0 ! Index of Z0a field
   integer            :: year                 ! year (0, ...) for nstep+1
   integer            :: mon                  ! month (1, ..., 12) for nstep+1
   integer            :: day                  ! day of month (1, ..., 31) for nstep+1
   integer            :: sec                  ! seconds into current date for nstep+1
   integer            :: mcdate               ! Current model date (yyyymmdd)
   character(len=*), parameter :: stream_name = 'prigentroughness'
   character(len=16), allocatable :: stream_varnames(:) ! array of stream field names
   !character(*), parameter :: subName = "('ch4finundatedstream::Init')"
   !-----------------------------------------------------------------------
   !if ( finundation_mtd /= finundation_mtd_h2osfc )then
      call this%InitAllocate( bounds )
      call control%ReadNML( bounds, NLFileName )

      allocate(stream_varnames(1))
      stream_varnames = (/"Z0a"/)  ! varname in cdf5_Z0a_Prigent-Globe-025x025-09262022.nc, in centimeter

      if (masterproc) then
         write(iulog,*) '  stream_varnames                  = ',stream_varnames
      end if

      if ( this%useStreams() )then
         call clm_domain_mct (bounds, dom_clm)

         call shr_strdata_create(sdat_rghn,name=stream_name,&
           pio_subsystem=pio_subsystem,               & 
           pio_iotype=shr_pio_getiotype(inst_name),   &
           mpicom=mpicom, compid=comp_id,             &
           gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,    &
           nxg=ldomain%ni, nyg=ldomain%nj,            &
           yearFirst=1997,                            &
           yearLast=1997,                             &
           yearAlign=1,                               &
           offset=0,                                  &
           domFilePath='',                            &
           !domFileName=trim(control%stream_fldFileName_ch4finundated), &
           domFileName=trim(control%stream_fldFileName_prigentroughness), &
           domTvarName='time',                        &
           domXvarName='LONGXY' ,                     &
           domYvarName='LATIXY' ,                     &  
           domAreaName='AREA',                        &
           domMaskName='LANDMASK',                    &
           filePath='',                               &
           !filename=(/trim(control%stream_fldFileName_ch4finundated)/),&
           filename=(/trim(control%stream_fldFileName_prigentroughness)/),&
           fldListFile=control%fldList,               &
           fldListModel=control%fldList,              &
           !fldListFile  = stream_varnames,                                    &
           !fldListModel = stream_varnames,                                    &
           fillalgo='none',                           &
           !mapalgo=control%ch4finundatedmapalgo,      &
           mapalgo=control%prigentroughnessmapalgo,    &
           calendar=get_calendar(),                   &
           taxmode='extend'                           )

         if (masterproc) then
            call shr_strdata_print(sdat_rghn,'CLM '//stream_name//' data')
         endif

         !if( finundation_mtd == finundation_mtd_ZWT_inversion )then
         !   index_ZWT0      = mct_avect_indexra(sdat%avs(1),'ZWT0')
         !   index_F0        = mct_avect_indexra(sdat%avs(1),'F0' )
         !   index_P3        = mct_avect_indexra(sdat%avs(1),'P3' )
         !else if( finundation_mtd == finundation_mtd_TWS_inversion )then
         !   index_FWS_TWS_A = mct_avect_indexra(sdat%avs(1),'FWS_TWS_A')
         !   index_FWS_TWS_B = mct_avect_indexra(sdat%avs(1),'FWS_TWS_B')
         !end if

         index_Z0a      = mct_avect_indexra(sdat_rghn%avs(1),'Z0a')

         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         !call get_curr_date(year, mon, day, sec)
         year = 1997
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat_rghn, mcdate, sec, mpicom, 'prigentrghn')

         ! Get the data
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
               this%prigent_rghn(g) = sdat_rghn%avs(1)%rAttr(index_Z0a,ig)
         end do
      end if
   !end if

  end subroutine Init

  !-----------------------------------------------------------------------
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type) :: this
    !
    ! !LOCAL VARIABLES:
    if ( trim(control%stream_fldFileName_prigentroughness) == '' )then
       UseStreams = .false.  ! this won't happen and UseStreams will always be true
    else
       UseStreams = .true. 
    end if
  end function UseStreams

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !use ch4varcon     , only: finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%prigent_rghn     (begg:endg))            ;  this%prigent_rghn     (:)   = nan

  end subroutine InitAllocate

  !==============================================================================
  subroutine CalcDragPartition(this, bounds, dpfct_rock)
    !
    ! !DESCRIPTION:
    ! Commented below by Danny M. Leung 31 Dec 2022
    ! Calculate the drag partition effect of friction velocity due to surface roughness following 
    ! Leung et al. (2022).  This module is used in the dust emission module DUSTMod.F90 for 
    ! calculating drag partitioning. The drag partition equation comes from Marticorena and 
    ! Bergametti (1995) with constants modified by Darmenova et al. (2009). Here it is assumed 
    ! that this equation is used only over arid/desertic regions, such that Catherine Prigent's
    ! roughness measurements represents mostly rocks. For more vegetated areas, the vegetation
    ! roughness and drag partitioning are calculated in the DustEmission subroutine. This 
    ! subroutine is used in the InitCold subroutine of DUSTMod.F90.
    !
    ! !USES:
    use PatchType               , only : patch
    use landunit_varcon         , only : istdlak
    use LandunitType            , only : lun
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type)             :: this
    type(bounds_type)              , intent(in)    :: bounds
    real(r8)                       , intent(inout) :: dpfct_rock(bounds%begp:)      ! [fraction] rock drag partition factor (roughness effect)
    !
    ! !LOCAL VARIABLES:
    integer  :: g, p, fp, l    ! Indices
    real(r8) :: z0s         ! smooth roughness length (m)

    ! constants
    real(r8), parameter :: D_p = 130e-6_r8           ! [m] Medium soil particle diameter, assuming a global constant of ~130 um following Leung et al. (2022)
    real(r8), parameter :: X = 10_r8                 ! [m] distance downwind of the roughness element (rock). Assume estimating roughness effect at a distance of 10 m following Leung et al. (2022)

    SHR_ASSERT_ALL_FL((ubound(dpfct_rock)        == (/bounds%endp/)), sourcefile, __LINE__)


    ! dmleung: this loop calculates the drag partition effect (or roughness effect) of rocks. We save the drag partition factor as a patch level quantity.
    z0s = 2_r8 * D_p / 30_r8 ! equation from Frank M. White (2006). Here we assume soil medium size is a global constant, and so is smooth roughness length.
    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       l = patch%landunit(p)
       if (lun%itype(l) /= istdlak) then
          dpfct_rock(p) = 1._r8 - ( log(this%prigent_rghn(g)*0.01_r8/z0s) / log(0.7_r8*(X/z0s)**0.8_r8) ) ! Calculating rock drag partition factor using Marticorena and Bergametti (1995). 0.01 is used to convert Z0a from centimeter to meter.
       end if
    end do

  end subroutine CalcDragPartition

  !==============================================================================

  subroutine ReadNML(this, bounds, NLFilename)
   !    
   ! Read the namelist data stream information.  
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use fileutils        , only : getavu, relavu
   !use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   !
   ! arguments
   implicit none
   class(streamcontrol_type) :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_prigentroughness = ' '
   character(len=CL)  :: stream_meshfile_prigentroughness = ' '
   character(len=CL)  :: prigentroughnessmapalgo = 'bilinear'
   character(len=*), parameter :: namelist_name = 'prigentroughness'    ! MUST agree with name in namelist and read
   character(len=*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(len=*), parameter :: subName = "('prigentroughness::ReadNML')"
   !character(len=*), parameter :: F00 = "('(ch4finundated_readnml) ',4a)"
   !-----------------------------------------------------------------------

   namelist /prigentroughness/ &               ! MUST agree with namelist_name above
        prigentroughnessmapalgo,  stream_fldFileName_prigentroughness, stream_meshfile_prigentroughness

   ! Default values for namelist

   ! Read ch4finundateddyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=prigentroughness,iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_fldFileName_prigentroughness, mpicom)
   call shr_mpi_bcast(stream_meshfile_prigentroughness , mpicom)
   call shr_mpi_bcast(prigentroughnessmapalgo            , mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  stream_fldFileName_prigentroughness = ',stream_fldFileName_prigentroughness
      write(iulog,*) '  stream_meshfile_prigentroughness    = ',stream_meshfile_prigentroughness
      write(iulog,*) '  prigentroughnessmapalgo             = ',prigentroughnessmapalgo
      write(iulog,*) ' '
   endif
   this%stream_fldFileName_prigentroughness = stream_fldFileName_prigentroughness
   this%stream_meshfile_prigentroughness    = stream_meshfile_prigentroughness
   this%prigentroughnessmapalgo             = prigentroughnessmapalgo
   this%fldList = "Z0a"
 end subroutine ReadNML

end module PrigentRoughnessStreamType
