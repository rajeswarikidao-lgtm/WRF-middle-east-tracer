!
! This module contains routines to initialize tracers, handle boundary conditions
! and other stuff related to tracers. As of WRFV3.2 it will users should compile 
! WRF-Chem for tracer runs, if they want full dispersion. Only when 
! WRF-Chem is compiled will turbulent and non-resolved convective transport be treated.
! When compiled with WRF_CHEM, tracer transport should work properly with nesting or
! supplied boundary conditions. Without compiling WRF-Chem, option 
! TRACER_TEST1 will partially work (no boundary conditions, 
! no turbulent transport, no subgrid-scale convection)
!
! Original version of the module is written by Georg Grell (Dec 2009). 
! Options for TRACER_TEST1 and TRACER_TEST2 supplied by Jeff Lee (NCAR)
! Current tracer options:
!
!   (1) TRACER_SMOKE: This needs the biomass burning module to also be active. 
!       It will then use smoke (CO emissions from fire) as tracer. One
!       variabe only.
!
!   (2) TRACER_TEST1 and TRACER_TEST2: 8 tracers, the only difference inbetween 
!       these options are p_tr17_3 and p_tr17_4, which are also filled with
!       CO emissions from fire for TRACER_TEST2. The other tracers are defined as:
!
!            tr17_1 : horizontal boundaries tracer
!            tr17_2 : horizontal boundaries tracer decaying with e-folding time of 1 day
!            tr17_3 : surface tracer (smoke for TRACER_TEST2)
!            tr17_4 : surface tracer (smoke for TRACER_TEST2)
!                     decaying with e-folding time of 1 day
!            tr17_5 : stratosphere tracer
!            tr17_6 : stratosphere tracer decaying with e-folding time of 1 day
!            tr17_7 : boundary layer tracer     
!            tr17_8 : boundary layer tracer decaying with e-folding time of 1 day



MODULE module_input_tracer
USE module_input_tracer_data
#if ( WRF_CHEM == 1 )
USE module_state_description, only:tracer_smoke,tracer_test1,tracer_test2,param_first_scalar,p_tr17_1,p_tr17_2,p_tr17_3,p_tr17_4,p_tr17_5,p_tr17_6,p_tr17_7,p_tr17_8
! ===== added by Rajeswari (tracer_opt=11 isotope tracers) : START =====
USE module_state_description, only:p_isotope1,p_isotope2
! ===== added by Rajeswari : END =====
#else
USE module_state_description, only:tracer_test1,tracer_test2,param_first_scalar,p_tr17_1,p_tr17_2,p_tr17_3,p_tr17_4,p_tr17_5,p_tr17_6,p_tr17_7,p_tr17_8
#endif
CONTAINS
   SUBROUTINE initialize_tracer (chem,chem_in_opt,         &
                                       tracer_opt,num_chem,&
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
      INTEGER,      INTENT(IN   )    :: chem_in_opt,tracer_opt,num_chem
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme,num_chem ), INTENT(INOUT) :: chem
#if ( WRF_CHEM == 1 ) 
      if(chem_in_opt == 1 )return
#endif
      if     (tracer_opt == TRACER_TEST1)then
       chem(:,:,:,:)=.0
#if ( WRF_CHEM == 1 )
      else if(tracer_opt == TRACER_TEST2)then
       chem(:,:,:,:)=.0
      else if(tracer_opt == TRACER_SMOKE)then
       chem(:,:,:,:)=.08
#endif
      endif
   END SUBROUTINE initialize_tracer
#if (EM_CORE == 1 ) 
   SUBROUTINE flow_dep_bdy_tracer  (  chem,                                       &
                               chem_bxs,chem_btxs,                                  &
                               chem_bxe,chem_btxe,                                  &
                               chem_bys,chem_btys,                                  &
                               chem_bye,chem_btye,                                  &
                               dt,                                              &
                               spec_bdy_width,z,                                &
                               have_bcs_chem,                        & 
                               u, v, tracer_opt, alt, & 
                               t,pb,p,t0,p1000mb,rcp,ph,phb,g, &
                               spec_zone, ic,           &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

!  This subroutine sets zero gradient conditions for outflow and a set profile value
!  for inflow in the boundary specified region. Note that field must be unstaggered.
!  The velocities, u and v, will only be used to check their sign (coupled vels OK)
!  spec_zone is the width of the outer specified b.c.s that are set here.
!  (JD August 2000)

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: tracer_opt
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: spec_zone,spec_bdy_width,ic
      REAL,         INTENT(IN   )    :: dt


      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) :: chem
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width), INTENT(IN   ) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width), INTENT(IN   ) :: chem_bys, chem_bye, chem_btys, chem_btye
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: z
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: alt
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: u
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: v
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,         &
          INTENT(IN   ) ::                                           &
                               ph,phb,t,pb,p
   real, INTENT (IN) :: g,rcp,t0,p1000mb

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: i_bdy_method
      real tempfac,convfac
      logical, optional    :: have_bcs_chem

      ibs = ids
      ibe = ide-1
      itf = min(ite,ide-1)
      jbs = jds
      jbe = jde-1
      jtf = min(jte,jde-1)
      ktf = kde-1

! i_bdy_method determines which "bdy_chem_value" routine to use
!   1=smoke, CO background
      i_bdy_method = 0
        if (tracer_opt == TRACER_TEST1 ) then
          i_bdy_method = 2
        end if   
#if ( WRF_CHEM == 1 )
        if (tracer_opt == TRACER_TEST2 ) then
          i_bdy_method = 2
        end if   
        if (tracer_opt == TRACER_SMOKE ) then
          i_bdy_method = 1
        end if
      if (have_bcs_chem) i_bdy_method =6
#endif
      if (ic .lt. param_first_scalar) i_bdy_method = 0

      IF (jts - jbs .lt. spec_zone) THEN
! Y-start boundary
        DO j = jts, min(jtf,jbs+spec_zone-1)
          b_dist = j - jbs
          DO k = kts, ktf
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbs+spec_zone)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
#if ( WRF_CHEM == 1 )
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bys(i,k,1),chem_btys(i,k,1),dt,ic)
#endif
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 
      IF (jbe - jtf .lt. spec_zone) THEN 
! Y-end boundary 
        DO j = max(jts,jbe-spec_zone+1), jtf 
          b_dist = jbe - j 
          DO k = kts, ktf 
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j+1) .gt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbe-spec_zone)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
#if ( WRF_CHEM == 1 )
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bye(i,k,1),chem_btye(i,k,1),dt,ic)
#endif
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (its - ibs .lt. spec_zone) THEN
! X-start boundary
        DO i = its, min(itf,ibs+spec_zone-1)
          b_dist = i - ibs
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(ibs+spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
#if ( WRF_CHEM == 1 )
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bxs(j,k,1),chem_btxs(j,k,1),dt,ic)   
#endif
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (ibe - itf .lt. spec_zone) THEN
! X-end boundary
        DO i = max(its,ibe-spec_zone+1), itf
          b_dist = ibe - i
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i+1,k,j) .gt. 0.)THEN
                chem(i,k,j) = chem(ibe-spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
#if ( WRF_CHEM == 1 )
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bxe(j,k,1),chem_btxe(j,k,1),dt,ic)
#endif
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

   END SUBROUTINE flow_dep_bdy_tracer
#else
#if ( WRF_CHEM == 1 )
   SUBROUTINE flow_dep_bdy_tracer  (  chem, chem_b,chem_bt,dt,                    &
                               spec_bdy_width,z,                                &
                               ijds, ijde,have_bcs_chem,                        & 
                               u, v, tracer_opt, alt, & 
                               t,pb,p,t0,p1000mb,rcp,ph,phb,g, &
                               spec_zone, ic,           &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

!  This subroutine sets zero gradient conditions for outflow and a set profile value
!  for inflow in the boundary specified region. Note that field must be unstaggered.
!  The velocities, u and v, will only be used to check their sign (coupled vels OK)
!  spec_zone is the width of the outer specified b.c.s that are set here.
!  (JD August 2000)

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: tracer_opt
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: ijds,ijde
      INTEGER,      INTENT(IN   )    :: spec_zone,spec_bdy_width,ic
      REAL,         INTENT(IN   )    :: dt


      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) :: chem
      REAL,  DIMENSION( ijds:ijde , kds:kde , spec_bdy_width, 4 ), INTENT(IN   ) :: chem_b
      REAL,  DIMENSION( ijds:ijde , kds:kde , spec_bdy_width, 4 ), INTENT(IN   ) :: chem_bt
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: z
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: alt
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: u
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: v
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,         &
          INTENT(IN   ) ::                                           &
                               ph,phb,t,pb,p
   real, INTENT (IN) :: g,rcp,t0,p1000mb

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: i_bdy_method
      real tempfac,convfac
      real       :: tracer_bv_def
      logical    :: have_bcs_chem

      tracer_bv_def = conmin
      ibs = ids
      ibe = ide-1
      itf = min(ite,ide-1)
      jbs = jds
      jbe = jde-1
      jtf = min(jte,jde-1)
      ktf = kde-1

      i_bdy_method = 0
        if (config_flags%tracer_opt == TRACER_SMOKE ) then
          i_bdy_method = 1
        end if
      if (have_bcs_chem) i_bdy_method =6
      if (ic .lt. param_first_scalar) i_bdy_method = 0

!----------------------------------------------------------------------
      IF (jts - jbs .lt. spec_zone) THEN
! Y-start boundary
        DO j = jts, min(jtf,jbs+spec_zone-1)
          b_dist = j - jbs
          DO k = kts, ktf
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbs+spec_zone)
              ELSE
                if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_b(i,k,1,P_YSB),chem_bt(i,k,1,P_YSB),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 
      IF (jbe - jtf .lt. spec_zone) THEN 
! Y-end boundary 
        DO j = max(jts,jbe-spec_zone+1), jtf 
          b_dist = jbe - j 
          DO k = kts, ktf 
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j+1) .gt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbe-spec_zone)
              ELSE
                if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_b(i,k,1,P_YEB),chem_bt(i,k,1,P_YEB),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (its - ibs .lt. spec_zone) THEN
! X-start boundary
        DO i = its, min(itf,ibs+spec_zone-1)
          b_dist = i - ibs
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(ibs+spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_b(j,k,1,P_XSB),chem_bt(j,k,1,P_XSB),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (ibe - itf .lt. spec_zone) THEN
! X-end boundary
        DO i = max(its,ibe-spec_zone+1), itf
          b_dist = ibe - i
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i+1,k,j) .gt. 0.)THEN
                chem(i,k,j) = chem(ibe-spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_b(j,k,1,P_XEB),chem_bt(j,k,1,P_XEB),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

   END SUBROUTINE flow_dep_bdy_tracer
#endif
#endif
   SUBROUTINE set_tracer(dtstep,ktau,pbl_h,tracer,t,tracer_opt,num_tracer,&
                         z,ht,ids,ide, jds,jde, kds,kde,                  & 
                               ims,ime, jms,jme, kms,kme,                 & 
                               its,ite, jts,jte, kts,kte                  )
      INTEGER,      INTENT(IN   )    :: ktau,tracer_opt,num_tracer
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer ), INTENT(INOUT) :: tracer
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme ), INTENT(IN) :: t,z
      REAL,  DIMENSION(ims:ime,jms:jme ), INTENT(IN) :: PBL_H,HT
      REAL,  INTENT(IN) :: dtstep
      INTEGER:: count_trop,count_pbl
!
! this is for tracer options tracer_test1 and tracer_test2
!
    factor_decay = 1./(86400./dtstep)

!-- decay, every time step (ktau), whole domain

    tracer(its:ite,kts:kte,jts:jte,p_tr17_2) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_2) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_4) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_4) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_6) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_6) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_8) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_8) * (1. - factor_decay)
 IF (ktau .ge. 2) THEN
    
!-- every time step, every grid point, restore some tracer

!(1)level 1 restore to 1.0
    if(tracer_opt == TRACER_TEST1)then   
       tracer(its:ite,kts,jts:jte,p_tr17_3)     = 1.0
       tracer(its:ite,kts,jts:jte,p_tr17_4)     = 1.0
    endif
       
    do i= its,ite
    do j= jts,jte
 
!(2)every level above tropopause (t minimum), restore to 1.0

!-- get levels of tropopause (count_trop)

       count_trop = minloc(t(i,kts:kte,j),1)

       tracer(i,count_trop:kte,j,p_tr17_5) = 1.0
       tracer(i,count_trop:kte,j,p_tr17_6) = 1.0

!(3)every level below pblh, restore to 1.0

!-- get levels in pbl (count_pbl)

       count_pbl = 0

       do k=kts,kte
          if ( (z(i,k,j)-ht(i,j)) .le. pbl_h(i,j) ) then
             count_pbl = count_pbl + 1
          endif
       end do

       if (count_pbl .ge. 1) then
          tracer(i,kts:count_pbl,j,p_tr17_7) = 1.0
          tracer(i,kts:count_pbl,j,p_tr17_8) = 1.0
       endif

    end do   ! j
    end do   ! i

 ENDIF   ! ktau  
   END SUBROUTINE set_tracer
! ===== added by Rajeswari (tracer_opt=11 isotope tracers) : START =====
!-----------------------------------------------------------------------
! set_isotope_tracer : called every time step from chem_driver when
! tracer_opt = 11. Tracer unit: Bq per kg of dry air.
! For each isotope:
!  (1) decay over the tile: A = A*exp(-ln2*dt/(T_half*3600)); T_half<=0 -> none
!  (2) release at cell (i,k,j):
!        pulse_start_h <= t < pulse_end_h  -> Q = pulse_rate      [Bq/h]
!        else t >= cont_start_h            -> Q = cont_rate       [Bq/h]
!        else                              -> Q = 0
!      dA [Bq/kg] = Q*dt/3600 / M_air,
!      M_air [kg] = (dx*dy/(msftx*msfty)) * dz8w / alt   (dry air in the cell)
!      Only the tile that owns the cell adds it (MPI/OpenMP safe).
! t = hours since simulation start (curr_secs from chem_driver).
!-----------------------------------------------------------------------
   SUBROUTINE set_isotope_tracer(dtstep, curr_secs, tracer, num_tracer,        &
                 alt, dz8w, dx, dy, msftx, msfty,                            &
                 hl1, i1, j1, k1, prate1, pst1, pend1, crate1, cst1,         &
                 hl2, i2, j2, k2, prate2, pst2, pend2, crate2, cst2,         &
                 ims,ime, jms,jme, kms,kme,                                  &
                 its,ite, jts,jte, kts,kte )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: num_tracer
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme
      INTEGER, INTENT(IN) :: its,ite, jts,jte, kts,kte
      REAL,    INTENT(IN) :: dtstep, curr_secs, dx, dy
      REAL,    INTENT(IN) :: hl1, prate1, pst1, pend1, crate1, cst1
      REAL,    INTENT(IN) :: hl2, prate2, pst2, pend2, crate2, cst2
      INTEGER, INTENT(IN) :: i1, j1, k1, i2, j2, k2
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer), INTENT(INOUT) :: tracer
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),            INTENT(IN)    :: alt, dz8w
      REAL, DIMENSION(ims:ime,jms:jme),                    INTENT(IN)    :: msftx, msfty
      REAL :: t_hours

      t_hours = curr_secs / 3600.
#if ( WRF_CHEM == 1 )
      CALL decay_and_release(p_isotope1, hl1, i1, j1, k1, prate1, pst1, pend1, crate1, cst1)
      CALL decay_and_release(p_isotope2, hl2, i2, j2, k2, prate2, pst2, pend2, crate2, cst2)
#endif

   CONTAINS

      SUBROUTINE decay_and_release(ip, hl, isrc, jsrc, ksrc, prate, pst, pend, crate, cst)
         INTEGER, INTENT(IN) :: ip, isrc, jsrc, ksrc
         REAL,    INTENT(IN) :: hl, prate, pst, pend, crate, cst
         REAL :: rate, air_mass
         REAL*8 :: rn_fdec   ! Rajeswari: Stage A3 (replaces decay_factor)

         IF (ip .lt. param_first_scalar) RETURN      ! species not allocated

!-- (1) radioactive decay
! Rajeswari: Stage A3 - quantised, double-precision decay factor shared with the chem species (rn_decay_factor);
! Rajeswari: the earlier per-step single-precision factor rounded to 1.0 for Cs-137 (dt=30 s); hl<=0 = no decay
         rn_fdec = rn_decay_factor(hl, DBLE(curr_secs), DBLE(dtstep))
         IF (rn_fdec .lt. 1.0D0) THEN
            tracer(its:ite,kts:kte,jts:jte,ip) = REAL( DBLE(tracer(its:ite,kts:kte,jts:jte,ip)) * rn_fdec )
         ENDIF

!-- (2) point-source release
         rate = 0.
         IF (t_hours .ge. pst .and. t_hours .lt. pend) THEN
            rate = prate
         ELSE IF (t_hours .ge. cst) THEN
            rate = crate
         ENDIF

         IF (rate .gt. 0.) THEN
            IF (isrc .ge. its .and. isrc .le. ite .and.        &
                jsrc .ge. jts .and. jsrc .le. jte .and.        &
                ksrc .ge. kts .and. ksrc .le. kte) THEN
               air_mass = dx*dy/(msftx(isrc,jsrc)*msfty(isrc,jsrc)) * dz8w(isrc,ksrc,jsrc) / alt(isrc,ksrc,jsrc)
               tracer(isrc,ksrc,jsrc,ip) = tracer(isrc,ksrc,jsrc,ip) + rate*dtstep/3600. / air_mass
            ENDIF
         ENDIF
      END SUBROUTINE decay_and_release

   END SUBROUTINE set_isotope_tracer
! ===== added by Rajeswari : END =====
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE bdy_tracer_value ( trac, trac_b, trac_bt, dt,ic)
                                  
    IMPLICIT NONE

    REAL,    intent(OUT)  :: trac
    REAL,    intent(IN)   :: trac_b
    REAL,    intent(IN)   :: trac_bt
    REAL,    intent(IN)   :: dt
    INTEGER, intent(IN)   :: ic

    REAL                  :: epsilc = 1.E-12
!   CHARACTER (LEN=80) :: message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     if( ntracer .GT. numtracer) then
!       message = ' Input_tracer_profile: wrong number of tracers'
!       return
!       CALL WRF_ERROR_FATAL ( message )
!     endif
     
      trac=max(epsilc,trac_b + trac_bt * dt)

      RETURN
  END SUBROUTINE bdy_tracer_value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ===== added by Rajeswari (radionuclide j->coarse coagulation transfer, Stage A2b) : START =====
! Rajeswari: moves tagged activity from the j mode to the coarse mode: iodaj->iodac and csaj->csac
! Rajeswari: amount moved = f*A_j with f = 1-exp(-L*dt*factor), L = rn_lossjc from coagrate_jc [s-1]
! Rajeswari: host aerosol number and moments are not changed (radionuclide mass is negligible)
   SUBROUTINE rn_coag_transfer(chem, num_chem, rn_lossjc, dtstep, factor,   &
                 ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte )
      USE module_state_description, only: p_iodaj, p_iodac, p_csaj, p_csac, param_first_scalar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: num_chem
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),          INTENT(IN)    :: rn_lossjc
      REAL, INTENT(IN) :: dtstep, factor
      INTEGER :: i, k, j
      REAL*8  :: x, f, da

! Rajeswari: do nothing if the radionuclide species are not part of this chem_opt
      IF ( p_iodaj .LT. param_first_scalar .OR. p_iodac .LT. param_first_scalar .OR. &
           p_csaj  .LT. param_first_scalar .OR. p_csac  .LT. param_first_scalar ) RETURN

      DO j = jts, jte
        DO k = kts, kte
          DO i = its, ite
            x = DBLE(rn_lossjc(i,k,j))*DBLE(dtstep)*DBLE(factor)
            IF ( x .LE. 0.0D0 ) CYCLE
! Rajeswari: fraction moved in this step; series form keeps full accuracy for very small x
            IF ( x .LT. 1.0D-5 ) THEN
              f = x*(1.0D0 - 0.5D0*x)
            ELSE
              f = 1.0D0 - EXP(-x)
            ENDIF
            da = f*DBLE(chem(i,k,j,p_iodaj))
            chem(i,k,j,p_iodaj) = REAL( DBLE(chem(i,k,j,p_iodaj)) - da )
            chem(i,k,j,p_iodac) = REAL( DBLE(chem(i,k,j,p_iodac)) + da )
            da = f*DBLE(chem(i,k,j,p_csaj))
            chem(i,k,j,p_csaj)  = REAL( DBLE(chem(i,k,j,p_csaj))  - da )
            chem(i,k,j,p_csac)  = REAL( DBLE(chem(i,k,j,p_csac))  + da )
          END DO
        END DO
      END DO
   END SUBROUTINE rn_coag_transfer
! ===== added by Rajeswari : END =====
! ===== added by Rajeswari (radioactive decay, Stage A3) : START =====
! Rajeswari: decay factor for one model step that ends at elapsed time t_end [s] and has length dt [s]
! Rajeswari: lambda = ln2/(hl*3600); decay is applied in quanta of DELTA = 1e-4 e-foldings:
! Rajeswari:   n(t) = floor(lambda*t/DELTA),  factor = exp(-DELTA*(n(t_end) - n(t_end-dt)))
! Rajeswari: why: in single precision exp(-lambda*dt) rounds to exactly 1.0 for Cs-137 with dt = 30 s,
! Rajeswari: so a per-step factor never decays Cs-137; with quanta the total decay is exact at each
! Rajeswari: quantum and never lags the exact solution by more than 0.01 %; it depends only on elapsed
! Rajeswari: time, so it needs no saved state and is restart-safe. hl <= 0 means no decay (factor = 1)
   REAL*8 FUNCTION rn_decay_factor(hl, t_end, dt)
      IMPLICIT NONE
      REAL,   INTENT(IN) :: hl            ! half-life [h]
      REAL*8, INTENT(IN) :: t_end, dt     ! elapsed time at the end of the step, step length [s]
      REAL*8, PARAMETER  :: delta = 1.0D-4
      REAL*8 :: lam, n1, n0
      rn_decay_factor = 1.0D0
      IF ( hl .LE. 0. ) RETURN
      lam = LOG(2.0D0) / ( DBLE(hl)*3600.0D0 )
      n1  = AINT( lam*MAX(t_end,     0.0D0)/delta )
      n0  = AINT( lam*MAX(t_end-dt,  0.0D0)/delta )
      IF ( n1 .GT. n0 ) rn_decay_factor = EXP( -delta*(n1 - n0) )
   END FUNCTION rn_decay_factor

! Rajeswari: radioactive decay of the chem_opt=100 radionuclide species, called every model step
! Rajeswari: I-131 (iodg, iodaj, iodac) uses hl_i, Cs-137 (csaj, csac) uses hl_cs [h]
   SUBROUTINE rn_decay(chem, num_chem, hl_i, hl_cs, t_end, dt,         &
                 ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte )
      USE module_state_description, only: p_iodg, p_iodaj, p_iodac, p_csaj, p_csac, param_first_scalar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: num_chem
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
      REAL,   INTENT(IN) :: hl_i, hl_cs
      REAL*8, INTENT(IN) :: t_end, dt
      REAL*8 :: fi, fc

! Rajeswari: do nothing if the radionuclide species are not part of this chem_opt
      IF ( p_iodg  .LT. param_first_scalar .OR. p_iodaj .LT. param_first_scalar .OR. &
           p_iodac .LT. param_first_scalar .OR. p_csaj  .LT. param_first_scalar .OR. &
           p_csac  .LT. param_first_scalar ) RETURN

      fi = rn_decay_factor(hl_i,  t_end, dt)
      fc = rn_decay_factor(hl_cs, t_end, dt)
! Rajeswari: multiply in double precision so the 1e-4 decrement is applied without rounding loss
      IF ( fi .LT. 1.0D0 ) THEN
         chem(its:ite,kts:kte,jts:jte,p_iodg)  = REAL( DBLE(chem(its:ite,kts:kte,jts:jte,p_iodg))  * fi )
         chem(its:ite,kts:kte,jts:jte,p_iodaj) = REAL( DBLE(chem(its:ite,kts:kte,jts:jte,p_iodaj)) * fi )
         chem(its:ite,kts:kte,jts:jte,p_iodac) = REAL( DBLE(chem(its:ite,kts:kte,jts:jte,p_iodac)) * fi )
      ENDIF
      IF ( fc .LT. 1.0D0 ) THEN
         chem(its:ite,kts:kte,jts:jte,p_csaj)  = REAL( DBLE(chem(its:ite,kts:kte,jts:jte,p_csaj))  * fc )
         chem(its:ite,kts:kte,jts:jte,p_csac)  = REAL( DBLE(chem(its:ite,kts:kte,jts:jte,p_csac))  * fc )
      ENDIF
   END SUBROUTINE rn_decay
! ===== added by Rajeswari : END =====
! ===== added by Rajeswari (radionuclide species, Stage C1) : START =====
! Rajeswari: point-source release into the chem_opt=100 radionuclide species, called every model step from chem_driver
! Rajeswari: identical source cell, time profile and Bq -> Bq/kg conversion as decay_and_release (tracer_opt=11):
! Rajeswari:   Q [Bq/h] = pulse rate for pulse_start <= t < pulse_end, else cont rate for t >= cont_start, else 0
! Rajeswari:   dA [Bq per kg dry air] = Q*dt/3600 / M_air,  M_air = dx*dy/(msftx*msfty) * dz8w / alt
! Rajeswari: I-131 (iso1 settings): fraction f_part into iodaj, (1 - f_part) into iodg
! Rajeswari: Cs-137 (iso2 settings): all into csaj; the coarse tags fill only through the j->c transfer
   SUBROUTINE rn_release(chem, num_chem, dtstep, curr_secs, alt, dz8w, dx, dy, msftx, msfty, f_part, &
                 i1, j1, k1, prate1, pst1, pend1, crate1, cst1,                                    &
                 i2, j2, k2, prate2, pst2, pend2, crate2, cst2,                                    &
                 ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte )
      USE module_state_description, only: p_iodg, p_iodaj, p_csaj, param_first_scalar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: num_chem
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),          INTENT(IN)    :: alt, dz8w
      REAL, DIMENSION(ims:ime,jms:jme),                  INTENT(IN)    :: msftx, msfty
      REAL,    INTENT(IN) :: dtstep, curr_secs, dx, dy, f_part
      REAL,    INTENT(IN) :: prate1, pst1, pend1, crate1, cst1
      REAL,    INTENT(IN) :: prate2, pst2, pend2, crate2, cst2
      INTEGER, INTENT(IN) :: i1, j1, k1, i2, j2, k2
      REAL :: t_hours, fp, da

! Rajeswari: do nothing if the radionuclide species are not part of this chem_opt
      IF ( p_iodg .LT. param_first_scalar .OR. p_iodaj .LT. param_first_scalar .OR. &
           p_csaj .LT. param_first_scalar ) RETURN

      t_hours = curr_secs / 3600.
      fp = MIN( MAX(f_part, 0.), 1. )      ! Rajeswari: particulate fraction limited to [0,1]

! Rajeswari: I-131
      da = rn_point_release(i1, j1, k1, prate1, pst1, pend1, crate1, cst1)
      IF ( da .GT. 0. ) THEN
         chem(i1,k1,j1,p_iodg)  = chem(i1,k1,j1,p_iodg)  + (1. - fp) * da
         chem(i1,k1,j1,p_iodaj) = chem(i1,k1,j1,p_iodaj) + fp * da
      ENDIF
! Rajeswari: Cs-137
      da = rn_point_release(i2, j2, k2, prate2, pst2, pend2, crate2, cst2)
      IF ( da .GT. 0. ) THEN
         chem(i2,k2,j2,p_csaj) = chem(i2,k2,j2,p_csaj) + da
      ENDIF

   CONTAINS

! Rajeswari: activity [Bq/kg] added to the source cell in this step; 0 when there is no release
! Rajeswari: or when the cell is not on this tile (only the owning tile adds it - MPI/OpenMP safe)
      REAL FUNCTION rn_point_release(isrc, jsrc, ksrc, prate, pst, pend, crate, cst)
         INTEGER, INTENT(IN) :: isrc, jsrc, ksrc
         REAL,    INTENT(IN) :: prate, pst, pend, crate, cst
         REAL :: rate, air_mass
         rn_point_release = 0.
         rate = 0.
         IF (t_hours .ge. pst .and. t_hours .lt. pend) THEN
            rate = prate
         ELSE IF (t_hours .ge. cst) THEN
            rate = crate
         ENDIF
         IF (rate .gt. 0.) THEN
            IF (isrc .ge. its .and. isrc .le. ite .and.        &
                jsrc .ge. jts .and. jsrc .le. jte .and.        &
                ksrc .ge. kts .and. ksrc .le. kte) THEN
               air_mass = dx*dy/(msftx(isrc,jsrc)*msfty(isrc,jsrc)) * dz8w(isrc,ksrc,jsrc) / alt(isrc,ksrc,jsrc)
               rn_point_release = rate*dtstep/3600. / air_mass
            ENDIF
         ENDIF
      END FUNCTION rn_point_release
   END SUBROUTINE rn_release
! ===== added by Rajeswari : END =====
! ===== added by Rajeswari (radionuclide species, Stage E2: large-scale wet scavenging) : START =====
! Rajeswari: wet removal of the radionuclide species (and optionally the host aerosol) by large-scale
! Rajeswari: (microphysics) precipitation, called every model step from chem_driver after rn_coag_transfer.
! Rajeswari: Column algorithm, top to bottom, per species:
! Rajeswari:   1. re-evaporation: where the rain flux F = precr decreases downward, the same fraction of
! Rajeswari:      the activity carried by the falling rain is returned to the air of that layer
! Rajeswari:   2. removal dC = C [1 - exp(-L dt)] in the layer; the removed activity joins the falling rain
! Rajeswari:   3. what reaches the ground is the wet deposition (RN_WDLS_*, Bq m-2)
! Rajeswari: Removal rate L (rn_wdep_opt = 1 base, 2 NAME reference):
! Rajeswari:   cloud layers (qc+qi > cld_thresh): in-cloud rainout, L = eps * base
! Rajeswari:      base = qlsink (opt 1; cloud water -> precipitation rate / qc from the microphysics)
! Rajeswari:      base = A * P**B with A = 8.4e-5 s-1, B = 0.79 (opt 2; NAME)
! Rajeswari:      eps = 1.0 accumulation, 0.5 coarse (as in the Grell convective scavenging), 0 Aitken
! Rajeswari:   other layers with rain falling through (P > 1e-4 mm/h): below-cloud washout
! Rajeswari:      opt 1: L = rn_wdep_scale * <A(d) P**B(d)> averaged over the lognormal mode (Wang, Zhang and
! Rajeswari:             Moran 2014, GMD 7, 799, rain, Table 8), mass- or number-weighted, with the MADE
! Rajeswari:             geometric mean diameters rn_dg* (m)
! Rajeswari:      opt 2: L = rn_wdep_scale * A * P**B (NAME, all particle sizes)
! Rajeswari:   iodg (treated as I2): in cloud L = f_aq * base (Henry's law partitioning into cloud water);
! Rajeswari:      below cloud L = H' P / H_col (Henry-limited washout, HYSPLIT form), H' = K_H R T
! Rajeswari: P = local rain rate from precr [kg m-2 s-1 = mm s-1]; test mode rn_wd_test_rain > 0 prescribes
! Rajeswari: a uniform P in every layer, below-cloud washout only, no re-evaporation.
   SUBROUTINE rn_wetscav_ls(chem, num_chem, numgas, moist, num_moist, precr, qlsink, t_phy, alt, dz8w, &
                 dgnuc, dgacc, dgcor, dtstep, config_flags, id,                                     &
                 wd_iodg, wd_iodaj, wd_iodac, wd_csaj, wd_csac,                                     &
                 ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte )
      USE module_configure, only: grid_config_rec_type
      USE module_scalar_tables, only: chem_dname_table
      USE module_state_description, only: p_iodg, p_iodaj, p_iodac, p_csaj, p_csac, p_qc, p_qi,  &
                                          p_nu0, p_ac0, p_corn, p_antha, p_seas, p_soila,         &
                                          param_first_scalar
      IMPLICIT NONE
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
      INTEGER, INTENT(IN) :: num_chem, numgas, num_moist, id
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),  INTENT(INOUT) :: chem
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_moist), INTENT(IN)    :: moist
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),           INTENT(IN)    :: precr, qlsink, t_phy, alt, dz8w, &
                                                                           dgnuc, dgacc, dgcor
      REAL, DIMENSION(ims:ime,jms:jme),                   INTENT(INOUT) :: wd_iodg, wd_iodaj, wd_iodac, &
                                                                           wd_csaj, wd_csac
      REAL, INTENT(IN) :: dtstep

! Rajeswari: species classes: 1 accumulation mass, 2 coarse mass, 3 iodine gas (Henry), 4 Aitken mass,
! Rajeswari:                  5 Aitken number, 6 accumulation number, 7 coarse number
      INTEGER, PARAMETER :: ncls = 7
      REAL*8,  PARAMETER :: pmin = 1.0D-4          ! Rajeswari: smallest rain rate treated [mm h-1]
      REAL*8,  PARAMETER :: rgas_atm = 0.082057D0  ! Rajeswari: gas constant [L atm mol-1 K-1]
      REAL*8,  PARAMETER :: rgas_si  = 8.31446D0   ! Rajeswari: gas constant [J mol-1 K-1]
      REAL*8,  PARAMETER :: a_name = 8.4D-5, b_name = 0.79D0   ! Rajeswari: NAME rain coefficients (Met Office)
      REAL*8,  PARAMETER :: eps_j = 1.0D0, eps_c = 0.5D0      ! Rajeswari: in-cloud efficiency, accumulation / coarse
      REAL*8,  PARAMETER :: cld_thresh = 1.0D-6               ! Rajeswari: cloud threshold qc+qi [kg kg-1]
      REAL*8,  PARAMETER :: kh_i2 = 2.84D0, dhr_i2 = 4300.0D0 ! Rajeswari: I2 Henry constant [M atm-1], dlnH/d(1/T) [K]
                                                              ! Rajeswari: (Eguchi et al. 1973, Sander compilation)
      INTEGER :: nsp, n, nv, i, j, k, ltr
      INTEGER, DIMENSION(num_chem) :: isp, icl, idg
      LOGICAL :: test, reevap, cloudy, anyrain
      REAL*8  :: p_mmh(kts:kte), wcar(num_chem), lamc(ncls)
      REAL*8  :: rmass, fe, rel, lam, f, dc, hcol, base, tk, heff, hdim, qliq, faq, qcl
      REAL*8  :: dgn, dga, dgc
      CHARACTER(LEN=32) :: nm

! Rajeswari: do nothing if the radionuclide species are not part of this chem_opt or the scheme is off
      IF ( config_flags%rn_wdep_opt .LE. 0 ) RETURN
      IF ( p_iodg  .LT. param_first_scalar .OR. p_iodaj .LT. param_first_scalar .OR. &
           p_iodac .LT. param_first_scalar .OR. p_csaj  .LT. param_first_scalar .OR. &
           p_csac  .LT. param_first_scalar ) RETURN

      test   = config_flags%rn_wd_test_rain .GT. 0.
      reevap = .NOT. test

! Rajeswari: list of species to scavenge: the five radionuclide species (idg = output index 1..5), then the host aerosol
      nsp = 0
      CALL add_sp(p_iodaj, 1, 2) ; CALL add_sp(p_csaj, 1, 4)
      CALL add_sp(p_iodac, 2, 3) ; CALL add_sp(p_csac, 2, 5)
      CALL add_sp(p_iodg, 3, 1)
      IF ( config_flags%rn_wdep_host_opt .EQ. 1 ) THEN
         DO nv = numgas+1, num_chem
            IF ( nv.EQ.p_iodg .OR. nv.EQ.p_iodaj .OR. nv.EQ.p_iodac .OR. nv.EQ.p_csaj .OR. nv.EQ.p_csac ) CYCLE
            IF ( nv .EQ. p_nu0 ) THEN
               CALL add_sp(nv, 5, 0)
            ELSE IF ( nv .EQ. p_ac0 ) THEN
               CALL add_sp(nv, 6, 0)
            ELSE IF ( nv .EQ. p_corn ) THEN
               CALL add_sp(nv, 7, 0)
            ELSE IF ( nv.EQ.p_antha .OR. nv.EQ.p_seas .OR. nv.EQ.p_soila ) THEN
               CALL add_sp(nv, 2, 0)
            ELSE
               nm  = chem_dname_table(id, nv)
               ltr = LEN_TRIM(nm)
               IF ( ltr .GT. 0 ) THEN
                  IF ( nm(ltr:ltr) .EQ. 'j' .OR. nm(ltr:ltr) .EQ. 'J' ) CALL add_sp(nv, 1, 0)
                  IF ( nm(ltr:ltr) .EQ. 'i' .OR. nm(ltr:ltr) .EQ. 'I' ) CALL add_sp(nv, 4, 0)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      DO j = jts, jte
      DO i = its, ite
! Rajeswari: rain rate in each layer [mm h-1]
         anyrain = .FALSE.
         hcol    = 0.0D0
         DO k = kts, kte
            IF ( test ) THEN
               p_mmh(k) = DBLE(config_flags%rn_wd_test_rain)
            ELSE
               p_mmh(k) = MAX( 0.0D0, DBLE(precr(i,k,j))*3600.0D0 )
            ENDIF
            IF ( p_mmh(k) .GT. pmin ) THEN
               anyrain = .TRUE.
               hcol    = hcol + DBLE(dz8w(i,k,j))
            ENDIF
            IF ( .NOT. test .AND. config_flags%rn_wdep_opt .EQ. 1 .AND. qlsink(i,k,j) .GT. 0. ) anyrain = .TRUE.
         ENDDO
         IF ( .NOT. anyrain ) CYCLE

         wcar(1:nsp) = 0.0D0
         DO k = kte, kts, -1
            rmass = DBLE(dz8w(i,k,j)) / DBLE(alt(i,k,j))     ! Rajeswari: dry-air mass of the layer [kg m-2]
! Rajeswari: 1. re-evaporation of falling rain between level k+1 and k
            IF ( reevap .AND. k .LT. kte ) THEN
               IF ( precr(i,k+1,j) .GT. 0. .AND. precr(i,k,j) .LT. precr(i,k+1,j) ) THEN
                  fe = DBLE(precr(i,k+1,j) - precr(i,k,j)) / DBLE(precr(i,k+1,j))
                  fe = MIN( 1.0D0, MAX( 0.0D0, fe ) )
                  DO n = 1, nsp
                     IF ( wcar(n) .LE. 0.0D0 ) CYCLE
                     rel     = fe*wcar(n)
                     wcar(n) = wcar(n) - rel
                     chem(i,k,j,isp(n)) = REAL( DBLE(chem(i,k,j,isp(n))) + rel/rmass )
                  ENDDO
               ENDIF
            ENDIF
! Rajeswari: 2. removal rates of the classes in this layer
            cloudy = .FALSE.
            qcl    = 0.0D0
            IF ( p_qc .GE. param_first_scalar ) qcl = MAX( 0.0D0, DBLE(moist(i,k,j,p_qc)) )
            IF ( .NOT. test ) THEN
               IF ( p_qi .GE. param_first_scalar ) THEN
                  cloudy = qcl + MAX( 0.0D0, DBLE(moist(i,k,j,p_qi)) ) .GT. cld_thresh
               ELSE
                  cloudy = qcl .GT. cld_thresh
               ENDIF
            ENDIF
            lamc(:) = 0.0D0
            tk = DBLE(t_phy(i,k,j))
            IF ( cloudy ) THEN
               IF ( config_flags%rn_wdep_opt .EQ. 1 ) THEN
                  base = MAX( 0.0D0, DBLE(qlsink(i,k,j)) )
               ELSE IF ( p_mmh(k) .GT. pmin ) THEN
                  base = a_name * p_mmh(k)**b_name
               ELSE
                  base = 0.0D0
               ENDIF
               lamc(1) = eps_j * base
               lamc(6) = lamc(1)
               lamc(2) = eps_c * base
               lamc(7) = lamc(2)
! Rajeswari: iodine gas: dissolved fraction in cloud water, f = H L / (H L + 1/(R T)), Henry's law with
! Rajeswari: H [mol m-3 Pa-1] = K_H [M atm-1] * 1000/101325 and L = liquid water volume fraction [m3 m-3]
               heff = kh_i2 * EXP( dhr_i2*(1.0D0/tk - 1.0D0/298.0D0) ) * 9.86923D-3
               qliq = qcl / DBLE(alt(i,k,j)) / 1000.0D0               ! Rajeswari: kg m-3 / 1000 kg m-3 = m3 m-3
               faq  = heff*qliq / ( heff*qliq + 1.0D0/(rgas_si*tk) )
               lamc(3) = faq * base
            ELSE IF ( p_mmh(k) .GT. pmin ) THEN
               IF ( config_flags%rn_wdep_opt .EQ. 1 ) THEN
                  dgn = DBLE(dgnuc(i,k,j)) ; IF ( dgn .LE. 0.0D0 ) dgn = 0.02D-6   ! Rajeswari: fallback before the
                  dga = DBLE(dgacc(i,k,j)) ; IF ( dga .LE. 0.0D0 ) dga = 0.15D-6   ! Rajeswari: first aerosol step
                  dgc = DBLE(dgcor(i,k,j)) ; IF ( dgc .LE. 0.0D0 ) dgc = 1.0D-6
                  lamc(1) = rn_lam_mode(dga*1.0D6, 2.0D0, p_mmh(k), 3)
                  lamc(2) = rn_lam_mode(dgc*1.0D6, 2.5D0, p_mmh(k), 3)
                  IF ( config_flags%rn_wdep_host_opt .EQ. 1 ) THEN
                     lamc(4) = rn_lam_mode(dgn*1.0D6, 1.7D0, p_mmh(k), 3)
                     lamc(5) = rn_lam_mode(dgn*1.0D6, 1.7D0, p_mmh(k), 0)
                     lamc(6) = rn_lam_mode(dga*1.0D6, 2.0D0, p_mmh(k), 0)
                     lamc(7) = rn_lam_mode(dgc*1.0D6, 2.5D0, p_mmh(k), 0)
                  ENDIF
               ELSE
                  lamc(1) = a_name * p_mmh(k)**b_name
                  lamc(2) = lamc(1) ; lamc(4) = lamc(1) ; lamc(5) = lamc(1) ; lamc(6) = lamc(1) ; lamc(7) = lamc(1)
               ENDIF
               lamc(1:2) = lamc(1:2) * DBLE(config_flags%rn_wdep_scale)
               lamc(4:7) = lamc(4:7) * DBLE(config_flags%rn_wdep_scale)
! Rajeswari: iodine gas below cloud: Henry-limited washout, L = H' P / H_col  (P in m s-1)
               heff = kh_i2 * EXP( dhr_i2*(1.0D0/tk - 1.0D0/298.0D0) )
               hdim = heff * rgas_atm * tk                                   ! Rajeswari: dimensionless H'
               IF ( hcol .GT. 0.0D0 ) lamc(3) = hdim * (p_mmh(k)/3.6D6) / hcol
            ENDIF
! Rajeswari: 3. apply the removal and add the removed activity to the falling rain
            DO n = 1, nsp
               lam = lamc(icl(n))
               IF ( lam .LE. 0.0D0 ) CYCLE
               f  = 1.0D0 - EXP( -lam*DBLE(dtstep) )
               dc = f*DBLE(chem(i,k,j,isp(n)))
               IF ( dc .LE. 0.0D0 ) CYCLE
               chem(i,k,j,isp(n)) = REAL( DBLE(chem(i,k,j,isp(n))) - dc )
               wcar(n) = wcar(n) + dc*rmass
            ENDDO
         ENDDO
! Rajeswari: activity reaching the ground = wet deposition [Bq m-2] (radionuclide species only)
         DO n = 1, nsp
            SELECT CASE ( idg(n) )
            CASE (1) ; wd_iodg(i,j)  = wd_iodg(i,j)  + REAL(wcar(n))
            CASE (2) ; wd_iodaj(i,j) = wd_iodaj(i,j) + REAL(wcar(n))
            CASE (3) ; wd_iodac(i,j) = wd_iodac(i,j) + REAL(wcar(n))
            CASE (4) ; wd_csaj(i,j)  = wd_csaj(i,j)  + REAL(wcar(n))
            CASE (5) ; wd_csac(i,j)  = wd_csac(i,j)  + REAL(wcar(n))
            END SELECT
         ENDDO
      ENDDO
      ENDDO

   CONTAINS

      SUBROUTINE add_sp(nvx, clx, dgx)
         INTEGER, INTENT(IN) :: nvx, clx, dgx
         IF ( nvx .LT. param_first_scalar ) RETURN
         nsp = nsp + 1
         isp(nsp) = nvx ; icl(nsp) = clx ; idg(nsp) = dgx
      END SUBROUTINE add_sp

   END SUBROUTINE rn_wetscav_ls

! Rajeswari: below-cloud scavenging coefficient by rain [s-1], Wang, Zhang and Moran (2014), GMD 7, 799-819,
! Rajeswari: Eqs. (5)-(7) with the rain constants of Table 8: L = A(d) R**B(d), d = wet diameter [um], R [mm h-1]
   REAL*8 FUNCTION rn_lam_wang(d_um, r_mmh)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: d_um, r_mmh
      REAL*8 :: d, r, x, la, bb
      d = MIN( 100.0D0, MAX( 1.0D-3, d_um ) )
      r = MIN( 100.0D0, MAX( 1.0D-4, r_mmh ) )
      x = LOG10(d)
      IF ( d .LE. 2.0D0 ) THEN
         la = -6.2609D0 + x*(0.6820D0 + x*(0.8676D0 + x*0.1282D0))
         bb =  0.7230D0 + 0.0303D0*x
      ELSE
         la = -14.707D0 + x*(51.043D0 + x*(-97.306D0 + x*(97.946D0 + x*(-53.923D0 + x*(15.311D0 + x*(-1.7510D0))))))
         bb = -0.6492D0 + x*(9.3483D0 + x*(-21.929D0 + x*(25.317D0 + x*(-15.395D0 + x*(4.7242D0 + x*(-0.5766D0))))))
      ENDIF
      rn_lam_wang = 10.0D0**la * r**bb
   END FUNCTION rn_lam_wang

! Rajeswari: rn_lam_wang averaged over a lognormal mode (number geometric mean diameter dg_um [um], width sg);
! Rajeswari: moment = 0 number-weighted, 3 mass-weighted (volume); 40-point midpoint rule over +-4.5 standard
! Rajeswari: deviations in ln d with Gaussian weights (error < 0.3 % against a dense integration)
   REAL*8 FUNCTION rn_lam_mode(dg_um, sg, r_mmh, moment)
      IMPLICIT NONE
      REAL*8,  INTENT(IN) :: dg_um, sg, r_mmh
      INTEGER, INTENT(IN) :: moment
      INTEGER, PARAMETER :: nq = 40
      REAL*8,  PARAMETER :: zmax = 4.5D0
      REAL*8 :: ls, mu, z, w, sw, s
      INTEGER :: q
      ls = LOG(sg)
      mu = LOG(dg_um) + DBLE(moment)*ls*ls
      s  = 0.0D0 ; sw = 0.0D0
      DO q = 1, nq
         z  = -zmax + (DBLE(q) - 0.5D0)*2.0D0*zmax/DBLE(nq)
         w  = EXP(-0.5D0*z*z)
         s  = s  + w*rn_lam_wang( EXP(mu + ls*z), r_mmh )
         sw = sw + w
      ENDDO
      rn_lam_mode = s/sw
   END FUNCTION rn_lam_mode
! ===== added by Rajeswari : END =====

! ===== added by Rajeswari (radionuclide species, Stage H: uptake of iodine gas on aerosol particles) : START =====
! Rajeswari: moves I-131 activity from the gas (iodg) to the particle tags on chemistry steps:
! Rajeswari:   rn_upt_opt = 1: coarse mode only -> iodac;  2: coarse -> iodac and accumulation -> iodaj
! Rajeswari: first-order rate per mode, k_m = N_m < k1(d) >, Fuchs-Sutugin rate of one particle
! Rajeswari:   k1(d) = 2 pi d D_g F(Kn, gamma),  F = (1+Kn) / (1 + (4/(3 gamma) + 0.377) Kn + 4/(3 gamma) Kn**2)
! Rajeswari:   Kn = 2 lambda / d,  lambda = 3 D_g / vbar  (Seinfeld and Pandis 2016, ch. 12; Jacob 2000)
! Rajeswari: kinetic limit k1 -> gamma vbar pi d**2 / 4; diffusion limit (gamma = 1, large d) k1 -> 2 pi d D_g
! Rajeswari: transfer over the chemistry step: dA = A_g [1 - exp(-k dt)], k = k_j + k_c, split k_j/k and k_c/k
! Rajeswari: host aerosol not changed (radionuclide mass is negligible); same gamma for both modes
   SUBROUTINE rn_uptake(chem, num_chem, t_phy, p_phy, alt, dz8w, dgacc, dgcor, dtstep, config_flags, &
                 kupt, upt_iodj, upt_iodc,                                                        &
                 ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte )
      USE module_configure, only: grid_config_rec_type
      USE module_state_description, only: p_iodg, p_iodaj, p_iodac, p_ac0, p_corn, param_first_scalar
      IMPLICIT NONE
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
      INTEGER, INTENT(IN) :: num_chem
      INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme, its,ite, jts,jte, kts,kte
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),          INTENT(IN)    :: t_phy, p_phy, alt, dz8w, dgacc, dgcor
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme),          INTENT(INOUT) :: kupt
      REAL, DIMENSION(ims:ime,jms:jme),                  INTENT(INOUT) :: upt_iodj, upt_iodc
      REAL, INTENT(IN) :: dtstep

      REAL*8, PARAMETER :: pi = 3.14159265358979D0
      REAL*8, PARAMETER :: rgas = 8.31446D0          ! Rajeswari: gas constant [J mol-1 K-1]
      REAL*8, PARAMETER :: mw_i2 = 0.25381D0         ! Rajeswari: molar mass of I2 [kg mol-1]
      REAL*8, PARAMETER :: sg_acc = 2.0D0, sg_cor = 2.5D0   ! Rajeswari: MADE geometric standard deviations
      INTEGER :: i, k, j, opt
      REAL*8  :: gam, tk, pa, vbar, dgas, nn, kj, kc, kt, x, f, da, rmass

      kupt(its:ite,kts:kte,jts:jte) = 0.
      opt = config_flags%rn_upt_opt
      gam = DBLE(config_flags%rn_upt_gamma)
      IF ( opt .LE. 0 .OR. gam .LE. 0.0D0 ) RETURN
      IF ( p_iodg .LT. param_first_scalar .OR. p_iodaj .LT. param_first_scalar .OR. &
           p_iodac .LT. param_first_scalar .OR. p_corn  .LT. param_first_scalar .OR. &
           p_ac0   .LT. param_first_scalar ) RETURN
      gam = MIN( 1.0D0, gam )

      DO j = jts, jte
      DO k = kts, kte
      DO i = its, ite
         tk   = DBLE(t_phy(i,k,j))
         pa   = DBLE(p_phy(i,k,j))
         IF ( tk .LE. 0.0D0 .OR. pa .LE. 0.0D0 ) CYCLE
         vbar = SQRT( 8.0D0*rgas*tk/(pi*mw_i2) )                                  ! Rajeswari: [m s-1]
! Rajeswari: I2 diffusivity [m2 s-1], rn_iodg_dvj [cm2 s-1] at 298 K and 1013.25 hPa, scaled as T**1.75 / p
         dgas = DBLE(config_flags%rn_iodg_dvj)*1.0D-4 * (tk/298.0D0)**1.75D0 * (101325.0D0/pa)
! Rajeswari: coarse mode (number [# kg-1] / alt [m3 kg-1] = # m-3; dry diameter [m])
         kc = 0.0D0
         nn = DBLE(chem(i,k,j,p_corn)) / DBLE(alt(i,k,j))
         IF ( nn .GT. 0.0D0 .AND. dgcor(i,k,j) .GT. 0. ) &
            kc = rn_kfs_mode(nn, DBLE(dgcor(i,k,j)), sg_cor, dgas, vbar, gam)
! Rajeswari: accumulation mode (wet diameter), option 2 only
         kj = 0.0D0
         IF ( opt .EQ. 2 ) THEN
            nn = DBLE(chem(i,k,j,p_ac0)) / DBLE(alt(i,k,j))
            IF ( nn .GT. 0.0D0 .AND. dgacc(i,k,j) .GT. 0. ) &
               kj = rn_kfs_mode(nn, DBLE(dgacc(i,k,j)), sg_acc, dgas, vbar, gam)
         ENDIF
         kt = kj + kc
         kupt(i,k,j) = REAL(kt)
         IF ( kt .LE. 0.0D0 ) CYCLE
! Rajeswari: fraction moved in this step; series form keeps full accuracy for very small k dt
         x = kt*DBLE(dtstep)
         IF ( x .LT. 1.0D-5 ) THEN
            f = x*(1.0D0 - 0.5D0*x)
         ELSE
            f = 1.0D0 - EXP(-x)
         ENDIF
         da = f*DBLE(chem(i,k,j,p_iodg))
         IF ( da .LE. 0.0D0 ) CYCLE
         chem(i,k,j,p_iodg)  = REAL( DBLE(chem(i,k,j,p_iodg))  - da )
         chem(i,k,j,p_iodaj) = REAL( DBLE(chem(i,k,j,p_iodaj)) + da*kj/kt )
         chem(i,k,j,p_iodac) = REAL( DBLE(chem(i,k,j,p_iodac)) + da*kc/kt )
! Rajeswari: column transfer diagnostics [Bq m-2]; dry-air mass of the layer = dz / alt [kg m-2]
         rmass = DBLE(dz8w(i,k,j)) / DBLE(alt(i,k,j))
         upt_iodj(i,j) = upt_iodj(i,j) + REAL( da*kj/kt*rmass )
         upt_iodc(i,j) = upt_iodc(i,j) + REAL( da*kc/kt*rmass )
      ENDDO
      ENDDO
      ENDDO
   END SUBROUTINE rn_uptake

! Rajeswari: uptake rate [s-1] of a gas on a lognormal mode: N [# m-3], number median diameter dg [m], width sg,
! Rajeswari: gas diffusivity dgas [m2 s-1], mean molecular speed vbar [m s-1], uptake coefficient gam.
! Rajeswari: k = N < k1(d) >, written as k = N dg**2 exp(2 ln**2 sg) < k1(d)/d**2 >_s, the average taken over the
! Rajeswari: surface-weighted distribution (lognormal, median ln dg + 2 ln**2 sg); 40-point midpoint rule over
! Rajeswari: +-4.5 standard deviations, as rn_lam_mode. Exact in the kinetic limit, where k1/d**2 is constant.
   REAL*8 FUNCTION rn_kfs_mode(nn, dg, sg, dgas, vbar, gam)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: nn, dg, sg, dgas, vbar, gam
      REAL*8,  PARAMETER :: pi = 3.14159265358979D0
      INTEGER, PARAMETER :: nq = 40
      REAL*8,  PARAMETER :: zmax = 4.5D0
      REAL*8  :: ls, mu, z, w, sw, s, d, lam, kn, fs
      INTEGER :: q
      rn_kfs_mode = 0.0D0
      IF ( nn .LE. 0.0D0 .OR. dg .LE. 0.0D0 .OR. gam .LE. 0.0D0 ) RETURN
      ls  = LOG( MAX(sg, 1.0D0 + 1.0D-6) )
      mu  = LOG(dg) + 2.0D0*ls*ls
      lam = 3.0D0*dgas/vbar
      s = 0.0D0 ; sw = 0.0D0
      DO q = 1, nq
         z  = -zmax + (DBLE(q) - 0.5D0)*2.0D0*zmax/DBLE(nq)
         w  = EXP(-0.5D0*z*z)
         d  = EXP(mu + ls*z)
         kn = 2.0D0*lam/d
         fs = (1.0D0 + kn) / (1.0D0 + (4.0D0/(3.0D0*gam) + 0.377D0)*kn + 4.0D0/(3.0D0*gam)*kn*kn)
         s  = s  + w*(2.0D0*pi*d*dgas*fs)/(d*d)
         sw = sw + w
      ENDDO
      rn_kfs_mode = nn * dg*dg*EXP(2.0D0*ls*ls) * s/sw
   END FUNCTION rn_kfs_mode
! ===== added by Rajeswari : END =====
END MODULE module_input_tracer
