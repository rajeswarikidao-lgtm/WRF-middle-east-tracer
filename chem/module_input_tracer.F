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
END MODULE module_input_tracer
