C=======================================================================
C                THIS IS A SIMPLIFIED SIB MODEL                              
C=======================================================================

CM THIS CODE OF YONGKANG XUES HAS BEEN MODIFIED TO RUN AS A SUBROUTINE
CM INSIDE AN ENKF RADIOMETRIC DATA ASSIMILATION FOR ESTIMATION OF SWE.
CM THE VEGIN1, VEGIN2 AND CNTROLS SUBROUTINES NO LONGER READ DATA FROM 
CM FILES, BUT EXTRACT IT FROM THESE ARRAYS; SOME VALUES THAT ARE 
CM STATIC HAVE BEEN HARDCODED INTO THESE SUBROUTINES.  THEREFORE, THE 
CM run.ctl FILE WHICH USED TO CONTAIN THE NAMES OF ALL OF THESE FILES
CM IS NOW OBSELETE.  FURTHERMORE, THE ALBEDO ADJUSTMENT BASED ON THE 
CM xadj VARIABLE IS ALSO OBSELETE AND NOT IMPLEMENTED HERE.  THE CODE
CM TO COMPUTE rnoffs AND rnoffb WAS ALSO REMOVED, AS IT WAS NOT BEING
CM USED IN THE PRESENT APPLICATION.

CM CONSTANT VEGETATION DATA IS CONTAINED IN THE VEGIN1ARG VARIABLE, 
CM VEGETATION DATA WHICH CHANGES EACH MONTH IN THE VEGIN2ARG VARIABLE.
CM N_U - NUMBER OF FORCING DATA AT EACH TIMESTEP
CM N_STEPS - NUMBER OF STEPS AT WHICH TO RUN THE MODEL
CM UARG - FORCING DATA
CM ALPHA - PARAMETERS DICTATING THE GROWTH OF THE SNOW GRAINS
CM T0 - INITIAL TIME DATA
CM X0 - INITIAL AUXILIARY DATA NEEDED TO INITIALIZE THE MODEL
CM Y0 - INITIAL SNOW STATE VARIABLES
CM Z - THE LATITUDE AND LONGITUDE OF THIS RUN
CM YOUT - THE OUTPUT ARRAY CONTAINING THE SNOW DATA

      subroutine ssib3(IVEG_TYPE,n_u,n_steps,uarg,alpha,t0, 
     &  x0,y0,z_in,yout,PIXEL,n_y,n_x,xout,zout,tend,replicate,rank,
     &  meas,n_a,f_veg,vcov_dat) ! ZenithAngle,Energy_out,Effective_PPT)

      COMMON /XRIB/RIB,temprib
      character*7 run
      INCLUDE 'comsib.in'
      include 'snow4.in'

cm  these lines to define sizes for arrays now used to initialize 
cm    ssib and call it as a subroutine
      integer n_u,n_steps,PIXEL,n_y,n_x,replicate,rank,meas,counter
      real uarg(n_u,n_steps)
      real alpha(n_a),t0(5),x0(12),y0(n_y),z_in(2),tend(5)
      real vcov_dat
cm    the dimension of f_veg is consistent with mvnrnd
      real yout(n_y,n_steps),xout(n_x,n_steps),zout(n_steps)
      real sweold,swenew,sag1,sag2
      real ALBEDOOUT(n_steps,2,2)
      real ZenithAngle(n_steps)
	  real Energy_out(8,n_steps)
	  real Effective_PPT(1,n_steps)
	  real f_veg(1,1)
c         real swenew_Vec(n_steps)
c         real sweold_Vec(n_steps)
c         real sag1_Vec(n_steps)
c         real sag2_Vec(n_steps)

cl ================ declaration of MEMLS parameters ====================
cl declaration for parameters in MEMLS3 and MEMLS4, NOV 16 2012
cl 4-layer memls parameters
         integer :: n_rlz4,n_rz4,n_yz4,n_ypz4,n_yrz4,n_zz4,n_zpz4,
     &  n_zrz4,n_cz4,n_xz4,n_az4,n_uz4,meas_switch4,n_bdrf4

      integer :: ierr4
      integer, dimension(1) :: iveg_type4
      integer, dimension(:), allocatable :: ctrl4
      real, dimension(1) :: tbs,tbs4
      real, dimension(:),allocatable :: freq4,theta_m4
      real, dimension(:,:),allocatable :: y4memls,xz4
      real, dimension(:,:,:),allocatable :: uz4
      real, dimension(1,1) :: f_veg_m4,albedo_m4
      real, dimension(1,2) :: bdrf4
      real, dimension(2,1) :: zz4,tbraw,tbraw4
      integer :: month_m4

cl 3-layer memls parameters
         integer :: n_rlz3,n_rz3,n_yz3,n_ypz3,n_yrz3,n_zz3,n_zpz3,
     &  n_zrz3,n_cz3,n_xz3,n_az3,n_uz3,meas_switch3,n_bdrf3

      integer :: ierr3
      integer, dimension(1) :: iveg_type3
      integer, dimension(:), allocatable :: ctrl3
      real, dimension(1) :: tbs3,CPEX
      real, dimension(:),allocatable :: freq3,theta_m3
      real, dimension(:,:),allocatable :: y3memls,xz3
      real, dimension(:,:,:),allocatable :: uz3
      real, dimension(1,1) :: f_veg_m3,albedo_m3
      real, dimension(1,2) :: bdrf3
      real, dimension(2,1) :: zz3,tbraw3
      integer :: month_m3
cl ================ MEMLS declaration done ====================



CM THE ctlpa AND nroot VARIABLE WILL BE HARDCODED IN (read from file before)
c     control stomatal resistance;
c     final stomatal resistance=ctlpa * stomatal resistance
      ctlpa=1.0
c     control root depth, nroot=1: root depth has no control: nroot not
c     =1: root depth is controled by rootp, which are read in vegin2.
      nroot=1
      f_veg(1,1)=0.0

cm  print out all inputs...
c      print *, iveg_type,n_u,n_steps,alpha,t0,x0,y0,z_in,PIXEL,
c     &  n_y,n_x,tend,replicate,rank,meas,n_a,f_veg,vcov_dat

CM SOME OF THESE FORMAT STATEMENTS WERE OBSELTE AND SO DELETED
    2 format(a12)
    3 format(a40)

CM INITIALIZE SOME VARIABLES

cm  put in a check that n_x is 13, not 12, since i used 12 for so long...
      if(n_x.ne.13)then
        print *, 'error! n_x=',n_x,'; n_x should be 13...'
      end if   

      CALL CNTROLS(MDAY,alpha,t0,x0,y0,z_in,pixel,n_a) 
      CALL CONSTS
  
      CUMERROR = 0
      ERRORCUM = 0
      xqlast1 = 0.
      xqlast2 = 0.
      xqlast3 = 0.
      xlastwet =0.
      xlastswe = 0.
      xlastcap1 = 0.
      
      sweold = 0.0
      swenew = 0.0


c     four lines changed by mike to allow time to be read in and used
c      iyr=1979
      iyr=YEAR
      mthst = MONTH
      ndyst = MDAY
      nhrst = TIME
      mthend = 240
      ndyend = 31
      nhrend = 24 
      nobs = 0
      XPG1 = 0.
      XPG2 = 0.
c these are for my own time keeping which replaces time in forcing file
      nyy=iyr
      nmm=mthst
      ndd=ndyst
      nhh=nhrst
c
C ** ASSIGN VEGETATION PARAMETERS FROM INPUT ARRAYS
      call vegin
      CALL VEGIN1(IVEG_TYPE) 
      CALL VEGIN2(IVEG_TYPE,nmm,f_veg,vcov_dat) 
  
c
cjyj--control parameter of snow model
      MDLSNO=0
      SD_CR=0.05
cm this part is changed by mike to allow snow density to be read in
cm      snden=3.75
      isnow=1
cm this park changed by mike
      SDEP=CAPAC(2)*snden
      call getinput 
      swe=CAPAC(2)    
cm this part is changed by mike, to allow dz(nd) to be read in

      snowdepth=dz(1)+dz(2)+dz(3)      
      
      if(mdlsno.eq.0.and.snowdepth.gt.sd_cr) then
        isnow=0
        call layern(tgs,0,nmm,ndd,nhh)
      endif
CS                                     10/13/98 
C ** CREATE THE COEFFICIENTS FOR SURFACE ALBEDO         
      XHC = 0.                                                          
      XHG = 0.                                                          
      XCI = 0.                                                         
      XCT = 0.                                                          
      XGI = 0.                                                         
      XGT = 0.                                                          
      XGS = 0.                                                          
      NNX = 0                                                           

      isnowhr=0
      iday=0

cm ********************************************
cm   START OF TIME LOOP
cm ********************************************
      DO ICTRL=1,N_STEPS 
cm    extract forcing data from uarg array

!      print *, ictrl

      icrash=1

      swdown=uarg(1,ICTRL)
      rainf=uarg(2,ICTRL)
      obsnow=uarg(3,ICTRL)
      dirdown=uarg(4,ICTRL)
      psur=uarg(5,ICTRL)
      tm=uarg(6,ICTRL)
      qair=uarg(7,ICTRL)
      windn=uarg(8,ICTRL)
      winde=uarg(9,ICTRL)
      afac=1.0


cl    dongyue constrain the soil tempc
      if(tgs.lt.273.16)then
          tgs=273.16
      end if

c     if(ictrl.eq.icrash) print *, 'U=',uarg(1:9,ictrl)

cm    check on whether it should be rain or snow...
      if(tm.gt.273.16.and.obsnow.gt.0.0)then
        rainf=obsnow
        obsnow=0.0
      end if
      
      if(tm.le.273.16.and.rainf.gt.0.0)then
        obsnow=rainf
        rainf=0.0
      end if

      if(obsnow.gt.0.0)then
        isnowhr=isnowhr+1
      end if

cm    conversions of forcing data
      UM=(windn**2+winde**2)**0.5 
      TPREC=(rainf+obsnow)*3600
cm    here surface pressure is converted from Pa to mbar      
      psur = psur / 100.0
      EM = (psur*qair)/0.6220

cm    update of time variables
      idkmax = 31
      lstmth = nmm
      do 3456, ikahan=0,19
      if ((nmm.eq.(4+12*ikahan)).or.
     & (nmm.eq.(9+12*ikahan)).or.
     & (nmm.eq.(11+12*ikahan)).or.
     & (nmm.eq.(6+12*ikahan))) then
      idkmax = 30
      endif
 3456 continue
      do 3457, ikahan=0,19
      if (nmm.eq.(2+12*ikahan)) then
      idkmax=28
      if ((nyy.eq.1980).or.(nyy.eq.1984).or.(nyy.eq.1988).
     & or.(nyy.eq.1992).or.(nyy.eq.1996).or.
     & (nyy.eq.2000).or.(nyy.eq.2004)) then
      idkmax = 29
      endif
      endif
 3457 continue
      if (nhh.eq.24) then
      nhh = 0
      ndd = ndd + 1
      if (ndd.eq.(idkmax+1)) then
      ndd = 1
      nmm = nmm + 1
      do 3458, ikahan=1,20
      if (nmm.eq.(1+12*ikahan)) then
      nmm = 1
      nyy = nyy + 1
      endif
 3458 continue
      endif
      endif
      nhh=nhh+1

       qsoil=0.
       wfsoil=0.
       solsoil=0.
       snroff=0.
       if (swdown.lt.0) swdown = 0.0

       if (tprec.le.0.) tprec = 0.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (lstmth.ne.nmm) then
        if ((nmm.eq.13).or.(nmm.eq.25).or.(nmm.eq.37).or.
     & (nmm.eq.49).or.(nmm.eq.61).or.(nmm.eq.73).or.
     & (nmm.eq.85).or.(nmm.eq.97).or.(nmm.eq.109).or.
     & (nmm.eq.121).or.(nmm.eq.133).or.(nmm.eq.145).or.
     & (nmm.eq.157).or.(nmm.eq.169).or.(nmm.eq.181).or.
     & (nmm.eq.193).or.(nmm.eq.205).or.(nmm.eq.217).or.
     & (nmm.eq.229)) rewind(2)
      CALL VEGIN2(IVEG_TYPE,nmm,f_veg,vcov_dat)
      endif

c      if (ictrl.eq.icrash) print *, 'ictrl=',ictrl,'after timekeeping'

cm     observation counter updated
       nobs = nobs + 1

      IF (MDLSNO.eq.0) THEN 
         iptype=2
cm       Mike is changing this to the SNTHERM cutoff temperature for snow
         If (TM.ge.273.16) iptype=1
         If (TPREC.le.0.) iptype=0
      END IF        
c
C ** START THE MAIN PROGRAM       
c
c      if (ictrl.eq.icrash) print *, 'ictrl=',ictrl,'before main prog.'

      RADN(3,2) = dirdown 
      UM = AMAX1(UM,0.25)                                               
      SWDOWN = AMAX1(SWDOWN,0.1)                                        
      IHOUR = NHH

c ** calculate solar zenith angle    
      CALL RADC2(SUNANG)
     
      ! Manu - is changing this so that the Zentih Angle is the same computed
      ! as a fucntion of slope and aspect from the disaggregation file
      SUNANG = ZenithAngle(ICTRL)
      SUNANG = AMAX1(0.01,SUNANG)  
c                                 
C ** CALCULATE THE CLOUD COVER USING AN EMPIRICAL EQUATION             
      CLOUD = (1160.*SUNANG - SWDOWN) / (963. * SUNANG)                 
      CLOUD = AMAX1(CLOUD,0.)                                           
      CLOUD = AMIN1(CLOUD,1.)                                           
      CLOUD = AMAX1(0.58,CLOUD)                                         

C ** SEPERATE THE SHORT WAVE RADIATION INTO DIFFERENT SPECTRUM         
      DIFRAT = 0.0604 / ( SUNANG-0.0223 ) + 0.0683                      
      IF ( DIFRAT .LT. 0. ) DIFRAT = 0.                                 
      IF ( DIFRAT .GT. 1. ) DIFRAT = 1.                 
      DIFRAT = DIFRAT + ( 1. - DIFRAT ) * CLOUD                         
      VNRAT = ( 580. - CLOUD*464. ) / ( ( 580. - CLOUD*499. )           
     &        + ( 580. - CLOUD*464. ) )                                 
C                           
      FRAC(1,1) = (1.-DIFRAT)*VNRAT
      FRAC(1,2) = DIFRAT*VNRAT
      FRAC(2,1) = (1.-DIFRAT)*(1.-VNRAT)
      FRAC(2,2) = DIFRAT*(1.-VNRAT)
      RADN(1,1) = FRAC(1,1)*SWDOWN
      RADN(1,2) = FRAC(1,2)*SWDOWN
      RADN(2,1) = FRAC(2,1)*SWDOWN
      RADN(2,2) = FRAC(2,2)*SWDOWN

Cm  initialize snow vars -mike
      IF  (MDLSNO.eq.0.and.ISNOW.EQ.0) THEN
         tsoil=TGS
         TGS=tssno(n)
         CAPAC(2)=swe
         SDEP=snowdepth
         ssss=sdep
      else
         ssss=capac(2)*snden
 
      END IF

c      if (ictrl.eq.icrash) print *, 'before radab call'

cm MAIN RADIATION CALLS
      CALL RADAB (ISNOW,MDLSNO,afac,ictrl,PIXEL,swenew,sweold,
     &     sag1,sag2)
  
c Manu set sweold = swe, this is needed in CLM albedo model
           sweold = swe

c      if (ictrl.eq.icrash) print *, 'after radab call'

      CALL ROOT1
      CALL STOMA1                                                       
      RSTUN = RST(1)                                 

C                                                                       
c     PPC = TPREC
c     PPL = TPREC-PPC
      PPL = TPREC                                                       
      PPC = TPREC-PPL 
 
C                                                                       
C ** WATER BALANCE CHES.   
c
      TOTWB = WWW(1) * POROS * ZDEPTH(1)                                
     &      + WWW(2) * POROS * ZDEPTH(2)                                
     &      + WWW(3) * POROS * ZDEPTH(3)                                
     &      + CAPAC(1) + CAPAC(2)                                       
      thelastmoist = www(1)*poros*zdepth(1) + www(2)*poros*zdepth(2)
     & + www(3)*poros*zdepth(3)
      thelastcapac1 = CAPAC(1)
      thelastcapac2 = CAPAC(2)

cm      if(meas.eq.2) then
cm        print *, 'www=',www,'capac=',capac,'totwb=',totwb
cm      end if

C                                                                       
c ** interception and runoff calculations
CS   Chanfe INTERC to INTERCS (****,***)                   ON 10/13/98 

c      if (ictrl.eq.icrash) print *, 'before intercs calls'

      CALL INTERCS (ISNOW,p0,CSOIL,dzsoil,CHISL)

c      if (ictrl.eq.icrash) print *, 'after intercs calls; ',
c     &  'isnow=',isnow     

      IF (MDLSNO.eq.0.and.ISNOW.EQ.0) THEN
         prcp=p0
         tkair=TM
c      if (ictrl.eq.icrash) print *, 'before getmet calls'
         CALL getmet (iptype,UM,nmm,nhh,ndd,ictrl,pixel)
                                                                       
c ** aerodynamic resistance and flux calculations
        solar=0.
        DO 1100 IVEG  = 2, 2
        DO 1100 IWAVE = 1, 2
        DO 1100 IRAD  = 1, 2
        solar=solar+RADFAC(IVEG,IWAVE,IRAD)*RADN(IWAVE,IRAD)
1100    CONTINUE

       CALL snow1st (dtt,TM,solsoil,ISNOW,nmm,ndd,nhh,ictrl,pixel,imike)
c      if (ictrl.eq.icrash) print *, 'before temrs2 calls'

c         if ( (ictrl.ge.5500).and.(ictrl.lt.5502) ) then
c         if  (ictrl.eq.5500) then
c         print*,ictrl,'SW1=',swe,'c(2)',capac(2),'SD',snowdepth,'DENS',snden
c      endif


       CALL TEMRS2 (MDLSNO,ISNOW,CHISL,tsoil,solsoil,meas,
     &   CSOIL,dzsoil,wfsoil,ictrl,pixel,y0,n_y,replicate,rank)


          call old 


c         if ( (ictrl.ge.5500).and.(ictrl.lt.5502) ) then
c      if  (ictrl.eq.5500) then
c         print*,ictrl,'SW2=',swe,'c(2)',capac(2),'SD',snowdepth,'DENS',snden
c      endif

      ELSE 
cm      if there is snowfall, find snowfall density (get_met)and 
cm        incorporate snowfall in the rest of the snowpack (newsnow)

c        if (ictrl.eq.icrash) print *, 'top of if-structure, isnow=1',
c     &    'tprec=',tprec,'iptype=',iptype,'gsize=',gsize

        if(TPREC.gt.0.and.iptype.ne.1)then
          prcp=TPREC/1000. 
          tkair=TM

c          if (ictrl.eq.icrash) print *, 'before getmet'
          call getmet(iptype,UM,nmm,ndd,nhh,ictrl,pixel)
c          if (ictrl.eq.icrash) print *, 'before newsnow'
cm        The 'getmet' subroutine sets 'iptype' to 0 if the depth of
cm        newsnowfall is less than some critical value.
          if(iptype.ne.0)  call newsnow(ISNOW,nmm,ndd,nhh,ictrl,pixel)
c          if (ictrl.eq.icrash) print *, 'after getmet, newsnow calls'
          endif

        if(gsize>0.) call graingrowth(ISNOW)

c        if (ictrl.eq.icrash) print *, 'before temrs1 calls'

        CALL TEMRS1 (MDLSNO,ISNOW,rank,replicate,pixel,ictrl,icrash)
 
c        if (ictrl.eq.icrash) print *, 'after temrs1 calls'

      END IF


cm UPDATE STATE VARIABLES AFTER N-R TEMPERATURE SOLUTION
      CALL UPDAT1 (MDLSNO,ISNOW,wfsoil,swe,snroff)

c      if (ictrl.eq.icrash) print *, 'after updat1 calls'

cm RELAYER SNOWPACK DEPENDING ON CHANGE IN DEPTH
      IF (MDLSNO.eq.0.and.ISNOW.EQ.0) THEN
         CAPAC(2)=swe
         swenew = swe

cm       mike is adding this on jan 31 05 to solve an issue
cm         that comes up when the snow melts completely away from a 
cm         pack with snowdepth > 0.05 m

      snowdepth=dzo(1)+dzo(2)+dzo(3) 

         If (snowdepth.lt.SD_CR) Then
           ISNOW=1
           call LAYER1 (CSOIL,TGS,dzsoil,h,w,snowdepth,
     &            swe,stemp,nd,gsize,gdia)
         Else
 
c          if (ictrl.le.5469) then 
           ISNOW=0
           call modnodenew
c          if (ictrl.eq.5500) print*,'1.MODcalled'
c          endif
   
         End if
      ELSE IF(MDLSNO.eq.0.and.ISNOW.GT.0) THEN
         If (capac(2)*snden.gt.SD_CR) Then
 
 
           swe=CAPAC(2)                           
           snowdepth=capac(2)*snden 
           ISNOW=0
   

   
           CALL LAYERN (TGS,1,nmm,nhh,ndd)
         Else
           ISNOW=1
cm    in order to output the snow, update snowdepth variable


           snowdepth = capac(2)*snden
          swe = capac(2)
   
cm    end of mikes changes
         End if
      END IF

c      if (ictrl.eq.icrash) print *, 'after layer calls'

      call set0

c      if (ictrl.eq.icrash) print *, 'before water balance checks'

      ROFF=ROFF+snroff 

c this next line was added by dr xue
c     rnoffs(nobs) = rnoffs(nobs) + snroff
C  

cm  WATER AND ENERGY BALANCE CHECKS

      ENDWB = WWW(1) * POROS * ZDEPTH(1)                                
     &      + WWW(2) * POROS * ZDEPTH(2)                                
     &      + WWW(3) * POROS * ZDEPTH(3)                                
     &      + CAPAC(1) + CAPAC(2) - PPL/1000. + ETMASS/1000. + ROFF     
      ERROR = TOTWB - ENDWB 

cm      if(meas.eq.2) then 
cm        print *, 'www=',www,'capac=',capac,'ppl=',ppl,'etmass=',etmass,
cm     &    'roff=',roff,'endwb=',endwb,'poros=',poros,'zdepth=',zdepth
cm      end if
      
      CUMERROR = CUMERROR + ERROR
      if ((nmm.eq.240).and.(ndd.eq.31).and.(nhh.eq.24))
     & write(6,*) 'CUMULATIVE ERROR = ',CUMERROR



c      if (ictrl.eq.49) print *, 'before energy balance checks, dtt=',
c     &  dtt

      CBAL = RADT(1) - CHF - (ECT+HC+ECI)/DTT                           
      GBAL = RADT(2) - SHF - (EGT+EGI+HG+EGS)/DTT                       
      ZLHS = RADT(1) + RADT(2) - CHF - SHF                              
      ZRHS = HFLUX + (ECT + ECI + EGT + EGI + EGS)/DTT
      RORRE = ZLHS-ZRHS
      ERRORCUM = ERRORCUM + RORRE
  
  
c Manu Extract the energy balance at the ground [snow level]
       Energy_out(1,ictrl)=GBAL
       Energy_out(2,ictrl)=RADT(2)
       Energy_out(3,ictrl)=SHF
           Energy_out(4,ictrl)=EGT
           Energy_out(5,ictrl)=EGI
           Energy_out(6,ictrl)=HG
           Energy_out(7,ictrl)=EGS
           Energy_out(8,ictrl)=ZLWU
c Extract the throughfall precipitation  [with units = meters/hr]
       Effective_PPT(1,ictrl)=p0 
      
   
 

c      if (ictrl.eq.49) print *, 'after energy balance calcs'

      IF(ABS (ZLHS - ZRHS) .GT. 0.010) THEN
         print *, 'warning: energy balance error greater than 0.01'

       IF(ABS (ZLHS - ZRHS) .GT. 4.000) THEN
         XPG1 = XPG1 + ZLHS                                             
         XPG2 = XPG2 + ZRHS                                             
         ectw=ect/dtt
         eciw=eci/dtt
         egtw=egt/dtt
         egiw=egi/dtt
         egsw=egs/dtt
 
           print *, ictrl,nmm,ndd,nhh,swdown,radn(3,2),zlwup,
     +               HFLUX,CHF,SHF,ectw,eciw,egtw,egiw,egsw,temprib,
     +               zlhs,zrhs
           print *, 'ssib3: energy balance error. snowdepth=',snowdepth,
     &    'bwo=',bwo,'tssno=',tssno,'tgs=',tgs,'gdia=',gdia,'capac(2)=',
     &    capac(2),'flo=',flo,'pixel=',pixel,'replicate=',replicate,
     &    'rank=',rank,'y0=',y0,'obsnow=',obsnow,'egi_old=',x0(9),
     &    'Tair=',uarg(6,ictrl),'ictrl=',ictrl,'Tair(i-1)=',
     &    uarg(6,ictrl-1),'x0=',x0,'zlat=',zlat,'zlong=',zlong,
     &    'LWDOWN=',uarg(4,ictrl),'radt(1)=',radt(1),'radt(2)=',radt(2)
      

         STOP

       END IF  
      END IF

c      if (ictrl.eq.49) print *, 'after balance checks'
c ** prepare for post-processing
c
c                                                                       
      wav    = sqrt(www(1)*www(2))                                      
      tenest = phsat * wav ** ( - bee)                                  
      tenes1 = phsat * (www(1) ** ( - bee) )                            
      tenes2 = phsat * (www(2) ** ( - bee) )                            
c
c ** store the output data for post-processing                         
c
c      if (ictrl.eq.49) print *, 'before x* calcs, hlat=',
c     & hlat

      xhc = hc / dtt + xhc                                              
      xhg = hg / dtt + xhg            
      xci = eci /hlat + xci                                             
      xct = ect /hlat + xct                                             
      xgi = egi /hlat + xgi                                             
      xgt = egt /hlat + xgt                                             
      xgs = egs /hlat + xgs 

cm CALCULATE OVERALL LATENT HEAT
      siblh = (((ECT + ECI + EGT + EGI + EGS)/HLAT
     & *(3150.19 - 2.378 * tm)) /3.6) 
cs  for france data albm=0.7, so Sun add folowing statement albm=0.4
cs on 03/29/99
c     albm=0.7
cs  sun add above statement on 03/29/99
      xmustar=0.7
      qmsensh=0.7
      albm=0.7
cs  sun add above statement on 03/29/99

c      if (ictrl.eq.49) print *, 'before output1 call...'



!cl ======================= DONGYUE ADD ITERATION ======================
!
!cl ------------ initialization -----------------
!! for 4 layer memls
!        n_rlz4=1
!        n_rz4=1
!        n_yz4=n_y_ssib !should be 14...
!        n_ypz4=1
!        n_yrz4=14
!        n_zz4=2 !number of observations
!        n_zpz4=2
!        n_zrz4=1
!        n_cz4=9
!        n_az4=8 !unused, actually...
!        n_xz4=13
!        n_uz4=9
!  !allocate & define ctrl 
!        allocate(ctrl4(n_cz4))
!        ctrl4=0
!        ctrl4(1)=0 !number of calculations to do, set in interfacez
!        ctrl4(2)=0 !number of snow layers, set in interfacez
!        ctrl4(3)=1 !atm_switch
!        ctrl4(4)=0 !veg_switch, set in interfacez
!  !these next numbers all pulled from old sizes.in file... careful!
!        ctrl4(5)=5 !number of auxiliary inputs for snow rtm... need to check!
!        ctrl4(6)=5 !number of snow inputs for snow rtm
!        ctrl4(7)=5 !number of inputs for canopy rtm
!        ctrl4(8)=4 !number inputs for atmospheric rtm
!        ctrl4(9)=1 !number of passive microwave channels
!  !allocate & set frequency and angle arrays
!  !allocate(freq(ctrl(9)),theta(ctrl(9)))
!        allocate(freq4(1),theta_m4(1))
!        freq4(1)=37.
!        theta_m4(1)=55.
!        month_m4=12 !hard-code to december
!  !note: just use local meas,rank as input to interfacez for debugging
!        ierr4=0 !dummy   
!        meas_switch4=1 !determines which measurements are used
!        iveg_type4(1)=veg_map
!        f_veg_m4(1,1)=0.
!        albedo_m4(1,1)=0. !dummy variable; not used in interfacez
!        n_bdrf4=1 !dummy variable to size bdrf dummy variable in interfacez
!        bdrf4=0. !dummy variable; not used in interfacez
!
!
!! for 3 layer memls
!        n_rlz3=1
!        n_rz3=1
!        n_yz3=n_y_ssib !should be 14...
!        n_ypz3=1
!        n_yrz3=14
!        n_zz3=2 !number of observations
!        n_zpz3=2
!        n_zrz3=1
!        n_cz3=9
!        n_az3=8 !unused, actually...
!        n_xz3=13
!        n_uz3=9
!  !allocate & define ctrl 
!        allocate(ctrl3(n_cz3))
!        ctrl3=0
!        ctrl3(1)=0 !number of calculations to do, set in interfacez
!        ctrl3(2)=0 !number of snow layers, set in interfacez
!        ctrl3(3)=1 !atm_switch
!        ctrl3(4)=0 !veg_switch, set in interfacez
!  !these next numbers all pulled from old sizes.in file... careful!
!        ctrl3(5)=5 !number of auxiliary inputs for snow rtm... need to check!
!        ctrl3(6)=5 !number of snow inputs for snow rtm
!        ctrl3(7)=5 !number of inputs for canopy rtm
!        ctrl3(8)=4 !number inputs for atmospheric rtm
!        ctrl3(9)=1 !number of passive microwave channels
!  !allocate & set frequency and angle arrays
!  !allocate(freq(ctrl(9)),theta(ctrl(9)))
!        allocate(freq3(1),theta_m3(1))
!        freq3(1)=37.
!        theta_m3(1)=55.
!        month_m3=12 !hard-code to december
!  !note: just use local meas,rank as input to interfacez for debugging
!        ierr3=0 !dummy   
!        meas_switch3=1 !determines which measurements are used
!        iveg_type3(1)=veg_map
!        f_veg_m3(1,1)=0.
!        albedo_m3(1,1)=0. !dummy variable; not used in interfacez
!        n_bdrf3=1 !dummy variable to size bdrf dummy variable in interfacez
!        bdrf3=0. !dummy variable; not used in interfacez
!cl ---------------- initialization done --------------------
!! determine whether go into the iteration
!      call ssib_layer(yout(1,ictrl-1),dz)
!      IF((ictrl.ne.1).and.((3600*obsnow).gt.2).and.
!     & ((yout(1,ictrl-1).ne.0.)))then
!
!!       4 layer memls       
!        allocate(y4memls(21,1),xz4(13,1),uz4(9,1,1))
!        y4memls=0.0
!        tbs4=0.0
!        xz4=0.0
!        uz4=0.0
!
!        y4memls(1:3,1)=dz(1:3)
!        y4memls(4,1)=3.6*obsnow
!        y4memls(5:7,1)=yout(2:4,ictrl-1)
!        y4memls(8,1)=yout(4,ictrl-1)
!        y4memls(9:11,1)=yout(5:7,ictrl-1)
!        y4memls(12,1)=yout(7,ictrl-1)
!        y4memls(13:15,1)=yout(8:10,ictrl-1)
!        y4memls(16,1)=yout(10,ictrl-1)
!        y4memls(17:19,1)=yout(11:13,ictrl-1)
!        y4memls(20,1)=yout(13,ictrl-1)
!        y4memls(21,1)=yout(14,ictrl-1)
!
!        xz4(:,1)=xout(:,ictrl-1)
!        uz4(:,1,1)=uarg(:,ictrl-1)
!
!! set water content to be 0
!        y4memls(13:16,1)=0.0
!
!
!! call interfacez4
!        call interfacez4(y4memls,zz4,n_rlz4,n_rz4,n_yz4,n_ypz4,n_yrz4,
!     & n_zz4,n_zpz4,n_zrz4,n_cz4,n_xz4,ctrl4,freq4,theta_m4,xz4,n_az4,
!     & n_u,uz4,month_m4,tbraw4,albedo_m4,ictrl,rank,iveg_type4,ierr4,
!     & f_veg_m4,bdrf4,n_bdrf4,meas_switch4,vcov_dat)
!
!        tbs4=tbs4+tbraw4(2,1)
!
!        !print*,'4 LAYER TB IS=====', TBS4
!        deallocate(y4memls,xz4,uz4)
!
!!       4 layer memls done 
!
!!       =========== begin run 3 layer memls ============       
!        allocate(y3memls(14,1),xz3(13,1),uz3(9,1,1))
!        y3memls=0.0
!        xz3=0.0
!        uz3=0.0
!        tbs3=0.0
!       
!
!        if(isnow.eq.0)then
!          y3memls(1,1)=sum(dzo)
!          y3memls(2,1)=bwo(1)
!          y3memls(3,1)=bwo(2)
!          y3memls(4,1)=bwo(3)
!          y3memls(5,1)=tssno(1)
!          y3memls(6,1)=tssno(2)
!          y3memls(7,1)=tssno(3)
!          y3memls(8,1)= flo(1)*bwo(1)/1000
!          y3memls(9,1)= flo(2)*bwo(2)/1000
!          y3memls(10,1)=flo(3)*bwo(3)/1000
!          y3memls(11,1)=gdia(1)
!          y3memls(12,1)=gdia(2)
!          y3memls(13,1)=gdia(3)
!          y3memls(14,1)=TGS
!        elseif(snowdepth.gt.0)then
!          y3memls(1,1)=snowdepth
!          y3memls(2,1)=capac(2)*1000/snowdepth
!          y3memls(3,1)=capac(2)*1000/snowdepth
!          y3memls(4,1)=capac(2)*1000/snowdepth
!          y3memls(5,1)=TGS
!          y3memls(6,1)=TGS
!          y3memls(7,1)=TGS
!          y3memls(8,1)=0.
!          y3memls(9,1)=0.
!          y3memls(10,1)=0.
!          y3memls(11,1)=gsize
!          y3memls(12,1)=gsize
!          y3memls(13,1)=gsize
!          y3memls(14,1)=TGS
!        else
!          y3memls(1:14,1)=0.0
!          y3memls(14,1)=TGS
!        endif
!
!! set water content to be 0
!        y3memls(8:10,1)=0.0
!
!
! 
!        xz3(1 ,1)=TD
!        xz3(2 ,1)=TC
!        xz3(3 ,1)=TA
!        xz3(4 ,1)=TM
!        xz3(5 ,1)=www(1)
!        xz3(6 ,1)=www(2)
!        xz3(7 ,1)=www(3)
!        xz3(8 ,1)=capac(1)
!        xz3(9 ,1)=egi
!        xz3(10,1)=gsize
!        xz3(11,1)=snden
!        xz3(12,1)=capac(2)
!        xz3(13,1)=roff
!
!        uz3(:,1,1)=uarg(:,ictrl)
!        
!! call interfacez3
!        call interfacez3(y3memls,zz3,n_rlz3,n_rz3,n_yz3,n_ypz3,n_yrz3,
!     & n_zz3,n_zpz3,n_zrz3,n_cz3,n_xz3,ctrl3,freq3,theta_m3,xz3,n_az3,
!     & n_u,uz3,month_m3,tbraw3,albedo_m3,ictrl,rank,iveg_type3,ierr3,
!     & f_veg_m3,bdrf3,n_bdrf3,meas_switch3,vcov_dat)
!
!        tbs3=tbs3+tbraw3(2,1)
!
!        !print*,'3 LAYER TB IS=====', TBS3
!        !print*, tbs3(1)-tbs4(1)
!        !print*,'ICTRL is', ictrl        
!        
!! iteration begins
!        counter=0
!        do while((tbs3(1)-tbs4(1)).gt.(0.005))
!          y3memls(11,1)=y3memls(11,1)+0.00001
!          call interfacez3(y3memls,zz3,n_rlz3,n_rz3,n_yz3,n_ypz3,n_yrz3,
!     & n_zz3,n_zpz3,n_zrz3,n_cz3,n_xz3,ctrl3,freq3,theta_m3,xz3,n_az3,
!     & n_u,uz3,month_m3,tbraw3,albedo_m3,ictrl,rank,iveg_type3,ierr3,
!     & f_veg_m3,bdrf3,n_bdrf3,meas_switch3,vcov_dat)
!          tbs3=0.0
!          tbs3=tbs3+tbraw3(2,1)
!          counter=counter+1
!        enddo
!
!        !print*,'itration times: ',counter
!! iteration ends 
!
!        if(isnow.eq.0)then
!          gdia(1)=y3memls(11,1)
!        elseif(snowdepth.gt.0)then
!          gsize=y3memls(11,1)
!        endif
!               
!        !print*,'================================'
!        
!        deallocate(y3memls,xz3,uz3)
!
!       ENDIF
!
!       deallocate(freq4,theta_m4,ctrl4,freq3,theta_m3,ctrl3)
!
!cl ======================= ITERATION DONE ==========================

      call output1(sunang,siblh,isnow,iday,n_steps,yout,n_y,xout,n_x,
     &  zout,tend,pixel)

c      if (ictrl.eq.49) print *, 'after output1 call...'

      END DO
      return
      end                                                               
cm **************************************************************
cm    END OF TIME LOOP
cm **************************************************************   

C=======================================================================
C                                                                       
       SUBROUTINE CNTROLS(MDAY,alpha,t0,x0,y0,z_in,ipixel,n_a) 
C                                                          1 AUGUST 1988
C=======================================================================
C                                                                       
C      INITIALISATION AND SWITCHES.                                   
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in'
      include 'snow4.in'
      dimension alpha(n_a),t0(5),x0(12),y0(14),z_in(2)
C
cm    READ(4,*)
cm    READ(4,*) DTT, ITRUNK, ILW
cm    READ(4,*) ZLAT, ZLONG, TIME,MONTH,MDAY,DAY,YEAR, NITER
cm    READ(4,*) ssisnow, snden
cm    READ(4,*) TC, TGS, TD, TA, TM, HT, QA
cm    READ(4,*) WWW(1),WWW(2),WWW(3),CAPAC(1),CAPAC(2)
cm    READ(4,*) DZ(1), DZ(2), DZ(3)
cm    READ(4,*) tssn(1),tssn(2),tssn(3)     
cm    READ(4,*) bw(1), bw(2), bw(3)
cm    READ(4,*) fl(1), fl(2), fl(3)    
cm    READ(4,*) gsize,gdia(1),gdia(2),gdia(3)
cm    READ(4,*) sg1,sg2,egi

CM INSTEAD OF READING THESE DATA FROM FILE, THEY ARE PASSED TO SSIB
CM   IN THE ALPHA,T0,X0,Y0,Z ARRAYS OR HARDCODED IN

CM    FIRST EXTRACT ABOVE ARRAYS TO COMMON VARS
      do i=1,3
        sg1(i)=alpha(i)
      end do
      sg2=alpha(4)
      sg3=alpha(5)
      sg4=alpha(6)

      time=t0(1)
      month=int(t0(2))
      mday=int(t0(3))
      day=t0(4)
      year=t0(5)
      TD=x0(1)
      TC=x0(2)
      TA=x0(3)
      TM=x0(4)
      WWW(1)=x0(5)
      WWW(2)=x0(6)
      WWW(3)=x0(7)
      CAPAC(1)=x0(8)
      egi=x0(9)
      gsize=x0(10)
      snden=x0(11)
      capac(2)=x0(12)
cm    note! there is no need to set roff to x0(13)...

      snowdepth=y0(1)
      bw(1)=y0(2)
      bw(2)=y0(3)
      bw(3)=y0(4)
      tssn(1)=y0(5)
      tssn(2)=y0(6)
      tssn(3)=y0(7)
cm      initialize fl as 0, then set only if bw>0
      fl(1)=0.
      fl(2)=0.
      fl(3)=0.
      if(bw(1).gt.0.) fl(1)=y0(8)/bw(1)*1000
      if(bw(2).gt.0.) fl(2)=y0(9)/bw(2)*1000
      if(bw(3).gt.0.) fl(3)=y0(10)/bw(3)*1000
      
      gdia(1)=y0(11)
      gdia(2)=y0(12)
      gdia(3)=y0(13)
      tgs=y0(14)
      zlat=z_in(1)
      zlong=z_in(2)

      tssno=0.
      flo=0.
      bwo=0.


CM DETERMINE DZ(1,2,3) FROM SNOWDEPTH USING LAYERING RULE 
      IF (snowdepth.le.0.05) THEN
        dz=0.0
      ELSE IF (snowdepth.gt.0.05.and.snowdepth.le.0.06) THEN     
        dz(1)=0.02
        dz(2)=0.02
        dz(3)=snowdepth- dz(1)- dz(2)
      ELSE IF ( snowdepth.gt.0.06.and.snowdepth.le.0.08) then
        dz(3)=0.02
        dz(2)=0.02
        dz(1)=snowdepth- dz(3)- dz(2)
      ELSE IF ( snowdepth.gt.0.08.and.snowdepth.le.0.62) then
        dz(3)=0.02
        dz(2)=(snowdepth- dz(3))*0.33333333 
        dz(1)=(snowdepth- dz(3))*0.66666667
      ELSE IF ( snowdepth.gt.0.62) then
        dz(3)=0.02
        dz(2)=0.20
        dz(1)=snowdepth- dz(3)- dz(2) 
      End IF 

cm  recompute capac(2) based on new snow variables, UNLESS, we have a
cm  snowpack that is greater than 0, but less than the critical depth.
cm  in that case, capac(2) is used to ke

      If(snowdepth.eq.0.)Then
        capac(2)=0.
      Else If (snowdepth.gt.0.05) Then
        capac(2)=(dz(1)*bw(1)+dz(2)*bw(2)+dz(3)*bw(3))/1000.
      Else
cm      set up one layer variables 

cm      this old way was to try to be fancy and preserve the densities
cm        determined by the update.  but i think it caused problems.
cm        the new way (directly below) resets the snow density to 4.3
cm        capac(2)=snowdepth*(bw(1)+bw(2)+bw(3))/3./1000.
cm        snden=snowdepth/capac(2)
          snden=4.3
          capac(2)=snowdepth/4.3
      End IF

CM ENTER HARDCODED VARIABLES
      DTT=3600.
c      ITRUNK=20
      ITRUNK=40
      ILW=1
      NITER=12
      SSISNOW=0.04
      HT=999.
      QA=999. 
 
C                  
      RETURN      
      END   
C=======================================================================
C                                                                       
      SUBROUTINE CONSTS                                                 
C                                                          1 AUGUST 1988
C=======================================================================
C                                                                       
C     INITIALIZATION OF PHYSICAL CONSTANTS                              
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in' 
C                                                                       
      PSUR     = 1000.                                                 
      CPAIR    = 1010.                                                  
      RHOAIR   = 1.225
cm    is this supposed to be the stefan-boltzman constant?
      STEFAN   = 5.669 * 10E-9                                          
      G        = 9.81                                                   
      VKC      = 0.41                                                   
      PIE      = 3.14159265                                             
      TIMCON   = PIE/86400.                                             
      CLAI     = 4.2 * 1000. * 0.2                                      
      CW       = 4.2 * 1000. * 1000.                                    
      TF       = 273.16                                                 
C-----------------------------------------------------------------------
C     N.B. : HLAT IS EXPRESSED IN J KG-1                                
C            SNOMEL IS EXPRESSED IN J M-1                               
C-----------------------------------------------------------------------
      HLAT     = ( 3150.19 - 2.378 * TM ) * 1000.                       
      SNOMEL   = 370518.5 * 1000.                                       
      PSY      = CPAIR / HLAT * PSUR / .622                             
C                                                                       
      RETURN                                                            
      END                                                               
C=====================================================================  
C                                                                       
      SUBROUTINE CROUT(A,N,M)                                           
C                                                         DECEM 1988    
C===================================================================== 
      DIMENSION A(N,M)                                                  
      DO 11 I=2,N                                                       
 11   A(1,I)=A(1,I)/A(1,1)                                              
      DO 3 K=2,N                                                        
      DO 2 IK=K,N                                                       
      DO 2 MI=2,K                                                       
 2    A(IK,K)=A(IK,K)-A(IK,MI-1)*A(MI-1,K)                              
      J1=K+1                                                            
      DO 3 J=J1,N                                                       
      DO 4 MJ=2,K                                                       
 4    A(K,J)=A(K,J)-A(K,MJ-1)*A(MJ-1,J)                                 
 3    A(K,J)=A(K,J)/A(K,K)                                              
      I1=N+1                                                            
      DO 5 I=I1,M                                                       
 5    A(1,I)=A(1,I)/A(1,1)                                              
      DO 7 JJ=I1,M                                                      
      DO 7 L=2,N                                                       
      DO 8 KL=2,L                                                       
 8    A(L,JJ)=A(L,JJ)-A(L,KL-1)*A(KL-1,JJ)                              
 7    A(L,JJ)=A(L,JJ)/A(L,L)                                            
      DO 10 JI=I1,M                                                     
      DO 10 K1=2,N                                                      
      K2=N-K1+2                                                         
      DO 10 K3=K2,N                                                     
      K4=N-K1+1                                                         
 10   A(K4,JI)=A(K4,JI)-A(K4,K3)*A(K3,JI)                               
      RETURN                                                            
      END                                                               
C====================================================================== 
C                                                                       
      SUBROUTINE DELRN ( RNCDTC, RNCDTG, RNGDTG, RNGDTC )               
C                                                                       
C====================================================================== 
C                                                                       
C     PARTIAL DERIVATIVES OF RADIATIVE AND SENSIBLE HEAT FLUXES         
C                                                                       
C---------------------------------------------------------------------- 
      include 'comsib.in' 
C                                                                       
      TC3 = TC * TC * TC                                                
      TG3 = TGS * TGS * TGS                                             
      FAC1 = ( 1. - ALBEDO(1,3,2) ) * ( 1.-THERMK ) * VCOVER(1)         
      FAC2 =   1. - ALBEDO(2,3,2)                                       
C                                                                       
      RNCDTC = - 2. * 4. * FAC1 * STEFAN * TC3                          
      RNCDTG = 4. * FAC1 * FAC2 * STEFAN * TG3                          
C                                                                       
      RNGDTG = - 4. * FAC2 * STEFAN * TG3                               
      RNGDTC = 4. * FAC1 * FAC2 * STEFAN * TC3                          
C                                                                       
      RETURN                                                            
      END                                                               
C====================================================================== 
C                                                                       
      SUBROUTINE DELHF ( HCDTC, HCDTG, HGDTG, HGDTC ) 
C                                                                       
C====================================================================== 
C                                                                       
C     PARTIAL DERIVATIVES OF SENSIBLE HEAT FLUXES                       
C                                                                       
C---------------------------------------------------------------------- 
      include 'comsib.in' 
C                                                                       

      RCP = RHOAIR * CPAIR                                              
      D1 = 1./RA + 1./RB + 1./RD                                        
      TA = ( TGS/RD + TC/RB + TM/RA ) / D1                              
C                                                                       

      HC = RCP * ( TC - TA ) / RB * DTT                                 
      HG = RCP * ( TGS - TA ) / RD * DTT                                
C---------------------------------------------------------------------- 
C     N.B. FLUXES EXPRESSED IN JOULES M-2                               
C---------------------------------------------------------------------- 
C                                                                       
      HCDTC = RCP / RB * ( 1./RA + 1./RD ) / D1                         
      HCDTG = - RCP / ( RB * RD ) / D1                                  
C                                                                       
      HGDTG = RCP / RD * ( 1./RA + 1./RB ) / D1                         
      HGDTC = - RCP / ( RD * RB ) / D1                                  
C                                                                       
      RETURN                                                            
      END                                                               
C====================================================================== 
C                                                                       
      SUBROUTINE DELEF ( ECDTC, ECDTG, EGDTG, EGDTC, DEADTC, DEADTG,    
     &        EC, EG, WC, WG, FC, FG, HR,MDLSNO,ISNOW)               
C                                                                       
C====================================================================== 
C                                                                       
C     PARTIAL DERIVATIVES OF LATENT HEAT FLUXES                         
C                                                                       
C---------------------------------------------------------------------- 
      include 'comsib.in' 
C                                                                       
C     RADD = 44.                                                        
      RCP = RHOAIR * CPAIR                                              
C---------------------------------------------------------------------- 
C     MODIFICATION FOR SOIL DRYNESS : HR = REL. HUMIDITY IN TOP LAYER   
C---------------------------------------------------------------------- 
C        
      HRR = HR                                                          
      IF ( FG .LT. .5 ) HRR = 1.                                        
C                                                                       
      RCC = RST(1)*FC + 2. * RB                                         
      COC = (1.-WC)/RCC + WC/(2.*RB)                                    
      RG = RST(2)*FG 
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
         RSURF=RSOIL
      ELSE                                         
         RSURF = RSOIL*FG
      END IF 
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HRR    
     &     + VCOVER(2)/(RSURF+RD+44.)*HRR                               
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)        
     &     + VCOVER(2)/(RSURF+RD+44.)                                   
      COG1 = COG1 + WG/RD * VCOVER(2)                                   
      COG2 = COG2 + WG/RD * VCOVER(2)                                   
C                                                                       
      D2 = 1./RA + COC + COG2                                           
      TOP = COC * ETC + COG1 * ETGS + EM/RA                             
      EA = TOP / D2                                                     
C                                                                       
      EC = ( ETC - EA ) * COC * RCP/PSY * DTT                           
C                                                                       
      EG = ( ETGS*COG1 - EA*COG2 ) * RCP/PSY * DTT                      
      DEADTC = GETC * COC / D2                                          

      DEADTG = GETGS * COG1 / D2                                        
C                                                                       
      ECDTC = ( GETC - DEADTC ) * COC * RCP / PSY                       
      ECDTG = - DEADTG * COC * RCP / PSY                                
C                                                                       
      EGDTG = ( GETGS*COG1 - DEADTG*COG2 ) * RCP / PSY                  
      EGDTC = - DEADTC * COG2 * RCP / PSY                               
C                                       
      RETURN                                                            
      END                                                               
C
C====================================================================   
C                                                                       
      SUBROUTINE FIT2(X1,Y1,DD,K,MM)                                    
C                                                       DECEM 1989      
C=====================================================================  
      PARAMETER ( NN = 30 )                                             
      DIMENSION X1(NN), Y1(11,NN), Y2(NN), SS(3,4),DD(11,3)     
      IOUT4 = 48                                                        
      DO 50 I = 1, 3                                                    
      DO 50 J = 1, 4                                                    
 50   SS(I,J) = 0.                                                      
      DO 100 I = 1, MM                                                  
      SS(1,2) = SS(1,2) + X1(I)                                         
      XX = X1(I) * X1(I)                                                
      SS(1,3) = SS(1,3) + XX                                            
      SS(2,3) = SS(2,3) + XX * X1(I)                                    
      SS(3,3) = SS(3,3) + XX * XX                                       
      SS(1,4) = SS(1,4) + Y1(K,I)                                       
      SS(2,4) = SS(2,4) + Y1(K,I) * X1(I)                               
      SS(3,4) = SS(3,4) + Y1(K,I) * XX                                  
 100  CONTINUE                                                          
      SS(1,1) = MM                                                      
      SS(2,1) = SS(1,2)                                                 
      SS(2,2) = SS(1,3)                                                 
      SS(3,1) = SS(2,2)                                                 
      SS(3,2) = SS(2,3)                                                 
      CALL CROUT(SS,3,4)                                                
      AA = SS(1,4)                                                      
      BB = SS(2,4)                                                      
      CC = SS(3,4)                                                      
      DO 200 I = 1, MM                                                  
      Y2(I) = AA + BB * X1(I) + CC * X1(I) * X1(I)                      
 200  CONTINUE                                                          
      DD(K,1) = AA                                                      
      DD(K,2) = BB                                                      
      DD(K,3) = CC                                                      
      XY2 = 0.                                                          
      DO 220 I = 1, MM                                                  
 220  XY2 = XY2 + (Y2(I) - Y1(K,I))**2                                  
      XY2 = SQRT(XY2/MM)                                                
      RETURN                                                            
      END       
C =====================================================================
c                                                                     
      SUBROUTINE INTERCS (ISNOW,p0,CSOIL,DZSOIL,CHISL)     
C                                                          1 AUGUST 1988
C=======================================================================
C                                                                       
C     CALCULATION OF (1) INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW 
C                    (2) SPECIFIC HEAT TERMS FIXED FOR TIME STEP        
C                                                                       
C     MODIFICATION 30 DEC 1985 : NON-UNIFORM PRECIPITATION             
C     ------------      CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY  
C                       RELATIONSHIP :-                                 
C                                                                       
C                                        F(X) = A*EXP(-B*X)+C           
C                                                                       
C                       THROUGHFALL, INTERCEPTION AND INFILTRATION      
C                       EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP      
C                       AND PROPORTION OF LARGE-SCALE PPN.              
C---------------------------------------------------------------------- 
      include 'comsib.in'
C                                                                       
      DIMENSION CAPACP(2), SNOWP(2), PCOEFS(2,2)                        
      DATA PCOEFS(1,1)/ 20. /, PCOEFS(1,2)/ .206E-8 /,                  
     &     PCOEFS(2,1)/ 0.0001 /, PCOEFS(2,2)/ 0.9999 /, BP /20. /      
C                                                                       
      AP = PCOEFS(2,1)                                                  
      CP = PCOEFS(2,2)                                                  
      TOTALP = PPC + PPL                                                
      IF(TOTALP.LT.1.E-8)GO TO 6000                                     
      AP = PPC/TOTALP * PCOEFS(1,1) + PPL/TOTALP * PCOEFS(2,1)          
      CP = PPC/TOTALP * PCOEFS(1,2) + PPL/TOTALP * PCOEFS(2,2)          
 6000 CONTINUE                                                          
C                                                                       
      ROFF = 0.                                                         
      THRU = 0.                                                         
      FPI  = 0.                                                         
C                                                                       
C---------------------------------------------------------------------- 
C     THERMAL CONDUCTIVITY OF THE SOIL, TAKING INTO ACCOUNT POROSITY    
C---------------------------------------------------------------------- 
C                                                                       
      THETA=WWW(1)*POROS                                                
      CHISL=( 9.8E-4+1.2E-3*THETA )/( 1.1-0.4*THETA )                   
      CHISL=CHISL*4.186E2                                               
C                                                                       
C---------------------------------------------------------------------- 
C     THERMAL DIFFUSIVITY AND HEAT CAPACITYOF THE SOIL                  
C---------------------------------------------------------------------- 
C                                                                       
      DIFSL=5.E-7                                                       
C                                                                       
      ROCS =CHISL/DIFSL                                                 
      D1   =SQRT(DIFSL*86400.0)                                         
      CSOIL=ROCS*D1/SQRT(PIE)/2.0 
C     YX2002 (test2)
      dzsoil=D1/SQRT(PIE)/2.0
      THALAS=0.
      OCEANS=0.
      POLAR=0.                                 
      CSOIL=CSOIL*(1.0-THALAS)+10.E10*OCEANS+POLAR*3.6*4.2E4            
C                                                                       
C                                                                       
      P0 = TOTALP * 0.001                                               

C                                                                       
C---------------------------------------------------------------------- 
C     INPUT PRECIPITATION IS GIVEN IN MM, CONVERTED TO M TO GIVE P0.    
C---------------------------------------------------------------------- 
C                                                                       
      DO 1000 IVEG = 1, 2                                               
C                                                                       
      SPWET1 = AMIN1 ( 0.05, CAPAC(IVEG))*CW                            
C                                                                       
      TS = TC                                                           
      SPECHT = ZLT(1) * CLAI                                            
      IF ( IVEG .EQ. 1 ) GO TO 1100                                     
      TS = TGS                                                          
      SPECHT = CSOIL                                                    
1100  CONTINUE                                                          
C                                                                       
      XSC = AMAX1(0., CAPAC(IVEG) - SATCAP(IVEG) )                      
      IF(IVEG.EQ.2 .AND. TS.LE.TF )GO TO 1170                           
      CAPAC(IVEG) = CAPAC(IVEG) - XSC                                   
      ROFF = ROFF + XSC                                                 
c     this line added new
c     rnoffs(nobs) = rnoffs(nobs) + XSC
1170  CONTINUE                                                          
      CAPACP(IVEG) = 0.                                                 
      SNOWP(IVEG) = 0.                                                  
C                                                                       
      IF( TS .GT. TF ) CAPACP(IVEG) = CAPAC(IVEG)                       
      IF( TS .LE. TF ) SNOWP(IVEG) = CAPAC(IVEG)                        
      CAPAC(IVEG) = CAPACP(IVEG)                                        
      SNOWW(IVEG) = SNOWP(IVEG)                                         
      ZLOAD = CAPAC(IVEG) + SNOWW(IVEG)                                 
C                                                                       
      FPI = ( 1.-EXP( - EXTK(IVEG,3,1) * ZLT(IVEG)/VCOVER(IVEG) ) )     
     &      * VCOVER(IVEG)                                              
      TTI = P0 * ( 1.-FPI )                                            
C                                                                       
C---------------------------------------------------------------------- 
C    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)            
C---------------------------------------------------------------------- 
C                                                                       
      XS = 1.                                                           
      IF ( P0 .LT. 1.E-9 ) GO TO 1150                                   
      ARG =  ( SATCAP(IVEG)-ZLOAD )/( P0*FPI*AP ) -CP/AP                
      IF ( ARG .LT. 1.E-9 ) GO TO 1150                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
1150  TEX = P0*FPI * ( AP/BP*( 1.- EXP( -BP*XS )) + CP*XS ) -           
     &      ( SATCAP(IVEG) - ZLOAD ) * XS                               
      TEX = AMAX1( TEX, 0. )                                            
C                                                                       
C---------------------------------------------------------------------- 
C    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION                    
C---------------------------------------------------------------------- 
C                                                                       
      THRU = TTI + TEX                                                  
      IF(IVEG.EQ.2.AND.TGS.LE.TF)THRU = 0.                              
C                                                                       
      PINF = P0 - THRU                                                  
      IF( TM .GT. TF ) CAPAC(IVEG) = CAPAC(IVEG) + PINF                 
      IF( TM .LE. TF ) SNOWW(IVEG) = SNOWW(IVEG) + PINF                 
C                                                                       
      IF( IVEG .EQ. 1 ) GO TO 1300                                      
      IF( TM .GT. TF ) GO TO 1200                                      
      SNOWW(IVEG) = SNOWP(IVEG) + P0                                    
      THRU = 0.                                                         
      GO TO 1300                                                        
C                                                                       

C---------------------------------------------------------------------- 
C    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )                  
C---------------------------------------------------------------------- 
C                                                                       
1200  EQUDEP = SATCO * DTT                                              
C                                                                       
      XS = 1.                                                           
      IF ( THRU .LT. 1.E-9 ) GO TO 1250                                 
      ARG = EQUDEP / ( THRU * AP ) -CP/AP                               
      IF ( ARG .LT. 1.E-9 ) GO TO 1250                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
1250  ROFFO = THRU * ( AP/BP * ( 1.-EXP( -BP*XS )) + CP*XS )            
     &       -EQUDEP*XS                                                 
      ROFFO = AMAX1 ( ROFFO, 0. )                                       
      ROFF = ROFF + ROFFO                                               
c this next line was added by dr xue
c     rnoffs(nobs) = rnoffs(nobs) + roffo
      WWW(1) = WWW(1) + (THRU - ROFFO) / ( POROS*ZDEPTH(1) ) 
1300  CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C    TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION                
C---------------------------------------------------------------------- 
C      
                                                                 
      DIFF = ( CAPAC(IVEG)+SNOWW(IVEG) - CAPACP(IVEG)-SNOWP(IVEG) )*CW 
      CCP = SPECHT + SPWET1                                             
      CCT = SPECHT + SPWET1 + DIFF                                      
C                                                                       
      TSD = ( TS * CCP + TM * DIFF ) / CCT                              
C                                                                       
      FREEZE = 0.                                                       
      IF ( TS .GT. TF .AND. TM .GT. TF ) GO TO 2000                     
      IF ( TS .LE. TF .AND. TM .LE. TF ) GO TO 2000                     
C                                                                       
      TTA = TS                                                          
      TTB = TM                                                          
      CCA = CCP                                                         
      CCB = DIFF                                                        
      IF ( TSD .GT. TF ) GO TO 2100                                     
C                                                                       
C---------------------------------------------------------------------- 
C    FREEZING OF WATER ON CANOPY OR GROUND                              
C---------------------------------------------------------------------- 
C                                                                       
      CCC = CAPACP(IVEG) * SNOMEL                                       
      IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW                        
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
C                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMIN1 ( CCC, FREEZE )) / SNOMEL                         
      IF(TSD .GT. TF)TSD = TF - 0.1                                     
C                                                                       
      
      GO TO 2000                                                        
C                                                                       
2100  CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C    MELTING OF SNOW ON CANOPY OR GROUND                                
C---------------------------------------------------------------------- 
C                                                                       
      CCC = - SNOWW(IVEG) * SNOMEL                                      
      IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW                      
C                                                                       
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
C                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMAX1( CCC, FREEZE )) / SNOMEL                          
      IF(TSD .LE. TF)TSD = TF - 0.1                                     
C                                                                       

2000  SNOWW(IVEG) = SNOWW(IVEG) + FREEZE                                
      CAPAC(IVEG) = CAPAC(IVEG) - FREEZE                                
C                                                                       
      IF( IVEG .EQ. 1 ) TC = TSD                                        
      IF( IVEG .EQ. 2 ) TGS = TSD                                       
      IF( SNOWW(IVEG) .LT. 0.0000001 ) GO TO 3000                       
c     modeified to force water into soil Xue Feb. 1994
c     ZMELT = 0.                                                        
c     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
c     IF ( TD .LE. TF ) ROFF = ROFF + CAPAC(IVEG)                       
      ZMELT = CAPAC(IVEG)                             
      CAPAC(IVEG) = 0.                                                  
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
C                                                                       
3000  CONTINUE                                                          
C                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
      SNOWW(IVEG) = 0.                                                  
C                                                                       
      P0 = THRU   
      

      IF (ISNOW.eq.0) go to 1001                                        
1000  CONTINUE                                                          

C                                                                       
C---------------------------------------------------------------------- 
C    CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.                  
C    N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN  
C    DEALING WITH VERY LATGE SNOWPACKS.                                 
C---------------------------------------------------------------------- 
C                                                                       
1001  CCX = ZLT(1) * CLAI + CAPAC(1) * CW                               
      SPWET = AMIN1 ( 0.05, CAPAC(2))*CW                                
      CG = (CSOIL + SPWET)                                              
C                                                                       
      RETURN                                                            
      END   
C  ================================================================= 
C=======================================================================RAD03700
C                                                                       RAD03710
      SUBROUTINE LONGRN( TRANC1, TRANC2, TRANC3)
C                                                          1 AUGUST 1988RAD03730
C=======================================================================RAD03740
C                                                                       RAD03750
C     CALCULATION OF DOWNWARD LONGWAVE.                                 RAD03760
C                                                                       RAD03770
C-----------------------------------------------------------------------RAD03780
      include 'comsib.in'                      
C                                                                       RAD03800
      DIMENSION TRANC1(2), TRANC2(2), TRANC3(2)
C                                                                       RAD03820
      IF(ILW .EQ. 1)GO TO 101                 
      IF(ILW .EQ. 2)GO TO 102                
      IF(ILW .EQ. 3)GO TO 103               
101   CONTINUE                             
C---------------------------------------------------------------------- RAD03870
C     DOWNWARD LONG-WAVE ASSUMED TO BE PROVIDED AS RADN(3,2)            RAD03880
C---------------------------------------------------------------------- RAD03890
      GO TO 200                           
C                                                                       RAD03910
102   CONTINUE                           
C---------------------------------------------------------------------- RAD03930
C     DOWNWARD LONG-WAVE FROM BRUNTS EQUATION, MONTEITH(1973), P37.    RAD03940
C---------------------------------------------------------------------- RAD03950
      ESKY = 0.53 + 0.06*SQRT(EM)       
      RADN(3,2)  =  ESKY*(1.+0.2*(CLOUD*CLOUD))*STEFAN*TM**4   
      GO TO 200                        
C                                                                       RAD03990
103   CONTINUE                                                
C---------------------------------------------------------------------- RAD04010
C     DOWNWARD LONG-WAVE FLUX CALCULATED AS RESIDUAL FROM MEASURED      RAD04020
C     NET RADIATION AND OUTGOING LONGWAVE RADIATION.                    RAD04030
C                                                                       RAD04040
C     CALCULATION OF ABSORBED FRACTIONS OF RADIATION ( EXPENDABLE )     RAD04050
C---------------------------------------------------------------------- RAD04060
      DO 2000 IWAVE = 1, 2                                   
C                                                                       RAD04080
      RAB(2,IWAVE,1) =  ( 1. - VCOVER(1) ) *                
     &  ( RADN(IWAVE,1) * ( 1. - ALBEDO(2,IWAVE,1) ) )     
      RAB(2,IWAVE,2) =  ( 1. - VCOVER(1) ) *              
     &    RADN(IWAVE,2) * ( 1. - ALBEDO(2,IWAVE,2) )     
C                                                                       RAD04130
      RAB(2,IWAVE,1) = RAB(2,IWAVE,1) + VCOVER(1) *     
     &  ( RADN(IWAVE,1) * ( TRANC1(IWAVE) * ( 1. - ALBEDO(2,IWAVE,1) ) +
     &    TRANC3(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) ) )               
      RAB(2,IWAVE,2) = RAB(2,IWAVE,2) + VCOVER(1) *                   
     &    RADN(IWAVE,2) * TRANC2(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) 
C                                                                       RAD04190
      RAB(1,IWAVE,1) =  VCOVER(1) *                                 
     &    RADN(IWAVE,1) * ( ( 1. - ALBEDO(1,IWAVE,1) ) -           
     &    TRANC1(IWAVE) * ( 1. - ALBEDO(2,IWAVE,1) ) -            
     &    TRANC3(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) )           
      RAB(1,IWAVE,2) =  VCOVER(1) *                             
     &    RADN(IWAVE,2) * ( ( 1. - ALBEDO(1,IWAVE,2) ) -       
     &    TRANC2(IWAVE) * ( 1. - ALBEDO(2,IWAVE,2) ) )        
2000  CONTINUE                                               
C                                                                       RAD04280
      SWAB = RAB(1,1,1) + RAB(1,1,2) + RAB(1,2,1) + RAB(1,2,2) + 
     &       RAB(2,1,1) + RAB(2,1,2) + RAB(2,2,1) + RAB(2,2,2)  
      SWUP = SWDOWN - SWAB                                     
      RADN(3,2) = RNETM - SWAB + ZLWUP                        
C                                                                       RAD04330
200   CONTINUE                                               
C                                                                       RAD04350
      RETURN                                                
      END                                                  
C====================================================================   
C                                                                       
       SUBROUTINE NEWTON(A1,Y,FINC,NOX,NONPOS,IWOLK,L)
C                                                                       
C====================================================================== 
C ** VERSION ACQUIRED FROM EROS 2/19/86.                                
C                                                                       
C=======================================================================
C                                                                       
C ** THE NEWTON RAPHSON ITERATIVE ROUTINE WILL BE USED TO GENERATE NEW  
C ** VALUES OF A1 IF DABSOLUTE VALUE OF Y IS GREATER THAN ERTOL;        
C ** A1 IS ESTIMATE, Y IS RESULTANT ERROR                               
C ** NOX IS EXIT CONDITION  (0=NO EXIT) OR (1 WHEN DABS(Y) LT ERTOL)    
C ** ERTOL IS THE DABSOLUTE VALUE OF Y NECESSARY TO OBTAIN AN EXIT      
C ** FINC IS INITIAL INCREMENT SIZE FOR SECOND ESTIMATE OF A1           
C ** NONPOS=0 IF QUANTITY TO BE MINIMIZED CAN BE LESS THAN ZERO;        
C ** NONPOS=1 IF QUANTITY CAN ONLY BE POSITIVE                          
C ** L IDENTIFIES WHICH QUANTITY IS BEING CALCULATED.                   
C                                                                       
C ** CONTROL VALUES: FINC,ERTOL,NOX,NONPOS,L:MUST BE SET BY USER        
C-----------------------------------------------------------------------
C                                                                       
CM  Mike modified this on June 12, 2006 so that if the objective value Y
CM   is identical in two successive iterations (if Y==Y1(L) ), we will use
CM   the derivative from the previous iteration.  The assumption is, if
CM   the objective value from the previous iteration is identical, the 
CM   derivative from that iteration wont be off by too much.

       DIMENSION   IWALK(3), NEX(3)                                     
       COMMON/TONNEW/  ZINC(3), A2(3), Y1(3)                            
       COMMON/NEWT/  ITER(3),PREV_DYDA
       DATA CONS/1.0/                                                   
C                                                                       
       ERTOL = 0.05 * FINC                                              
       IWALK(L) = IWOLK                                                 
       NEX(L)=NOX                                                       
C                                                                       
       IF ( ITER(L) .GE. 490 ) GO TO 160                                
       IF (ERTOL .LT. 0.00000001) ERTOL=0.000001                        
       IF (ABS(Y) .LE. ERTOL) GO TO 150                                 
       IF((ABS(Y-Y1(L))).LE.0.01*ERTOL .AND. IWALK(L).EQ.0 ) GO TO 8    
C                                                                     
       IF(ABS(Y1(L)).GT.ERTOL) GO TO 1                                  
       A2(L)=A1                                                         
       A1=A1-Y                                                          
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       ITER(L)=1                                                        
       IF (IWALK(L) .EQ. 3) GO TO 101                                   
       IWALK(L)=0                                                       
       GO TO 101                                                        
   1   ITER(L)=ITER(L)+1                                                
       IF(ITER(L) .EQ. 10) IWALK(L)=1                                   
       IF(IWALK(L) .NE. 0) GO TO 2                                      
       IF(ABS(Y) .GT. ERTOL) GO TO 3                                    
       NEX(L)=1                                                         
       GO TO 150                                                        
   3   if( abs(y-y1(l)).lt.1e-5 ) then
CM     Mike is adding this if statement to augment the original line of
CM       code, in order to handle the case where Y==Y1(L).  In that case,
CM       we will use the derivative from the previous iteration, then 
CM       continue as normal.
c         print*, 'Warning in SSiB/Newton: Objective value for this',
c     &     'iteration identical to that of previous iteration. Using',
c     &     'derivative from previous iteration instead of',
c     &     ' recalculating; prev_dyda=',prev_dyda
         a=a1-y/prev_dyda
c         print *, 'Done with new calculation'
       else
         A=A1-Y*(A1-A2(L))/(Y-Y1(L))                                   
       end if
       IF(ABS(A-A1).GT.(10.0*FINC))                                    
     &     A=A1+10.0*FINC*SIGN(CONS,(A-A1))                      
CM     Mike is adding this line to save the derivative at this iteration
       prev_dyda=(y-y1(l))/(a1-a2(l))

       A2(L)=A1                                                         
       A1=A                                                             
       Y1(L)=Y                                                          
       GO TO 101                                                        
   2   IF(IWALK(L).EQ.2)GO TO 4                                         
       IF(IWALK(L).EQ.3) GO TO 6                                        
       IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO  3                    
       ZINC(L)=(A1-A2(L))/4.0                                           
       A1=A2(L)+ZINC(L)                                                 
       IWALK(L)=2                                                       
       NEX(L)=0                                                         
       GO TO 101                                                        
   4   IF(SIGN(CONS,Y) .EQ.SIGN(CONS,Y1(L))) GO TO 5                    
       ZINC(L)=-ZINC(L)/4.0                                             
       A2(L)=A1                                                         
       A1=A1+ZINC(L)                                                    
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       GO TO 101                                                        
   5   A2(L)=A1                                                         
       A1=A1+ZINC(L)                                                    
       Y1(L)=Y                                                          
       NEX(L)=0                                                         
       GO TO 101                                                        
   6   IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO 7                     
       IWALK(L)=1                                                       
       GO TO 2                                                          
   7   A2(L) = A1                                                       
       A1 = A1+FINC                                                     
       Y1(L)=Y                                                          
       NEX(L) = 0                                                       
       GO TO 101                                                        
   8   A1 = A1 + FINC*2.0                                               
       NEX(L)=0                                                         
       GO TO 101                                                        
160    CONTINUE                                                         
c      WRITE(6,900) Y, A1                                               
C      STOP                                                             
 900   FORMAT ( 3X,' FAILURE TO CONVERGE AFTER 490 ITERATIONS',         
     & /, 3X,' Y = ',2G12.5,2X,I14)                                     
C                                                                       
 150   NEX(L) = 1                                                       
       ZINC(L)=0.0                                                      
       ITER(L) = 0                                                      
       IWALK(L)=0                                                       
       Y1(L)=0.0                                                        
       Y=0.0                                                            
       A2(L)=0.0                                                        
 101   CONTINUE                                                         
       IF(NONPOS.EQ.1.AND.A1.LT.0.0) A1=A2(L)/2.0                       
       NOX = NEX(L)                                                     
       IWOLK = IWALK(L)                                                 
C                                                                       
       RETURN                                                           
       END                                                              
C ===================================================================    RAD000
       SUBROUTINE RADAB(ISNOW,MDLSNO,afac,ictrl,ipixel,swenew,sweold,
     &        sag1,sag2)
	  
C                                                         1 AUGUST 1988 RAD00040
C=======================================================================RAD00050
C                                                                       RAD00060
C     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT       RAD00070
C     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY                     RAD00080
C                                                                       RAD00090
C-----------------------------------------------------------------------RAD00100
c     DIMENSION XQQ(2,2,2)                                              RAD00110
      include 'comsib.in'                                          
      include 'snow4.in'
!      real sag1, sag2
C                                                                       RAD00140
      DIMENSION TRANC1(2), TRANC2(2), TRANC3(2)                  
C                                                                       RAD00160
      F = SUNANG
      if (ictrl.eq.1) then
       sag1=0.0
      end if
C                                                                       RAD00180
C---------------------------------------------------------------------- RAD00190
C     CALCULATION OF MAXIMUM WATER STORAGE VALUES.                      RAD00200
C---------------------------------------------------------------------- RAD00210
C                                                                       RAD00220
      FMELT = 1.0                                              
      IF ( ABS(TF-TGS) .LT. 0.5 ) FMELT = 0.6                 
      SATCAP(1) =  ZLT(1) * 0.0001                           
      SATCAP(2) =  ZLT(2) * 0.0001                          
CS  Sun change following DEPCOV     10/13/98
      IF (ISNOW.eq.0) THEN 
         DEPCOV = AMAX1( 0., (snowdepth-Z1))
      ELSE
         DEPCOV = AMAX1( 0., (CAPAC(2)*snden-Z1))          
      END IF
CS----------------------------10/13/98 
      DEPCOV = AMIN1( DEPCOV, (Z2-Z1)*0.95 )                   
      SATCAP(1) = SATCAP(1) * ( 1. - DEPCOV / ( Z2 - Z1 ) )   
C                                                                       RAD00300
C---------------------------------------------------------------------- RAD00310
C                                                                       RAD00320
      DO 1000 IWAVE = 1, 2                                   
C                                                                       RAD00340
      DO 2000 IVDUM = 1, 2                                  
C                                                                       RAD00360
      IF ( IVDUM .EQ. 1 ) IVEG = 2                         
      IF ( IVDUM .EQ. 2 ) IVEG = 1                        
C---------------------------------------------------------------------- RAD00390
C     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO            RAD00400
C         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING  RAD00410
C         SNOW TRANSMITTANCE = 0.20, 0.54                               RAD00420
C                                                                       RAD00430
C---------------------------------------------------------------------- RAD00440
c      if (ictrl.eq.49) print *, 'before modification for snow...',
c     &  'iwave=',iwave,'ivdum=',ivdum

      SCOV = 0.                                                     
      IF( IVEG .EQ. 2 ) GO TO 100                                  
      IF( TC .LE. TF ) SCOV =  AMIN1( 0.5, CAPAC(1) / SATCAP(1) ) 
100   CONTINUE                                                   
      REFF1 = ( 1. - SCOV ) * REF(IVEG,IWAVE,1) + SCOV * ( 1.2 -
     &        IWAVE * 0.4 ) * FMELT                            
      REFF2 = ( 1. - SCOV ) * REF(IVEG,IWAVE,2) + SCOV * ( 1.2 -
     &        IWAVE * 0.4 ) * FMELT                                    
      TRAN1 = TRAN(IVEG,IWAVE,1) * ( 1. - SCOV )                      
     &        + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT )         
     &        * TRAN(IVEG,IWAVE,1)                                  
      TRAN2 = TRAN(IVEG,IWAVE,2) * ( 1. - SCOV )                      
     &        + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT ) * 0.9   
     &        * TRAN(IVEG,IWAVE,2)                                  
C                                                                       RAD00590
C---------------------------------------------------------------------- RAD00600
C                                                                       RAD00610
      SCAT = GREEN(IVEG)*( TRAN1 + REFF1 ) +( 1. - GREEN(IVEG) ) *    
     &       ( TRAN2 + REFF2)                                        
      CHIV = CHIL(IVEG)                                             
C                                                                       RAD00650
      IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01                       
      AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV                
      BB = 0.877 * ( 1. - 2. * AA )                                    
C                                                                     
      PROJ = AA + BB * F                                             
      EXTKB = ( AA + BB * F ) / F                                   
      ZMEW = 1. / BB * ( 1. - AA / BB * ALOG ( ( AA + BB ) / AA ) )
      ACSS = SCAT / 2. * PROJ / ( PROJ + F * BB )                     
      ACSS = ACSS * ( 1. - F * AA / ( PROJ + F * BB ) * ALOG ( ( PROJ  
     *       +   F * BB + F * AA ) / ( F * AA ) ) )                  
C                                                                       RAD00760
      EXTK( IVEG, IWAVE, 1 ) = PROJ / F * SQRT( 1.-SCAT )           
      EXTK( IVEG, IWAVE, 2 ) = 1. / ZMEW * SQRT( 1.-SCAT )         
      EXTK( IVEG, 3, 1 ) = AA + BB                                
      EXTK( IVEG, 3, 2 ) = 1./ZMEW                               
C                                                                       RAD00810
      UPSCAT = GREEN(IVEG) * TRAN1 + ( 1. - GREEN(IVEG) ) * TRAN2     
      UPSCAT = 0.5 * ( SCAT + ( SCAT - 2. * UPSCAT ) *               
     *         (( 1. - CHIV ) / 2. ) ** 2 )                         
C                                                                       RAD00850
      BETAO = ( 1. + ZMEW * EXTKB ) / ( SCAT * ZMEW * EXTKB ) * ACSS 
C                                                                       RAD00870
C---------------------------------------------------------------------- RAD00880
C                                                                       RAD00890
C     DICKINSONS VALUES                                                RAD00900
C                                                                       RAD00910
      BE = 1. - SCAT + UPSCAT                                     
      CE = UPSCAT                                                   
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )              
      IF ( ABS(BOT) .GT. 1.E-10) GO TO 200                       
      SCAT = SCAT* 0.98                                         
      BE = 1. - SCAT + UPSCAT                                  
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )         
200   CONTINUE                                               
      DE = SCAT * ZMEW * EXTKB * BETAO                      
      FE = SCAT * ZMEW * EXTKB * ( 1. - BETAO )            
C---------------------------------------------------------------------- RAD01020
C                                                                       RAD01030
      CCE = DE * BE - ZMEW * DE * EXTKB + CE * FE         
      FFE = BE * FE + ZMEW * FE * EXTKB + CE * DE        
C                                                                       RAD01060
 
      TORE = -CCE / BOT                                 
      SIGE = -FFE / BOT                                
C                                                                       RAD01090
      PSI = SQRT(BE**2 - CE**2)/ZMEW                  
C                                                                       RAD01110
C---------------------------------------------------------------------- RAD01120
C     REDUCTION IN EXPOSED HEIGHT OF UPPER STOREY AS SNOW ACCUMULATES   RAD01130
C                                                                       RAD01140
CS Sun Change following SDEP to SDEP=snowdepth on 10/13/98

c      if (ictrl.eq.49) print *, 'radab: reduction in height...',
c     &  'iwave=',iwave,'ivdum=',ivdum,'z2=',z2,'z1=',z1,'vcover=',
c     &   vcover

      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
          SDEP=snowdepth
      ELSE 
          SDEP = CAPAC(2) *snden
      END IF
CS 

      FAC = ( SDEP - Z1 ) / ( Z2 - Z1 )              
      FAC = AMAX1( 0., FAC )                        
      FAC = AMIN1( 0.99, FAC )                     
C                                                                       RAD01190

      ZAT = ZLT(IVEG) / VCOVER(IVEG)              
      IF ( IVEG .EQ. 1 ) ZAT = ZAT * (1.-FAC)    
C                                                                       RAD01220
      POWER1 = AMIN1( PSI*ZAT, 50. )            
      POWER2 = AMIN1( EXTKB*ZAT, 50. )         
      EPSI = EXP( - POWER1 )                  
      EK = EXP ( - POWER2 )                  
C                                                                       RAD01270
      ROSB = SOREF(IWAVE)                   
      ROSD = SOREF(IWAVE)                  
      IF ( IVEG .EQ. 2 ) GO TO 300        
      ROSB = ALBEDO(2,IWAVE,1)           
      ROSD = ALBEDO(2,IWAVE,2)          
300   CONTINUE                         
C                                                                       RAD01340

c      if (ictrl.eq.49) print *, 'radab: reduction in height...',
c     &  'rosd=',rosd

      GE = ROSB / ROSD                
C                                                                       RAD01360
C-----------------------------------------------------------------------RAD01370
C     CALCULATION OF DIFFUSE ALBEDOS                                    RAD01380
C-----------------------------------------------------------------------RAD01390
C                                                                       RAD01400
c      if (ictrl.eq.49) print *, 'radab: before calc. of diff alb...',
c     &  'iwave=',iwave,'ivdum=',ivdum

      F1 = BE - CE / ROSD                                             
      ZP = ZMEW * PSI                                                
C                                                                   
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -                     
     &      ( BE - ZP ) * ( F1 + ZP ) * EPSI                      
      ALPHA = CE * ( F1 - ZP ) / EPSI / DEN                      
      BETA = -CE * ( F1 + ZP ) * EPSI / DEN                     
      F1 = BE - CE * ROSD                                      
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI           
C                                                                       RAD01500
      GAMMA = ( F1 + ZP ) / EPSI / DEN                       
      DELTA = - ( F1 - ZP ) * EPSI / DEN                    
C                                                                       RAD01530
      ALBEDO(IVEG,IWAVE,2) =  ALPHA + BETA                 
C     XQQ(IVEG,IWAVE,2) = ALBEDO(IVEG, IWAVE, 2)                        RAD01550
C                                                                       RAD01560
      IF ( IVEG .EQ. 1 ) GO TO 400                        
      SCOV2 = 0.                                         
      IF ( TGS .LE. TF ) SCOV2 = AMIN1( 1., CAPAC(2) / 0.004 )         
      ALBEDO(2,IWAVE,2) =                                             
     & ROSD * ( 1. - VCOVER(2) ) + ALBEDO(2,IWAVE,2) * VCOVER(2)     
      ALBEDO(2,IWAVE,2) =                                           
     & ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,2) + SCOV2 *                
     & ( 1.2-IWAVE*0.4 ) *                                     
     & FMELT                                                  
cm
cm   below in direct albedo section is the call to my albedo function, 
cm    which will overwrite diffuse albedo, if there is snow
cm

400   CONTINUE                                               
C                                                                       RAD01670
      TRANC2(IWAVE) = GAMMA * EPSI + DELTA / EPSI           
C                                                                       RAD01690
C-----------------------------------------------------------------------RAD01700
C     CALCULATION OF DIRECT ALBEDOS                                     RAD01710
C-----------------------------------------------------------------------RAD01720
C                                                                       RAD01730
c      if (ictrl.eq.49) print *, 'radab: before calc. dir. albs...',
c     &  'iwave=',iwave,'ivdum=',ivdum

      F1 = BE - CE / ROSD                                  
      ZMK = ZMEW * EXTKB                                  
C                                                                       RAD01760
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -           
     &      ( BE - ZP ) * ( F1 + ZP ) * EPSI                           
      ALPHA = ( DE - TORE * ( BE + ZMK ) ) * ( F1 - ZP ) / EPSI -     
     &        ( BE - ZP ) * ( DE - CE*GE - TORE * ( F1 + ZMK ) ) * EK
      ALPHA = ALPHA / DEN                                           
      BETA = ( BE + ZP ) * (DE - CE*GE - TORE * ( F1 + ZMK ))* EK -
     &       ( DE - TORE * ( BE + ZMK ) ) * ( F1 + ZP ) * EPSI    
      BETA = BETA / DEN                                          
      F1 = BE - CE * ROSD                                       
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI            
      GAMMA = - SIGE * ( F1 + ZP ) / EPSI -                   
     &        ( FE + CE * GE * ROSD + SIGE * ( ZMK - F1 ) ) * EK       
      GAMMA = GAMMA / DEN                                             
      DELTA = ( CE * GE * ROSD + FE + SIGE * ( ZMK - F1 ) ) * EK     
     &        + SIGE * ( F1 - ZP ) * EPSI                           
      DELTA = DELTA / DEN                                          
C                                                                       RAD01930
      ALBEDO(IVEG,IWAVE,1) = TORE + ALPHA + BETA                  
C     XQQ(IVEG,IWAVE,1) = ALBEDO(IVEG, IWAVE, 1)                        RAD01950
C                                                                       RAD01960
C---------------------------------------------------------------------- RAD01970
C                                                                       RAD01980

      IF( IVEG .EQ. 1 ) GO TO 500                                      
      ALBEDO(2,IWAVE,1) = ROSB * ( 1. - VCOVER(2) )                   
     &                    + ALBEDO(2,IWAVE,1) * VCOVER(2)            
      ALBEDO(2,IWAVE,1) = ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,1) +      
     &                    SCOV2 * ( 1.2-IWAVE*0.4 ) * FMELT        
                                                                  
cm    here is the call to compute the albedo over the snow surface
cm
      if(ISNOW==0.and.IVEG==2)then
      
      ! Manu [Oct2010] changes the snow albedo calculation, I think the best one is the CLM
      ! parameterization but still need to perform further albedo assesments
        
c COMMENT/UNCOMMENT THE FOLLOWING FOR SNTHERM:
        ! call sntalb(gdia(3),SUNANG,ALBEDO,salbo,IWAVE,afac)
        
c COMMENT/UNCOMMENT THE FOLLOWING FOR CLM:
        ! call snowage(3600.0,tgs,swenew*1000.0,sweold*1000.0,sag1,sag2)		 
        !  call albsnowCLM(sunang,sag2,ALBEDO(2,1:2,1:2))
         sag1 = sag2

      end if

c      if (ictrl.eq.49) print *, 'radab: after call to sntherm...',
c     &  'iwave=',iwave,'ivdum=',ivdum,'albedo(2,2,2)=',albedo(2,2,2),
c     &  'gida(3)=',gdia(3),'afac=',afac
cm
500   CONTINUE                                                         
C                                                                       RAD02060
      TRANC1(IWAVE) = EK                                              
      TRANC3(IWAVE) = SIGE * EK + GAMMA * EPSI + DELTA / EPSI        
C                                                                       RAD02090
2000  CONTINUE                                                      
C                                                                       RAD02110
C---------------------------------------------------------------------- RAD02120
C     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES    RAD02130
C     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND              RAD02140
C---------------------------------------------------------------------- RAD02150
C                                                                       RAD02160
c      if (ictrl.eq.49) print *, 'before calc. of terms modify sw_in...',
c     &  'iwave=',iwave,'ivdum=',ivdum

      RADFAC(2,IWAVE,1) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,1) )  
     &       + VCOVER(1) * ( TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) ) 
     &       + TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )            
C                                                                       RAD02200
      RADFAC(2,IWAVE,2) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,2) )  
     &       + VCOVER(1) *  TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) )  
C                                                                       RAD02230
      RADFAC(1,IWAVE,1) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,1) )     
     &       - TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) )             
     &       - TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )          
C                                                                       RAD02270
      RADFAC(1,IWAVE,2) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,2) )       
     &       - TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )             
C                                                                       RAD02300
C     XQQ(1,IWAVE,1) = RADFAC(1,IWAVE,1)                                RAD02310
C     XQQ(1,IWAVE,2) = RADFAC(1,IWAVE,2)                                RAD02320
C     XQQ(2,IWAVE,1) = RADFAC(2,IWAVE,1)                                RAD02330
C     XQQ(2,IWAVE,2) = RADFAC(2,IWAVE,2)                                RAD02340
C                                                                       RAD02350
C                                                                       RAD02360
C---------------------------------------------------------------------- RAD02370
C     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB )                     RAD02380
C---------------------------------------------------------------------- RAD02390
C                                                                       RAD02400


      DO 3000 IRAD = 1, 2                                              
      SALB(IWAVE,IRAD) = ( 1.-VCOVER(1) ) * ALBEDO(2,IWAVE,IRAD) +    
     &                   VCOVER(1) * ALBEDO(1,IWAVE,IRAD)            
3000  CONTINUE                                                      

C                                                                       RAD02450
C---------------------------------------------------------------------- RAD02460
C     SAVING OF EXTINCTION COEFFICIENTS ( PAR ) FOR STOMAT CALCULATION  RAD02470
C---------------------------------------------------------------------- RAD02480
c      if (ictrl.eq.49) print *, 'before calc. of extinction coeff...',
c     &  'iwave=',iwave,'ivdum=',ivdum

      IF ( IWAVE .EQ. 2 ) GO TO 600                                    
      RADSAV(1) = 1. - VCOVER(1)                                      
     &          + VCOVER(1) * ( TRANC1(IWAVE) + TRANC3(IWAVE) )      
      RADSAV(2) = 1. - VCOVER(1) + VCOVER(1) * TRANC2(IWAVE)        
C     XQQ(1,1,1) = RADSAV(1)                                            RAD02530
C     XQQ(1,2,1) = RADSAV(2)                                            RAD02540
600   CONTINUE                                                     
C                                                                       RAD02560
1000  CONTINUE                                                    

CM THESE LINES (FROM HERE TO THE 730 MARK) ARE AN ALBEDO ADJUSTMENT THAT IS
CM NOW OBSELETE, THOUGH IT DIDNT APPEAR TO BE IN AFFECT WHEN I RECEIVED
CM THE CODE.
c
c     albedo adjustment ==============================================
c     if (SDEP.gt.0.) then
c        xadj=xadj1
c     else
c        xadj=0.0
c     end if
c     temp test
c     xadj=0.
c     if (xadj.eq.0.) go to 730
c     xx = radfac(1,1,2) + radsav(2)
c     xy = radfac(1,1,1) + radsav(1)
c     ssum = salb(1,1)*frac(1,1) + salb(1,2)*frac(1,2)+
c    &       salb(2,1)*frac(2,1) + salb(2,2)*frac(2,2)
c     for diffuse albedo
c     do 650 iwave = 1, 2
c     salb(iwave,2) = salb(iwave,2) + xadj * salb(iwave,2) / ssum
c     x0 = 1. - salb(iwave,2)
c     x1 = radfac(1,iwave,2) + radfac(2,iwave,2)
c     x2 = radfac(1,iwave,2) / x1
c     x3 = radfac(2,iwave,2) / x1
c     radfac(1,iwave,2) = x0 * x2
c     radfac(2,iwave,2) = x0 * x3
c     if (salb(iwave,2).gt.1..or.radfac(1,iwave,2).gt.1..or.
c    &    radfac(2,iwave,2).gt.1..or.salb(iwave,2).lt.0..or.
c    &    radfac(1,iwave,2).lt.0..or.radfac(2,iwave,2).lt.0.) then
c          write(6,640) nymdh,iwave,salb(iwave,2),radfac(1,iwave,2),
c     &                 radfac(2,iwave,2)
c     end if
c650  continue
c640  format(1x,'unrealistic value, dif',2i12,4e11.4)
c     for direct albedo
c     do 750 iwave = 1, 2
c     salb(iwave,1) = salb(iwave,1) + xadj * salb(iwave,1) / ssum
c     x0 = 1. - salb(iwave,1)
c     x1 = radfac(1,iwave,1) + radfac(2,iwave,1)
c     x2 = radfac(1,iwave,1) / x1
c     x3 = radfac(2,iwave,1) / x1
c     radfac(1,iwave,1) = x0 * x2
c     radfac(2,iwave,1) = x0 * x3
c     radsav(1) =  xy - radfac(1,1,1)
c     radsav(2) =  xx - radfac(1,1,2)
c     if (salb(iwave,1).gt.1..or.radfac(1,iwave,1).gt.1..or.
c    &    radfac(2,iwave,1).gt.1..or.salb(iwave,1).lt.0..or.
c    &    radfac(1,iwave,1).lt.0..or.radfac(2,iwave,1).lt.0.) then
c          write(6,740) nymdh,iwave,salb(iwave,1),radfac(1,iwave,1),
c     &                  radfac(2,iwave,1)
c     end if
c750  continue
c740  format(1x,'unrealistic value',2i12,4e11.4)
c730  continue
***************** end adjustment *******************************
c     temporary change

      sibswup = radn(1,1)*salb(1,1) + radn(1,2)*salb(1,2)
     &                     + radn(2,1)*salb(2,1) + radn(2,2)*salb(2,2)

      if ((swdown.gt.0.1).and.(sibswup.gt.0.1)) then
         sibalbedo = sibswup / swdown
         if (sibalbedo.gt.1.) then
            sibswup =  0.
            sibalbedo = -99.
         endif
      else
         sibswup = -9999. 
         sibalbedo = -99.
      endif

C                                                                       RAD02580
C     WRITE(48, 1010) F                                                 RAD02590
C     WRITE(48, 1010) RADFAC(1,1,1), RADFAC(1,2,1), RADFAC(1,1,2),      RAD02600
C    *   RADFAC(1,2,2), RADFAC(2,1,1), RADFAC(2,2,1), RADFAC(2,1,2),    RAD02610
C    *   RADFAC(2,2,2), SALB(1,1), SALB(1,2), SALB(2,1), SALB(2,2),     RAD02620
C    *   RADSAV(1), RADSAV(2)                                           RAD02630
C1010 FORMAT(5X, 4F11.4)                                                RAD02640
C---------------------------------------------------------------------- RAD02650
C                                                                       RAD02660
C     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND        RAD02670
C                                                                       RAD02680
C---------------------------------------------------------------------- RAD02690
C                                                                       RAD02700
      TC4 = TC * TC * TC * TC                                         
      TG4 = TGS * TGS * TGS * TGS                                    
C                                                                       RAD02730
      ZKAT = EXTK(1,3,2) * ZLT(1) / VCOVER(1)                       
      ZKAT = AMIN1( 50. , ZKAT )                                   
      ZKAT = AMAX1( 1.E-5, ZKAT )                                 
      THERMK = EXP(-ZKAT)                                        
C                                                                       RAD02780
      FAC1 =  VCOVER(1) * ( 1.-THERMK )                         
      FAC2 =  1.                                               
      CLOSS =  2. * FAC1 * STEFAN * TC4                       
      CLOSS =  CLOSS - FAC2 * FAC1 * STEFAN * TG4            
      GLOSS =  FAC2 * STEFAN * TG4                          
      GLOSS =  GLOSS - FAC1 * FAC2 * STEFAN * TC4          
C                                                                       RAD02850
      ZLWUP =  FAC1 * STEFAN * TC4 + (1. - FAC1 ) * FAC2 * STEFAN * TG4
      TGEFF = SQRT( SQRT ( ( ZLWUP / STEFAN ) ) )                     
C                                                                       RAD02880
      RADSAV(3) = EXTK(1,1,1)                                        
      RADSAV(4) = EXTK(1,1,2)                                       
      RADSAV(5) = EXTK(2,1,1)                                      
      RADSAV(6) = EXTK(2,1,2)                                     
      RADSAV(7) = THERMK                                         
      RADSAV(8) = EXTK(1,3,1)                                   
      RADSAV(9) = EXTK(2,3,1)                                  
      RADSAV(10)= CLOSS                                       
      RADSAV(11)= GLOSS                                      
      RADSAV(12)= TGEFF                                     

C-----------------------------------------------------------------------RAD03000
C                                                                       RAD03010
      CALL LONGRN( TRANC1, TRANC2, TRANC3)                         
C                                                                       RAD03030
C-----------------------------------------------------------------------RAD03040
C                                                                       RAD03050
      CALL RADUSE                                                 
C                                                                       RAD03070
C-----------------------------------------------------------------------RAD03080
C                                                                       RAD03090
      RETURN                                                     
      END                                                       
C
c                                                                       RAD03120
      SUBROUTINE RADUSE                                        
C                                                          1 AUGUST 1988RAD03160
C=======================================================================RAD03170
C                                                                       RAD03180
C     CALCULATION OF ABSORPTION OF RADIATION BY SURFACE                 RAD03190
C                                                                       RAD03200
C-----------------------------------------------------------------------RAD03210
      include 'comsib.in'                                     
C                                                                       RAD03230
      P1F         = RADSAV(1)                                
      P2F         = RADSAV(2)                               
      EXTK(1,1,1) = RADSAV(3)                              
      EXTK(1,1,2) = RADSAV(4)                             
      EXTK(2,1,1) = RADSAV(5)                            
      EXTK(2,1,2) = RADSAV(6)                           
      THERMK      = RADSAV(7)                          
      EXTK(1,3,1) = RADSAV(8)                         
      EXTK(2,3,1) = RADSAV(9)                        
      CLOSS       = RADSAV(10)                      
      GLOSS       = RADSAV(11)                     
      TGEFF       = RADSAV(12)                    
C                                                                       RAD03360
C---------------------------------------------------------------------- RAD03370
C     SUMMATION OF SHORT-WAVE RADIATION ABSORBED BY CANOPY AND GROUND   RAD03380
C---------------------------------------------------------------------- RAD03390
C                                                                       RAD03400
      RADT(1) = 0.                               
      RADT(2) = 0.                              
C                                                                       RAD03430
      DO 1000 IVEG  = 1, 2                     
      DO 1000 IWAVE = 1, 2                    
      DO 1000 IRAD  = 1, 2                   
C                                                                       RAD03470
      RADT(IVEG) = RADT(IVEG)+RADFAC(IVEG,IWAVE,IRAD)*RADN(IWAVE,IRAD)
C                                                                    
1000  CONTINUE                                                      
C                                                                       RAD03510
C                                                                       RAD03520
      RADT(1) = RADT(1) + RADN(3,2)*VCOVER(1)*(1.- THERMK)         
     &        - CLOSS                                             
      RADT(2) = RADT(2) + RADN(3,2)*( 1.-VCOVER(1)*(1-THERMK) )  
     &        - GLOSS                                           
C                                                                       RAD03570
C---------------------------------------------------------------------- RAD03580
C                                                                       RAD03590
      PAR(1) = RADN(1,1) + RADN(1,2) + 0.001                   
      PD(1) = ( RADN(1,1) + 0.001 ) / PAR(1)                  
      P1 = P1F * RADN(1,1) + 0.001                           
      P2 = P2F * RADN(1,2)                                  
      PAR(2) = P1 + P2                                     
      PD(2) = P1 / PAR(2)                                 
C                                                                       RAD03660
      RETURN                                             
      END                                               
C                                                                       RAD03690

C====================================================================   
C                                                                       
      SUBROUTINE RADC2(COSZ)                                            
C                                                                       
C=====================================================================  
C                                                                       
C ** SOLAR ZENITH ANGLE COMPUTATION; DOWNCOMING RADIATION AT BOTTOM.    
C                                                                       
      include 'comsib.in' 
C                                                                       
      DAYSPY = 365.                                                     
      IF ( AMOD( YEAR, 4. ) .EQ. 0. ) DAYSPY = 366.                     
C                                                                       
C ** JULIAN DAY AND TIME UPDATE; SKIP ON 1ST TIME STEP (INITIALIZED)  
C
      IF(ITER .EQ. 1)GO TO 10                                           
      TIME = TIME + DTT / 3600.                                         
      IF ( TIME .GE. 23.99 ) TIME = 0.0                                 
      DAY = DAY +  DTT / 86400.                                         
C                                                                       
   10 CONTINUE                                                          
C                                                                       
      IF ( DAY .GT. DAYSPY ) YEAR = YEAR + 1.                           
      IF ( DAY .GT. DAYSPY ) DAY = DAY - DAYSPY                         
C                                                                       
C ** SOLAR DECLINATION CALCULATION                                      
C                                                                       
      DECMAX = PIE * ( 23.5 / 180.)                                     
      SOLS   = ( 4141./24. ) + AMOD( YEAR+3., 4. ) * 0.25               
C                                                                       
      SEASON = ( DAY - SOLS ) / 365.2                                   
      DEC    = DECMAX * COS ( 2. * PIE * SEASON )                       
C                                                                       
      RFD  = PIE / 180.                                                 
      SIND = SIN( DEC )                                                 
      COSD = COS( DEC )                                                 
      HAC  = -TAN( ZLAT * RFD )*TAN( DEC )                              
      HAC  = AMIN1(HAC,1.0)                                             
      HAC  = AMAX1(HAC,-1.0)                                            
C                                                                       
C **  H IS THE HALF-DAY LENGTH (IN RADIANS)                             
C                                                                       
      H   = ACOS(HAC)                                                  
      DAWN= -H                                                          
      DUSK= +H                                                          
      SR  = 12.-(H/(15.*RFD))                                           
      SS  = 12.+(H/(15.*RFD))                                           
C                                                                       
C     WRITE(6,53)SR,SS                                                  
C  53 FORMAT(1X,' SUNRISE IS AT ',F8.2,' HOURS AND SUNSET IS AT ',F8.2, 
C    &' HOURS .')                                                       
C                                                                       
C ** CALCULATION OF SOLAR ANGLE (SUNANG)                                
C                                                                       
      COSHR = COS( - PIE + (TIME + 0.5*DTT/3600.) / 24. * 2. * PIE )    
ccc   test (12/2/92)
c     COSHR = COS( - PIE + (TIME + 1.5*DTT/3600.) / 24. * 2. * PIE )    
c     COSHR = COS( - PIE + (TIME) / 24. * 2. * PIE )    
C                                                                       
      SUNANG = SIN( ZLAT*RFD ) * SIND + COS ( ZLAT*RFD ) * COSD * COSHR 
                                                                        
      SUNANG = ( SUNANG + SQRT( SUNANG * SUNANG ) ) / 2.                
C$    SUNANG = SUNANG + 0.01745 * ( 1. - SUNANG )                       
      COSZ= SUNANG + 0.01745 * ( 1. - SUNANG )                          
C                                                                       
C ** CALCULATION OF INCOMING SHORT-WAVE RADIATION : EMPIRICAL           
C ** BYPASSED FOR SIBDUM AND RUTHE (OFF-LINE) TESTS.                    
C                                                                       
C$    SWI = SUNANG*1360.0*( 0.23 + 0.57 * 0.8 )                         
C$                                                                      
C$    DIRR = 1.0634 - 0.0902 / ( 0.0150 + SUNANG )                      
C$    IF ( DIRR .GT. 0.94 ) DIRR = 0.94                                 
C$    IF ( DIRR .LT. 0.02 ) DIRR = 0.02                                 
C$    DIFF = 1. - DIRR                                                 
C$    RADN(1,1) = DIRR * 0.5 * SWI                                      
C$    RADN(1,2) = DIFF * 0.5 * SWI                                      
C$    RADN(2,1) = 0.25 * SWI                                            
C$    RADN(2,2) = 0.25 * SWI                                            
C$    RADN(3,2) = 200.0                                                 
C                                                                       
      RETURN                                                            
      END                                                               
C=====================================================================  
C                                                                       
      SUBROUTINE RASIT5(TRIB,CTNI,CUNI,FTT,FVV,rank,pixel,replicate,
     &   ictrl,iter_num)
C                                                                       
C=======================================================================
C     CUU AND CTT ARE LINEAR  (A SIMPLIFIED VERSION, XUE ET AL. 1991)  
C                                                                      
      include 'comsib.in'
      COMMON / XRIB/ RIB,temprib                                        
cm    debug information
      integer rank,pixel,replicate,ictrl,iter_num
      FS(X) = 66.85 * X                                              
      FT(X) = 0.904 * X                                          
      FV(X) = 0.315 * X                                              

C                                                                      
C                                                                      
C     CU AND CT ARE THE FRICTION AND HEAT TRANSFER COEFFICIENTS.    
C     CUN AND CTN ARE THE NEUTRAL FRICTION AND HEAT TRANSFER         
C     COEFFICIENTS.                                                     
C                                                                       
C                                                                       
      G2= 0.75                                                          
      G3= 0.75                                                          
      Z22 = Z2                                                          
      ZL = Z2 + 11.785 * Z0                                             
cm      this line from ratko in 05/29/06 email
cm      if(zwind.le.xdd.or.zl.le.xdd) xdd=min(zwind,zl)-0.1d0
cm      assume xdd==d, then we have:
      if(zwind.le.d.or.zl.le.d) d=min(zwind,zl)-0.1d0

      Z2 = D + Z0                                                       

      CUNI = ALOG((ZWIND-D)/Z0)/VKC                                     
      IF (ZL.LT.ZWIND) THEN                                             
         XCT1 = ALOG((ZWIND-D)/(ZL-D))                                  
         XCT2 = ALOG((ZL-D)/(Z2-D))                                     
         XCTU2 = ALOG((ZL-D)/(Z22-D))                                   
         CTNI = (XCT1 + G3 * XCT2) / VKC                                
      ELSE                                                              
         XCT2 =  ALOG((ZWIND-D)/(Z2-D))                                 
         XCTU2 =  ALOG((ZWIND-D)/(Z22-D))                               
         CTNI = G3 * XCT2 /VKC                                          
      END IF                                                            

C        NEUTRAL VALUES OF USTAR AND VENTMF                             
C                                                                       
         USTARN=UM/CUNI                                                 
         VENTN =RHOAIR   /CTNI*USTARN                                   
      IF (ZL.LT.ZWIND) THEN                                             
         U2 = UM - 1. / VKC * USTARN * (XCT1 + G2 *                     
     *           XCTU2)                                                 
      ELSE                                                              
         U2 = UM - 1. / VKC * USTARN * G2 * XCTU2                      
      END IF                                                            
cm    this line added based on ratkos suggestion, 5/30/06
      U2=AMAX1(U2,0.01)

C                                                                       
C     STABILITY BRANCH BASED ON BULK RICHARDSON NUMBER.                 
C                                                                       
      THVGM= TRIB-TM                                                    
      IF (TA.EQ.0.) THVGM = 0.                                          
      RIB  = -THVGM*G*(ZWIND-D) / (TM*(UM-U2)**2)                       
      RIB  = MAX(-10. E0 ,RIB)                                          
      RIB  = MIN( .1643 E0 ,RIB)                                        

C                                                                       
C     NON-NEUTRON CORRECTION  (SEE XUE ET AL(1991))                     
      IF(RIB.LT.0.0)THEN                                                
         GRIB = +RIB                                                    
         GRZL = +RIB*(ZL-D)/(ZWIND-D)                                   
         GRZ2 = +RIB*(Z2-D)/(ZWIND-D)                                   
         FVV =  FV(GRIB)                                                
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = FT(GRIB) + (G3-1.) * FT(GRZL) - G3 * FT(GRZ2)        
         ELSE                                                           
             FTT = G3*(FT(GRIB) - FT(GRZ2))                             
         END IF                                                         
         CUI = CUNI + FVV
         CTI = CTNI + FTT
      ELSE                                                              
         RZL = RIB/(ZWIND-D)*(ZL-D)                                     
         RZ2 = RIB/(ZWIND-D)*(Z2-D)                                     
         FVV = FS(RIB)                                                  
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = FS(RIB) + (G3-1) * FS(RZL) - G3 * FS(RZ2)            
         ELSE                                                           
             FTT = G3 * (FS(RIB) - FS(RZ2))                             
         END IF                                                         
 312     continue
         CUI = CUNI + FVV
cm       in a certain bizarre combination of roughness length, wind speed,
cm       and air temperature, combined with a certain value of TRIB set by 
cm       the n-r loop, we get cui=0., and then divide by it a few lines 
cm       down.  therefore i'm making sure cui doesn't get set to 0.
cm       see journal: june 21, 2006
         if(abs(cui).lt.1.E-5) cui=0.01
         CTI = CTNI + FTT
      ENDIF                                                             
 310  CONTINUE                                                          
C                                                                       
      USTAR =UM/CUI                                                     
      RAF = CTI / USTAR                                                 
      IF (RAF.LT.0.80) RAF = 0.80                                       
C                                                                       
      RA  = RAF                                                         
      XTEM1 = RAF                                                       
C                                                                       
      UEST  = USTAR                                                     
      DRAG = RHOAIR * UEST*UEST                                         
      Z2 = Z22                                                          
 1010 FORMAT(1X,'RIB,CTI,CUI,CTN,CUN',I10,7E10.3)                       
 1011 FORMAT(1X,'RIB,RAF,USTAR,UM,U2',7E10.3)                           
      RETURN                                                            
      END                                                               
C=======================================================================
C                                                                       
       SUBROUTINE RBRD1                                                 
C                                                          DECMB    1988
C=======================================================================
C                                                                       
C      CALCULATION OF RB AND RD AS FUNCTIONS OF U2 AND TEMPERATURES     
C                                                                       
C-----------------------------------------------------------------------
       include 'comsib.in' 
C                                                                       
       RB  = 1.0/(SQRT(U2)/RBC+ZLT(1)*.004)                             
C                                                                       
       TGTA = TGS- TA                                                   
       TEMDIF = ( TGTA + SQRT(TGTA*TGTA) ) / 2. + 0.1                   
       FIH = SQRT( 1. + 9. * G *TEMDIF * Z2 / TGS / ( U2*U2) )          
       RD  = RDC / U2 / FIH                                             
C                                                                       
       RETURN                                                           
       END                                                              
C=======================================================================
C                                                                       
      SUBROUTINE ROOT1                                                  
C                                                          1 DEC 1988   
C=======================================================================
C                                                                       
C    CALCULATION OF SOIL MOISTURE POTENTIALS IN ROOT ZONE OF EACH       
C    VEGETATION LAYER AND SUMMED SOIL+ROOT RESISTANCE                   
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in'
C                                                                       
      DO 1000 IL = 1, 3                                                 
      PHSOIL(IL) = PHSAT * AMAX1( 0.05, WWW(IL) ) ** ( - BEE )          
 1000 CONTINUE                                                          
C                                                                       
C-----------------------------------------------------------------------
C     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE USED FOR SOURCE      
C-----------------------------------------------------------------------
C                                                                       
C                                                                       
      PHROOT(1) = PHSOIL(1)-0.01                                        
C                                                                       
      DO 1200 I = 2 ,3                                                  
 1200 PHROOT(1) = AMAX1( PHROOT(1), PHSOIL(I) )                         
      PHROOT(2) = PHROOT(1)                                             
C                                                                       
C                                                                       
      RETURN                                                            
      END  
c                                                             
C=======================================================================
c
      subroutine run2                                                   
c                                                                               
c=======================================================================        
c    calculation of interflow, infiltration excess and loss to                  
c    groundwater .  all losses are assigned to variable 'roff' .                
c----------------------------------------------------------------------         
      include 'comsib.in' 
c                                                                               
      dimension temw(3), temwp(3), temwpp(3),  
     &          aaa(2) , bbb(2)  , ccc(2)   , qqq(2)                    
c                                                                               
      do 1000 i = 1, 3                                                  
c                                                                               
      temw(i)   = amax1( 0.03, www(i) )                                 
      temwp(i)  = temw(i) ** ( -bee )                                   
      temwpp(i) = amin1( 1., temw(i)) ** ( 2.*bee+ 3. )                 
1000  continue                                                          
c                                                                               
c-----------------------------------------------------------------------        
c                                                                               
c    calculation of gravitationally driven drainage from w(3) : taken           
c    as an integral of time varying conductivity.addition of liston             
c    baseflow term to original q3g to insure flow in                            
c    dry season. modified liston baseflow constant scaled                       
c    by available water.                                                        
c                                                                               
c     q3g (q3) : equation (62) , SE-86                                          
c                                                                               
c-----------------------------------------------------------------------        
c                                                                               
      pows = 2.*bee+2.                                                  
      q3g = temw(3)**(-pows) + satco/zdepth(3)/poros*slope*pows*dtt     
      q3g = q3g ** ( 1. / pows )                                        
      q3g = - ( 1. / q3g - www(3) ) * poros * zdepth(3) / dtt           
      q3g = amax1( 0., q3g )                                            
      q3g = amin1( q3g, www(3)*poros*zdepth(3)/dtt )                    
c                                                                               
      q3g = q3g + 0.002*poros*zdepth(3)*0.5 / 86400. * www(3)           
c                                                                               
c----------------------------------------------------------------------         
c                                                                               
c    calculation of inter-layer exchanges of water due to gravitation           
c    and hydraulic gradient. the values of w(x) + dw(x) are used to             
c    calculate the potential gradients between layers.                          
c    modified calculation of mean conductivities follows ME-82 ), 
c    reduces recharge flux to top layer.                      
c                                                                               
c      dpdw           : estimated derivative of soil moisture potential         
c                       with respect to soil wetness. assumption of             
c                       gravitational drainage used to estimate likely          
c                       minimum wetness over the time step.                     
c                                                                               
c      qqq  (q     )  : equation (61) , SE-86                                   
c             i,i+1                                                             
c            -                                                                  
c      avk  (k     )  : equation (4.14) , ME-82                                 
c             i,i+1                                                             
c                                                                               
c----------------------------------------------------------------------         
c                                                                               
      wmax = amax1( www(1), www(2), www(3), 0.05 )                     
      wmax = amin1( wmax, 1. )                                        
      pmax = wmax**(-bee)                                            
      wmin = (pmax-2./( phsat*(zdepth(1)+2.*zdepth(2)+zdepth(3))))     
     &        **(-1./bee)                                              
      wmin = amin1( www(1), www(2), www(3), wmin )                     
      wmin = amax1( wmin, 0.02 )                                       
      pmin = wmin**(-bee)                                              
      dpdw = phsat*( pmax-pmin )/( wmax-wmin )                         
c                                                                      
      do 2000 i = 1, 2                                                 
c                                                                               
      rsame = 0.                                                       
      avk  = temwp(i)*temwpp(i) - temwp(i+1)*temwpp(i+1)               
      div  = temwp(i+1) - temwp(i)                                     
      if ( abs(div) .lt. 1.e-6 ) rsame = 1.                            
      avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame )              
      avkmin = satco * amin1( temwpp(i), temwpp(i+1) )                 
      avkmax = satco * amax1( temwpp(i), temwpp(i+1) )*1.01            
      avk = amax1( avk, avkmin )                                       
      avk = amin1( avk, avkmax )                                       
c                                                                               
c-----------------------------------------------------------------------        
c     conductivities and base flow reduced when temperature drops below         
c     freezing.                                                                 
c-----------------------------------------------------------------------        
c                                                                               
      tsnow = amin1 ( tf-0.01, tgs ) 
      areas = amin1(0.999,13.2*snoww(2))
      tgg = tsnow*areas + tgs*(1.-areas) 
      ts    = tgg*(2-i) + td*(i-1)                                     
      props = ( ts-(tf-10.) ) / 10.                                    
      props = amax1( 0.05, amin1( 1.0, props ) )                       
      avk  = avk * props                                               
      q3g  = q3g * props                                               
c                                                                               
c-----------------------------------------------------------------------        
c     backward implicit calculation of flows between soil layers.               
c-----------------------------------------------------------------------        
c                                                                               
      dpdwdz = dpdw * 2./( zdepth(i) + zdepth(i+1) )                  
      aaa(i) = 1. + avk*dpdwdz*( 1./zdepth(i)+1./zdepth(i+1) )        
     &            *dtt/poros                                          
      bbb(i) =-avk *   dpdwdz * 1./zdepth(2)*dtt/poros                
      ccc(i) = avk * ( dpdwdz * ( www(i)-www(i+1) ) + 1. +            
     &           (i-1)*dpdwdz*q3g*1./zdepth(3)*dtt/poros )            
2000  continue                                                        
c                                                                               
      denom  = ( aaa(1)*aaa(2) - bbb(1)*bbb(2) )                       
      rdenom = 0.                                                      
      if ( abs(denom) .lt. 1.e-6 ) rdenom = 1.                         
      rdenom = ( 1.-rdenom)/( denom + rdenom )                         
      qqq(1)   = ( aaa(2)*ccc(1) - bbb(1)*ccc(2) ) * rdenom            
      qqq(2)   = ( aaa(1)*ccc(2) - bbb(2)*ccc(1) ) * rdenom           
c                                                                               
c-----------------------------------------------------------------------        
c     update wetness of each soil moisture layer due to layer interflow         
c        and base flow.                                                         
c-----------------------------------------------------------------------        
c                                                                               
      www(3) = www(3) - q3g*dtt/(poros*zdepth(3))
      roff = roff + q3g * dtt                                          
c                                                                               
      do 3000 i = 1, 2                                                
c                                                                               
      qmax   =  www(i)   * (poros*zdepth(i)  /dtt)                     
      qmin   = -www(i+1) * (poros*zdepth(i+1)/dtt)                     
      qqq(i) = amin1( qqq(i),qmax)                                     
      qqq(i) = amax1( qqq(i),qmin)                                     
      www(i)   =   www(i)   - qqq(i)/(poros*zdepth(i)  /dtt)           
      www(i+1) =   www(i+1) + qqq(i)/(poros*zdepth(i+1)/dtt)           
3000  continue    
c this next line added by dr xue
cm    rnoffb(nobs) = rnoffb(nobs) + q3g*dtt
c                                                                               
      do 4000 i = 1, 3                                                 
      excess = amax1(0.,(www(i) - 1.))                                 
      www(i) = www(i) - excess                                         
      roff   = roff   + excess * poros*zdepth(i)                       
c this next was added by dr xue
c       *** LOAD IN as root-drainage for PILPS
c
c     if (i.lt.2) then
c       rnoffs(nobs)= rnoffs(nobs)+ excess*POROS*ZDEPTH(I)
c     else
c       rnoffb(nobs)= rnoffb(nobs)+ excess*POROS*ZDEPTH(I)
c     endif
 4000  continue
c                                                                               
c-----------------------------------------------------------------------        
c     prevent negative values of www(i)                                         
c-----------------------------------------------------------------------        
c                                                                               
      do 4002 i = 1,2                                                  
      deficit   = amax1 (0.,(1.e-12 - www(i)))                         
      www (i)   = www(i) + deficit                                     
      www (i+1) = www(i+1) - deficit * zdepth(i) / zdepth (i+1)        
 4002 continue
      www(3)    = amax1 (www(3),1.e-12)                                
c                                                                               
      return                                                           
      end                                                              
C                                                                       
C=======================================================================
C                                                                       
      SUBROUTINE DETERM ( DM, DET )                                     
C                                                                       
C====================================================================== 
C                                                                       
C     CALCULATION OF DETERMINANT  OF 3*3 MATRIX                         
C                                                                       
C---------------------------------------------------------------------- 
      DIMENSION DM(3,3)                                                 
C                                                                       
      DET = 0.                                                          
      DET = DM(1,1) * ( DM(2,2) * DM(3,3) - DM(2,3) * DM(3,2) )         
      DET = DM(1,2) * ( DM(2,3) * DM(3,1) - DM(2,1) * DM(3,3) ) + DET   
      DET = DM(1,3) * ( DM(2,1) * DM(3,2) - DM(2,2) * DM(3,1) ) + DET   
C                                                                       
      RETURN                                                            
      END                                                               
C=======================================================================
C                                                                       
      SUBROUTINE SNOWM (MDLSNO,ISNOW,wfsoil,swe)                        
C                                                          1 AUGUST 1988
C=======================================================================
C                                                                       
C     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES          
C     N.B. THIS VERSION DEALS WITH REFREEZING OF WATER                  
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in' 
C                                                                       
      DO 1000 IVEG = 1, 2                                               
C     
CS    Sun Add following part for snow melting and water flux to soil (wfsoil)
CS is  greater zero                                  10/13/98
      IF (ISNOW.eq.0.and.IVEG.eq.2) THEN
          ZMELT=wfsoil
          WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) ) 
          CAPAC(2)=SWE
          GO TO 1000 
      END IF 
CS                                                   10/13/98
      CCT = CCX                                                         
      TS = TC                                                           
      DTS = DTC                                                         
      FLUX = CHF                                                        
      IF ( IVEG .EQ. 1 ) GO TO 100                                      
      CCT = CG                                                          
      TS = TGS                                                          
      DTS = DTG                                                         
      FLUX = CCT * DTG / DTT                                            
                                                                        
100   CONTINUE                                                          
C                                                                       
      TTA = TS - DTS                                                    
      TTB = TS                                                          
      SNOWW(IVEG) = 0.                                                  
      IF ( TTA .LE. TF ) SNOWW(IVEG) = CAPAC(IVEG)                      
      CAPAC(IVEG) = CAPAC(IVEG) - SNOWW(IVEG)                           
      IF ( TTA .GT. TF .AND. TTB .GT. TF ) GO TO 200                    
      IF ( TTA .LE. TF .AND. TTB .LE. TF ) GO TO 200                    
C                                                                       
      DTF = TF - TTA                                                    
      DTIME1 = CCT * DTF /  FLUX                                        
      HF = FLUX*(DTT-DTIME1)                                            
      FCAP = - CAPAC(IVEG)  * SNOMEL                                    
      SPWET = AMIN1( 5. , SNOWW(IVEG) )                                 
      IF ( DTS .GT. 0. ) FCAP =  SPWET * SNOMEL                         
      DTIME2 = FCAP / FLUX                                              
      DTF2 =   FLUX * (DTT-DTIME1-DTIME2)/CCT                           
      TN = TF + DTF2                                                    
      TS = TF - 0.1                                                     
      IF (ABS(HF) .GE.ABS(FCAP) ) TS = TN                               
      CHANGE = HF                                                       
      IF (ABS(CHANGE) .GE.ABS(FCAP) ) CHANGE = FCAP                     
C                                                                       
      CHANGE = CHANGE / SNOMEL                                          
      SNOWW(IVEG) = SNOWW(IVEG) - CHANGE                                
      CAPAC(IVEG) = CAPAC(IVEG) + CHANGE                                
C                                                                       
      IF ( IVEG .EQ. 1 ) TC = TS                                        
      IF ( IVEG .EQ. 2 ) TGS = TS                                       
      IF ( SNOWW(IVEG) .LT. 0.00001 ) GO TO 200                         
c     ZMELT = 0.                                                        
c     modeified to force water into soil Xue Feb. 1994
      ZMELT = CAPAC(IVEG)                             
c     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
c     IF ( TD .LE. TF ) ROFF = ROFF + CAPAC(IVEG)                      
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
      CAPAC(IVEG) = 0. 
200   CONTINUE                                                          
C                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
C                                                                       
1000  CONTINUE                                                          
CS Sun changes  following statatement which is alwayes functioned 
CS    in Xues code            10/13/98 
      IF (ISNOW.ne.0.or.MDLSNO.ne.0) THEN           
         FLUXEF = SHF - CCT*DTG/DTT                    
         TD = TD + FLUXEF / ( CG * 2. * SQRT ( PIE*365. ) ) * DTT 
      END IF 
CS   10/13/98      
C                                                                       
      RETURN                                                            
      END                                                               
C=======================================================================
C                                                                       
      SUBROUTINE STOMA1                                                 
C                                                         19 DECEMB 1988
C=======================================================================
C                                                                       
C     CALCULATION OF PAR-LIMITED STOMATAL RESISTANCE                    
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in' 
C                                                                       
      DO 1000 IVEG = 1, 2                                               
C                                                                       
      AT = ZLT(IVEG) / VCOVER(IVEG)                                     
C                                                                       
      IF (SUNANG .LE. 0.02) THEN                                        
         XABC = RSTPAR(IVEG,1) / RSTPAR(IVEG,2) + RSTPAR(IVEG,3)        
         RST(IVEG) = 0.5 / XABC * AT                                    
         IF (RST(IVEG) .LT. 0.) RST(IVEG) = 0.00001                     
         GO TO 1010                                                     
      END IF                                                            
C                                                                       
      GAMMA = ( RSTPAR(IVEG,1) + RSTPAR(IVEG,2) * RSTPAR(IVEG,3) ) /    
     *          RSTPAR(IVEG,3)                                          
C                                                                       
      POWER1 = AMIN1( 50., AT * EXTK(IVEG,1,1) )                        
      POWER2 = AMIN1( 50., AT * EXTK(IVEG,1,2) )                        
C                                                                       
C-----------------------------------------------------------------------
C     ROSS INCLINATION FUNCTION                                         
C-----------------------------------------------------------------------
C                                                                       
      AA = 0.5 - 0.633 * CHIL(IVEG)- 0.33 * CHIL(IVEG)* CHIL(IVEG)      
      BB = 0.877 * ( 1. - 2. * AA )                                     
C                                                                       
C-----------------------------------------------------------------------
C     COMBINED ESTIMATE OF K-PAR USING WEIGHTS FOR DIFFERENT COMPONENTS 
C-----------------------------------------------------------------------
C                                                                       
      ZAT = ALOG( ( EXP(-POWER1) + 1. )/2. ) * PD(IVEG)                 
     &      / ( POWER1/AT )                                             
      ZAT = ZAT + ALOG( ( EXP(-POWER2) + 1. )/2. )                      
     & * ( 1. - PD(IVEG) ) / ( POWER2/AT )                              
C                                                                       
      POW1 = AMIN1( 50., (POWER1*ZAT/AT) )                              
      POW2 = AMIN1( 50., (POWER2*ZAT/AT) )                              
C                                                                       
      ZK = 1. / ZAT * ALOG( PD(IVEG) * EXP ( POW1 )                     
     &      + ( 1. - PD(IVEG) ) * EXP ( POW2 ) )                        
C                                                                       
C                                                                       
      POW = AMIN1( 50., ZK*AT )                                         
      EKAT = EXP ( POW )                                                
C                                                                       
      AVFLUX = PAR(IVEG) * ( PD(IVEG) / SUNANG * ( AA + BB * SUNANG )   
     &      + ( 1. - PD(IVEG) )*( BB / 3. + AA * 1.5                    
     &      + BB / 4. * PIE ))                                          
C                                                                       
      RHO4 = GAMMA / AVFLUX                                             
C                                                                       
      RST(IVEG) = RSTPAR(IVEG,2)/GAMMA * ALOG(( RHO4 * EKAT + 1. ) /    
     *              ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) - ALOG (( RHO4 + 1. / EKAT ) /              
     *              ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) / ( ZK * RSTPAR(IVEG,3) )                   
C                                                                       
C---------------------------------------------------------------------- 
C     MODIFICATIONS FOR GREEN FRACTION : RST UPRIGHT                    
C---------------------------------------------------------------------- 
C                                                                       
1010  RST(IVEG) = 1. / ( RST(IVEG) * GREEN(IVEG) + 0.0000001)           
1000  CONTINUE                                                          
C                                                                       
      rst(1) = rst(1) * ctlpa
      RETURN                                                            
      END                                                               
C=======================================================================
C                                                                       
      SUBROUTINE TEMRS1(MDLSNO,ISNOW,rank,replicate,pixel,ictrl,icrash)
C                                                          19 JULY 1989 
C=======================================================================
C     A SIMPLIFIED VERSION (XUE ET AL. 1991)                            
C     CORE ROUTINE: CALCULATION OF CANOPY AND GROUND TEMPERATURE        
C     INCREMENTS OVER TIME STEP, FLUXES DERIVED.                        
C-----------------------------------------------------------------------
C                                                                       
C     SUBROUTINES IN THIS BLOCK : TEMRS1                                
C     -------------------------   DELRN                                 
C                                 DELHF                                 
C                                 DELEF                                 
C                                 STRES1                                
      include 'comsib.in' 

cm    for debugging purposes...
      integer rank,replicate,pixel

      COMMON / XRIB/ RIB ,temprib                                       
      COMMON/TONNEW/ ZINC(3), A2(3), Y1(3)                              
      COMMON/NEWT/ ITEX(3)                                              
      DIMENSION RSTM(2)

C                                                                       
C---------------------------------------------------------------------- 
C     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE     
C     GE(X) IS D E(X) / D ( TEMP )                                      
C---------------------------------------------------------------------- 
C                                                                       
      E(X) = EXP( 21.18123 - 5418. / X ) / .622                         
      GE(X) = EXP( 21.18123 - 5418. / X ) * 5418.                       
     $        / (X*X) / .622                                            
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: top of subroutine,hlat=',
c     &  hlat, 'satcap(1:2)=',satcap(1:2)

      ETC   = E(TC)                                                     
      ETGS  = E(TGS)                                                    
      GETC  = GE(TC)                                                    
      GETGS = GE(TGS)                                                   
      HLAT     = ( 3150.19 - 2.378 * TM ) * 1000.                       
      PSY      = CPAIR / HLAT * PSUR / .622                             
      RCP = RHOAIR * CPAIR                                              
C     RADD = 44.                                                        
C                                                                       
      WC = AMIN1( 1., CAPAC(1)/SATCAP(1) )                              
      WG = AMIN1( 1., CAPAC(2)/SATCAP(2) )                              
C                                                                       
C---------------------------------------------------------------------- 
C      RSOIL FUNCTION FROM FIT TO CAMILLO AND GURNEY (1984) DATA.       
C      WETNESS OF UPPER 0.5 CM OF SOIL CALCULATED FROM APPROXIMATION    
C      TO MILLY FLOW EQUATION WITH REDUCED (1/50 ) CONDUCTIVITY IN      
C      TOP LAYER.                                                       
C---------------------------------------------------------------------- 
C                                                                       
c     WT = WWW(1) + 0.75 * ZDEPTH(1) / ( ZDEPTH(1) + ZDEPTH(2) )        
c    &     * (WWW(1) - (WWW(2)**2)/WWW(1) ) / 2. * 50.                  
c     FAC = AMIN1( WT, 0.99 )                                           
c     FAC = AMAX1( FAC, WWW(1) * 0.1 )                            
C
c------------------------------------------------------------
c *** soil resistance calculation alteration Y.K. Xue Feb. 1994**
c------------------------------------------------------------
c                                  
c      if (ictrl.eq.icrash) print *, 'temrs1: before soil res. calcs'

      FAC = AMIN1( www(1), 0.99 )                                       
      FAC = AMAX1( FAC, 0.02 )                            
      RSOIL =  101840. * (1. - FAC ** 0.0027) 
c 
c------------------------------------------------------------                                                                       
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: tgs divide, tgs=',
c     & TGS

      PSIT = PHSAT * FAC ** (- BEE )                                    
      ARGG = AMAX1(-10.,(PSIT*G/461.5/TGS))                             
      HR = EXP(ARGG)                                                    
C                                                                       
C---------------------------------------------------------------------- 
C     ALTERATION OF AERODYNAMIC TRANSFER PROPERTIES IN CASE OF SNOW     
C     ACCUMULATION.                                                     
C---------------------------------------------------------------------- 
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: adjust aero props for snow'

      RESD = D                                                          
      RESZ0 = Z0                                                        
      RESRDC = RDC                                                      
      RESRBC = RBC                                                      
      RESV2 = VCOVER(2)                                                 
C                                                                       
      IF ( TGS .GT. TF ) GO TO 100                                      
C                                                                       
      SDEP = CAPAC(2) *snden
      SDEP = AMIN1( SDEP, (Z2*0.95) )                                   
      D = Z2 - ( Z2-D ) / Z2 * ( Z2 - SDEP )                            
      Z0 = Z0 / ( Z2-RESD ) * ( Z2-D )                                  
      RDC = RDC * ( Z2-SDEP ) / Z2                                      
      RBC = RBC * Z2 / ( Z2-SDEP )                                      
      VCOVER(2) = 1.                                                    
      WG = AMIN1( 1., CAPAC(2) / 0.004 )                                
      RST(2) = RSOIL                                                    
C                                                                       
100   CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C      CALCULATION OF EA, TA, RA, RB, RD AND SOIL MOISTURE STRESS       
C      FOR THE BEGINNING OF THE TIME STEP                               
C                                                                       
C---------------------------------------------------------------------- 
      IFIRST = 1                                                        
      ICOUNT = 0                                                        
      IONCE = 1                                                         
C                                                                       
      TGEN = TGS                                                        
      TCEN = TC                                                         
      FC = 1.                                                           
      FG = 1.                                                           
      TA = TGS                                                          
CC    TA = TM                                                           
      TRIB = TM                                                         
      EA = EM                                                           
      HT = 0.                                                           
      IONCE = 0                                                         
C                                                                       

1000  CONTINUE                                                          
      ICOUNT = ICOUNT + 1                                               
      CALL RASIT5(TRIB,CTNI,CUNI,FTT,FVV,rank,pixel,replicate,ictrl,0)
      IF ( IFIRST .EQ. 1 ) CALL RBRD1                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
      TA = ( TGS/RD + TC/RB + TM/RA ) / D1                              
      HT = ( TA - TM ) * RCP / RA                                       
      RCC = RST(1)*FC + 2. * RB                                         
      COC = (1.-WC)/RCC + WC/(2.*RB)                                    
      RG = RST(2)*FG
cdk adding if condition to match other code
c     if (mdlsno.eq.0.and.isnow.eq.0) then
c     rsurf = rsoil
c     else
      RSURF = RSOIL*FG
c     endif
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HR     
     &       + VCOVER(2)/(RSURF+RD+44.)*HR                              
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)        
     &       + VCOVER(2)/(RSURF+RD+44.)                                 
      COG1 = COG1 + WG/RD * VCOVER(2)                                   
      COG2 = COG2 + WG/RD * VCOVER(2)                                   
      D2 = 1./RA + COC + COG2                                           
      TOP = COC * ETC + COG1 * ETGS + EM / RA                           
      EA = TOP / D2                                                     
      DROP = AMAX1( 0., (E(TA)-EA) )                                    
C                                                                       
C---------------------------------------------------------------------- 
c
      CALL STRES1 ( IFIRST , RSTM)                                      
C---------------------------------------------------------------------- 
C                                                                       
      IFIRST = 0                                                        
      ERIB = EA                                                         
      TRIB = TA                                                         
CCC                                                                     
      IF ( ICOUNT .LE. 4 ) GO TO 1000                                   
C                                                                       
C---------------------------------------------------------------------- 
C
      CALL DELRN ( RNCDTC, RNCDTG, RNGDTG, RNGDTC )                     
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     DEW CALCULATION : DEW CONDITION IS SET AT BEGINNING OF TIME STEP. 
C     IF SURFACE CHANGES STATE DURING TIME STEP, LATENT HEAT FLUX IS    
C     SET TO ZERO.                                                      
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      IF ( EA .GT. ETC ) FC = 0.                                        
      IF ( EA .GT. ETGS) FG = 0.                                        
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     WET FRACTION EXHAUSTION TEST : IF CAPAC(X) IS EXHAUSTED IN        
C     A TIME STEP, INTERCEPTION LOSS IS LIMITED TO CAPAC(X).            
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C---------------------------------------------------------------------- 
C     START OF NON-NEUTRAL RESISTANCE CALCULATION LOOP                  
C---------------------------------------------------------------------- 
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: before non-neutral calcs'

      I = 0                                                             
C                                                                      
C    ----- INITIALIZE NEWTON-RAPHSON ITERATIVE ROUTINE FOR RASIT 3,5,8  
                    NOX = 0                                             
                 NONPOS = 1                                             
                  IWALK = 0                                             
                     LX = 2                                             
                   FINC = 1.                                            
                   ITEX(LX) = 0.                                        
                   ZINC(LX) = 0.                                        
                   A2(LX)   = 0.                                       
                   Y1(LX)   = 0.                                        
2000  CONTINUE                                                          

      CALL RASIT5(TRIB,CTNI,CUNI,FTT,FVV,rank,pixel,replicate,ictrl,i)

C---------------------------------------------------------------------- 
C                                                                       
      CALL DELHF ( HCDTC, HCDTG, HGDTG, HGDTC ) 
C                                                                       

      CALL DELEF ( ECDTC, ECDTG, EGDTG, EGDTC, DEADTC, DEADTG, EC, EG , 
     &             WC, WG, FC, FG, HR,MDLSNO,ISNOW )                    

C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     CALCULATION OF COEFFICIENTS OF TEMPERATURE TENDENCY EQUATIONS     
C        C - CANOPY                                                     
C        G - GROUND                                                     
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      CCODTC = CCX / DTT - RNCDTC + HCDTC + ECDTC                       
      CCODTG = - RNCDTG + HCDTG + ECDTG                                 
      CCORHS = RADT(1) - ( HC + EC ) / DTT                              
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      GCODTG = CG / DTT + TIMCON*CG*2. - RNGDTG + HGDTG + EGDTG         
      GCODTC = - RNGDTC + HGDTC + EGDTC                                 
      GCORHS = RADT(2) - TIMCON*CG*2. * ( TGS -TD ) - ( HG + EG ) / DTT 
C                                                                       
      DENOM = CCODTC * GCODTG - CCODTG * GCODTC                         
C                                                                       
      DTC = ( CCORHS * GCODTG - CCODTG * GCORHS ) / DENOM               
      DTG = ( CCODTC * GCORHS - CCORHS * GCODTC ) / DENOM               

C                                                                       
C---------------------------------------------------------------------- 
C     CHECK IF INTERCEPTION LOSS TERM HAS EXCEEDED CANOPY STORAGE       
C---------------------------------------------------------------------- 
C                                                                       
C                                                                       
      ECPOT = ( (ETC - EA) + (GETC - DEADTC)*DTC - DEADTG*DTG )         
      ECI = ECPOT * WC /(2.*RB) * RCP/PSY * DTT                         
      ECIDIF=AMAX1(0.0,(ECI-CAPAC(1)*1.E3*HLAT))                        
      ECI   =AMIN1(ECI,(    CAPAC(1)*1.E3*HLAT))                        
C                                                                       
      EGPOT = ( (ETGS - EA) + (GETGS - DEADTG)*DTG - DEADTC*DTC )       
      EGI = EGPOT * VCOVER(2) * WG/RD * RCP/PSY * DTT                   
      EGIDIF=AMAX1(0.0,(EGI-CAPAC(2)*1.E3*HLAT))                        
      EGI   =AMIN1(EGI,(    CAPAC(2)*1.E3*HLAT))
	  													  						                          
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      TGEN = TGS + DTG                                                  
      TCEN = TC + DTC                                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
      TAEN = ( TGEN / RD + TCEN / RB + TM / RA ) / D1                   
C                                                                       
      HEND = ( TAEN - TM ) * RCP / RA + (ECIDIF + EGIDIF)/DTT           
      Y= TRIB - TAEN                                                    
      I = I + 1                                                         
      HT   = HEND                                                       
      IF ( I .GT. ITRUNK ) GO TO 200                                    

C                                                                       

      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX)


      IF(NOX.NE.1)GO TO 2000                                            
C                                                                       
200   CONTINUE                                                          
C     IQIN = IQIN + I                                                   
C     IF (I.GT.10) IQIN1 = IQIN1 + 1                                    
1010  FORMAT(1X,I3,1X,'TR1B,Y,RA,RIB,EGDF',7E11.4)                      
1011  FORMAT(1X,'HEND,HT,Y,TA,TC,TG,ECDF',8E11.4)                       
1012  FORMAT(5X,I10,I5)                                                 
c      WRITE(6,1014)  RIB                                                
1014  FORMAT(5X,F12.5)                                                  
C                                                                       
C---------------------------------------------------------------------- 
C     EXIT FROM NON-NEUTRAL CALCULATION                                 
C                                                                       
C     EVAPOTRANSPIRATION FLUXES CALCULATED FIRST ( J M-2 )              
C---------------------------------------------------------------------- 
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: after non-neutral calcs'

      HRR = HR                                                          
      IF ( FG .LT. .5 ) HRR = 1.                                        
      if (ISNOW.eq.0) then
      RSURF = RSOIL*FG                                                  
      else
      RSURF = RSOIL*FG
      endif
C                                                                       

      COCT = (1.-WC)/RCC                                                
      COGT = VCOVER(2) * (1.-WG)/( RG + RD )                            
      COGS1 = (1.-VCOVER(2)) / ( RD + RSURF ) * HRR                     
     &        + VCOVER(2) / ( RD + RSURF + 44.) * HRR                   

      COGS2 = COGS1 / HRR                                               
C                                                                       
      ECT = ECPOT * COCT * RCP/PSY * DTT                                
C                                                                       
      EGT = EGPOT * COGT * RCP/PSY * DTT                                
      EGS = (ETGS + GETGS*DTG ) * COGS1                                 
     &      - ( EA + DEADTG*DTG + DEADTC*DTC ) * COGS2                  
      EGS = EGS * RCP/PSY * DTT                                         
      EGSMAX = WWW(1) / 2. * ZDEPTH(1) * POROS * HLAT * 1000.           
      EGIADD = AMAX1( 0., EGS - EGSMAX )                                
      EGS = AMIN1 ( EGS, EGSMAX )                                       
      EGIDIF = EGIDIF + EGIADD                                          
C                                                                       
C---------------------------------------------------------------------- 
C     SENSIBLE HEAT FLUX CALCULATED WITH LATENT HEAT FLUX CORRECTION    
C---------------------------------------------------------------------- 
      HC = HC + (HCDTC*DTC + HCDTG*DTG)*DTT + ECIDIF                    
      HG = HG + (HGDTC*DTC + HGDTG*DTG)*DTT + EGIDIF                    
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     TEST OF DEW CONDITION. LATENT HEAT FLUXES SET TO ZERO IF SIGN     
C     OF FLUX CHANGES OVER TIME STEP.EXCESS ENERGY DONATED TO SENSIBLE  
C     HEAT FLUX.                                                        
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      ECF = SIGN( 1., ECPOT )                                           
      EGF = SIGN( 1., EGPOT )                                           
      DEWC = FC * 2. - 1.                                               
      DEWG = FG * 2. - 1.                                               
C                                                                       
      IF(DEWC*ECF.GT.0.0) GO TO 300                                     
      HC = HC + ECI + ECT                                               
      ECI = 0.                                                          
      ECT = 0.                                                          
300   IF(DEWG*EGF.GT.0.0) GO TO 400                                     
      HG = HG + EGS + EGI + EGT                                         
      EGS = 0.                                                          
      EGI = 0.                                                          
      EGT = 0.                                                          
400   CONTINUE                                                          
C                                                                       
      EC = ECI + ECT                                                    
      EG = EGT + EGS + EGI                                              

C                                                                       
C---------------------------------------------------------------------- 
C     ADJUSTMENT OF TEMPERATURES AND VAPOR PRESSURE , CALCULATION OF    
C     SENSIBLE HEAT FLUXES.                                             
C---------------------------------------------------------------------- 
C                                                                       
c      if (ictrl.eq.icrash) print *, 'temrs1: before sensible calcs'

      TC  = TCEN                                                        
      TGS = TGEN                                                        
      TA  = TAEN                                                        
      EA = EA + DEADTC*DTC + DEADTG*DTG                                 
C                                                                       
      RADT(1) = RADT(1) + RNCDTC*DTC + RNCDTG*DTG                       
      RADT(2) = RADT(2) + RNGDTC*DTC + RNGDTG*DTG 
c
c ** simulated net all-wave radiation **
c
c     sibnet(nmm,ndd,nhh) = RADT(1) + RADT(2)                   
C                                                                       
      CHF = CCX / DTT * DTC                                             
      SHF = CG / DTT * DTG + TIMCON*CG*2. * ( TGS - TD )                
C                                                                       
      ZLWUP = ZLWUP - RNCDTC * DTC / 2.                                 
     &              - RNGDTG * DTG * (1.-VCOVER(1)*(1.-THERMK) )        
C                                                                       
      IF ( TGS .GT. TF ) GO TO 500                                      
      EGS = EG - EGI                                                    
      EGT = 0.                                                          
500   CONTINUE                                                          
C                                                                       
      VCOVER(2) = RESV2                                                 
      D = RESD                                                          
      Z0 = RESZ0                                                        
      RDC = RESRDC                                                      
      RBC = RESRBC                                                      
C                                                                       
      RETURN                                                            
      END       
C 
C=======================================================================
CS Sun Change original SUBROUTINE TEMRS1 into  TEMRS2  on 10/13/98
CS                            
      SUBROUTINE TEMRS2(MDLSNO,ISNOW,CHISL,tsoil,solsoil,meas,
     &  CSOIL,dzsoil,wfsoil,ictrl,pixel,y0,n_y,replicate,rank) 
CS                                        10/13/98  
c                                     
C                                                          19 JULY 1989 
C=======================================================================
C     A SIMPLIFIED VERSION (XUE ET AL. 1991)                            
C     CORE ROUTINE: CALCULATION OF CANOPY AND GROUND TEMPERATURE        
C     INCREMENTS OVER TIME STEP, FLUXES DERIVED.                        
C-----------------------------------------------------------------------
C                                                                       
C     SUBROUTINES IN THIS BLOCK : TEMRS1                                
C     -------------------------   DELRN                                 
C                                 DELHF                                 
C                                 DELEF                                 
C                                 STRES1    
      include 'comsib.in'
CS  Sun add following statements          10/13/98
      include 'snow4.in'
c sun Adds Local variables               
      real work(nd),work1(nd),delth(nd)
      data delth/nd*0.0/
      integer n_y,replicate,rank,meas
      real y0(n_y)
CS                                       10/13/98 
      COMMON / XRIB/ RIB ,temprib                                       
      COMMON/TONNEW/ ZINC(3), A2(3), Y1(3)                              
      COMMON/NEWT/ ITEX(3)                                              
      DIMENSION RSTM(2)

cm add definition for pixel as an integer
      integer pixel
C                                                                       
C---------------------------------------------------------------------- 
C     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE     
C     GE(X) IS D E(X) / D ( TEMP )                                      
C---------------------------------------------------------------------- 
C

      E(X) = EXP( 21.18123 - 5418. / X ) / .622                         
      GE(X) = EXP( 21.18123 - 5418. / X ) * 5418.                       
     $        / (X*X) / .622                                            
C                                                                       
      ETC   = E(TC)                                                     
      ETGS  = E(TGS)                                                    
      GETC  = GE(TC)                                                    
      GETGS = GE(TGS)                                                   
      HLAT     = ( 3150.19 - 2.378 * TM ) * 1000.                       
      PSY      = CPAIR / HLAT * PSUR / .622                             
      RCP = RHOAIR * CPAIR                                              
C     RADD = 44. 
      WC = AMIN1( 1., CAPAC(1)/SATCAP(1) )  
CS  SUN CHANGE foolowing  statement to one new 10/13/98 
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
        WG=1.0
      ELSE 
        WG = AMIN1( 1., CAPAC(2)/SATCAP(2) )
      END IF                              

CS    on 10/13/98                                                                       
C---------------------------------------------------------------------- 
C      RSOIL FUNCTION FROM FIT TO CAMILLO AND GURNEY (1984) DATA.       
C      WETNESS OF UPPER 0.5 CM OF SOIL CALCULATED FROM APPROXIMATION    
C      TO MILLY FLOW EQUATION WITH REDUCED (1/50 ) CONDUCTIVITY IN      
C      TOP LAYER.                                                       
C---------------------------------------------------------------------- 
C                                                                       
c     WT = WWW(1) + 0.75 * ZDEPTH(1) / ( ZDEPTH(1) + ZDEPTH(2) )        
c    &     * (WWW(1) - (WWW(2)**2)/WWW(1) ) / 2. * 50.                  
c     FAC = AMIN1( WT, 0.99 )                                           
c     FAC = AMAX1( FAC, WWW(1) * 0.1 )                            
C
c------------------------------------------------------------
c *** soil resistance calculation alteration Y.K. Xue Feb. 1994**
c------------------------------------------------------------
c                                  
      FAC = AMIN1( www(1), 0.99 )                                       
      FAC = AMAX1( FAC, 0.02 )  
CS Sun fixed following RSOIL equation as equal to     10/13/98
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN                         
        RSOIL=10000000000.
      ELSE 
        RSOIL =  101840. * (1. - FAC ** 0.0027)
      END IF 
CS                                                   10/13/98 
c------------------------------------------------------------                                                                       
C                                                                       
      PSIT = PHSAT * FAC ** (- BEE )                                    
      ARGG = AMAX1(-10.,(PSIT*G/461.5/TGS)) 
      HR = EXP(ARGG)                                                   
C                                                                       
C---------------------------------------------------------------------- 
C     ALTERATION OF AERODYNAMIC TRANSFER PROPERTIES IN CASE OF SNOW     
C     ACCUMULATION.                                                     
C---------------------------------------------------------------------- 
C                                                                       
      RESD = D                                                          
      RESZ0 = Z0                                                        
      RESRDC = RDC                                                      
      RESRBC = RBC                                                      
      RESV2 = VCOVER(2)                                                 
C                        


                                               
      IF ( TGS .GT. TF ) GO TO 100                                      
CS Sun Change following statement into another one: SDEP=snowdepth  10/13/98  
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN    
         SDEP=snowdepth
      ELSE 
         SDEP = CAPAC(2) *snden
      END IF
	  
	    
CS                                    10/13/98 
      SDEP = AMIN1( SDEP, (Z2*0.95) )                                   
      D = Z2 - ( Z2-D ) / Z2 * ( Z2 - SDEP )                            
      Z0 = Z0 / ( Z2-RESD ) * ( Z2-D )                                  
      RDC = RDC * ( Z2-SDEP ) / Z2                                      
      RBC = RBC * Z2 / ( Z2-SDEP )                                      
      VCOVER(2) = 1. 
CS  Sun change the WG to WG=1.0       10/13/98                
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN      
        WG=1.0
      ELSE
        WG = AMIN1( 1., CAPAC(2) / 0.004 )   
      END IF   
CS                                    10/13/98                          
      RST(2) = RSOIL                                                    
C                                                                       
100   CONTINUE                                                          
C                                                                             
C---------------------------------------------------------------------- 
C                                                                       
C      CALCULATION OF EA, TA, RA, RB, RD AND SOIL MOISTURE STRESS       
C      FOR THE BEGINNING OF THE TIME STEP                               
C                                                                       
C---------------------------------------------------------------------- 
      IFIRST = 1                                                        
      ICOUNT = 0                                                        
      IONCE = 1                                                         
C                                                                       
      TGEN = TGS                                                        
      TCEN = TC                                                         
      FC = 1.                                                           
      FG = 1.                                                           
      TA = TGS                                                          
CC    TA = TM                                                           
      TRIB = TM                                                         
      EA = EM                                                           
      HT = 0.                                                           
      IONCE = 0                                                         
C                                                                       
1000  CONTINUE                                                          
      ICOUNT = ICOUNT + 1                                               
      CALL RASIT5(TRIB,CTNI,CUNI,FTT,FVV,rank,pixel,replicate,ictrl,0)
      IF ( IFIRST .EQ. 1 ) CALL RBRD1                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
      TA = ( TGS/RD + TC/RB + TM/RA ) / D1                              
      HT = ( TA - TM ) * RCP / RA                                       
      RCC = RST(1)*FC + 2. * RB                                         
      COC = (1.-WC)/RCC + WC/(2.*RB)                                    
      RG = RST(2)*FG 
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
           RSURF = RSOIL
      ELSE 
           RSURF = RSOIL*FG
      END IF
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HR     
     &       + VCOVER(2)/(RSURF+RD+44.)*HR                              
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)        
     &       + VCOVER(2)/(RSURF+RD+44.)                                 
      COG1 = COG1 + WG/RD * VCOVER(2)                                   
      COG2 = COG2 + WG/RD * VCOVER(2)                                   
      D2 = 1./RA + COC + COG2                                           
      TOP = COC * ETC + COG1 * ETGS + EM / RA                           
      EA = TOP / D2                                                     
      DROP = AMAX1( 0., (E(TA)-EA) )                                    
C                                                                       
C---------------------------------------------------------------------- 
c
      CALL STRES1 ( IFIRST , RSTM)                                      
C---------------------------------------------------------------------- 
C                                                                       
      IFIRST = 0                                                        
      ERIB = EA                                                         
      TRIB = TA                                                         
CCC                                                                     
      IF ( ICOUNT .LE. 4 ) GO TO 1000                                   
C                                                                       
C---------------------------------------------------------------------- 
C
      CALL DELRN ( RNCDTC, RNCDTG, RNGDTG, RNGDTC )                     
C                                                                             
C----------------------------------------------------------------------       
C                                                                       
C     DEW CALCULATION : DEW CONDITION IS SET AT BEGINNING OF TIME STEP. 
C     IF SURFACE CHANGES STATE DURING TIME STEP, LATENT HEAT FLUX IS    
C     SET TO ZERO.                                                      
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      IF ( EA .GT. ETC ) FC = 0.                                        
      IF ( EA .GT. ETGS) FG = 0.                                        
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     WET FRACTION EXHAUSTION TEST : IF CAPAC(X) IS EXHAUSTED IN        
C     A TIME STEP, INTERCEPTION LOSS IS LIMITED TO CAPAC(X).            
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C---------------------------------------------------------------------- 
C     START OF NON-NEUTRAL RESISTANCE CALCULATION LOOP                  
C---------------------------------------------------------------------- 
C                                                                       
      II = 0                                                            
C                                                                      
C    ----- INITIALIZE NEWTON-RAPHSON ITERATIVE ROUTINE FOR RASIT 3,5,8  
                    NOX = 0                                             
                 NONPOS = 1                                             
                  IWALK = 0                                             
                     LX = 2                                             
                   FINC = 1.                                            
                   ITEX(LX) = 0.                                        
                   ZINC(LX) = 0.                                        
                   A2(LX)   = 0.                                       
                   Y1(LX)   = 0.   
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
**********************************************************
c Next loop, we calculate the thermal conductivities   ****
c and specific heat                                    ****
c**********************************************************
         CALL  TPROPTY(CHISL,dzsoil)     

cm write qk and dsol values to file
cm       write(23,*) dsol(1),dsol(2),dsol(3),qk(1),qk(2),qk(3)

c*************************************************************
c Next we calculate the balances of energy and water    ******
c*************************************************************
c
         tssn(n+1) = tkair
c
c*************************************************************
         icount = 0
         do i=1,n
           work(i)   = tssno(i)
           work1(i)  = dliqvol(i)
         end do
         hx    = 0.0
         NK=n
      ELSE 
         NK=1
      END IF     
         RADDWN=solsoil
         RADDWN=RADDWN+dsol(1)+dsol(2)
         RNG = RADT(2) - RADDWN
         RADT(2)=RNG             

      do 57 ik = NK , 1 , -1
ccccc Next calculate snow layers temperatures and densities      
        IF (ISNOW.ne.0) go to 2000    
         If((sso(ik).lt.1d0.and.porosity(ik).gt.0d0))then
            udum0 = dzo(ik)*(porosity(ik) -work1(ik))
			
	
			
            if(udum0.lt.0.0)then
              print *, ' udum0 is WRONG in thermal.f,obsnow=',obsnow
              Stop
            end if
			
			
		
			
            if(wf(ik+1).gt.udum0)then
               uuu=udum0
               snroff = snroff + (wf(ik+1)-udum0)
c       write (100,*)'O udum0 wf(ik+1) snroff',udum0,wf(ik+1),snroff
      hroff=hroff+(wf(ik+1)-udum0)*cl*rhowater*(tssn(ik+1)-273.16)
            else
            uuu=wf(ik+1)
            endif
            dhp(ik+1)=(uuu*cl*rhowater*(tssn(ik+1)-273.16))/dtt
            w(ik)=wo(ik)+ uuu



            bwo(ik)=rhowater*w(ik)/dzo(ik)
            cto(ik)=(bwo(ik)/920.0)*1.9e+6
            dmlto(ik)=w(ik)*dlm*rhowater
            if (ho(ik).lt.-dmlto(ik)) then
			

			
               fio(ik)=1.0
               flo(ik)=0.0 
               tssno(ik)=( ho(ik)+dmlto(ik))/(cto(ik)*dzo(ik))+273.16
	
	   		   
c	  print*,'ho(',ik,')',ho(ik),'dmlto',dmlto(ik),ictrl,tssno(ik)
	  
				   			   
            else 
               tssno(ik)=273.16
               fio(ik)=-ho(ik)/dmlto(ik) 
               flo(ik)=1.0-fio(ik)
			   
			   
c	  print*,'ho(',ik,')',ho(ik),'dmlto',dmlto(ik),ictrl,tssno(ik)
			   
            end if 
            blo(ik)=bwo(ik)*flo(ik)
            bio(ik)=bwo(ik)*fio(ik) 
            dliqvol(ik)=blo(ik)/rhowater
            dicevol(ik)=bio(ik)/dice
         Else
             w(ik)=wo(ik)
            snroff = snroff +wf(ik+1)
c           write (100,*)'H wf(ik+1) snroff',wf(ik+1),snroff
      hroff=hroff+wf(ik+1)*cl*rhowater*(tssn(ik+1)-273.16)
            dhp(ik+1) = 0.0
         End if


cs  Sun add. It is important because tssno(n) is changed here on 1/25/99 .
         TGS=tssno(NK)
cs  0n 1/25/99
c*************************************************************
         If (ik.lt.Nk) Then
c Next:  ik < n 
c
           if(ik.ne.1) then
             b1 = dsol(ik) + delth(ik)
     &          + qk(ik+1)*(tssn(ik+1)-tssno(ik)) + qk(ik)*work(ik-1)

           else
             b1 = dsol(ik) + delth(ik)
     &          + qk(ik+1)*(tssn(ik+1)-tssno(ik)) + qk(ik)*tsoil
           endif
c
             b2 = - qk(ik)
c Important: delth(ik) must be initialized after using.
           delth(ik) = 0.0
c
         End if
c
         dmlt(ik)=w(ik)*dlm*rhowater
         If (ik.lt.NK.and.ik.ge.1) Then 
            fff = -( ho(ik) + (b1+b2*273.16)*dtt ) 
     &          / ( rhowater*dlm*w(ik) )
c when snow temperature equals 273.16
cc
            if(fff.gt.0.0.and.fff.le.1.0) then
                 ICASE=2
                 fi(ik)=fff
            else if (fff.gt.1.0) then
                 ICASE=1
            else if (fff.le.0.0) then
                 ICASE=3
            end if
        End if
        If (ik.lt.NK) go to 3000   
c
CS     Sun add  above paragraph        on 10/13/98            
2000  CONTINUE      



C                                                                       
      CALL RASIT5(TRIB,CTNI,CUNI,FTT,FVV,rank,pixel,replicate,ictrl,ii)
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      CALL DELHF ( HCDTC, HCDTG, HGDTG, HGDTC )   
C                                                                       
      CALL DELEF ( ECDTC, ECDTG, EGDTG, EGDTC, DEADTC, DEADTG, EC, EG , 
     &             WC, WG, FC, FG, HR,MDLSNO,ISNOW ) 
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     CALCULATION OF COEFFICIENTS OF TEMPERATURE TENDENCY EQUATIONS     
C        C - CANOPY                                                     
C        G - GROUND                                                     
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       


      CCODTC = CCX / DTT - RNCDTC + HCDTC + ECDTC                       
      CCODTG = - RNCDTG + HCDTG + ECDTG                                 
      CCORHS = RADT(1) - ( HC + EC ) / DTT                              
C                                                                       
C---------------------------------------------------------------------- 
CS Sun Change following original GCODCG into new one     10/13/98
      IF (ISNOW.eq.0) THEN 
          GCODTG= cto(n)*dzo(n)/DTT - RNGDTG + HGDTG + EGDTG + qk(n)
      ELSE 
          GCODTG = CG / DTT + TIMCON*CG*2. - RNGDTG + HGDTG + EGDTG 
      END IF 
      GCODTC = - RNGDTC + HGDTC + EGDTC  
CS  From NOW ON WE REALLY GET INTO SNOW PART    !!!!. ON 10/13/98 

      IF (MDLSNO.ne.0.or.ISNOW.ne.0) THEN
        GCORHS = RADT(2)-TIMCON*CG*2.*( TGS -TD )-( HG + EG )/ DTT
      ELSE  
         fi(n)=1.0
         GCORHS1 = ho(n)/DTT+RNG - ( HG + EG ) / DTT +dhp(n+1)
cs   &         - qk(n)*(tssno(n)-tssno(n-1))-cto(n)*dzo(n)*(tssno(n)-273.16)/DTT
     &   - qk(n)*(TGS-tssno(n-1))-cto(n)*dzo(n)*(tssno(n)-273.16)/DTT
         GCORHS = GCORHS1+ rhowater*dlm*w(n)*fi(n)/DTT
      END IF 
C                                                                       
      DENOM = CCODTC * GCODTG - CCODTG * GCODTC                         
C                                                                       

      DTC = ( CCORHS * GCODTG - CCODTG * GCORHS ) / DENOM               
      DTG = ( CCODTC * GCORHS - CCORHS * GCODTC ) / DENOM 
	  

CS  Sun add following part here for inserting snow routing on 10/13/98 
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN 
        If ((TGS+DTG).le.273.16) Then
          TGSNEW=(TGS+DTG)
          ICASE=1
cs Sun debug on 1998/12/14  end        
          h(NK)=( TGSNEW-273.16)*cto(n)*dzo(n)-fi(NK)*w(NK)*dlm*rhowater
          Dh_DTT_DTG=(h(NK)-ho(NK))/DTT
          tonm1=tssno(NK-1)
          qkn=qk(n) 
        Else
          DTG=273.16-TGS
          DTC= (CCORHS - CCODTG*DTG)/CCODTC
      fi(NK)=(GCODTC*DTC+GCODTG*DTG-GCORHS1)/(rhowater*dlm*w(n))*DTT
          if (fi(NK).ge.0.0.and.fi(NK).le.1.0) then 
             h(NK)=-fi(n)*w(n)*dlm*rhowater
             Dh_DTT_DTG=(h(NK)-ho(NK))/DTT
             ICASE=2
             tonm1=tssno(NK-1)
             qkn=qk(n) 
          else if (fi(NK).lt.0.)then
             h(NK)= -fi(NK)*w(NK)*dlm*rhowater
             Dh_DTT_DTG=(h(NK)-ho(NK))/DTT
             ICASE=3
             tonm1=tssno(NK-1)
             qkn=qk(n) 
             fff=fi(NK)
             fi(NK)=0.0
          end if
        End if
      END IF
      


C---------------------------------------------------------------------- 
C     CHECK IF INTERCEPTION LOSS TERM HAS EXCEEDED CANOPY STORAGE       
C---------------------------------------------------------------------- 
C                                                                       
      ECPOT = ( (ETC - EA) + (GETC - DEADTC)*DTC - DEADTG*DTG )         
      ECI = ECPOT * WC /(2.*RB) * RCP/PSY * DTT                         
      ECIDIF=AMAX1(0.0,(ECI-CAPAC(1)*1.E3*HLAT))                        
      ECI   =AMIN1(ECI,(    CAPAC(1)*1.E3*HLAT))                        
C                                                                       
      EGPOT = ( (ETGS - EA) + (GETGS - DEADTG)*DTG - DEADTC*DTC )  
	  					        
      EGI = EGPOT * VCOVER(2) * WG/RD * RCP/PSY * DTT

      EGIDIF=AMAX1(0.0,(EGI-CAPAC(2)*1.E3*HLAT))                        
      EGI   =AMIN1(EGI,(    CAPAC(2)*1.E3*HLAT))

	  

cm    i think sometimes too much snow evaporates
cm    add a check to make sure less snow evaporates than top layer
cm      EGIDIF=AMAX1(0.0,(EGI-0.9*w(3)*1.E3*HLAT))                        
cm      EGI   =AMIN1(EGI,(    0.9*w(3)*1.E3*HLAT))

C---------------------------------------------------------------------- 
C                                                                       
      TGEN = TGS + DTG                                                  
      TCEN = TC + DTC                                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
      TAEN = ( TGEN / RD + TCEN / RB + TM / RA ) / D1                   
C                                                                       
      HEND = ( TAEN - TM ) * RCP / RA + (ECIDIF + EGIDIF)/DTT           
      Y= TRIB - TAEN                                                    
      II = II + 1                                                       
      HT   = HEND                                                       
      IF ( II .GT. ITRUNK ) GO TO 200                                  
C     

                                                                  
      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX)                      
      IF(NOX.NE.1)GO TO 2000
200   CONTINUE




CS  Sun add following part here for inserting snow routing on 10/13/98 
3000  IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN 
         IF (ICASE.eq.1.and.ik.eq.NK) THEN
            tssn(NK)=TGS+DTG
         END IF  
         If (ik.eq.NK) then 


            SNOFAC = HLAT / (HLAT + SNOMEL /1000.) 
            egidw =  EGI*SNOFAC /HLAT/1000.
            w(n)=w(n)-egidw
            swe=swe-egidw
            dzold=dzo(n)
            dzo(n)=dzo(n)-egidw*rhowater/bwo(n)
cs sun: following way to correct h(n) may lead to unballance of energy.
            ho(n)=ho(n)*dzo(n)/dzold
            capac(2)=swe 
			
			
			
         End if 
		 
   

      CALL snresult(ik,ICASE,qsoil,wfsoil,tsoil,
     *    b1,b2,fff,delth,ictrl,pixel,y0,n_y,replicate,rank)

      END IF
 57   continue  
 
		 
		     
     
CS sun add following parts on 12/5/98   start 
clwp  11/17/2000, Li add following sentence to recalculate the snowdepth
c      SNOWDEPTH=DZO(1)+DZO(2)+DZO(3)
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN 

			
        swe=w(1)+w(2)+w(3)
        capac(2)=swe
        if((dz(1)+dz(2)+dz(3)).ne.0.0) then
        snden=(bw(1)*dz(1)+bw(2)*dz(2)+bw(3)*dz(3))
     *               /(dz(1)+dz(2)+dz(3))
        snden=1000./snden
        endif
			
c		endif
		
	
		
      ENDIF
	  
		  

CS sun add above  parts on 12/5/98  end 
C     IQIN = IQIN + I                                                   
C     IF (I.GT.10) IQIN1 = IQIN1 + 1                                    
1010  FORMAT(1X,I3,1X,'TR1B,Y,RA,RIB,EGDF',7E11.4)                      
1011  FORMAT(1X,'HEND,HT,Y,TA,TC,TG,ECDF',8E11.4)                       
1012  FORMAT(5X,I10,I5)                                                 
c      WRITE(6,1014)  RIB                                                
1014  FORMAT(5X,F12.5)                                                  
C                                                                       
C---------------------------------------------------------------------- 
C     EXIT FROM NON-NEUTRAL CALCULATION                                 
C                                                                       
C     EVAPOTRANSPIRATION FLUXES CALCULATED FIRST ( J M-2 )              
C---------------------------------------------------------------------- 
C
      HRR = HR
      IF ( FG .LT. .5 ) HRR = 1.
cs SUn change RSURF = RSOIL*FG into followings: 02/03/99 start
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
           RSURF = RSOIL
      ELSE 
           RSURF = RSOIL*FG
      END IF
cs sun 03/02/99  end   
C                                                                        
      COCT = (1.-WC)/RCC                                                
      COGT = VCOVER(2) * (1.-WG)/( RG + RD )                            
      COGS1 = (1.-VCOVER(2)) / ( RD + RSURF ) * HRR                     
     &        + VCOVER(2) / ( RD + RSURF + 44.) * HRR                   
      COGS2 = COGS1 / HRR                                               
C                                                                       
      ECT = ECPOT * COCT * RCP/PSY * DTT                                
C                                                                       
      EGT = EGPOT * COGT * RCP/PSY * DTT 
      EGS = (ETGS + GETGS*DTG ) * COGS1                                 
     &      - ( EA + DEADTG*DTG + DEADTC*DTC ) * COGS2                  
      EGS = EGS * RCP/PSY * DTT 
CS  Sun add following statement on 10/13/98
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) EGS=0.0
CS                                 10/13/98                                        
      EGSMAX = WWW(1) / 2. * ZDEPTH(1) * POROS * HLAT * 1000.           
      EGIADD = AMAX1( 0., EGS - EGSMAX )                                
      EGS = AMIN1 ( EGS, EGSMAX )                                       
      EGIDIF = EGIDIF + EGIADD                                          
C                                                                       
C---------------------------------------------------------------------- 
C     SENSIBLE HEAT FLUX CALCULATED WITH LATENT HEAT FLUX CORRECTION    
C---------------------------------------------------------------------- 
      HC = HC + (HCDTC*DTC + HCDTG*DTG)*DTT + ECIDIF                    
      HG = HG + (HGDTC*DTC + HGDTG*DTG)*DTT + EGIDIF 
C---------------------------------------------------------------------- 
C                                                                       
C     TEST OF DEW CONDITION. LATENT HEAT FLUXES SET TO ZERO IF SIGN     
C     OF FLUX CHANGES OVER TIME STEP.EXCESS ENERGY DONATED TO SENSIBLE  
C     HEAT FLUX.                                                        
C                                                                       
C---------------------------------------------------------------------- 
C 
cs Sun add following one statement IF (ISNOW.eq.0.and.MDLSNO.eq.0) go to 
cs 401  CONTINUE  to skip folloing statements from
CS  ECF = SIGN( 1., ECPOT ) to 400   CONTINUE 
      IF (ISNOW.eq.0.and.MDLSNO.eq.0) go to 401
      ECF = SIGN( 1., ECPOT )                                           
      EGF = SIGN( 1., EGPOT )                                           
      DEWC = FC * 2. - 1.                                               
      DEWG = FG * 2. - 1.                                               
C                                                                       
      IF(DEWC*ECF.GT.0.0) GO TO 300
      HC = HC + ECI + ECT                                               
      ECI = 0.                                                          
      ECT = 0.                                                          
300   IF(DEWG*EGF.GT.0.0) GO TO 400                                     
      HG = HG + EGS + EGI + EGT                                         
      EGS = 0.                                                          
      EGI = 0.                                                          
      EGT = 0.                                                          
400   CONTINUE
401   CONTINUE
C                                                                       
      EC = ECI + ECT                                                    
      EG = EGT + EGS + EGI      
C                                                                       
C---------------------------------------------------------------------- 
C     ADJUSTMENT OF TEMPERATURES AND VAPOR PRESSURE , CALCULATION OF    
C     SENSIBLE HEAT FLUXES.                                             
C---------------------------------------------------------------------- 
C                                             
cs sun add following new statement 02/04/99
      TGSOLD=TGS
cs sun end                               
      TC  = TCEN                                                        
      TGS = TGEN 
CS Sun add following statement: 10/13/98
      IF (ISNOW.eq.0) tssn(n)=TGS
CS                              10/13/98  
      TA  = TAEN                                                        
      EA = EA + DEADTC*DTC + DEADTG*DTG                                 
C                                                                       
      RADT(1) = RADT(1) + RNCDTC*DTC + RNCDTG*DTG                       
      RADT(2) = RADT(2) + RNGDTC*DTC + RNGDTG*DTG 
c
c ** simulated net all-wave radiation **
c
c     sibnet(nmm,ndd,nhh) = RADT(1) + RADT(2)                   
C                                                                       
      CHF = CCX / DTT * DTC
cs  sun change the original statement:   on 12/14/98 
cs  SHF = CG / DTT * DTG + TIMCON*CG*2. * ( TGS - TD)
cs into following part where CG / DTT * DTG is rplaced by Dh_DTT_DTG
      IF (ISNOW.eq.0) THEN 
        SHF=  Dh_DTT_DTG - dhp(n+1)+ qkn*(TGSOLD-tonm1)
     &  +qkn*DTG
      ELSE                                             
        SHF = CG / DTT * DTG + TIMCON*CG*2. * ( TGS - TD )                
      END IF
C                                                                       
      ZLWUP = ZLWUP - RNCDTC * DTC / 2.                                 
     &              - RNGDTG * DTG * (1.-VCOVER(1)*(1.-THERMK) )        
C                                                                       
      IF ( TGS .GT. TF ) GO TO 500                                      
      EGS = EG - EGI                                                    
      EGT = 0.                                                          
500   CONTINUE                                                          
C                                                                       
      VCOVER(2) = RESV2                                                 
      D = RESD                                                          
      Z0 = RESZ0                                                        
      RDC = RESRDC                                                      
      RBC = RESRBC 
CS   Sun add next paragrapg to get soil surface temperature TGS  10/13/98
      IF (MDLSNO.eq.0.and.ISNOW.eq.0) THEN
         TGS=tsoil
         atmp= (qsoil+solsoil)/csoil
         btmp=2.*3.1416/86400. 
         ctmp=csoil*btmp/csoil/(365.*3.1416)**0.5
         TGS=(tsoil+atmp*DTT+btmp*DTT*TD/(1.+ctmp*DTT))/
     &       (1.+btmp*DTT*(1.-ctmp*DTT/(1.+ctmp*DTT)))
         TD=(ctmp*DTT*TGS+TD)/(1.+ctmp*DTT)  

cm  trying to fix weird yampa valley pixel ... mike 14 nov 2014
cm         tgs=273.16

      END IF
      RETURN                                                            
      END
C

C====================================================================== 
C                                                                       
      SUBROUTINE STRES1( IFIRST ,RSTM)                                  
C                                                                       
C====================================================================== 
C                                                                       
C     CALCULATION OF ADJUSTMENT TO LIGHT DEPENDENT STOMATAL RESISTANCE  
C     BY TEMPERATURE, HUMIDITY AND STRESS FACTORS                       
C     SIMPLIFIED SEE XUE ET AL(1991)                                    
C                                                                       
C         RSTFAC(IVEG,1) = FD                                           
C         RSTFAC(IVEG,2) = FP                                           
C         RSTFAC(IVEG,3) = FT                                           
C         RSTFAC(IVEG,4) = FTPD                                         
C                                                                       
C---------------------------------------------------------------------- 
      include 'comsib.in' 
      DIMENSION RSTM(2), DEP(3)                                  
C---------------------------------------------------------------------- 
C     HUMIDITY, TEMPERATURE AND TRANSPIRATION FACTORS                   
C---------------------------------------------------------------------- 
C                                                                       
      DO 1000 IVEG = 1, 2                                               
C                                                                       
      TV = TC                                                           
      ETV = ETC                                                         
      RAIR = RB * 2.                                                    
      IF ( IVEG .EQ. 1 ) GO TO 100                                      
      TV = TGS                                                          
      ETV = ETGS                                                        
      RAIR = RD                                                         
100   CONTINUE                                                          
C                                                                       
      TV = AMIN1 ( ( TU(IVEG) - 0.1 ), TV )                             
      TV = AMAX1 ( ( TLL(IVEG) + 0.1 ), TV )                            
C                                                                       
      IF( IFIRST .EQ. 0 ) GO TO 200                                     
      RSTM(IVEG) = RST(IVEG)                                            
      D2 = ( TU(IVEG) - TOPT(IVEG) ) / ( TOPT(IVEG) - TLL(IVEG) )       
      D1 = 1. /(( TOPT(IVEG) - TLL(IVEG) )*                             
     $        EXP( ALOG( TU(IVEG) - TOPT(IVEG))*D2))                    
      RSTFAC(IVEG,3) = D1*( TV-TLL(IVEG)) * EXP(ALOG(TU(IVEG)-TV)*D2)   
C                                                                       
      IF (RSTFAC(IVEG,3).LT.0.) RSTFAC(IVEG,3) = 0.                     
      IF (RSTFAC(IVEG,3).GT.1.) RSTFAC(IVEG,3) = 1.                     
C                                                                       
C                                                                       
C---------------------------------------------------------------------- 
C      SIMPLIFIED CALCULATION OF LEAF WATER POTENTIAL FACTOR , FP       
C---------------------------------------------------------------------- 
C                                                                       
      if (nroot.eq.1) then
      XROT = ROOTD(1)                                                   
      DO 7400 I = 1, 3                                                  
 7400 DEP(I) = 0.                                                       
      DO 7500 I = 1, 3                                                  
      DEP(I) = MIN(ZDEPTH(I), XROT)                                     
      XROT = XROT - ZDEPTH(I)                                           
      IF (XROT.LE.0.) GO TO 7410                                        
 7500 CONTINUE                                                          
 7410 CONTINUE                                                          
      XDR = (PHSOIL(1) * DEP(1) + PHSOIL(2) * DEP(2)                    
     *      +PHSOIL(3) * DEP(3)) /ROOTD(1)                              
      else
      XDR = PHSOIL(1) * rootp(1) + PHSOIL(2) * rootp(2)
     *      +PHSOIL(3) * rootp(3)
      end if
      XDR = - XDR                                                       
      IF (XDR .LE. 0.001) XDR = 0.001                                   
      XDR = ALOG (XDR)                                                  
      EXPONENT = AMAX1(-86.0, (- PH1(IVEG) * (PH2(IVEG) - XDR)) )
      RSTFAC(IVEG,2) = 1. - EXP(EXPONENT)
c     RSTFAC(IVEG,2) = 1. - EXP(- PH1(IVEG) * (PH2(IVEG) - XDR))       
      IF (RSTFAC(IVEG,2).GT.1.) RSTFAC(IVEG,2) = 1.                     
      IF (RSTFAC(IVEG,2).LT.0.) RSTFAC(IVEG,2) = 0.                     
C                                                                       
200   RST(IVEG) = RSTM(IVEG)                                            
C                                                                       
      EPOT = ETV - EA                                                   
      EPOT = AMAX1(0.0001,(ETV-EA)) 
c
c
c               ***** PJS mod 10/9/92 ***** 
c ***** based on Verma FIFE-87 function for C4 grasses *****
c                                    
      rstfac(iveg,1) = 1./ ( 1 + defac(iveg)*drop )                          
c
      IF (RSTFAC(IVEG,1).LT.0.) RSTFAC(IVEG,1) = 0.                     
      IF (RSTFAC(IVEG,1).GT.1.) RSTFAC(IVEG,1) = 1.                     
C                                                                       
C                                                                       
C---------------------------------------------------------------------- 
C     VALUE OF FP FOUND                                                 
C---------------------------------------------------------------------- 
C                                                                       
300   FTPD = RSTFAC(IVEG,1) * RSTFAC(IVEG,2) * RSTFAC(IVEG,3)           
      RSTFAC(IVEG,4) = AMAX1( FTPD, 0.00001 )                           
C---------------------------------------------------------------------- 
C                                                                       
                                                                        
      RST(IVEG) = RST(IVEG) / RSTFAC(IVEG,4) / VCOVER(IVEG)             
C                                                                       
      RST(IVEG) = AMIN1( RST(IVEG), 100000. )                           
C                                                                       
1000  CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                          
c      
cs ********************************************************************
      SUBROUTINE UPDAT1 (MDLSNO,ISNOW,wfsoil,swe,snroff) 
CS ********************************************************************
C                                                                       
C     UPDATING OF SOIL MOISTURE STORES AND INTERCEPTION CAPACITY        
C                                                                       
C-----------------------------------------------------------------------
      include 'comsib.in' 
      DIMENSION EF(3)                                                   
C                                                                       
C---------------------------------------------------------------------- 
C     EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY       
C     ( HLAT*1000.) LOSS IS IN M M-2                                    
C     MASS TERMS ARE IN KG M-2 DT-1                                     
C---------------------------------------------------------------------- 
C                                                                       
      SNOFAC = HLAT / ( HLAT + SNOMEL /1000. )                          
      FACKS = 1.                                                        
      IF ( (TC-DTC) .LE. TF ) FACKS = SNOFAC                            
      IF ( (ECT+ECI) .GT. 0.) GO TO 100                                 
      ECI = ECT + ECI                                                   
      ECT = 0.                                                          
      FACKS = 1. / FACKS                                                
100   CAPAC(1)=CAPAC(1) - ECI*FACKS/HLAT/1000.                          
C                                                                       
      ECMASS = ( ECT + ECI * FACKS ) / HLAT                             
c     ECMASS = ( ECT + ECI ) / HLAT
C     
cs    Sun add following statement  IF (ISNOW.EQ.0)  go to 201 on 12/5/98 
c this next line fixed the water balance
      IF (ISNOW.EQ.0)  FACKS = SNOFAC
      IF (ISNOW.EQ.0)  go to 201   
      FACKS = 1.                                                        
      IF ( (TGS-DTG) .LE. TF ) FACKS = SNOFAC
      IF ( (EGT+EGI) .GT. 0. ) GO TO 200                                
      EGI = EGT + EGI                                                   
      EGT = 0.                                                          
      FACKS = 1. / FACKS  
200   CAPAC(2)=CAPAC(2) - EGI*FACKS/HLAT/1000.
C                                                                       
201   EGMASS = ( EGT + EGS + EGI * FACKS ) / HLAT                       
c201   EGMASS = ( EGT + EGS + EGI ) /HLAT
C                                                                       
      ETMASS = ECMASS + EGMASS                                          
C                                                                       
      HFLUX = ( HC + HG ) / DTT                                         
C                                                                       
C---------------------------------------------------------------------- 
C      DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE            
C---------------------------------------------------------------------- 
C                                                                       
      DO 1000 IVEG = 1, 2                                               
      IF ( CAPAC(IVEG) .GT. 0.000001 ) GO TO 300                        
      WWW(1) = WWW(1) + CAPAC(IVEG) / ( POROS*ZDEPTH(1) )               
      CAPAC(IVEG) = 0.
300   CONTINUE                                                          
1000  CONTINUE   
C---------------------------------------------------------------------- 
C     SNOWMELT / REFREEZE CALCULATION                                   
C---------------------------------------------------------------------- 
CS  Sun Change following  CALL SNOWM to  SNOWM (ISNOW,wfsoil,swe) 
CS   10/13/98  
      CALL SNOWM (MDLSNO,ISNOW,wfsoil,swe)   
CS 10/13/98 
C---------------------------------------------------------------------- 
C     BARE SOIL EVAPORATION LOSS                                        
C---------------------------------------------------------------------- 
C                                                                       
      WWW(1) = WWW(1) - EGS / HLAT / 1000. / ( POROS * ZDEPTH(1) )      
C                                                                       
C---------------------------------------------------------------------- 
C   EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE                     
C---------------------------------------------------------------------- 
C                                                                       
      DO 2000 IVEG = 1, 2                                               
C                                                                       
      IF ( IVEG .EQ. 1 ) ABSOIL = ECT / HLAT / 1000.                    
      IF ( IVEG .EQ. 2 ) ABSOIL = EGT / HLAT / 1000.                    
C                                                                       
      EF(2) = 0.                                                        
      EF(3) = 0.                                                        
      TOTDEP = ZDEPTH(1)                                                
C                                                                       
      DO 3000 IL = 2, 3                                                 
      TOTDEP = TOTDEP + ZDEPTH(IL)                                      
C                                                                       
C     DIV = AMAX1 ( 1., ( PHSOIL(IL) - PHL(IVEG) ) )                    
C                                                                       
      IF ( ROOTD(IVEG) .LT. TOTDEP ) GO TO 400                          
C                                                                       
      EF(IL) = ZDEPTH(IL) / ROOTD(IVEG)                                 
      GO TO 500                                                         
C                                                                       
400   CONTINUE                                                          
      EF(IL) = ROOTD(IVEG) - TOTDEP + ZDEPTH(IL)                        
      EF(IL) = EF(IL) / ROOTD(IVEG)                                     
      GO TO 600                                                         
C                                                                       
500   CONTINUE                                                          
3000  CONTINUE                                                          
C                                                                       
600   EFT = EF(2) + EF(3)                                               
      EF(2) = EF(2) / EFT                                               
      EF(3) = EF(3) / EFT                                               
C                                                                       
      DO 4000 IL = 2, 3                                                 
      WWW(IL) = WWW(IL) - ABSOIL * EF(IL) / ( POROS * ZDEPTH(IL) )
4000  CONTINUE                                                          
C                                                                       
2000  CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C     CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO         
C     GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .       
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      DO 5000 IL = 1, 2                                                 
      IF ( WWW(IL) .GT. 0. ) GO TO 700                                  
      WWW(IL+1) = WWW(IL+1) + WWW(IL) * ZDEPTH(IL)/ZDEPTH(IL+1)         
      WWW(IL) = 0.                                                      
700   CONTINUE                                                          
5000  CONTINUE
C                                                                       
c     IF ( TD .LT. TF ) GO TO 800                                       
C 
      call run2
C 
800   CONTINUE                                                          
C                                                                       
      IF (WWW(1) .GT.1.) WWW(2) = WWW(2) + (WWW(1)-1.) * ZDEPTH(1)/     
     &                             ZDEPTH(2)                            
      IF (WWW(1) .GT.1.) WWW(1) = 1.                                    
      IF (WWW(2) .GT.1.) WWW(3) = WWW(3) + (WWW(2)-1.) * ZDEPTH(2)/     
     &                             ZDEPTH(3)                            
      IF (WWW(2) .GT.1.) WWW(2) = 1.                                    
      IF (WWW(3) .GT.1.) ROFF   = ROFF + (WWW(3)-1.)*POROS*ZDEPTH(3)    
c this line added new
c     if (www(3).gt.1) rnoffb(nobs) = rnoffb(nobs) + (www(3)-1.)*
c    & poros*zdepth(3)
      IF (WWW(3) .GT.1.) WWW(3) = 1.                                    
C                       
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VEGIN1(IVEG_TYPE)
C                                                     17 NOVEMBER 1990
C=======================================================================
C
C     READ VEGETATION PHYSIOLOGY
C
C    SURFACE PARAMETERS ARE READ IN SAME ORDER AS IN GCM
C    SUBROUTINE SIBINP. ONLY EXCEPTION IS THAT 1-D VERSION READS IN
C    SITE SPECIFIC PARAMETERS CORB1 ... ZMET .
C
C     VARIABLES THAT ENTER THROUGH COMSIB:
C        SUBSCRIPTS (IV, IW, IL) :
C              IV = VEGETATION STORY; 1 = TOP AND 2 = BOTTOM
C              IW = RADIATION WAVELENGTH; 1 = VISIBLE, 2 = NEAR
C                   INFRARED AND 3 = THERMAL INFRARED
C              IL = VEGETATION STATE; 1 = LIVE (GREEN) AND
C                   2 = DEAD (STEMS AND TRUNK)
C
C   TRAN(IV,IW,IL): LEAF TRANSMITTANCE
C   REF (IV,IW,IL): LEAF REFLECTANCE
C   RSTPAR(IV,IW) : PAR-DEPENDENT LEAF STOMATAL RESISTANCE COEFFICIENTS
C                          A =(J/M**3) B = 2(W/M**2) C = 3(S/M)
C   SOREF(IW)     : SOIL REFLECTANCE
C   CHIL(IV)      : LEAF ANGLE DISTRIBUTION FACTOR
C   TOPT(IV)      : OPTIMUM TEMPERATURE FOR STOMATAL FUNCTIONING
C   TL(IV)        : LOWER TEMPERATURE LIMIT FOR STOMATAL FUNCTIONING
C   TU(IV)        : UPPER TEMPERATURE LIMIT FOR STOMATAL FUNCTIONING
C   DEFAC(IV)     : VAPOR PRESSURE DEFICIT PARAMETER
C   PH1(IV)       : coefficient for soil water effect - mike
C   PH2(IV)       : soil water potential for wilting point - mike
C   ROOTD(IV)     : ROOTING DEPTH
C   BEE           : SOIL WETNESS EXPONENT
C   PHSAT         : SOIL TENSION AT SATURATION
C   SATCO         : HYDRAULIC CONDUCTIVITY AT SATURATION
C   POROS         : SOIL POROSITY
C   ZDEPTH        : DEPTH OF 3 SOIL MOISTURE LAYERS
C   Z0            : ROUGHNESS LENGTH
C   D             : ZERO PLANE DISPLACEMENT
C   ZLT(IV)       : LEAF AREA INDEX
C   GREEN(IV)     : GREEN LEAF FRACTION
C   VCOVER(IV)    : VEGETATION COVER FRACTION
C
C     VARIABLES ( SPECIFIC TO SIB 1-D VERSION ONLY ) FROM COMSIB
C
C      ZWIND  : REFERENCE HEIGHT FOR WIND MEASUREMENT
C      ZMET   : REFERENCE HEIGHT FOR TEMPERATURE, HUMIDITY MEASUREMENT
C        THE ABOVE ARE GENERATED FROM SIBX + MOMOPT OUTPUT
C
C
C-----------------------------------------------------------------------

cm These data are now passed to ssib subroutine in a single array, vegin1.
cm   This was done so that the model might be initialized without reading
cm   from a file each time, since the model will be called many thousands
cm   of times on each of several processors. -Mike, Aug 2005

      include 'comsib.in'
cm      real data_in(58)

cm      READ(1,*)
cm     READ(1,*) (TRAN(1,IW,1), IW=1,3), (TRAN(1,IW,2), IW=1,3)
cm    READ(1,*) (TRAN(2,IW,1), IW=1,3), (TRAN(2,IW,2), IW=1,3)
cm      READ(1,*) (REF (1,IW,1), IW=1,3), (REF (1,IW,2), IW=1,3)
cm      READ(1,*) (REF (2,IW,1), IW=1,3), (REF (2,IW,2), IW=1,3)
cm      READ (1,*) (RSTPAR(1,IWAVE), IWAVE=1,3),
cm     & (RSTPAR(2,IWAVE), IWAVE=1,3)
cm      READ (1,*) (SOREF(IWAVE), IWAVE=1,3),(CHIL(IV), IV=1,2)
cm      READ (1,*) (TOPT(IV), IV=1,2),
cm     & (TLL(IV),  IV=1,2), (TU(IV),   IV=1,2)
cm      READ (1,*) (DEFAC(IV),IV=1,2), (PH1(IV),  IV=1,2),
cm     & (PH2(IV),  IV=1,2)
cm      READ (1,*) (ROOTD(IV),IV=1,2),BEE,PHSAT
cm      READ (1,*) SATCO, POROS, SLOPE
cm      READ (1,*) (ZDEPTH(IDEP), IDEP=1,3)
cm      READ(1,*) ZWIND

c      j=1
c      do iv=1,2
c        do iw=1,3
c          do il=1,2
c            TRAN(IV,IW,IL)=data_in(j)
c            j=j+1
c          end do
c        end do
c      end do
c      do iv=1,2
c        do iw=1,3
c          do il=1,2
c            REF(IV,IW,IL)=data_in(j)
c            j=j+1
c          end do
c        end do
c      end do
c      do iv=1,2
c        do iw=1,3
c          RSTPAR(IV,IW)=data_in(j)
c          j=j+1
c        end do
c      end do
c      do iw=1,3
c        SOREF(IW)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        CHIL(IV)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        TOPT(IV)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        TLL(IV)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        TU(IV)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        DEFAC(IV)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        ph1(iv)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        ph2(iv)=data_in(j)
c        j=j+1
c      end do
c      do iv=1,2
c        rootd(iv)=data_in(j)
c        j=j+1
c      end do
c      BEE=data_in(j)
c      j=j+1
c      PHSAT=data_in(j)
c      j=j+1
c      SATCO=data_in(j)
c      j=j+1
c      POROS=data_in(j)
c      j=j+1
c      SLOPE=data_in(j)
c      j=j+1
c      do idep=1,3
c        ZDEPTH(IDEP)=data_in(j)
c        j=j+1
c      end do 
c      ZWIND=data_in(j)

      do iv=1,2
        do iw=1,3
          do il=1,2
            TRAN(IV,IW,IL)=TRAN0(IVEG_TYPE,IV,IW,IL)
          end do
        end do
      end do
      do iv=1,2
        do iw=1,3
          do il=1,2
            REF(IV,IW,IL)=REF0(IVEG_TYPE,IV,IW,IL)
          end do
        end do
      end do
      do iv=1,2
        do iw=1,3
          RSTPAR(IV,IW)=RSTPAR0(IVEG_TYPE,IV,IW)
        end do
      end do
      do iw=1,3
        SOREF(IW)=SOREF0(IVEG_TYPE,IW)
      end do
      do iv=1,2
        CHIL(IV)=CHIL0(IVEG_TYPE,IV)
      end do
      do iv=1,2
        TOPT(IV)=TOPT0(IVEG_TYPE,IV)
      end do
      do iv=1,2
        TLL(IV)=TL0(IVEG_TYPE,IV)
      end do
      do iv=1,2
        TU(IV)=TU0(IVEG_TYPE,IV)
      end do
      do iv=1,2
        DEFAC(IV)=DEFAC0(IVEG_TYPE,IV)
      end do
      do iv=1,2
        ph1(iv)=PH10(IVEG_TYPE,IV)
      end do
      do iv=1,2
        ph2(iv)=PH20(IVEG_TYPE,IV)
      end do
      do iv=1,2
        rootd(iv)=ROOTD0(IVEG_TYPE,IV)
      end do
      BEE=BEE0(IVEG_TYPE)
      PHSAT=PHSAT0(IVEG_TYPE)
      SATCO=SATCO0(IVEG_TYPE)
      POROS=POROS0(IVEG_TYPE)
cm      SLOPE=SOPE0(IVEG_TYPE)
      SLOPE=0.
      do idep=1,3
        ZDEPTH(IDEP)=DEPTH0(IVEG_TYPE,IDEP)
      end do 
cm    NOTE: ZWIND defined in vegin.f
cm      ZWIND=data_in(j)

      RETURN
      END
c
C========================================================================
C
      SUBROUTINE VEGIN2(IVEG_TYPE,mymonth,f_veg,vcov_dat)
C                                                     12 OCTOBER 1993
C-----------------------------------------------------------------------
      include 'comsib.in'

      dimension f_veg(1,1)

cm These data are now passed to ssib subroutine in a single array, vegin1.
cm   This was done so that the model might be initialized without reading
cm   from a file each time, since the model will be called many thousands
cm   of times on each of several processors. -Mike, Aug 2005

cm The vcov_dat variable represents external data for the vegetation cover
cm   fraction that is assigned to the upper storey VCOVER variable instead
cm   of the default values from the SSiB input file. - Mike, Aug 2006 

cm      real data_in(12,12)
        real LAIlocal 
C
cm      READ(2,*)
cm      READ(2,*) (ZLT(IV), IV=1,2),(GREEN(IV), IV=1,2),Z2,Z1
cm      READ(2,*) (VCOVER(IV), IV=1,2),Z0,D
cm      READ(2,*) RBC,  RDC
c

c      ZLT(1)=data_in(mymonth,1)
c      ZLT(2)=data_in(mymonth,2)
c      GREEN(1)=data_in(mymonth,3)
c      GREEN(2)=data_in(mymonth,4)
c      Z2=data_in(mymonth,5)
c      Z1=data_in(mymonth,6)
c      VCOVER(1)=data_in(mymonth,7)
c      VCOVER(2)=data_in(mymonth,8)
c      Z0=data_in(mymonth,9)
c      D=data_in(mymonth,10)
c      RBC=data_in(mymonth,11)
c      RDC=data_in(mymonth,12)

c     multiply leaf area index values by f_veg
      ZLT(1)=min(ZLT0(IVEG_TYPE,mymonth,1)*exp(f_veg(1,1)),8.25)

c Manu FEbruary 2011 change this to be ZLT=ZLT*VCOVER because it seems that ZLT
c                    is considered to be the area averaged LAI
c leave this commented for now, and read in original ssib vegin files...
c see Veg analysis on my blog for more details...
      !print*,'readZLT',ZLT0(IVEG_TYPE,mymonth,1) 
      !ZLT(1)=min(ZLT0(IVEG_TYPE,mymonth,1)*exp(f_veg(1,1))
      !&           VCOVER(1)
  
      
      ZLT(2)=min(ZLT0(IVEG_TYPE,mymonth,2)*exp(f_veg(1,1)),8.25)
      GREEN(1)=GREEN0(IVEG_TYPE,mymonth,1)
      GREEN(2)=GREEN0(IVEG_TYPE,mymonth,2)
      Z2=Z20(IVEG_TYPE,mymonth)
      Z1=Z10(IVEG_TYPE,mymonth)
c      VCOVER(1)=VCOVER0(IVEG_TYPE,mymonth,1)
      VCOVER(1)=vcov_dat/100.
      VCOVER(2)=VCOVER0(IVEG_TYPE,mymonth,2)
      Z0=Z00(IVEG_TYPE,mymonth)
      D=D0(IVEG_TYPE,mymonth)
      RBC=RBC0(IVEG_TYPE,mymonth)
      RDC=RDC0(IVEG_TYPE,mymonth)

      RETURN
      END
c
c ========================================================================
c **********      Snowmodel is inserted in the following.     ************
c ========================================================================
c ========================================================================
      subroutine getinput
      include 'snow4.in' 
c.......................
c       ssisnow commented by mike to allow it to be read in
c       ssisnow  = 0.04
       n        = 3
       flmin    = 0.03
       flmax    = 0.10
       dzmin    = 0.002
       womin    = 0.0004
       cl       = 4212.7
       dlm      = 3.335d5
       rhowater = 1000.0
       dice     = 920.0
       return 
       end
c***********************************************************************
c GETMET reads in met data from file, and generates required met
c parameters (OLD)
cm  now it doesn not read anything, it performs a few basic snow calculations
c***********************************************************************
      subroutine getmet(iptype,UM,nmm,ndd,nhh,ictrl,pixel)
      include 'snow4.in'
      integer pixel

c     assign snowfall density, bifall 
      IF (prcp.gt.1.E-9) THEN
        if(iptype.eq.2)then
cm        iptype == 2 corresponds to snowfall
          prcps=prcp
          prcpw=0.0
cm        compute bifall and dfall
          if(tkair.gt.258.16)then
          Ab=1.4*(278.15-tkair)**(-1.15)+0.008*UM**1.7
          bifall=500*(1.-0.951*exp(-Ab))
          else
          Ab=0.008*UM**1.7
          bifall=500.*(1-0.904*exp(-Ab))
          end if
cm TEST//
cm          bifall=100.
cm // TEST

cm        new snow grain diameter


        IF((ictrl.gt.849).and.(ictrl.lt.969))then
      
          if(bifall.gt.1000.)then
            ddfall=0.0
          elseif(bifall.gt.400.)then
            ddfall=1.0*(sg3+1.1E-13*400**sg4)
          else
            ddfall=1.0*(sg3+1.1E-13*bifall**sg4)
          end if

        ELSE       

          if(bifall.gt.1000.)then
            ddfall=0.0
          elseif(bifall.gt.400.)then
            ddfall=sg3+1.1E-13*400**sg4
          else
            ddfall=sg3+1.1E-13*bifall**sg4
          end if
     
        ENDIF

        else if(iptype.eq.1)then
cm        rainfall       
          prcpw=prcp
          prcps=0.0
          flfall=1.0
          fifall=0.
          blfall=1000.0
          bifall=0.0
          ddfall=0.0
        endif
      ELSE 
cm      no precipitation at all  
        prcpw=0.0
        prcps=0.0
        iptype = 0
        return
      END IF 

CM THIS PART ADDED 9/15/05 BASED ON SUB2.F OFFLINE SNOW CODE
       IF(IPTYPE.NE.1)THEN
         flfall=0.
         fifall=1.
         blfall=0.
       END IF

       return
       end
c ***************************************************************************
c    LAYERN divides the one layer information into n layers information after
c    snowdepth greater  than sd_cr( criticla depth to determine whether snow
c    layer is seperated from soil body or not.
cm   LAYERN is called in two places - when the model initialises as well as 
cm     any time when the snow depth increases from z<sd_cr to z>sd_cr.  
cm     the 'swtch' variable is set to 0 for initializing and 1 for 
cm     transitioning past the critical depth
c  **************************************************************************
       subroutine layern (tg,swtch,im,id,ih)
cm     swtch variable added by mike oct 03
       integer swtch
       include 'snow4.in'
       include 'comsib.in'
c         n=3
          IF(snowdepth.gt.0.05.and.snowdepth.le.0.06) THEN     
             dzo(1)=0.02
             dzo(2)=0.02
             dzo(3)=snowdepth- dzo(1)- dzo(2)
          ELSE IF ( snowdepth.gt.0.06.and.snowdepth.le.0.08) then
             dzo(3)=0.02
             dzo(2)=0.02
             dzo(1)=snowdepth- dzo(3)- dzo(2)
          ELSE IF ( snowdepth.gt.0.08.and.snowdepth.le.0.62) then
             dzo(3)=0.02
             dzo(2)=(snowdepth- dzo(3))*0.33333333 
             dzo(1)=(snowdepth- dzo(3))*0.66666667
          ELSE IF ( snowdepth.gt.0.62) then
             dzo(3)=0.02
             dzo(2)=0.20
             dzo(1)=snowdepth- dzo(3)- dzo(2)        
         End IF 
cm       changed by mike to allow initial conditions to be read in at model start
         do  777 i=1,n
             if (swtch.eq.0)then
cm             when we are initializing the model, use the temp and density read
cm                in cntrols subroutine	  
                 tssno(i)=tssn(i)
                 bwo(i)=bw(i)
             else if (swtch.eq.1)then
cm         when the depth is transitioning past critical point, use the 1-layer
cm                model results
                 if(tg.lt.273.16)then
                   tssno(i)=tg
                 else
                   tssno(i)=273.16
                 end if
                 bwo(i)=swe*rhowater/snowdepth  
                 gdia(i)=gsize
             end if
 777     continue
c*********************************************************************
c  Next we will calculate the initial variables for time step going on
c*********************************************************************
          do 666 i=1,n
cm      mike is adding this line to make initialization yield non-zero wo values
          if (dzo(i).eq.0.) dzo(i)=dz(i)

          wo(i)=(bwo(i)*dzo(i))/rhowater
          cto(i)=(bwo(i)/920.0)*1.9e+6
          bto(i)=bwo(i)
          IF (tssno(i).eq.273.16)then
cm              changed by mike to allow fl to be read in initially
                if (swtch.eq.0.)then
cm                when we are initializing the model, use the water content read
cm                 in cntrols subroutine	  
                   flo(i)=fl(i)
                else if (swtch.eq.1.)then
cm              when the depth is transitioning past critical point, assume the 
cm                  moisture is equal to the minimum moisture parameter
                    flo(i)= flmin
                end if                              
              fio(i)=1.0- flo(i)
              ho(i)=(-1.0)*wo(i)*fio(i)*dlm*rhowater
              blo(i)=bwo(i)*flo(i)
              bio(i)=bwo(i)*fio(i)
              dliqvol(i)=blo(i)/rhowater
              dicevol(i) = bio(i)/dice
          ELSE IF(tssno(i).lt.273.16) THEN
              flo(i)=0.0
              fio(i)=1.0
              dmlto(i)=wo(i)*dlm*rhowater
              ho(i)=(tssno(i)-273.16)*cto(i)*dzo(i)-dmlto(i)
              blo(i)=0.0
              bio(i)=bwo(i)
              dliqvol(i)=0.0
              dicevol(i) = bio(i)/dice
          ENDIF
666     continue
        return
        end
c ************************************************************************
c    LAYER1 combine n layers informations into one layer information after
c    snowdepth less than sd_cr
c************************************************************************
      subroutine layer1 (csoil,tgs,dzsoil,h,w,snowdepth,
     *           swe,stemp,nd,gsize,gdia)
       parameter (dice=917.0, rhowater=1000.0,dlm=3.335d5)
       dimension h(nd),w(nd),gdia(nd)
       swe=w(1)+w(2)+w(3)
C      YX2002 (test2)
       snh=h(1)+h(2)+h(3)+csoil*(tgs-273.16)
c      snh=h(1)+h(2)+h(3)+csoil*(tgs-273.16)*dzsoil
       dmlto=swe*dlm*rhowater
       scv=1.9e+6*(swe/snowdepth)/dice
       if (snh.gt.0.0) then
C      YX2002 (test2)
          stemp=snh/(swe*4.18*10**6.+csoil)+273.16
c         stemp=snh/(swe*4.18*10**6.+csoil*dzsoil)+273.16
       else if (snh.gt.-dmlto) then
          stemp=273.16
       else 
C      YX2002 (test2)
          stemp=(snh+dmlto)/(scv*snowdepth+csoil)+273.16
c         stemp=(snh+dmlto)/(scv*snowdepth+csoil*dzsoil)+273.16
       end if 
cm     compute GSIZE from mass-weighted average
       gsize=(w(1)*gdia(1)+w(2)*gdia(2)+w(3)*gdia(3))/
     &        (w(1)+w(2)+w(3))
       return
       end
c **********************************************************:
      subroutine modnodesun
      include 'snow4.in'
      IF (snowdepth.le.0.06) goto 100
      IF (snowdepth.le.0.08) goto 200
      IF (snowdepth.le.0.62) goto 300
      IF (snowdepth.gt.0.62) goto 400
 100  CONTINUE
      DZ2=0.02
      DZ3=DZ2
      DZ1=snowdepth-( DZ2+DZ3)
      DDZ1=DZ1-dzo(1)
      IF (DDZ1.GT.0.0) THEN
          DDZ1=MIN(DDZ1,dzo(2))
          CALL COMBO (DDZ1,dzo(1),dzo(2),wo(1),wo(2),ho(1),ho(2),
     &    tssno(1),tssno(2),bwo(1),bwo(2),bio(1),bio(2),blo(1),blo(2),
     &    bto(1),bto(2),fio(1),fio(2),flo(1),flo(2),cto(1),cto(2),
     &    dliqvol(1),dliqvol(2),dicevol(1),dicevol(2))
      ELSE
          DDZ1=-DDZ1
           CALL COMBO (DDZ1,dzo(2),dzo(1),wo(2),wo(1),ho(2),ho(1),
     &    tssno(2),tssno(1),bwo(2),bwo(1),bio(2),bio(1),blo(2),blo(1),
     &    bto(2),bto(1),fio(2),fio(1),flo(2),flo(1),cto(2),cto(1),
     &    dliqvol(2),dliqvol(1),dicevol(2),dicevol(1))
      END IF
      SUN23=dzo(2)+dzo(3)
      DZ2=0.5*SUN23
      DZ3= DZ2
      DDZ2=DZ2-dzo(2)
      IF (DDZ2.GT.0.0) THEN
          CALL COMBO (DDZ2,dzo(2),dzo(3),wo(2),wo(3),ho(2),ho(3),
     &    tssno(2),tssno(3),bwo(2),bwo(3),bio(2),bio(3),blo(2),blo(3),
     &    bto(2),bto(3),fio(2),fio(3),flo(2),flo(3),cto(2),cto(3),
     &    dliqvol(2),dliqvol(3),dicevol(2),dicevol(3))
      ELSE
          DDZ2=-DDZ2
          CALL COMBO (DDZ2,dzo(3),dzo(2),wo(3),wo(2),ho(3),ho(2),
     &    tssno(3),tssno(2),bwo(3),bwo(2),bio(3),bio(2),blo(3),blo(2),
     &    bto(3),bto(2),fio(3),fio(2),flo(3),flo(2),cto(3),cto(2),
     &    dliqvol(3),dliqvol(2),dicevol(3),dicevol(2))
      END IF
      GO to 500
 200  CONTINUE
      DZ1=0.02
      DDZ1=DZ1-dzo(1)
      IF (DDZ1.GT.0.0) THEN
          DDZ1=MIN(DDZ1,dzo(2))
          CALL COMBO (DDZ1,dzo(1),dzo(2),wo(1),wo(2),ho(1),ho(2),
     &    tssno(1),tssno(2),bwo(1),bwo(2),bio(1),bio(2),blo(1),blo(2),
     &    bto(1),bto(2),fio(1),fio(2),flo(1),flo(2),cto(1),cto(2),
     &    dliqvol(1),dliqvol(2),dicevol(1),dicevol(2))
      ELSE
          DDZ1=-DDZ1
           CALL COMBO (DDZ1,dzo(2),dzo(1),wo(2),wo(1),ho(2),ho(1),
     &    tssno(2),tssno(1),bwo(2),bwo(1),bio(2),bio(1),blo(2),blo(1),
     &    bto(2),bto(1),fio(2),fio(1),flo(2),flo(1),cto(2),cto(1),
     &    dliqvol(2),dliqvol(1),dicevol(2),dicevol(1))
      END IF
      SUN23=dzo(2)+dzo(3)
      DZ2=0.02
      DDZ2=DZ2-dzo(2)
      IF (DDZ2.GT.0.0) THEN
          CALL COMBO (DDZ2,dzo(2),dzo(3),wo(2),wo(3),ho(2),ho(3),
     &    tssno(2),tssno(3),bwo(2),bwo(3),bio(2),bio(3),blo(2),blo(3),
     &    bto(2),bto(3),fio(2),fio(3),flo(2),flo(3),cto(2),cto(3),
     &    dliqvol(2),dliqvol(3),dicevol(2),dicevol(3))
      ELSE
          DDZ2=-DDZ2
          CALL COMBO (DDZ2,dzo(3),dzo(2),wo(3),wo(2),ho(3),ho(2),
     &    tssno(3),tssno(2),bwo(3),bwo(2),bio(3),bio(2),blo(3),blo(2),
     &    bto(3),bto(2),fio(3),fio(2),flo(3),flo(2),cto(3),cto(2),
     &    dliqvol(3),dliqvol(2),dicevol(3),dicevol(2))
      END IF
      GO to 500
 300  CONTINUE
      DZ1=0.02
      DDZ1=DZ1-dzo(1)
      IF (DDZ1.GT.0.0) THEN
          DDZ1=MIN(DDZ1,dzo(2))
          CALL COMBO (DDZ1,dzo(1),dzo(2),wo(1),wo(2),ho(1),ho(2),
     &    tssno(1),tssno(2),bwo(1),bwo(2),bio(1),bio(2),blo(1),blo(2),
     &    bto(1),bto(2),fio(1),fio(2),flo(1),flo(2),cto(1),cto(2),
     &    dliqvol(1),dliqvol(2),dicevol(1),dicevol(2))
      ELSE
          DDZ1=-DDZ1
           CALL COMBO (DDZ1,dzo(2),dzo(1),wo(2),wo(1),ho(2),ho(1),
     &    tssno(2),tssno(1),bwo(2),bwo(1),bio(2),bio(1),blo(2),blo(1),
     &    bto(2),bto(1),fio(2),fio(1),flo(2),flo(1),cto(2),cto(1),
     &    dliqvol(2),dliqvol(1),dicevol(2),dicevol(1))
      END IF
      SUN23=dzo(2)+dzo(3)
      DZ2=0.33333333*SUN23
      DDZ2=DZ2-dzo(2)
      IF (DDZ2.GT.0.0) THEN
          CALL COMBO (DDZ2,dzo(2),dzo(3),wo(2),wo(3),ho(2),ho(3),
     &    tssno(2),tssno(3),bwo(2),bwo(3),bio(2),bio(3),blo(2),blo(3),
     &    bto(2),bto(3),fio(2),fio(3),flo(2),flo(3),cto(2),cto(3),
     &    dliqvol(2),dliqvol(3),dicevol(2),dicevol(3))
      ELSE
          DDZ2=-DDZ2
          CALL COMBO (DDZ2,dzo(3),dzo(2),wo(3),wo(2),ho(3),ho(2),
     &    tssno(3),tssno(2),bwo(3),bwo(2),bio(3),bio(2),blo(3),blo(2),
     &    bto(3),bto(2),fio(3),fio(2),flo(3),flo(2),cto(3),cto(2),
     &    dliqvol(3),dliqvol(2),dicevol(3),dicevol(2))
      END IF
      GO to 500
 400  CONTINUE
      DZ1=0.02
      DDZ1=DZ1-dzo(1)
      IF (DDZ1.GT.0.0) THEN
          DDZ1=MIN(DDZ1,dzo(2))
          CALL COMBO (DDZ1,dzo(1),dzo(2),wo(1),wo(2),ho(1),ho(2),
     &    tssno(1),tssno(2),bwo(1),bwo(2),bio(1),bio(2),blo(1),blo(2),
     &    bto(1),bto(2),fio(1),fio(2),flo(1),flo(2),cto(1),cto(2),
     &    dliqvol(1),dliqvol(2),dicevol(1),dicevol(2))
      ELSE
          DDZ1=-DDZ1
           CALL COMBO (DDZ1,dzo(2),dzo(1),wo(2),wo(1),ho(2),ho(1),
     &    tssno(2),tssno(1),bwo(2),bwo(1),bio(2),bio(1),blo(2),blo(1),
     &    bto(2),bto(1),fio(2),fio(1),flo(2),flo(1),cto(2),cto(1),
     &    dliqvol(2),dliqvol(1),dicevol(2),dicevol(1))
      END IF
      SUN23=dzo(2)+dzo(3)
      DZ2=0.2
      DDZ2=DZ2-dzo(2)
      IF (DDZ2.GT.0.0) THEN
          CALL COMBO (DDZ2,dzo(2),dzo(3),wo(2),wo(3),ho(2),ho(3),
     &    tssno(2),tssno(3),bwo(2),bwo(3),bio(2),bio(3),blo(2),blo(3),
     &    bto(2),bto(3),fio(2),fio(3),flo(2),flo(3),cto(2),cto(3),
     &    dliqvol(2),dliqvol(3),dicevol(2),dicevol(3))
      ELSE
          DDZ2=-DDZ2
          CALL COMBO (DDZ2,dzo(3),dzo(2),wo(3),wo(2),ho(3),ho(2),
     &    tssno(3),tssno(2),bwo(3),bwo(2),bio(3),bio(2),blo(3),blo(2),
     &    bto(3),bto(2),fio(3),fio(2),flo(3),flo(2),cto(3),cto(2),
     &    dliqvol(3),dliqvol(2),dicevol(3),dicevol(2))
      END IF
 500  Return
      End
C **************************************************************************
      SUBROUTINE COMBO (DDZ2,dzp,dzm,wp,wm,hp,hm,
     &    tp,tm,bwp,bwm,bip,bim,blp,blm,
     &    btp,btm,fip,fim,flp,flm,ctp,ctm,
     &    dliqvolp,dliqvolm,dicevolp,dicevolm)
c
      dice=920.
      rhowater=1000.
      dlm = 3.335d5
      RATIO= DDZ2/dzm
      dzp=dzp + RATIO*dzm
      wp = wp + RATIO*wm
      hp = hp + RATIO*hm
      bwp= wp*rhowater/dzp
      btp= bwp
      ctp= (1.9e6)*(bwp/920.0)
      dmlt=wp*rhowater*dlm
      if(hp.ge.(-1.0)*dmlt)then
         tp=273.16
         fip=(-1.0)*hp/dmlt
         flp=1.0-fip
         blp=bwp*flp
         bip=bwp*fip
         dliqvolp = blp/rhowater
         dicevolp = bip/dice
      else
c
         flp=0.0
         fip=1.0
         tp=(hp+dmlt)/(ctp*dzp)+273.16
         bip=bwp
         blp=0.0
         dliqvolp = 0.0
         dicevolp = bip/dice
      endif
c
      dzm=dzm - RATIO*dzm
      wm = wm - RATIO*wm
      hm = hm - RATIO*hm
      Return
      End
c***********************************************************************
      subroutine old
      include 'snow4.in'
c
      do 20 i=1,n
            tssno(i)=tssn(i)
            bwo(i)=bw(i)
            blo(i)=bl(i)
            bio(i)=bi(i)
            ho(i)=h(i)
            flo(i)=fl(i)
            fio(i)=fi(i)
            wo(i)=w(i)
            dzo(i)=dz(i)
            sso(i)=ss(i)
            cto(i)=ct(i)
            bto(i)=bt(i)
            dmlto(i)=dmlt(i)
 20   continue
      return
      end
c *********************************************************************
      subroutine set0
      include 'snow4.in'
c
      do 100 i=n+1,nd
            tssno(i)=0.0
            bwo(i)=0.0
            blo(i)=0.0
            bio(i)=0.0
            ho(i)=0.0
            flo(i)=0.0
            fio(i)=0.0
            wo(i)=0.0
            dzo(i)=0.0
            sso(i)=0.0
            cto(i)=0.0
            bto(i)=0.0
            dmlto(i)=0.0
100   continue
c.....................
      do 200 i=1,nd
           wf(i)=0.0
           dhp(i)=0.0
200   continue
c
      return
      end
c ************************************************************************        
c SNOWALB calculate the snow albedo based on the scheme from Loth and Graf
c ************************************************************************
       SUBROUTINE snowalb (snowpile,snowdepth,TGS,rhowater,dtt,bifall,
     *                 SUNANG,CLOUD,soilalbd,snowalbd)
c
       parameter(a=1.0,b=1.0,c=-1.3,na=3)
       cosz=SUNANG
       dfall=rhowater*snowpile/bifall
       IF (snowdepth.le.0.001.and.dfall.le.0.001) THEN
            snowalbd=soilalbd
       ELSE  
         IF (snowpile.gt.0.) snowalbd=snowalbd+
     &                      0.1*100.*rhowater*snowpile/bifall
         If (snowdepth.lt.0.25) Then 
            if (TGS.lt.273.16) then
                snowalbd=snowalbd-0.0061*dtt/86400.
            else 
                snowalbd=snowalbd-0.071*dtt/86400. 
            end if
         Else
            snowalbd=0.5+(snowalbd-0.5)*exp(-0.24*dtt/86400.0)
         End if
         sinz=min(cosz,0.866)
         f1=CLOUD*CLOUD
         f2=exp(1.0-(1.0-sinz*sinz))
         f=a*f1+b*f2+c*f1*f2
         snowalbd=snowalbd+(1.0-snowalbd)*snowalbd**na*f
         snowalbd=min(snowalbd,0.85)
       END IF  
       return
       end  
c********************************************************************
      subroutine snow1st (dtt,TM,solsoil,ISNOW,nmm,ndd,nhh,ictrl,pixel,
     &  imike)
      include 'snow4.in'
      integer pixel
c
      tkair  = TM
      prcp   = prcpw+prcps
      snroff = 0.0
      hroff  = 0.0
      dksatsnow=0.01
c  rain
c......

      if(prcpw.gt.0.0)then
           wf(n+1)=amin1(prcpw, dksatsnow*dtt)
           dhp(n+1)=(wf(n+1)/dtt)*cl*rhowater*(tkair-273.16)
           snroff =snroff+(prcpw-wf(n+1) )
           hroff=hroff+(prcpw-wf(n+1))*cl*rhowater*(tkair-273.16)
      else if(prcps.gt.0.0)then 
c  snow, add new nodes
c.....................
           wf(n+1)=0.0
           dhp(n+1)=0.0
           call newsnow(ISNOW,nmm,ndd,nhh,ictrl,pixel)
      endif

c*********************************
c     Grain growth rate for snow
c*********************************
cm    write(22,*) nmm,ndd,nhh,tssno(1),tssno(2),tssno(3)
      call graingrowth(ISNOW)

c*********************************
c     Compaction rate for snow
c*********************************
      do 277 i=n,1,-1
      dicevol(i) = bio(i)/dice
      dliqvol(i) = blo(i)/rhowater
      porosity(i)=1.0-dicevol(i)
      porosity(i)=amin1(porosity(i),1.0)
      porosity(i)=amax1(porosity(i),0.0)
      so(i)=ssisnow
      if(porosity(i).ne.0.0) so(i)=dliqvol(i)/porosity(i)
      sso(i)=(so(i)-ssisnow)/(1.0-ssisnow)
277   continue
      overburden=0.0
      do 377 i=n,1,-1
      overburden=overburden+ wo(i)*rhowater
      call compact(bio(i),tssno(i),blo(i),overburden,pdzdtc(i),
     &     sso(i),dice)
 377  continue

c***********************************************
c    Calculate some variables after new snowfall
c***********************************************
       do 390 i = 1,n
         if((sso(i).lt.1.0.and.porosity(i).gt.0.0))then
           dzot=dzo(i)*(1d0+pdzdtc(i)*dtt)
           dzo(i)=amax1(dzot,dzmin)
c
            if(wo(i).gt.womin)then
               bwo(i)=(wo(i)*rhowater)/dzo(i)
               if (bwo(i).gt.920.0) then
                   bwo(i)=920.0
                   dzo(i)=(wo(i)*rhowater)/bwo(i)
               end if 
            endif
c
            blo(i)=bwo(i)*flo(i)
            bio(i)=bwo(i)*fio(i)
            bto(i)=bwo(i)
        end if
c
        dicevol(i) = bio(i)/dice
        dliqvol(i) = blo(i)/rhowater
        dummy = dliqvol(i) + dicevol(i) 
        if(dummy.gt.1.0)then
          dliqvol(i) = 1.0 - dicevol(i)
          blo(i) = dliqvol(i)*rhowater
          bwo(i) = blo(i) + bio(i)
          dzo(i)=(wo(i)*rhowater)/bwo(i)
        endif
        cto(i)=(bwo(i)/920.0)*1.9e+6 
c
         porosity(i)=1.0-dicevol(i)
         if(porosity(i) .gt. 1.0)porosity(i)=1.0
         if(porosity(i) .lt. 0.0)porosity(i)=0.0
c
         if(porosity(i).gt.0.0)then
           so(i)=dliqvol(i)/porosity(i)
         else
           so(i)=ssisnow
         endif
c
        if(so(i).gt.ssisnow)then
          sso(i)=(so(i)-ssisnow)/(1.0-ssisnow)
        else
          sso(i)=0.0
        endif
cccccc  dmass is for using to calculate dsol in sdsol.f
      dmass(i)=bto(i)*dzo(i)
 390  continue

      snowdepth=dzo(1)+dzo(2)+dzo(3)

c*********************************************
c     Optical parameters and solar extinction
c*********************************************
      IF (solar .gt. 0d0 ) THEN    
        call sdsol(dsol,dmass,n,solar,solsoil,gsize,gdia,isnow)
      ELSE 
        do 112 i=1,n
           dsol(i)=0d0
112     continue
        solsoil=0.0 
      END IF       
c      write(*,*) dsol
      return
      end
c******************************************************************
c     COMPACT calculates the  natural compaction rate of snow cover
c******************************************************************
      subroutine compact(bi,t,bl,overburden,pdzdt,ss,dice)
      data c2,c3,c4,c5/23d-3,2.777d-6,0.04,2.0/
      data dm/150/
      data eta0/0.9d6/
      if(bi  .ge. dice .or. ss .ge. 1.)return
      ddz1=-c3*exp(-c4*(273.15-t))
      if(bi .gt. dm)  ddz1=ddz1*exp(-46.0d-3*(bi-dm))
      if(bl .gt. 0.01)ddz1=ddz1*c5
ccccc compaction due to overburden
      ddz2=-overburden*exp(-0.08*(273.15-t)-c2*bi)/eta0
ccccc compaction occurring during melt has been taken into account in thermal.f
      ddz3=0d0
      pdzdt=ddz1+ddz2+ddz3
      return
      end
c***********************************************************************
      subroutine newsnow(ISNOW,im,id,ih,ictrl,pixel)
      include 'comsib.in'
      include 'snow4.in'
      integer pixel
ccccc calculate rate of change in element thickness due to snow falling
ccccc Precip has just started or previous top node is full. Initiate a new node.
c**********************************************************************

      dzfall=obsnow*3600./1000.*rhowater/bifall

cm    compute new grain diameters
      if(ISNOW>0)then
cm    update one-layer model grain diameter, GSIZE, then return
        
cm       if (ictrl.eq.3252) print *, 'before gsize update...'
        gsize=(ddfall*dzfall*bifall+gsize*swe*rhowater)/
     &        (dzfall*bifall+swe*rhowater)
C        write(*,*) 'no combd called'

cm       if (ictrl.eq.3252) print *, 'after gsize update...'

        return
      else
cm      update three-layer model grain diameters, gdia(nd)
        call combd
      endif

      dzo(n)=dzo(n)+dzfall
      wo(n)=wo(n)+prcp
      bwo(n)=(wo(n)*rhowater)/dzo(n)
      cto(n)=1.9e+6*(bwo(n)/920.0)
      dum=(tkair-273.16)*cto(n)*dzfall
     &        -(1.0-flfall)*(blfall+bifall)*dlm*dzfall
      ho(n)=ho(n)+dum
      dmlto(n)=wo(n)*rhowater*dlm
      if (ho(n). ge. -dmlto(n)) then
        tssno(n)=273.16
        fio(n)=-ho(n)/dmlto(n)
        flo(n)=1.0-fio(n)
        blo(n)=bwo(n)*flo(n)
        bio(n)=bwo(n)*fio(n)
        dliqvol(n)=blo(n)/rhowater
        dicevol(n)=bio(n)/dice
      else
ccccc when snow temperature is below 273.16
        fio(n)=1.0
        flo(n)=0.0
        bio(n)=bwo(n)
        blo(n)=0.0
        dliqvol(n)=0.0
        dicevol(n)=bio(n)/dice
        wf(n)=0.0
        tssno(n)=(ho(n)+dmlto(n))/(cto(n)*dzo(n))+273.16
      end if
cm    add call to modnodnew, which re-layeres snowpack
      if (ISNOW==0) then
        
c       if (ictrl.le.5469) then
        call modnodenew
c       endif

      endif
      return
      end

c*******************************************************************
c     SDSOL computes absorption of solar radiation within snow cover
c*******************************************************************
      subroutine sdsol(dsol,dmass,n,solar,solsoil,gsize,gdia,isnow)
      parameter(nd = 4)
      integer n
      real dsol(nd),dmass(nd),fext(nd),gdia(nd)
c     gsize commented by mike to allow it to be read in
c     gsize   = 5.d-4  
      bext    = 400.0
      cv      = 3.795d-3
      depth   = 30
      do i=1,n
         fext(i) = 0.0
      enddo
c
      tmass = 0.0
      do 10 i=1,n
         j=n+1-i
         tmass=tmass+dmass(j)
         if(tmass.gt.depth) goto 30
cm       compute fext based on grain diameter: gsize or gdia(nd)
         if(isnow>0)then
           fext(j)=exp(-cv*dmass(j)/sqrt(gsize))
         else
           fext(j)=exp(-cv*dmass(j)/sqrt(gdia(j)))
         endif
         if(j .eq. n) fext(n)=exp(-bext*2d-3)*fext(n)
 10   continue
c*************
30    tsolt = solar
      do 20 i=1,n
      j=n+1-i
         if(tsolt .le. 0d0)then
            dsol(j)=0d0
            tsolb=0.0
         else
            tsolb=tsolt*fext(j)
            dsol(j)=tsolt-tsolb
            tsolt=tsolb
         end if
 20   continue
      solsoil = tsolb
      return
      end
c*******************************************************************
c     TPROPTY calculate the thermal conductivities and specific heat
c*******************************************************************
      subroutine tpropty(thksoil,dzsoil)
      include 'snow4.in' 
ccccc this is thermal conductivity for snow from R.Jordan(1991)(2.4)
      do 37 i=1,n
         thkice=2.290d0
         thkair=2.30d-2
         thk(i) = thkair+(7.75d-5 *bwo(i)+ 1.105d-6*
     &            bwo(i)*bwo(i))*(thkice -thkair)+0.1
 37   continue
ccccc calculate the ratio of thermal conductivity 
ccccc at the ineterface between two layers(2.7)
      do 47 i=2,n
      qk(i)=2.0*thk(i)*thk(i-1)/(thk(i)*dzo(i-1)+thk(i-1)*dzo(i))
47    continue
C     YX2002 (test2) but do nothing at this stage
      qk(1)= 2.0*thk(1)*thksoil/(thk(1)*dzsoil+thksoil*dzo(1))
      return
      end
c**************************************************
      subroutine snresult(i,ICASE,qsoil,wfsoil,tsoil,
     *    b1,b2,fff,delth,ictrl,pixel,y0,n_y,replicate,rank)
       include 'snow4.in'
       include 'comsib.in'
       dimension delth(nd)
       data bwe/200.0/
       integer pixel,ictrl,replicate,rank
       integer n_y
       real y0(n_y)
       
CS 
       hx=0.0 
       IF (ICASE.EQ.1)   THEN
             fi(i)=1.0
             fl(i)=0.0
             dz(i)=dzo(i)
             bw(i)=(w(i)*rhowater)/dz(i)
            if((w(i)/dz(i)).lt.0.05.or.(w(i)/dz(i))
     &       .gt.(dice/1000.0))then
              bw(i) = bwo(i)
              dz(i) = (w(i)*rhowater)/bw(i)
             endif
             bi(i)=bw(i)
             bl(i)=0.0
             bt(i)=bw(i)
             wf(i)=0.0
             if (i.eq.1) wfsoil=0.0
             dliqvol(i)=0.0
             dicevol(i)=bi(i)/dice
             ct(i)=(bw(i)/920.0)*1.9e+6
             if  (i.eq.n)  then 
               h(i)=ct(i)*dz(i)*(tssn(i)-273.16)-rhowater*dlm*w(n)*fi(n)
             else 
               tssn(i) = ( ho(i) + ct(i)*dz(i)*273.16 + b1*dtt
     &            + rhowater*dlm*w(i) )
     &            / ( ct(i)*dz(i) - b2*dtt )
               h(i) = ho(i) + (b1+b2*tssn(i))*dtt
             end if
             if(tssn(i).gt.273.16) then
c               print *, 'tssn(i)>273.16, delta=',tssn(i)-273.16
               if(tssn(i).gt.273.16+5E-5)then
c               print *, 'y0=',y0,'tssn=',tssn(1:3),'b1=',b1,'b2=',b2,
c     &           'i=',i,'ct(i)=',ct(i),'ho(i)=',ho(i),'dz(i)=',dz(i),
c     &           'w(i)=',w(i),'fi(n)=',fi(n),'itctrl=',ictrl,
c     &           'replicate=',replicate,'pixel=',pixel,'rank=',rank,
c     &           'obsnow=',obsnow
               Stop' Snow Temp. Wrong in thermal.f'
               end if
             end if
c*****************************************************************         
        ELSE IF   (ICASE.EQ.2)  THEN
c when snow temperature equals 273.16
             fl(i)=1.0-fi(i)
             tssn(i)=273.16
             wf(i)=0.0
c
          If(bwo(i).ge.bwe) Then
             if(fl(i).gt.flmin)then
                wf(i) = w(i)-(fi(i)/(1.0-flmin))*w(i)
                w(i)  = (fi(i)/(1.0-flmin))*w(i)
                dum   = wf(i)
                fl(i)=flmin
                fi(i)=1.0-fl(i)
             endif
          Else
c...................................
             flm = flmin+(flmax-flmin)*((bwe-bwo(i))/bwe)
             if(fl(i).gt.flm)then
               wf(i) = w(i)-(fi(i)/(1.0-flm))*w(i)
               w(i)  = (fi(i)/(1.0-flm))*w(i)
               dum   = wf(i)
               fl(i)=flm
               fi(i)=1.0-fl(i)
             endif
          Endif
c.................................................
c
          If( wf(i).gt.0.0) Then
             if(i.ne.1)then
               wf(i)=amin1(dum, dksatsnow*dtt)
               snroff = snroff + (dum - wf(i))
c              write (100,*)'B dum wf(i) snroff',dum,wf(i),snroff
               hroff=hroff+(dum-wf(i))*cl*rhowater*(tssn(i)-273.16)
             else
ctest2
                if(www(1).ge.1.0) then
                    snroff = snroff + wf(i)
                    wfsoil=0.0
                else
                   slwet=www(1)*poros*zdepth(1)
                   www(1)=(slwet+wf(i))/(poros*zdepth(1))
                   if(www(1).gt.1.0) then
                      snroff = snroff + (www(1)-1.0)*poros*zdepth(1)
                      www(1)=1.0
                   endif
                wfsoil=0.0
                endif
ctest2
c               snroff = snroff + wf(i)
c               wfsoil=0.0
c               write (100,*)'C wf(i) snroff',wf(i),snroff
                hroff=hroff + wf(i)*cl*rhowater*(tssn(i)-273.16)
             endif
          Endif
cccccc next concerning compaction occurring during melt
             xnodalmelt=bio(i)*dzo(i)-w(i)*rhowater*fi(i)
          If(xnodalmelt.gt.0.0.and.bio(i)*dzo(i).gt.0.0
     &       .and.(bio(i).lt.250.0.or.(i.eq.n.and.
     &        bio(i).lt.400.0))) Then
              ddz3=-xnodalmelt/(bio(i)*dzo(i))
              dz(i)=dzo(i)*(1.0+ddz3)
          Else
              dz(i)=dzo(i)
          Endif
              bw(i)=(w(i)*rhowater)/dz(i)
c.............................................
          If((w(i)/dz(i)).lt.0.05.or.(w(i)/dz(i))
     &        .gt.(dice/1000.0)) Then
            bw(i) = bwo(i)
            dz(i) = (w(i)*rhowater)/bw(i)
          Endif

          bi(i)=bw(i)*fi(i)
          bl(i)=bw(i)-bi(i)
          bt(i)=bw(i)
          ct(i)=(bw(i)/920.0)*1.9e+6
          dliqvol(i)=bl(i)/rhowater
          dicevol(i)=bi(i)/dice
          h(i)=(-1.0)*w(i)*fi(i)*dlm*rhowater
cc***********************************************
        ELSE IF (ICASE.EQ.3) THEN
c            i=n
c           else if(fff.le.0.0) then
cccccc next calculate ponding condition.
            fl(i) = 1.0
            fi(i) = 0.0
c           dz(i) = w(i)
            wf(i) = w(i)
            dum=  wf(i)
            dz(i) = 10e-15
            w(i) = 10e-15
            bw(i) =rhowater
            bl(i)=bw(i)
            bi(i)=0.0
            dliqvol(i) = 1.0
            dicevol(i) = 0.0
            ct(i)=(bw(i)/920.0)*1.9e+6
            tssn(i) = 273.16
            h(i) = 0.0
c 
            If (i.eq.n) Then
               if (i.eq.1) then
                  hx=(-1.0)*w(i)*fff*dlm*rhowater/dtt 
                  snroff=wf(1)+snroff 
                  wfsoil=0.0
               else 
                  wf(i)=amin1(dum, dksatsnow*dtt)
                  snroff = snroff + (dum - wf(i))
                  hroff=hroff+(dum-wf(i))*cl*rhowater*(tssn(i)-273.16)
                  delth(n-1) = (-1.0)*w(i)*fff*dlm*rhowater/dtt 
               end if 
            Else 
              if(i.eq.1)then
                 hx         = ho(i)/dtt + b1+b2*tssn(i)
ctest2
                if(www(1).ge.1.0) then
                    snroff = snroff + wf(i)
                    wfsoil=0.0
                else
                   slwet=www(1)*poros*zdepth(1)
                   www(1)=(slwet+wf(i))/(poros*zdepth(1))
                   if(www(1).gt.1.0) then
                      snroff = snroff + (www(1)-1.0)*poros*zdepth(1)
                      www(1)=1.0
                   endif
                wfsoil=0.0
                endif
ctest2
c                 snroff = snroff + wf(i)
c                 wfsoil= 0.0
              else
                 wf(i)=amin1(dum, dksatsnow*dtt)
                 snroff = snroff + (dum - wf(i))
                 hroff=hroff+(dum-wf(i))*cl*rhowater*(tssn(i)-273.16)
                delth(i-1) = ho(i)/dtt + b1+b2*tssn(i)
              endif
            End if
c
        END IF
cS  Calculate the heat flux into  the soil: qsoil     on 10/13/98.
cS  qsoil : downward is positive [ W/m**2]
        if (i.eq.1) qsoil =  qk(1)*(tssn(1) - tsoil) + hx 
cs                                                       10/13/98 
        return
        end 
c------------------------------------------------------------------------
c
      subroutine output1(cosz,siblh,isnow,iday,n_steps,yout,n_y,
     &  xout,n_x,zout,tend,pixel)
c
c------------------------------------------------------------------------
c  this subroutine stores states to output array
c    by mike
c------------------------------------------------------------------------

      include 'comsib.in'
      include 'snow4.in'

      integer n_y,n_x,n_steps,pixel
      dimension yout(n_y,n_steps), xout(n_x,n_steps),
     &  zout(n_steps),tend(5)
      
c     formerly configured to output daily, now outputs hourly (but 
c       still uses iday index... ought to be redone with ihour 
c       or something

        iday=iday+1
        if(isnow.eq.0)then
          yout(1,iday)=sum(dzo)
          yout(2,iday)=bwo(1)
          yout(3,iday)=bwo(2)
          yout(4,iday)=bwo(3)
          yout(5,iday)=tssno(1)
          yout(6,iday)=tssno(2)
          yout(7,iday)=tssno(3)
          yout(8,iday)= flo(1)*bwo(1)/1000
          yout(9,iday)= flo(2)*bwo(2)/1000
          yout(10,iday)=flo(3)*bwo(3)/1000
          yout(11,iday)=gdia(1)
          yout(12,iday)=gdia(2)
          yout(13,iday)=gdia(3)
          yout(14,iday)=TGS
        elseif(snowdepth.gt.0)then
          yout(1,iday)=snowdepth
          yout(2,iday)=capac(2)*1000/snowdepth
          yout(3,iday)=capac(2)*1000/snowdepth
          yout(4,iday)=capac(2)*1000/snowdepth
          yout(5,iday)=TGS
          yout(6,iday)=TGS
          yout(7,iday)=TGS
          yout(8,iday)=0.
          yout(9,iday)=0.
          yout(10,iday)=0.
          yout(11,iday)=gsize
          yout(12,iday)=gsize
          yout(13,iday)=gsize
          yout(14,iday)=TGS
        else
          yout(1:14,iday)=0.
          yout(14,iday)=TGS
        end if

cm        xout(1 ,iday)=TD !deep soil temp
cm        xout(2 ,iday)=TC !canopy temp
cm        xout(3 ,iday)=TA !air temp
cm        xout(4 ,iday)=TM !temperature at reference height
cm        xout(5 ,iday)=www(1) !soil saturation in layer 1
cm        xout(6 ,iday)=www(2) !soil saturation in layer 2
cm        xout(7 ,iday)=www(3) !soil saturation in layer 3
cm        xout(8 ,iday)=capac(1) !canopy storage
cm        xout(9 ,iday)=egi !snow evaporation    
cm        xout(10,iday)=gsize !one layer grain size
cm        xout(11,iday)=snden !one layer snow density
cm        xout(12,iday)=capac(2) !snow water equivalent (used for one-layer) 

        xout(1 ,iday)=TD 
        xout(2 ,iday)=TC 
        xout(3 ,iday)=TA 
        xout(4 ,iday)=TM 
        xout(5 ,iday)=www(1) 
        xout(6 ,iday)=www(2) 
        xout(7 ,iday)=www(3) 
        xout(8 ,iday)=capac(1) 
        xout(9 ,iday)=egi 
        xout(10,iday)=gsize 
        xout(11,iday)=snden 
        xout(12,iday)=capac(2)
cm      the output runoff is a new addition...
        xout(13,iday)=roff

c record albedo

        if(sibswup.eq.-9999.)then
          zout(iday)=-9999.
        else
          zout(iday)=sibswup/swdown
        end if

c if this is the final simulation time, record the tf vector
        if(iday.eq.n_steps)then
cm          tend(1)=real(nhh) !Hour
cm          tend(2)=real(nmm) !Month
cm          tend(3)=real(ndd) !Day of the Month
cm          tend(4)=real(day) !Julian Day
cm          tend(5)=real(nyy) !Year
          tend(1)=real(nhh) 
          tend(2)=real(nmm) 
          tend(3)=real(ndd) 
          tend(4)=real(day) 
          tend(5)=real(nyy) 
        end if
      return
      end
c
c============================================================================
c
       subroutine errors
c
c----------------------------------------------------------------------------
c ** this subroutine is to calculate the RMS error
c ** between SSiB computed and observed values
c----------------------------------------------------------------------------
c
       include 'comsib.in'
c
c ** initialise counts and set up loop controls (month and day)
c
       iswup = 0
       inet = 0
       iustar = 0
       isht = 0
       ilh = 0
       ievap = 0
       itg = 0
       ishf = 0
       ndays = nobs/24
       km = mthst
       kd = ndyst
c
       do 5000, i=1,ndays
cs sun add following part on 02/24/99 start
      kdmax = 31
      do 2416, ikahan=0,19
      if ((km.eq.(4+12*ikahan)).or.
     & (km.eq.(9+12*ikahan)).or.
     & (km.eq.(11+12*ikahan)).or.
     & (km.eq.(6+12*ikahan))) then
      kdmax = 30
      endif
 2416 continue
      do 2417, ikahan=0,19
      if (km.eq.(2+12*ikahan)) then
      kdmax=28
      if ((km.eq.14).or.(km.eq.62).or.(km.eq.110).
     & or.(km.eq.158).or.(km.eq.206)) then
      kdmax = 29
      endif
      endif
 2417 continue
c ** increment loop controls (month and day)
c
         kd = kd+1
         if (kd.gt.kdmax) then
           km = km + 1
c          if(km.gt.12) km=km-12
           kd = 1
         endif
c
 5000  continue
      return
      end
c----------------------------------------
cjyj   this is latest version modified by jyj
         subroutine modnode
c*****
      include 'snow4.in'

      if(dzo(n) .lt. dzmin .or.wo(n).lt.womin
     &   .or.bwo(n).eq.rhowater.or.flo(n).eq.1.0) then
         call combojyj(3,2,0.0)
         ppp=amax1(2./3. , 1.-(0.02-0.001)/dzo(2))
         call combojyj(2,3,ppp)
      endif
c*****
      if(dzo(3).gt.0.02) then
         ppp=1./3.
         ppp=0.02/dzo(3)
         call combojyj(3,2,ppp)
      endif        
c*****
      if(snowdepth.lt.0.62) then
             if((dzo(2)/dzo(1)-0.5).gt.1d-1)then
             ppp=(1+dzo(1)/dzo(2))*(1./3.)
             call combojyj(2,1,ppp)
         else 
             ppp=(1+dzo(2)/dzo(1))*(2./3.)
             call combojyj(1,2,ppp)
         endif
      else
         if(dzo(2).gt.0.20)then
             ppp=2./3.
             call combojyj(2,1,ppp)
         else 
             ppp=1.0-(0.38*0.5-dzo(2))/dzo(1)
             call combojyj(1,2,ppp)
         endif
      endif
c*****
      return
      end
c======================================================
      subroutine combojyj(ndiv,nnew,ratio)
      include 'snow4.in'
c......................
      dzo(nnew)=dzo(nnew)+(1-ratio)*dzo(ndiv)
       ho(nnew)= ho(nnew)+(1-ratio)* ho(ndiv)
       wo(nnew)= wo(nnew)+(1-ratio)* wo(ndiv)
      bwo(nnew)= wo(nnew)*rhowater/dzo(nnew)
      ctit=1.9e6*bwo(nnew)/920.0
      hmlt=wo(nnew)*rhowater*dlm
      if(ho(nnew).ge.(-1.0*hmlt)) then
         tssno  (nnew) = 273.16
         fio    (nnew) = (-1.0)*ho(nnew)/hmlt
         flo    (nnew) = 1.0-fio(nnew)
         blo    (nnew) = bwo(nnew)*flo(nnew)
         bio    (nnew) = bwo(nnew)*fio(nnew)
         dliqvol(nnew) = blo(nnew)/rhowater
      else
c
         flo  (nnew)=0.0
         fio  (nnew)=1.0
         tssno(nnew)=(ho(nnew)+hmlt)/(ctit*dzo(nnew))+273.16
         bio  (nnew)=bwo(nnew)
         blo  (nnew)=0.0
         bto  (nnew)=bwo(nnew)
         dliqvol(nnew) = 0.0
      endif

      dzo(ndiv)=ratio*dzo(ndiv)
       ho(ndiv)=ratio* ho(ndiv)
       wo(ndiv)=ratio* wo(ndiv)

       return
       end
c **********************************************************:
      subroutine modnodenew
      include 'snow4.in'
clwp  10/30/2000, for the adjustment of layers 2,3
       IF (snowdepth.le.0.06) then
         DZ1=0.02
         DZ2=DZ1
         DZ3=snowdepth-( DZ2+DZ1)
       ELSE
         DZ3=0.02
       ENDIF
c     to get the expected change of top layer of snow
         DDZ3=DZ3-dzo(3)
c     to get the expected change of top layer of snow
      IF (DDZ3.GT.0.0) THEN
          DDZ3=MIN(DDZ3,dzo(2))
          CALL COMBO (DDZ3,dzo(3),dzo(2),wo(3),wo(2),ho(3),ho(2),
     &    tssno(3),tssno(2),bwo(3),bwo(2),bio(3),bio(2),blo(3),blo(2),
     &    bto(3),bto(2),fio(3),fio(2),flo(3),flo(2),cto(3),cto(2),
     &    dliqvol(3),dliqvol(2),dicevol(3),dicevol(2))
      ELSE
          DDZ3=-DDZ3
           CALL COMBO (DDZ3,dzo(2),dzo(3),wo(2),wo(3),ho(2),ho(3),
     &    tssno(2),tssno(3),bwo(2),bwo(3),bio(2),bio(3),blo(2),blo(3),
     &    bto(2),bto(3),fio(2),fio(3),flo(2),flo(3),cto(2),cto(3),
     &    dliqvol(2),dliqvol(3),dicevol(2),dicevol(3))
      END IF
clwp  10/30/2000, for the adjustment of layers 1,2
      SUM12=dzo(1)+dzo(2)
      IF (snowdepth.le.0.06) THEN
      DZ2=0.5*SUM12
        ELSE IF (snowdepth.gt.0.06.and.snowdepth.le.0.08) THEN
        DZ2=0.02
          ELSE IF (snowdepth.gt.0.08.and.snowdepth.le.0.62) THEN
          DZ2=0.33333333*SUM12
            ELSE IF (snowdepth.gt.0.62) THEN
            DZ2=0.20
      ENDIF
c     to get the expected change of middle layer of snow
         DDZ2=DZ2-dzo(2)
c     to get the expected change of middle layer of snow
      IF (DDZ2.GT.0.0) THEN
          CALL COMBO (DDZ2,dzo(2),dzo(1),wo(2),wo(1),ho(2),ho(1),
     &    tssno(2),tssno(1),bwo(2),bwo(1),bio(2),bio(1),blo(2),blo(1),
     &    bto(2),bto(1),fio(2),fio(1),flo(2),flo(1),cto(2),cto(1),
     &    dliqvol(2),dliqvol(1),dicevol(2),dicevol(1))
      ELSE
          DDZ2=-DDZ2
          CALL COMBO (DDZ2,dzo(1),dzo(2),wo(1),wo(2),ho(1),ho(2),
     &    tssno(1),tssno(2),bwo(1),bwo(2),bio(1),bio(2),blo(1),blo(2),
     &    bto(1),bto(2),fio(1),fio(2),flo(1),flo(2),cto(1),cto(2),
     &    dliqvol(1),dliqvol(2),dicevol(1),dicevol(2))
      END IF
      Return
      End

c*******************************************************************
c     COMBD computes new grain diameters for each layer after 
c        depth update
c        by Mike, 3/16/04
c*******************************************************************
      subroutine COMBD
      include 'snow4.in'
      include 'comsib.in'
      real gdiao(3),gdia_uc(1),gdia_c(1),gdia_uco(1),gdia_co(1)
      real ext_dep_zn1,ext_dep_zn2,ext_dep_zn3,ext_dep_z11,
     &      ext_dep_z21,ext_dep_z31,ext_dep_z12,ext_dep_z22,ext_dep_z32,
     &      ext_dep_z13,ext_dep_z23,ext_dep_z33
      znew=dzfall+dzo(3)+dzo(2)+dzo(1)
c     determine snow depths after udpate     
      IF (znew.le.0.06) THEN
          DZ1=0.02
          DZ2=DZ1
          DZ3=znew-(DZ2+DZ1)
      ELSE IF (znew.gt.0.06.and.znew.le.0.08) THEN
          DZ3=0.02
          DZ2=0.02
          DZ1=znew-(DZ3+DZ2)
      ELSE IF (znew.gt.0.08.and.znew.le.0.62) THEN
          DZ3=0.02
          DZ2=0.33333333*(znew-DZ3)
          DZ1=znew-(DZ3+DZ2)
      ELSE IF (znew.gt.0.62) THEN
          DZ3=0.02
          DZ2=0.20
          DZ1=znew-(DZ3+DZ2)
      ENDIF
c     compute grain diameters after update
c
c     1. Assignment
c     1.1 Assign z3+
c     case 3.1
      if(dzfall>=DZ3) then
        zn3=DZ3
        z33=0.
        z23=0.
        z13=0.
c     case 3.2
      elseif (dzfall+dzo(3)>=DZ3) then
        zn3=dzfall
        z33=DZ3-zn3
        z23=0.
        z13=0.
c     case 3.3
      else
        zn3=dzfall
        z33=dzo(3)
        z23=DZ3-zn3-z33
        z13=0.
      endif
c     1.2 assign z2+
      if(dzfall>=DZ3+DZ2) then
        zn2=DZ2
        z32=0.
        z22=0.
        z12=0.
      elseif(dzfall+dzo(3)>=DZ3+DZ2) then
        zn2=dzfall-zn3
        z32=DZ2-zn2
        z22=0.
        z12=0.
      elseif(dzfall+dzo(3)+dzo(2).ge.DZ3+DZ2.and.z33.lt.dzo(3)) then
        zn2=0.
        z32=dzo(3)-z33
        z22=DZ2-z32
        z12=0.
      elseif(dzfall+dzo(3)+dzo(2).ge.DZ3+DZ2.and.z33.eq.dzo(3)) then
        zn2=0.
        z32=0.
        z22=DZ2
        z12=0.
      else
        zn2=0.
        z32=0.
        z22=DZ3-dzfall-dzo(3)
        z12=DZ2-z22
      endif
c     1.3 assign z1+
      if (DZ1<=dzo(1)) then
        z11=DZ1
        z21=0.
        z31=0.
        zn1=0.
      elseif (DZ1<=dzo(1)+dzo(2)) then
        z11=dzo(1)
        z21=DZ1-z11
        z31=0.
        zn1=0.
      elseif (DZ1<=dzo(1)+dzo(2)+dzo(3)) then
        z11=dzo(1)
        z21=dzo(2)
        z31=DZ1-z11-z21
        zn1=0.
      else
        z11=dzo(1)
        z21=dzo(2)
        z31=dzo(3)
        zn1=DZ1-dzo(1)-dzo(2)-dzo(3)
      endif
c     2. compute new grain diameters
      gdiao(3)=gdia(3)
      gdiao(2)=gdia(2)
      gdiao(1)=gdia(1)
      gdia_uco(1)=gdia_uc(1)
      gdia_co(1)=gdia_c(1)     


c     compute extiction depths of each contributed layears,Dongyue, 11 apr,2012 
      ext_dep_zn1=exp(-5*zn1/cos(55/180*3.1415926)) 
      ext_dep_z11=exp(-5*z11/cos(55/180*3.1415926))
      ext_dep_z21=exp(-5*z21/cos(55/180*3.1415926))
      ext_dep_z31=exp(-5*z31/cos(55/180*3.1415926))
      ext_dep_zn2=exp(-5*zn2/cos(55/180*3.1415926))
      ext_dep_z12=exp(-5*z12/cos(55/180*3.1415926))
      ext_dep_z22=exp(-5*z22/cos(55/180*3.1415926))
      ext_dep_z32=exp(-5*z32/cos(55/180*3.1415926))
      ext_dep_zn3=exp(-5*zn3/cos(55/180*3.1415926))
      ext_dep_z23=exp(-5*z23/cos(55/180*3.1415926))
      ext_dep_z33=exp(-5*z33/cos(55/180*3.1415926))

      gdia(3)=((bifall*zn3*ddfall+bio(3)*z33*gdiao(3)+
     &         bio(2)*z23*gdiao(2))/(bifall*zn3+bio(3)*z33
     &          +bio(2)*z23))
      gdia(2)=((bifall*zn2*ddfall+bio(3)*z32*gdiao(3)+bio(2)*z22
     &           *gdiao(2)+bio(1)*z12*gdiao(1))/(bifall*zn2+bio(3)
     &           *z32+bio(2)*z22+bio(1)*z12))    
      gdia(1)=((bifall*zn1*ddfall+bio(3)*z31*gdiao(3)+bio(2)
     &         *z21*gdiao(2)+bio(1)*z11*gdiao(1))/(bifall*zn1
     &         +bio(3)*z31+bio(2)*z21+bio(1)*z11))

        

!cl    dongyue track bottom layer grain size change 16 June,12 
!       IF (gdia(1).LT.gdiao(1)) THEN
!          cpex=0.6123*(log(1+(rainf+obsnow)*8640)/log(3.6))
!          if(cpex.gt.1.0)then
!           cpex=1.
!          endif
!       ENDIF

      return
      end

c*******************************************************************
c     graingrowth computes growth of grain diameter
c        by Mike, 3/22/04
c*******************************************************************
      subroutine graingrowth(ISNOW)

      include 'snow4.in'
      include 'comsib.in'
      real stheta(3)

c     define a few physical and theoretical constants, parameters
      uvlim=1E-6
      de0=0.9E-4
      c1l=5.276E8              
      c1i=8.047E9
      Lvl=2.505E6
      Lvi=2.838E6
      Rw=461.296

      IF(ISNOW>0)THEN
c       evolve grain diameter for one layer
c       compute evaporation in units of kgm^-2s^-1
        Uvai=egi/Lvi/dtt
c       compute rate of grain growth
        if(abs(Uvai)<1E-6)then
          gsizdt=sg1(3)/gsize*abs(Uvai)
        else
          gsizdt=sg1(3)/gsize*1E-6
        endif
c       compute new grain diamter and return to MAIN
        gsize=gsize+gsizdt*dtt
      return
      ELSE
c       compute stheta in each layer
        DO i=1,3
          stheta(i)=flo(i)*bio(i)/1000
        END DO

c       evolve grain diameter for three layers      
        if(stheta(3)<=0.02)then
          Uvai=egi/Lvi/dtt
        else
          Uvai=egi/Lvl/dtt
        endif

        DO i=1,3
c         calculate df for each layer, as product of de and ckt
          dee(i)=de0*(1000/psur)*(tssno(i)/273.16)**6
          if(stheta(i)<=0.02)then
            ckt(i)=c1i/(tssno(i)**2)*(Lvi/Rw/tssno(i)-1)*
     &               exp(-Lvi/Rw/tssno(i))
          else
            ckt(i)=c1l/(tssno(i)**2)*(Lvl/Rw/tssno(i)-1)*
     &               exp(-Lvl/Rw/tssno(i))
          endif
          df(i)=dee(i)*ckt(i)       
        END DO
c       calculate rate of change of grain diameter in each layer

        DO i=1,3
          if(stheta(i)<=1E-4)then
            call vapflux(i)
            if(abs(uvbar)<uvlim)then
              dddt(i)=sg1(i)/gdia(i)*abs(uvbar)
            else
              dddt(i)=sg1(i)/gdia(i)*uvlim
            endif
          elseif(stheta(i)<0.09)then
            dddt(i)=sg2/gdia(i)*(stheta(i)+0.05)
          else
            dddt(i)=sg2/gdia(i)*0.14
          endif
          gdia(i)=gdia(i)+dddt(i)*dtt
        END DO
      ENDIF      

      return
      end

c*******************************************************************
c     VAPFLUX computes vapor flux through the snowpack
c        by Mike, 3/22/04
c*******************************************************************

      subroutine VAPFLUX(i)
      include 'snow4.in'
      include 'comsib.in'
      
c     compute flux through the upper node boundary      
      if(i==3)then
        Uvu=Uvai
      else
        Uvu=2*df(i+1)*df(i)*(tssno(i+1)-tssno(i))/
     &         (dzo(i+1)*df(i)+dzo(i)*df(i+1))
      endif

c     compute flux through the lower boundary      
      if(i==1)then
        Uvl=0
      else
        Uvl=2*df(i)*df(i-1)*(tssno(i)-tssno(i-1))/
     &         (dzo(i)*df(i-1)+dzo(i-1)*df(i))
      endif
c     compute average flux 

      if(i==1)then
        uvbar=abs(Uvu)
      else
        uvbar=0.5*(abs(Uvu)+abs(Uvl))
      endif

      return
      end

c************************************************************************
c     SNTALB utilizes the albedo routine used in SNTHERM, getmet routine
c        by Mike, 4/8/04
c************************************************************************

      subroutine SNTALB(Dn,cosz,ALBEDO,salbo,IWAVE,albfac)

c       Dn is the grain diameter in the surface layer - gdia(3) elsewhere
c        cosz is the cosine of the zenith angle - SUNANG elsewhere in SAST
c        These formulae are described in Marks and Dozier, 1992:
c        Climate and Energy Exchange at the Snow Surface in the Alpine Region of 
c        the Sierra Nevada 2. Snow Cover Energy Balance.  WRR 28(11): 3043-3054

      dimension ALBEDO(2,3,2)
      dimension salbo(2,2)      

c     compute three values used in algorithm
      sqrtr=sqrt(0.50*Dn)
      dum=0.57735026918963
    !  tcosz=max(cosz,0d0)

c     save copies of original albedos
      salbo(IWAVE,1)=ALBEDO(2,IWAVE,1)
      salbo(IWAVE,2)=ALBEDO(2,IWAVE,2)     

c     compute visible or nir albedo, based on IWAVE
      select case (IWAVE)
        case(1)
c         Estimate visible albedos for diffuse radiation incident on 
c           horizontal surface. getmet, line 156
          ALBEDO(2,IWAVE,2)=albfac*(1d0-sqrtr*2d0)
c         Estimate visible albedos for direct radiation at given zenith angle.
c           getmet, line 168
          ALBEDO(2,IWAVE,1)=albfac*(1d0-sqrtr*2d0+
     &                          sqrtr*1.575d0*(dum-tcosz))     
        case(2)
c         Estimate NIR albedos for diffuse radiation incident on 
c           horizontal surface. getmet, line 156
          ALBEDO(2,IWAVE,2)=albfac*(0.85447*exp(-21.23*sqrtr))
c         Estimate visible albedos for direct radiation at given zenith angle.
c           getmet, line 168
          ALBEDO(2,IWAVE,1)=albfac*(0.85447*exp(-21.23*sqrtr)+
     &                       (2.4d0*sqrtr+0.12d0)*(dum-tcosz))

        end select

      return
      end

C **************************************************
