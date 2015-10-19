! This code has a weak probe beam. Beam 2
!*************************************************************************
! Code for solving 3D bead dynamics in the presence of an ion trap
! UNSCALED variables throughout.
!  KAPP2= amplitude decay= 10^5-10^6 Hz 
! WK= wavenumber sin kx= sin WK*X k= 5.9*10^6 
! XM= mass of the bead= 9.2*10^-18Kg for r=100nm bead
!                     = 9.2*10^-21Kg for r=10nm bead
! E1= Trap Input pumping amplitude 
! A= field coupling. U(x)=hbar A |a|^2 cos ^2 (kx-phas1)
! DETUN1= trap detuning= about zero MHz
! use only real variables for speed
! STATE(10)= array with phase space variables
! STATE(1)=re(alpha1); STATE(2)=Imag (alpha1) Trap beam
! STATE(3)=re(alpha2); STATE(4)=Imag (alpha2) Weak probe beam
! STATE(5)=dX/dt=P_x/m; STATE(6)=X centre of mass variables
! STATE(7)=dY/dt=P_y/m; STATE(8)=y centre of mass variables
! STATE(9)=dZ/dt=P_z/m; STATE(10)=Z centre of mass variables
!  
! PHAS1= phase shift of trap field, between trap and probe beams 
! GammaM=mechanical damping due to gas
! GAMMN= optical cool rate , numerical
! nav= used to average over oscillations to get gamma
! nbreak= break time period into nbreak chunks in RK
! OFFSET=DC offset field March 19 2014
!*********************************************
! NB no intent(in) or intent(out) used- add these later!!
!************************************************************************()

IMPLICIT NONE      
integer  ::ll,kk,it,in,ii,jj,NWELLS,NPERIOD,NTOT
integer:: mid,nstoch,NPM,ierr
integer  ::nav,nbreak,nrkutta

! time propagation
double precision::pi,pi2,DT,TDEL,TIN,DET,TEMpeq,F0,freq,Total,gammbig,test
double precision:: TPERIOD,ASQ,Atrap,TT,DETUN1,DETUN2,OFF
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,Press,TEMP
double precision::RHO,WK,PIN,PIN2,QX,QY,QZ,Q
double precision::RTRAP,V0,trapfreq,omega,Om,omsec,energy

! omega=trap drive; omegam=mech freq.
! calculated input parameters
double precision:: E1,E2,PHAS1,tem1,AOUTr,AOUTi
double precision:: KAPP2,XM,A,OMEGAM,AMP,Vamp,XV
double precision:: x0,Wk0,check,GAMMAM,OMTRAPsq,VOPT,W2
double precision:: XT0,YT0,ZT0,XX,YY,ZZ,XP0,YP0,ZP0,Avcos,Scatter,RR,pointpcle

! read parameter file with input values
include 'PROBE2_OMIT.h'

! arrays for the FTs
DOUBLE PRECISION, DIMENSION(NPERIOD):: OMSTOR
DOUBLE PRECISION, DIMENSION(NPERIOD):: SR1,SI1,SR2,SI2,SR3,SI3
DOUBLE PRECISION, DIMENSION(NPERIOD)::GR,GI,GR3,GI3
DOUBLE PRECISION, DIMENSION(NPERIOD)::SPEC1,SPEC2,SPEC3
DOUBLE PRECISION, DIMENSION(NPERIOD)::SP1,SP2,SP3
DOUBLE PRECISION:: SXXR(NPERIOD),SXXI(NPERIOD)

! arrays 
DOUBLE PRECISION, DIMENSION(NTOT):: STATEeq,STATE,STATEout
DOUBLE PRECISION, DIMENSION(NPERIOD):: XT,YT,ZT,PX,GAMMAT,XMEAN,Signal,Time

! initialise random number generator
call random_seed()
pi=dacos(-1.d0)
pi2=2.d0*pi
PHAS1=0.d0
DETUN1=DET*pi2

! set the offset between traps 
! this setting looks like the cooling curves!
OFF=0.d-5

! stochastic realisations
nstoch=1
NPM=NPERIOD

! averaging for cooling calculations
nav=100

! break each mechanical freq oscillation into nbreak
! chunks for plotting
nbreak=10

! subdivide further into nrkutta steps for propagation
nrkutta=4000
 
OPEN(8,file="NUMGAMMA_OMIT.dat",status="unknown")
OPEN(12,file="TRAJECTORY_OMIT.dat",status="unknown")
OPEN(14,file="TRAJECTORXYZ_OMIT.dat",status="unknown")
OPEN(16,file="TRANSMISSION_OMIT.dat",status="unknown")
OPEN(18,file="SPECTRA_OMIT.dat",status="unknown")

write(6,*)'Charge number, Voltage '
write(6,100) Q/1.6*1d19,V0

! open loop over detuning
do 10 ii=-7,9

    DETUN2=50.d3+(ii-1)*5d3
    ! DETUN2=abs(DET)
    DETUN2=DETUN2*pi2

    ! Zero FT functions
    do jj=1,NPERIOD
        OMSTOR(jj)=0.d0
        SR1(jj)=0.d0
        SI1(jj)=0.d0
        SR2(jj)=0.d0
        SI2(jj)=0.d0
        SR3(jj)=0.d0
        SI3(jj)=0.d0
        GR(jj)=0.d0
        GI(jj)=0.d0
        SPEC1(jj)=0.d0
        SPEC2(jj)=0.d0
        SPEC3(jj)=0.d0
        SP1(jj)=0.d0
        SP2(jj)=0.d0
        SP3(jj)=0.d0
    enddo
       
    ! initialise by finding an equilibrium state for photon field and particle
    ! for that input power; calculate relevant parameters
    CALL EQUIL(E1,E2,A,DETUN1,DETUN2,GammaM,W2,XM,Kapp2,OMtrapsq,STATEeq,WK0,OMEGAM,omega,AMP,Vamp)
        
    !INITIALISE POSITIONS AND MOMENTA
    write(6,*)' Energy in secular motion'

    ! try to guess initial energy of particle
    ! assume we start NWELLS  out from ion trap centre and accelerate back
    NWELLS=100
    XT0=NWELLS/WK*pi
    Qx=2*OMTRAPsq/omega/omega
    omsec=OMTRAPsq/sqrt(2.d0)/omega
    Energy=XM*omsec**2*XT0**2
    write(6,100) energy
    write(6,*) 'secular frequency omegas and f'
    write(6,100) omsec,omsec/2/pi
    write(6,*)'Max speed from secular motion'
    XV=sqrt(2*energy/XM)
    write(6,100) XV
  
    ! set all x,y,z intial values to zero
    do kk=5,NTOT
        STATEeq(kk)=0.d0
    enddo

    !  reset position in x:
    STATEeq(4)=xt0
    write(6,*) 'x,vx,vy,vz'
    write(6,200)STATEeq(4),STATEeq(3),STATEeq(5),STATEeq(7)
    write(6,*)'trap beam detuning, kappa/2 and omegam='
    write(6,200) detun1,kapp2,omegam

    ! period of each oscillation at equilibrium
    OM=max(omegam,abs(detun1))
    write(6,*)'maximum frequency scale max(omM,omega),period'
    TPERIOD=2.d0*pi/OM
    write(6,200)OM,TPERIOD

    ! plot/analyse data every time interval TDEL
    TDEL=TPERIOD/nbreak
    ! TDEL is *further* subdivided into nrkutta intervals later

    ! now time evolve these initial conditions with Runge Kutta propagator
    ! call time propagator to evolve initial state for a time interval
    ! =TDEL in each pass. So true time =TIN +(it-1)*TDEL
    ! first evolve the Runge Kutta with high damping to settle the particle
    ! in the chosen well
    ! assume pressure= 10 mbar
    Tin=0.d0
    gammbig=1.d4
    ! fill state arrays with estimates
    do 20 ll=1,NTOT
        STATE(ll)=STATEeq(ll)
    20 enddo

    do it=1,20
        ! evolve for a time TDEL*20 so offset and gravity equilibrate
        CALL RKINI(nrkutta,STATE,STATEout,TDEL,Tin,DETUN1,DETUN2,XM,OMTRAPsq,W2,GAMMBIG,A,E1,E2,OMEGA,kapp2,OFF)
        Tin=Tin+Tdel
    enddo   

    ! now you have a cold particle- now run with appropriate noise
    ! reset the position so as to add one quarter of the well half width
    ! State(4)=State(4)-5.32d-7/6

    Tin=0.d0
    Do 30 it=1,Nperiod
        ! store bead coordinates as a function of time
        PX(it)=STATE(3)
        XT(it)=STATE(4)
        YT(it)=STATE(6)
        ZT(it)=STATE(8)

        ! if bead jumps out, put it back at well=NWELLS!!
        test=STATE(4)*WK/pi-NWELLS
        if(abs(test).gt.0.4d0)then 
            XT0=NWELLS/WK*pi
            write(6,*)'jumped out, well=',STATE(4)*WK/pi,XT0
            STATE(4)=XT0
            STATE(3)=0.d0
        endif

        ! OUTPUT   
        AOUTr=STATE(1)
        AOUTi=STATE(2)
        ASQ=sqrt(AOUTr**2+ AOUTi**2)

        ! Fill parameters for the FT:trap
        SR3(it)=ASQ
        SI3(it)=0.d0
        SR1(it)=AOUTr
        SI1(it)=0.d0
        TIME(it)=TIN
        XP0=XT(it)
        YP0=YT(it)
        ZP0=ZT(it)
        SR2(it)=XP0
        SI2(it)=0.
                    
        ! work out the scattering rate of light averaged over the nanosphere
        RR=R0
        Scatter=Avcos(WK,RR,XP0)

        ! point particle
        ! Scatter=cos(WK*XP0)**2
        ! Below the envelope has no 2*y^2/w^2 because factor of two is in W2 already
        tem1=Scatter*EXP(-2*(YP0**2+ZP0**2)/W2)
        pointpcle=cos(WK*XP0)**2*EXP(-2*(YP0**2+ZP0**2)/W2)
        write(12,200)Tin,ASQ
        write(14,200)Tin,WK*XP0/pi,WK*YP0/pi,WK*ZP0/pi
     

        ! Evolve  positions in time with basic 4th order Runge Kutta.
        ! nrkutta time increments
         CALL RK(nrkutta,STATE,STATEout,TDEL,Tin,DETUN1,DETUN2,XM,OMTRAPsq,W2,GAMMAM,A,E1,E2,OMEGA,kapp2,OFF)
            
        ! update the time
        Tin=Tin+TDEL
                     
        ! overwrite initial state before going round the loop 30 again
        do 31 ll=1,NTOT
            STATE(ll)=Stateout(ll)
        31 enddo

        ! close time loop
    30 enddo

    ! Work out Fourier transforms
    ! Calculate spectrum of  output square term by a discrete Fourier transform
    call sfft(SR1,SI1,NPM,NPM,NPM,1,ierr)
    F0=2.*pi/(TDEL*(nperiod-1))
    do kk=1,NPM
        ! if need be sum over stochastic realisations  
        Spec1(kk)=sqrt(SR1(kk)**2+SI1(kk)**2)
    enddo

    ! Calculate spectrum of position
    ! by a discrete Fourier transform
    F0=2.*pi/(TDEL*(nperiod-1))
    call sfft(SR2,SI2,NPM,NPM,NPM,1,ierr)
     
    do kk=1,NPM
        ! if need be sum over stochastic realisations  
        Spec2(kk)=sqrt(SR2(kk)**2+SI2(kk)**2)
    enddo

    ! Calculate spectrum of complex FT output field
    call sfft(SR3,SI3,NPM,NPM,NPM,1,ierr)
    do kk=1,NPM
        ! if need be sum over stochastic realisations  
        Spec3(kk)=sqrt(SR3(kk)**2+SI3(kk)**2)
    enddo

    ! Reorder FTs from 1-N to -N/2, N/2
    mid=NPM/2
    do kk=1,NPM
        if(kk.le.mid)then
            in=kk-1
            else
            in=-(NPM+1-kk)
        endif 
        SR1(in+mid+1)=pi2*SPEC1(kk)
        SR2(in+mid+1)=pi2*SPEC2(kk)
        SR3(in+mid+1)=pi2*SPEC3(kk)
    enddo

    ! remove  spike at freq=0!!
    SR3(mid)=0.5*(SR3(mid-2)+SR3(mid+2))
    SR2(mid)=0.5*(SR2(mid-2)+SR2(mid+2))
    SR1(mid)=0.5*(SR1(mid-2)+SR1(mid+2))
    SR3(mid+1)=SR3(mid)
    SR2(mid+1)= SR2(mid)
    SR1(mid+1)= SR1(mid)

    do in=1,NPM
        jj=in-mid-1
        freq= jj*F0/pi2
        if((freq.gt.0d0).and.(freq.lt.2e5))then
            write(18,200) freq,(SR1(in))/E1,(SR2(in)),(SR3(in))/E1
        endif
        omstor(in)=freq
    enddo

    ! estimate the transmission based on integrated FT of the probe
    ! use the complex FT
    call NORM(NPM,TOTAL,SR3,OMSTOR)
    write(6,*)'detun2,total', Detun2/2/pi,Total
    write(16,200) Detun2/2/pi,Total

10 enddo

100   FORMAT(2E14.6,1x,4(E14.6))
200  FORMAT(5E16.8)

STOP
END

!*************************************************************************
!*************************************************************************
Function Avcos(WK,RR,XP0)
! Average the scattering rates over finite volume of nanosphere 
! Avcos=0.5-3/[16(kr)^3] cos 2kx_0 *(2kr*cos 2kr-sin 2kr)
!*************************************************************************
!*************************************************************************

IMPLICIT NONE 
integer  ::ii,jj
double precision:: pi,WK,RR,XP0
double precision:: WKR,WK0,tem1,tem2,tem3,Avcos

pi=dacos(-1.d0)
WKR=WK*RR
WK0=WK*XP0  

tem1=cos(2*WK0)
tem2=cos(2*WKR)
tem3=sin(2*WKR)
AVcos=0.5d0-3.d0/(16.*WKR**3)*tem1*(2*WKR*tem2+tem3)

! point particle limit
! Avcos=0.5+0.5*tem1

return
end

!*************************************************************************
!*************************************************************************
subroutine iontrap(OMTRAPsq,omega,Qx,Qy,Qz,XT0,YT0,ZT0,Tin,XX,YY,ZZ,XV)
! work out ion trap trajectories analytically 
! work out ion trap trajectories analytically 
!**************************************************************************
!**************************************************************************
IMPLICIT NONE 
integer  ::ii,jj
double precision::OMTRAPsq,omega,secular
double precision:: XV,coeff
double precision:: Qx,Qy,Qz,XT0,YT0,ZT0,Tin,XX,YY,ZZ

Qx=2*OMTRAPsq/omega/omega
coeff=Qx*omega/2/sqrt(2.d0)
secular=coeff*Tin

! first term is secular motion second is micromotion
XX=XT0*cos(secular)+ Qx*XT0/2*cos(secular)*cos(omega*Tin)

! velocity
XV=-coeff*XT0*sin(secular)-coeff*Qx*XT0/2*sin(secular)*cos(omega*Tin)
XV=XV-Qx*XT0/2*cos(secular)*omega*sin(omega*Tin)
Qy=Qx
YY=YT0*cos(secular)+ Qy*YT0/2*cos(secular)*cos(omega*Tin)
Qz=-2*Qx
ZZ=ZT0*cos(secular)+ Qz*ZT0/2*cos(secular)*cos(omega*Tin)

return
end

!*************************************************************************
!*************************************************************************
 SUBROUTINE EQUIL(E1,E2,A,DETUN1,DETUN2,GammaM,W2,XM,Kapp2,OMtrapsq,STATEeq,WK0,OMEGAM,omega,AMP,Vamp)
! subroutine below is provided by user and specifies initial state
! and calculates  relevant parameters 
!*************************************************************************
!*************************************************************************
                
IMPLICIT NONE 
integer  ::ii,m,jj,NTOT,NPERIOD
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,Press,TEMP
double precision::RHO,WK,PIN,PIN2,Q
double precision::RTRAP,V0,trapfreq,omega,DET

include 'PROBE2_OMIT.h'
double precision:: DETUN1,DETUN2,E1,E2,PHAS1,W2,W2M,welln,NWELLS
double precision:: KAPP,ENERGY,A,XM,KAPP2
double precision::pi,pi2
double precision::Polaris,VOL,OMOPT,OMTRAPsq,QMICRO,AMP,Vamp,Astab
double precision:: X0,WK0,OMEGAM,Gammam,omsec,omsecf,Beta
double precision::C1,C2,C3,C4,COOL1,ASQ1,ASQ2,XWELL,XEQ,WKXEQ     
DOUBLE PRECISION, DIMENSION(NTOT):: STATEeq

! zero eq. initial values
STATEeq=0.d0
PI=dacos(-1.d0)
XM=RHO*4.*pi/3.*R0**3
write(6,*)'Mass='
write(6,100) XM
Polaris=4.*pi*EPSI0*(EPSR-1.)/(EPSR+2.)*R0**3
write(6,*)'Polarisability='
write(6,100) Polaris
OMOPT=C*WK

! now waist is waist radius
W2=waist**2
VOL=XL*Pi*W2
A=OMOPT*POLARIS/2./VOL/EPSI0
write(6,*)'A/(2pi)='
write(6,100) A/pi/2.
KAPP2=pi*c/finesse/XL/2.d0
write(6,*)'kappa/2='
write(6,100) Kapp2

! trap beam equilibrium
E1=KAPP2*PIN/2./OMOPT/hbar
E1=sqrt(E1)

! probe beam equilibrium
E2=KAPP2*PIN2/2./OMOPT/hbar
E2=sqrt(E2)
write(6,*)'trap and probe beam amplitudes E1,E2='
write(6,100) E1,E2

! Peter 1.d-4 mBar => 0.125Hz
! now take usual expression eg Levitated review by Li Geraci etc
! 1 bar= 10^ 5 pascal; Press is in mbar = 10^ 2 pascal
! gamma=16 P/(pi*v*rho*R)
! v=speed of air=500 /s
GAMMAM=1600.*press/pi
GAMMAM=GAMMAM/500/RHO/R0
write(6,*)'mechanical damping* 2pi'
write(6,100)GAMMAM
OMTRAPsq=2.*Q*V0/XM/RTRAP**2
omega=trapfreq*2*pi
QMICRO=OMtrapsq/omega**2
write(6,*)'micromotion q/2='
write(6,100) QMICRO
Energy=3./2.d0*boltz*temp
AMP=2*energy/XM/OMtrapsq
AMP=sqrt(Amp)
write(6,*)'thermal amplitude of motion(m)in well numbers='
write(6,100) WK*AMP/2/pi
write(6,*) 'speed equivalent to thermal motion'
Vamp=sqrt(2.*Energy/XM)
write(6,100)vamp
WK0=0.d0

! We calculate equilibrium photon fields and x_0
! real part of photon field for trap
! this version uses 2*PHAS=pi/2 for equilibria
C1=KAPP2**2+(DETUN1+A*cos(WK0-phas1)**2)**2
STATEeq(1)=E1*(DETUN1+A*cos(WK0-phas1)**2)/C1

! Imaginary part of photon field 1
STATEeq(2)=-E1*KAPP2/C1

! |alpha1|^2
ASQ1=E1*E1/C1
write(6,*)'trap beam photon number in cavity='
write(6,100) ASQ1

! Set KE to a fraction of the optical trap depth.
Energy=hbar*A*ASQ1
Vamp=sqrt(2.*Energy/XM)
write(6,*)'velocity max equiv to opt.trap depth'
write(6,100) Vamp

! stability a parameter
! first work out 2/waist**2
W2M=2.d0/W2 

! trapped
! Astab=8*hbar*A*ASQ1*W2M/XM/omega**2
! untrapped
Astab=4*hbar*A*ASQ1*W2M/XM/omega**2
write(6,*)'Astability,q'
write(6,*) Astab,Qmicro*2
omsec=omega*Qmicro/sqrt(2.d0)/2/pi
Beta=sqrt(Astab+2*Qmicro**2)
Omsecf=0.5*Beta*omega/2/pi
write(6,*)'Astability,q/2,Beta'
write(6,*) Astab,Qmicro,Beta
write(6,*)' y secular frequency: uncorrected by field,corrected'
write(6,*) omsec,omsecf
Beta=sqrt(Astab+2*1.7**2*Qmicro**2)
Omsecf=0.5*Beta*omega/2/pi
write(6,*)' z secular frequency: uncorrected by field,corrected'
write(6,*) omsec*1.7,omsecf
write(6,*)'trapped'
Beta=sqrt(2*Astab+2*Qmicro**2)
Omsecf=0.5*Beta*omega/2/pi
write(6,*)'y TRAPPED  secular frequency:corrected'
write(6,*) omsec,omsecf

! equilibrium oscillation frequency
OMEGAM=2.*WK**2*A*hbar/XM
OMEGAM=OMEGAM*ASQ1*cos(2.*WK0)
OMEGAM=sqrt(OMEGAM)
write(6,*)'Mechanical frequency f in Hz'
write(6,100) omegam/2/pi
write(6,*)'Mechanical frequency omega M,detuning, kappa2'
write(6,100) omegam,detun1,kapp2

! analytical optomechanical cooling rate USING ONLY TRAP BEAM
! Formula given by Linear response theory or Perturbation theory
! given in 2012 PRA Rapid
COOL1=0.d0
C1=4*ASQ1*Kapp2*WK**2*A**2*hbar/omegaM/xm

! use the time averaged rate - average over 2pi/omega
! Gammopt=2*kappa2*A**2*hbar/omega/xm [k X_eq]**2
! Xeq= damped micromotion=omtrap**2/omegam**2 *ASQ1*X_well
! ************************************************
! X_well= position of well at which particle localises
! here assume it is well number welln; then rescale
! eg well 1000 gives (1000/welln)^2 times the cooling.
welln=100
WKXEQ=OMTRAPsq/OMEGAM**2*pi*welln
write(6,*)'Trap freq^2 omsq; omsq/omegam^2'
write(6,100)OMTRAPsq,OMTRAPsq/OMEGAM**2
C1=C1*(WKXEQ)**2
C2=DETUN1+A*(cos(WK0))**2
write(6,*)'shifted detuning=det+A'
write(6,100)C2/2.d0/pi

! here we neglect opto shift!
C3=(C2+omegam)**2+kapp2**2
C3=1.d0/C3
C4=(C2-omegam)**2+kapp2**2
C4=1.d0/C4
COOL1=-C1*(C3-C4)

! Halve the cooling rate because of the period average
! < cos omega t**2>=0.5
COOL1=COOL1*0.5
write(6,*)'cooling rate for well=',welln
write(6,*)'optomechanical cooling rate, mechanical damping,kXeq'
write(6,100)COOL1,GAMMAM,WKXEq
write(6,*)'resolved sideband cooling formula'
C3=OMTRAPsq/OMEGAM**2
C2=pi**2*A/Kapp2*welln**2
C2=C2*C3**2*OMEGAM

100   format(4D16.8)

RETURN
END


!*************************************************************************
!*************************************************************************
SUBROUTINE RK(nstep,STATE,STATEout,TDEL,TIME,DETUN1,DETUN2,XM,OMTRAPsq,W2,GAMMAM,A,E1,E2,OMEGA,kapp2,OFF)
! 4th order  Runge Kutta propagator to evolve for a time interval =TIME
!*************************************************************************
!*************************************************************************

IMPLICIT NONE 
integer  ::ll,kk,it,jj,nstep,NTOT,NPERIOD,idum
double precision:: cnoise,XNAV
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,Press,TEMP
double precision::RHO,WK,PIN,PIN2,Q
double precision::RTRAP,V0,trapfreq,DET

include 'PROBE2_OMIT.h'

! input parameters
double precision:: DETUN1,DETUN2,E1,E2,PHAS1,Omech
double precision:: KAPP2,omega,A,XM,ASQ
double precision:: OMTRAPsq,W2,W2M,GAMMAM,coeff,GRAVITY,OFF,OFFSET

! time propagation
double precision::pi,DT,DEL
double precision::TT,TIME,TDEL
double precision::etare,etaim,sigma,SNR,Xnoise,GAM
DOUBLE PRECISION, DIMENSION(NTOT):: STATE, STATEout
DOUBLE PRECISION, DIMENSION(NTOT):: STATE0,XKR

! seed for random number
! idum=1543
! coupling coefficient that gives frequency of mech osc
coeff=A*hbar*WK/XM
!GRAVITY=9.8d0
GRAVITY=0.d0
OFFSET=OFF
W2M=2.d0/W2

! estimate local/instantaneous mechanical frequency
! total number of photons
ASQ=STATEout(1)**2+STATEout(2)**2
Omech=2.*WK**2*A*hbar/XM

! ignore *cos(2.*WK0)
Omech=Omech*ASQ
Omech=sqrt(omech)

! work out noise scaling term
cnoise=hbar*omech/XM
cnoise=sqrt(cnoise)

! number of quanta
XNAV=BOLTZ*TEMP/hbar/omech
cnoise=sqrt(2.d0*XNAV+1.d0)*cnoise

! SNR=1.d7 makes particle jump out of well in old trap.1.d5 usual
! in this version, allow 1 millisecond to cool and establish equilibrium before turning on gravity and noise
SNR=sqrt(2.D0*KAPP2)
Xnoise=sqrt(gammam*2.d0)*cnoise*2.

! for our gaussian noise set variance to be sqrt(DT)
DT=TDEL/nstep
sigma=sqrt(DT)

do 1 ll=1,Ntot
    STATE0(ll)=STATE(ll)
    STATEout(ll)=STATE(ll)
1 enddo

TT=TIME
do 10 it=1,nstep
    ! XK1=F(tt,STATE(tt))*Dt
    DEL=0.d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 11 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/6.d0
        STATE(ll)=STATE0(ll)+XKR(ll)/2.d0
    11 enddo

    ! Xk2= F(tt+Dt/2,STATE(tt)+xk1/2)*Dt
    DEL=0.5d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 12 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/3.d0
        STATE(ll)=STATE0(ll)+XKR(ll)/2.d0
    12 enddo
    
    ! XK3=F(tt+Dt/2), STATE(tt) +xk2/2)*Dt
    DEL=0.5d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 13 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/3.d0
        STATE(ll)=STATE0(ll)+XKR(ll)
    13 enddo

    ! XK4=F(tt+Dt, STATE(tt)+xk3)*Dt
    DEL=1.d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 14 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/6.d0
    14 enddo 

    call Gasdev(etare,etaim,sigma,idum)
    ! write(6,*)etare,etaim,SNR
    ! add noise to trap optical field
    ! add noise to x and px
    !STATEout(1)=STATEout(1)+etaim*SNR
    !STATEout(3)=STATEout(3)+etare*xnoise
    ! call Gasdev(etare,etaim,sigma,idum)
    ! STATEout(3)=STATEout(3)+etare*SNR

    ! now Stateout contains psi(tt+DT). Reset state0
    do 16 ll=1,Ntot
        STATE0(ll)= STATEout(ll)
    16 enddo
    TT=TT+DT    
10 ENDDO

RETURN
END

!*************************************************************************
!*************************************************************************
SUBROUTINE FUNC(DT,STATE,XKR,TT,DEL,DETUN1,DETUN2,Coeff,OMTRAPsq,GAMMAM,kapp2,E1,E2,A,W2M,omega,GRAVITY,OFFSET)
! works out the derivatives d\alpha/dt, d\alpha*/dt , d^2x/dt^2,dx/dt,
! also multiplies by DT
!*************************************************************************
!*************************************************************************

IMPLICIT NONE
integer  ::ll,kk,it,jj,Nstep,NTOT,NPERIOD
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,Press,TEMP
double precision::RHO,WK,PIN,PIN2,Q
double precision::RTRAP,V0,trapfreq,omega,DET

include 'PROBE2_OMIT.h'

! input parameters
double precision:: DETUN1,DETUN2,E1,E2,RAT,W2M,OMTRAPsq,VOPT,GAMMAM
double precision:: KAPP2,A,XM

! time propagation
double precision::pi,pi2,DT,DS1,DS2,ASQ1,ASQ2,WKX,OFFSET,GRAVITY,Dprobe
double precision::TIME,TT,DEL,XX,velox,YY,Veloy,ZZ,Veloz
double precision:: Wide,AWIDE,Coeff,VION,VFIELD,VYZ,coswk,coswk2
double precision::ALPR1,ALPI1,ALPR2,ALPI2
DOUBLE PRECISION, DIMENSION(NTOT):: STATE
DOUBLE PRECISION, DIMENSION(NTOT):: XKR,DX

! in case of explicit time driving: TIME=TT+DT*DEL
! take C1 out into parameter file
! Vtrap=2*Q*V0/RTRAP**2/XM
! W2=(fullwaist/2)**2
! W2M=2.d0/W2
! Coeff=A*HBAR*WK/XM
! OFFSET is  of the ion trap relative to optical beam
! if trap is lower by 100 wells OFF=+/-532d-9*100=5.32d-5 m
! Define: Vion,W2,Vtrap,Wide,PHAS1,GammaM,Coeff
! write(6,100)W2,A,PHAS1,DETUN1,VOPT
! put in offset between ion trap and optical beam

TIME=TT+DT*DEL

! IN THIS VERSION try sin drive to catch fast
!VION=OMTRAPsq*cos(omega*time)
VION=OMTRAPsq*sin(omega*time)

Dprobe=DETUN2*time

ALPR1=STATE(1)
ALPI1=STATE(2)
Velox=STATE(3)
XX=STATE(4)
WKX=WK*XX
COSWK=COS(WKX)
cosWK2=COSWK*COSWK
Veloy=STATE(5)
YY=STATE(6)
Veloz=STATE(7)
ZZ=STATE(8)

!Wide=exp(-(YY*YY+ZZ*ZZ)*W2M)
wide=1.d0  
ASQ1=ALPR1*ALPR1+ALPI1*ALPI1
AWIDE=WIDE*ASQ1

! optical shift in detuning.In 3D depends on y,z
! VOPT=A*COS(WKX-PHAS1)**2*WIDE
VOPT=A*coswk2*WIDE

! Vfield=Coeff*sin(2.*WKX)*AWIDE
! add optical fields without PHAS1 in this model
Vfield=Coeff*sin(2.*WKX)*AWIDE
VYZ=-2*coeff/WK*coswk2*AWIDE*W2M

! both fields have same detuning
DS1=DETUN1+VOPT

! write(6,*)DETUN1,VOPT
! optical field, real and imaginary
! trap drive with iE1       

DX(1)=-DS1*ALPI1-KAPP2*ALPR1-E1-E2*cos(Dprobe)
DX(2)=DS1*ALPR1-KAPP2*ALPI1+E2*sin(Dprobe)

! x,y,z coords
DX(3)=-Vfield-GammaM*Velox-VION*XX
DX(4)=Velox

! with gravity and offset
DX(5)=-GammaM*Veloy - VION*(YY+OFFSET)+VYZ*YY-GRAVITY
DX(6)=Veloy

! modify for Peter trap coeff=1.7 not 2
DX(7)=-GammaM*Veloz+ 1.7*VION*ZZ+VYZ*ZZ
DX(8)=Veloz

! now multiply by *DT
do 30 ll=1,NTOT
    XKR(ll)=DX(ll)*DT
30 ENDDO

100 format(5E14.6)

RETURN
END

!*************************************************************************
!*************************************************************************
Subroutine Gasdev(etare,etaim,sigma,idum)
! Gaussian deviates
!*************************************************************************
!*************************************************************************

IMPLICIT NONE
!variable decs
double precision  v1,v2,R,Fac
double precision  sigma
double precision  etare,etaim
integer  idum
double precision x,y,ran1
1 call random_number(x)
call random_number(y)

! random number between -1,+1
v1=2.d0*(x)-1.d0
v2=2.d0*(y)-1.d0

! in unit circle
R=v1*v1+v2*v2

if (R.GE.1.0.OR.R.EQ.0.0) then
    goto 1
endif

Fac=sqrt(-2.d0*sigma*sigma*log(R)/R)
etare=v1*fac/dsqrt(2.d0)
etaim=v2*fac/dsqrt(2.d0)

end subroutine Gasdev

!*************************************************************************
!*************************************************************************
SUBROUTINE sfft(a, b, ntot, n, nspan, isn, ierr)
!*************************************************************************
!*************************************************************************
!  MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE USING
!  MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
!  BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
!  MODIFIED BY A. H. MORRIS, NSWC/DL, DAHLGREN VA
!  ARRAYS A AND B ORIGINALLY HOLD THE REAL AND IMAGINARY COMPONENTS OF
!  THE DATA, AND RETURN THE REAL AND IMAGINARY COMPONENTS OF THE
!  RESULTING FOURIER COEFFICIENTS.
!  MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN ARRAY ELEMENT
!  SUCCESSOR FUNCTION, WITHOUT LIMIT ON THE NUMBER OF IMPLIED MULTIPLE
!  SUBSCRIPTS.
!  THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
!  THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
!  NTOT IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
!  N IS THE DIMENSION OF THE CURRENT VARIABLE.
!  NSPAN/N IS THE SPACING OF CONSECUTIVE DATA VALUES
!  WHILE INDEXING THE CURRENT VARIABLE.
!  THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX EXPONENTIAL,
!  AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
!  A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3) IS COMPUTED BY
!  CALL SFFT(A, B, N1*N2*N3, N1, N1, 1, IERR)
!  CALL SFFT(A, B, N1*N2*N3, N2, N1*N2, 1, IERR)
!  CALL SFFT(A, B, N1*N2*N3, N3, N1*N2*N3, 1, IERR)
!  FOR A SINGLE-VARIATE TRANSFORM,
!  NTOT = N = NSPAN = (NUMBER OF COMPLEX DATA VALUES), E.G.
!  CALL SFFT(A, B, N, N, N, 1, IERR)
!  THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX ARRAY A,
!  THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO GIVE THE CORRECT INDEXING
!  INCREMENT AND A(2) USED TO PASS THE INITIAL ADDRESS FOR THE SEQUENCE
!  OF IMAGINARY VALUES, E.G.
!  CALL SFFT(A, A(2), NTOT, N, NSPAN, 2, IERR)
!  ARRAYS NFAC(MAXN), NP(MAXP), AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF)
!  ARE USED FOR TEMPORARY STORAGE.
!  MAXN MUST BE >= THE NUMBER OF FACTORS OF N
!  MAXF MUST BE >= THE MAXIMUM PRIME FACTOR OF N.
!  MAXP MUST BE > THE NUMBER OF PRIME FACTORS OF N.
!  IN ADDITION, MAXN IS ASSUMED TO BE ODD.
!  IF THE SQUARE-FREE PORTION K OF N HAS TWO OR MORE PRIME FACTORS,
!  THEN MAXP MUST BE >= K-1.
!  IERR IS A VARIABLE. IERR IS SET TO 0 IF NO INPUT ERRORS ARE
!  DETECTED. OTHERWISE, IERR IS ASSIGNED ONE OF THE VALUES
!  IERR=1    N IS LESS THAN 1
!  IERR=2    N HAS MORE THAN MAXN FACTORS
!  IERR=3    N HAS A PRIME FACTOR GREATER THAN MAXF OR THE SQUARE-FREE
!  PORTION OF N IS GREATER THAN MAXP+1

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)
REAL(dp), INTENT(IN OUT)  :: a(*)
REAL(dp), INTENT(IN OUT)  :: b(*)
INTEGER, INTENT(IN)       :: ntot
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: nspan
INTEGER, INTENT(IN)       :: isn
INTEGER, INTENT(OUT)      :: ierr

!  ARRAY STORAGE IN NFAC FOR A MAXIMUM OF 15 FACTORS OF N.
!  IF N HAS MORE THAN ONE SQUARE-FREE FACTOR, THE PRODUCT OF THE
!  SQUARE-FREE FACTORS MUST BE <= 210
INTEGER   :: nfac(15), np(209)

!  ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23
REAL(dp)  :: at(23), ck(23), bt(23), sk(23)

INTEGER   :: i, inc, j, jc, jf, jj, k, k1, k2, k3, k4, kk, ks, kspan, kspnn,  &
kt, l, m, max, maxf, maxn, maxp, nn, nt, num
REAL(dp)  :: aa, aj, ajm, ajp, ak, akm, akp, bb, bj, bjm, bjp, bk, bkm, bkp,  &
c1, c2, c3, c72, cd, rad, radf, s1, s2, s3, s72, s120, sd, u, v

!  EQUIVALENCE (i,ii)
!  THE FOLLOWING CONSTANTS SHOULD AGREE WITH THE ARRAY DIMENSIONS.
maxn = 15
maxf = 23
maxp = 209
!  SET THE FOLLOWING CONSTANTS
!     RAD = 2.0*PI
!     S72 = SIN(RAD/5.0)
!     C72 = COS(RAD/5.0)
!     S120 = SQRT(0.75)
rad = 6.2831853071796_dp
s72 = .951056516295154_dp
c72 = .309016994374947_dp
s120 = .86602540378444_dp

ierr = 0
IF(n-1 < 0) THEN
GO TO  1000
ELSE IF (n-1 == 0) THEN
GO TO   960
END IF
inc = isn
IF(isn >= 0) GO TO 10
s72 = -s72
s120 = -s120
rad = -rad
inc = -inc
10 nt = inc*ntot
ks = inc*nspan
kspan = ks
nn = nt-inc
jc = ks/n
radf = rad*jc*0.5
i = 0
jf = 0

!  DETERMINE THE FACTORS OF N
m = 0
k = n
MAX = maxn/2
GO TO 20

15 IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = 4
k = l
20 l = k/16
IF(k == l*16) GO TO 15
j = 3
jj = 9
GO TO 30

25 IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = j
k = k/jj
30 IF(MOD(k,jj) == 0) GO TO 25
j = j+2
jj = j**2
IF(j <= maxf .AND. jj <= k) GO TO 30
IF(k > 4) GO TO 40
kt = m
nfac(m+1) = k
IF(k /= 1) m = m+1
GO TO 80

40 l = k/4
IF(k /= l*4) GO TO 50
IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = 2
k = l
kt = m
IF(k == 1) GO TO 85
50 kt = m
IF(k-1 > maxp) GO TO 1002
num = maxn-kt-kt
j = 2
60 IF(MOD(k,j) /= 0) GO TO 70
m = m+1
nfac(m) = j
num = num-1
k = k/j
IF(k == 1) GO TO 80
IF(num <= 0) GO TO 1001
70 l = (j+1)/2
j = l+l+1
IF(j <= maxf) GO TO 60
GO TO 1002

80 IF(kt == 0) GO TO 100
85 j = kt
90 m = m+1
nfac(m) = nfac(j)
j = j-1
IF(j /= 0) GO TO 90

!  COMPUTE FOURIER TRANSFORM
100 sd = radf/kspan
cd = 2.0*SIN(sd)**2
sd = SIN(sd+sd)
kk = 1
i = i+1
IF(nfac(i) /= 2) GO TO 400

!  TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
kspan = kspan/2
k1 = kspan+2
210 k2 = kk+kspan
ak = a(k2)
bk = b(k2)
a(k2) = a(kk)-ak
b(k2) = b(kk)-bk
a(kk) = a(kk)+ak
b(kk) = b(kk)+bk
kk = k2+kspan
IF(kk <= nn) GO TO 210
kk = kk-nn
IF(kk <= jc) GO TO 210
IF(kk > kspan) GO TO 800
220 c1 = 1.0-cd
s1 = sd
230 k2 = kk+kspan
ak = a(kk)-a(k2)
bk = b(kk)-b(k2)
a(kk) = a(kk)+a(k2)
b(kk) = b(kk)+b(k2)
a(k2) = c1*ak-s1*bk
b(k2) = s1*ak+c1*bk
kk = k2+kspan
IF(kk < nt) GO TO 230
k2 = kk-nt
c1 = -c1
kk = k1-k2
IF(kk > k2) GO TO 230
u = sd*s1+cd*c1
v = sd*c1-cd*s1
ak = c1-u
s1 = s1+v

!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN ONE MAY SUBSTITUTE
!     C1 = AK
c1 = 1.5-0.5*(ak*ak+s1*s1)
s1 = c1*s1
c1 = c1*ak
kk = kk+jc
IF(kk < k2) GO TO 230
k1 = k1+inc+inc
kk = (k1-kspan)/2+jc
IF(kk <= jc+jc) GO TO 220
GO TO 100

!  TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
320 k1 = kk+kspan
k2 = k1+kspan
ak = a(kk)
bk = b(kk)
aj = a(k1)+a(k2)
bj = b(k1)+b(k2)
a(kk) = ak+aj
b(kk) = bk+bj
ak = -0.5*aj+ak
bk = -0.5*bj+bk
aj = (a(k1)-a(k2))*s120
bj = (b(k1)-b(k2))*s120
a(k1) = ak-bj
b(k1) = bk+aj
a(k2) = ak+bj
b(k2) = bk-aj
kk = k2+kspan
IF(kk < nn) GO TO 320
kk = kk-nn
IF(kk <= kspan) GO TO 320
GO TO 700

!  TRANSFORM FOR FACTOR OF 4
400 IF(nfac(i) /= 4) GO TO 600
kspnn = kspan
kspan = kspan/4
410 c1 = 1.0
s1 = 0.0
420 k1 = kk+kspan
k2 = k1+kspan
k3 = k2+kspan
akp = a(kk)+a(k2)
akm = a(kk)-a(k2)
ajp = a(k1)+a(k3)
ajm = a(k1)-a(k3)
a(kk) = akp+ajp
ajp = akp-ajp
bkp = b(kk)+b(k2)
bkm = b(kk)-b(k2)
bjp = b(k1)+b(k3)
bjm = b(k1)-b(k3)
b(kk) = bkp+bjp
bjp = bkp-bjp
IF(isn < 0) GO TO 450
akp = akm-bjm
akm = akm+bjm
bkp = bkm+ajm
bkm = bkm-ajm
IF(s1 == 0.0) GO TO 460
430 a(k1) = akp*c1-bkp*s1
b(k1) = akp*s1+bkp*c1
a(k2) = ajp*c2-bjp*s2
b(k2) = ajp*s2+bjp*c2
a(k3) = akm*c3-bkm*s3
b(k3) = akm*s3+bkm*c3
kk = k3+kspan
IF(kk <= nt) GO TO 420

440 u = sd*s1+cd*c1
v = sd*c1-cd*s1
c2 = c1-u
s1 = s1+v
!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN ONE MAY SUBSTITUTE
!     C1 = C2
c1 = 1.5-0.5*(c2*c2+s1*s1)
s1 = c1*s1
c1 = c1*c2
c2 = c1*c1-s1*s1
s2 = 2.0*c1*s1
c3 = c2*c1-s2*s1
s3 = c2*s1+s2*c1
kk = kk-nt+jc
IF(kk <= kspan) GO TO 420
kk = kk-kspan+inc
IF(kk <= jc) GO TO 410
IF(kspan == jc) GO TO 800
GO TO 100

450 akp = akm+bjm
akm = akm-bjm
bkp = bkm-ajm
bkm = bkm+ajm
IF(s1 /= 0.0) GO TO 430
460 a(k1) = akp
b(k1) = bkp
a(k2) = ajp
b(k2) = bjp
a(k3) = akm
b(k3) = bkm
kk = k3+kspan
IF(kk <= nt) GO TO 420
GO TO 440

!  TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
510 c2 = c72**2-s72**2
s2 = 2.0*c72*s72
520 k1 = kk+kspan
k2 = k1+kspan
k3 = k2+kspan
k4 = k3+kspan
akp = a(k1)+a(k4)
akm = a(k1)-a(k4)
bkp = b(k1)+b(k4)
bkm = b(k1)-b(k4)
ajp = a(k2)+a(k3)
ajm = a(k2)-a(k3)
bjp = b(k2)+b(k3)
bjm = b(k2)-b(k3)
aa = a(kk)
bb = b(kk)
a(kk) = aa+akp+ajp
b(kk) = bb+bkp+bjp
ak = akp*c72+ajp*c2+aa
bk = bkp*c72+bjp*c2+bb
aj = akm*s72+ajm*s2
bj = bkm*s72+bjm*s2
a(k1) = ak-bj
a(k4) = ak+bj
b(k1) = bk+aj
b(k4) = bk-aj
ak = akp*c2+ajp*c72+aa
bk = bkp*c2+bjp*c72+bb
aj = akm*s2-ajm*s72
bj = bkm*s2-bjm*s72
a(k2) = ak-bj
a(k3) = ak+bj
b(k2) = bk+aj
b(k3) = bk-aj
kk = k4+kspan
IF(kk < nn) GO TO 520
kk = kk-nn
IF(kk <= kspan) GO TO 520
GO TO 700

!  TRANSFORM FOR ODD FACTORS
600 k = nfac(i)
kspnn = kspan
kspan = kspan/k
IF(k == 3) GO TO 320
IF(k == 5) GO TO 510
IF(k == jf) GO TO 640
jf = k
s1 = rad/k
c1 = COS(s1)
s1 = SIN(s1)
ck(jf) = 1.0
sk(jf) = 0.0
j = 1
630 ck(j) = ck(k)*c1+sk(k)*s1
sk(j) = ck(k)*s1-sk(k)*c1
k = k-1
ck(k) = ck(j)
sk(k) = -sk(j)
j = j+1
IF(j < k) GO TO 630
640 k1 = kk
k2 = kk+kspnn
aa = a(kk)
bb = b(kk)
ak = aa
bk = bb
j = 1
k1 = k1+kspan
650 k2 = k2-kspan
j = j+1
at(j) = a(k1)+a(k2)
ak = at(j)+ak
bt(j) = b(k1)+b(k2)
bk = bt(j)+bk
j = j+1
at(j) = a(k1)-a(k2)
bt(j) = b(k1)-b(k2)
k1 = k1+kspan
IF(k1 < k2) GO TO 650
a(kk) = ak
b(kk) = bk
k1 = kk
k2 = kk+kspnn
j = 1
660 k1 = k1+kspan
k2 = k2-kspan
jj = j
ak = aa
bk = bb
aj = 0.0
bj = 0.0
k = 1
670 k = k+1
ak = at(k)*ck(jj)+ak
bk = bt(k)*ck(jj)+bk
k = k+1
aj = at(k)*sk(jj)+aj
bj = bt(k)*sk(jj)+bj
jj = jj+j
IF(jj > jf) jj = jj-jf
IF(k < jf) GO TO 670
k = jf-j
a(k1) = ak-bj
b(k1) = bk+aj
a(k2) = ak+bj
b(k2) = bk-aj
j = j+1
IF(j < k) GO TO 660
kk = kk+kspnn
IF(kk <= nn) GO TO 640
kk = kk-nn
IF(kk <= kspan) GO TO 640

!  MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
700 IF(i == m) GO TO 800
kk = jc+1
710 c2 = 1.0-cd
s1 = sd
720 c1 = c2
s2 = s1
kk = kk+kspan
730 ak = a(kk)
a(kk) = c2*ak-s2*b(kk)
b(kk) = s2*ak+c2*b(kk)
kk = kk+kspnn
IF(kk <= nt) GO TO 730
ak = s1*s2
s2 = s1*c2+c1*s2
c2 = c1*c2-ak
kk = kk-nt+kspan
IF(kk <= kspnn) GO TO 730
u = sd*s1+cd*c1
v = sd*c1-cd*s1
c2 = c1-u
s1 = s1+v
!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN THEY MAY BE DELETED.
c1 = 1.5-0.5*(c2*c2+s1*s1)
s1 = c1*s1
c2 = c1*c2
kk = kk-kspnn+jc
IF(kk <= kspan) GO TO 720
kk = kk-kspan+jc+inc
IF(kk <= jc+jc) GO TO 710
GO TO 100

!  PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES
!  PERMUTATION FOR SQUARE FACTORS OF N
800 np(1) = ks
IF(kt == 0) GO TO 890
k = kt+kt+1
IF(m < k) k = k-1
j = 1
np(k+1) = jc
810 np(j+1) = np(j)/nfac(j)
np(k) = np(k+1)*nfac(j)
j = j+1
k = k-1
IF(j < k) GO TO 810
k3 = np(k+1)
kspan = np(2)
kk = jc+1
k2 = kspan+1
j = 1
IF(n /= ntot) GO TO 850

!  PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
820 ak = a(kk)
a(kk) = a(k2)
a(k2) = ak
bk = b(kk)
b(kk) = b(k2)
b(k2) = bk
kk = kk+inc
k2 = kspan+k2
IF(k2 < ks) GO TO 820
830 k2 = k2-np(j)
j = j+1
k2 = np(j+1)+k2
IF(k2 > np(j)) GO TO 830
j = 1
840 IF(kk < k2) GO TO 820
kk = kk+inc
k2 = kspan+k2
IF(k2 < ks) GO TO 840
IF(kk < ks) GO TO 830
jc = k3
GO TO 890

!  PERMUTATION FOR MULTIVARIATE TRANSFORM
850 k = kk+jc
860 ak = a(kk)
a(kk) = a(k2)
a(k2) = ak
bk = b(kk)
b(kk) = b(k2)
b(k2) = bk
kk = kk+inc
k2 = k2+inc
IF(kk < k) GO TO 860
kk = kk+ks-jc
k2 = k2+ks-jc
IF(kk < nt) GO TO 850
k2 = k2-nt+kspan
kk = kk-nt+jc
IF(k2 < ks) GO TO 850
870 k2 = k2-np(j)
j = j+1
k2 = np(j+1)+k2
IF(k2 > np(j)) GO TO 870
j = 1
880 IF(kk < k2) GO TO 850
kk = kk+jc
k2 = kspan+k2
IF(k2 < ks) GO TO 880
IF(kk < ks) GO TO 870
jc = k3
890 IF(2*kt+1 >= m) RETURN
kspnn = np(kt+1)

!  PERMUTATION FOR SQUARE-FREE FACTORS OF N
j = m-kt
nfac(j+1) = 1
900 nfac(j) = nfac(j)*nfac(j+1)
j = j-1
IF(j /= kt) GO TO 900
kt = kt+1
nn = nfac(kt)-1
jj = 0
j = 0
GO TO 906

902 jj = jj-k2
k2 = kk
k = k+1
kk = nfac(k)
904 jj = kk+jj
IF(jj >= k2) GO TO 902
np(j) = jj
906 k2 = nfac(kt)
k = kt+1
kk = nfac(k)
j = j+1
IF(j <= nn) GO TO 904

!  DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
j = 0
GO TO 914

910 k = kk
kk = np(k)
np(k) = -kk
IF(kk /= j) GO TO 910
k3 = kk
914 j = j+1
kk = np(j)
IF(kk < 0) GO TO 914
IF(kk /= j) GO TO 910
np(j) = -j
IF(j /= nn) GO TO 914
maxf = inc*maxf
!  REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
GO TO 950

924 j = j-1
IF(np(j) < 0) GO TO 924
jj = jc
926 kspan = jj
IF(jj > maxf) kspan = maxf
jj = jj-kspan
k = np(j)
kk = jc*k+i+jj
k1 = kk+kspan
k2 = 0
928 k2 = k2+1
at(k2) = a(k1)
bt(k2) = b(k1)
k1 = k1-inc
IF(k1 /= kk) GO TO 928
932 k1 = kk+kspan
k2 = k1-jc*(k+np(k))
k = -np(k)
936 a(k1) = a(k2)
b(k1) = b(k2)
k1 = k1-inc
k2 = k2-inc
IF(k1 /= kk) GO TO 936
kk = k2
IF(k /= j) GO TO 932
k1 = kk+kspan
k2 = 0
940 k2 = k2+1
a(k1) = at(k2)
b(k1) = bt(k2)
k1 = k1-inc
IF(k1 /= kk) GO TO 940
IF(jj /= 0) GO TO 926
IF(j /= 1) GO TO 924

950 j = k3+1
nt = nt-kspnn
i = nt-inc+1
IF(nt >= 0) GO TO 924
960 RETURN

!  ERROR FINISH - THERE IS AN INPUT ERROR
1000 ierr = 1
RETURN
1001 ierr = 2
RETURN
1002 ierr = 3
RETURN
END SUBROUTINE sfft

!*************************************************************************
!*************************************************************************
SUBROUTINE NORM(NPTS,Total,SXXR,OMSTOR)
!*************************************************************************
!*************************************************************************
                
IMPLICIT NONE 
integer  ::ii,jj,NTOT,NPERIOD,NPTS
double precision::pi,pi2,DEL,XRE,XIM,OM,Total,F1,F2
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,DET,Press,TEMP
double precision::RHO,WK,PIN,PIN2,Q
double precision::RTRAP,V0,trapfreq

include 'PROBE2_OMIT.h'

DOUBLE PRECISION, DIMENSION(NPTS):: OMSTOR
DOUBLE PRECISION:: SXXR(NPTS), SXXI(NPTS)

pi2=2.d0*dacos(-1.d0)

! integrate between two frequencies
F1=-2.5d5
F2=-2.d5

! integrate the  spectrum of bead
! so far only trapezoidal rule- improve later
XRE=0.d0
XIM=0.d0
DEL=ABS(OMSTOR(2)-OMSTOR(1))
do ii=1,NPTS-1
    if((OMstor(ii).gt.F1).and.(OMstor(ii).lt.F2))then 
        XRE=XRE+0.5d0*( SXXR(ii)+SXXR(ii+1))
    endif
enddo

Total=XRE*DEL/abs(F1-F2)

100 format(2I3,1x,6E14.6)

return
end 
 
!*************************************************************************
!*************************************************************************
SUBROUTINE RKINI(nstep,STATE,STATEout,TDEL,TIME,DETUN1,DETUN2,XM,OMTRAPsq,W2,GAMMAM,A,E1,E2,OMEGA,kapp2,OFF)
! 4th order  Runge Kutta propagator to evolve for a time interval =TIME
! this version has no noise and works with high pressure to capture
! the particle
!*************************************************************************
!*************************************************************************
IMPLICIT NONE 
integer  ::ll,kk,it,jj,nstep,NTOT,NPERIOD,idum
double precision:: cnoise,XNAV
double precision::R0,EPSR,EPSI0,C,hbar,BOLTZ
double precision::waist,XL,Finesse,Press,TEMP
double precision::RHO,WK,PIN,PIN2,Q
double precision::RTRAP,V0,trapfreq,DET

include 'PROBE2_OMIT.h'
! input parameters

double precision:: DETUN1,DETUN2,E1,E2,PHAS1,Omech
double precision:: KAPP2,omega,A,XM,ASQ
double precision:: OMTRAPsq,W2,W2M,GAMMAM,coeff,GRAVITY,OFFSET,OFF

! time propagation
double precision::pi,DT,DEL
double precision::TT,TIME,TDEL
double precision::etare,etaim,sigma,SNR,Xnoise,GAM,GRAV
DOUBLE PRECISION, DIMENSION(NTOT):: STATE, STATEout
DOUBLE PRECISION, DIMENSION(NTOT):: STATE0,XKR

! coupling coefficient that gives frequency of mech osc
coeff=A*hbar*WK/XM
!GRAV=9.8d0
GRAV=0.d0
OFFSET=OFF
W2M=2.d0/W2

! in this version, allow 1 millisecond to cool and establish equilibrium before turning on gravity and noise
Xnoise=sqrt(gammam*2.d0)*cnoise*2.
GAM=gammam

! ramp up gravity and offset slowly so system equilibrates
gravity=GRAV*time/TDEL/20
offset=OFF*time/TDEL/20

! for our gaussian noise set variance to be sqrt(DT)
DT=TDEL/nstep
sigma=sqrt(DT)

do 1 ll=1,Ntot
    STATE0(ll)=STATE(ll)
    STATEout(ll)=STATE(ll)
1 enddo

TT=TIME
do 10 it=1,nstep
    ! XK1=F(tt,STATE(tt))*Dt
    DEL=0.d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 11 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/6.d0
        STATE(ll)=STATE0(ll)+XKR(ll)/2.d0
    11 enddo

    ! Xk2= F(tt+Dt/2,STATE(tt)+xk1/2)*Dt
    DEL=0.5d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 12 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/3.d0
        STATE(ll)=STATE0(ll)+XKR(ll)/2.d0
    12 enddo

    ! XK3=F(tt+Dt/2), STATE(tt) +xk2/2)*Dt
    DEL=0.5d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 13 ll=1,Ntot
        STATEout(ll)=STATEout(ll)+XKR(ll)/3.d0
        STATE(ll)=STATE0(ll)+XKR(ll)
    13 enddo

    ! XK4=F(tt+Dt, STATE(tt)+xk3)*Dt
    DEL=1.d0
    CALL FUNC(DT,STATE,XKR,TT,DEL,DETUN1,Detun2,Coeff,OMTRAPsq,GAM,kapp2,E1,E2,A,W2M,&
    omega,GRAVITY,offset)

    ! update Stateout with runge Kutta increment
    do 14 ll=1,Ntot
    STATEout(ll)=STATEout(ll)+XKR(ll)/6.d0
    14 enddo 

    ! now Stateout contains psi(tt+DT). Reset state0
    do 16 ll=1,Ntot
        STATE0(ll)= STATEout(ll)
    16 enddo

    TT=TT+DT
10 ENDDO

RETURN
END