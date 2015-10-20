! input parameters
! NPERIOD=number of optical trap periods
! NTOT=number of equations so 8=>1 optical mode +3D
! R0= sphere radius
! XL=cavity length 
! Pin=input power in Watts
! Press =air pressure in millibars
! XL=cavity length, waist
! rho=sphere density
! WK=2*pi/lambda=k
! DET=detuning in KHz
! Q= sphere charge
! Press= gas pressure in millibar
! waist is  waist radius
! RR0 and V0 = ion trap parameters
PARAMETER(NPERIOD=10000,NTOT=8,R0=200*1.d-9,&
EPSR=1.45d0**2,Epsi0=8.854d-12,c=3.d8,hbar=1.05d-34,BOLTZ=1.4d-23,&
! old mirrors
! waist=140.*1.d-6,XL=3.7d-2,Finesse=15000.d0,DET=-0.288d6,&

! new mirrors
!waist=60.d0*1.d-6,XL=1.3d-2,Finesse=0.4d5,DET=-50.d3,&
waist=60.d0*1.d-6,XL=1.3d-2,Finesse=0.4d5,&
!Press=1.e-4,rho=2198.,WK=5.9d6,Pin=2.15d-3,PIN2=0.01*Pin,&
Press=1.e-4,rho=2198.,WK=5.9d6,&
Q=1.d0*1.6d-19,RTRAP=1.d-3/sqrt(2.d0),V0=600.d0,trapfreq=1500.,TEMP=300.,mechfreq=50000.d0,G=1d8)