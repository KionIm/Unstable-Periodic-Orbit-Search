module orbit
   implicit none
   save
   integer(8), parameter :: MM = 21
   integer, parameter :: n = 4*(MM+1)*(MM+1) -5 + 1		! dimension of system, 
                              !   including unknown period T
   integer            :: mgmres = 300	! max GMRES iterations
   integer            :: nits = 5000	! max Newton iterations
   double precision   :: rel_err = 1d-8	! relative error |F|/|x|

   double precision :: del = -1d0	! These rarely need changing
   double precision :: mndl = 1d-20	!   for any problem
   double precision :: mxdl = 1d+20
   !double precision :: gtol = 1d-4
   double precision :: gtol = 1d-3 !Hookstep
   double precision :: epsJ = 1d-7 !1d-7

   double precision :: tol	
   integer :: info
   integer :: ndts	! number of timesteps to take in period T
   logical :: fixT,fixLamda	! Fix T for equilibrium, rather than PO solution

   REAL(8):: VAR(4*(MM+1)*(MM+1))

end module orbit

module definition
   use orbit
   implicit none
   save
   REAL(8),PARAMETER::PI=3.1415926535897932385D0
   REAL(8),PARAMETER::SQRT3=1.7320508075688772935D0
   INTEGER(8),PARAMETER::NM=MM+1
   INTEGER(8),PARAMETER::LM=(MM+1)*(MM+1)
   !INTEGER(8),PARAMETER::JM=16,IM=32
   !INTEGER(8),PARAMETER::JM=18,IM=36
   INTEGER(8),PARAMETER::JM=32,IM=64
   REAL(8)::W1(JM*IM),W2(JM*IM*2),WORK(LM*4,3),DVAR(LM*4)
   INTEGER(8)::II
   REAL(8)::TIM
   REAL(8)::GAM
   REAL(8)::T(IM*3/2),R(((MM+1)*(2*NM-MM-1)+1)/4*3+(2*NM-MM)*(MM+1)/2+MM+1)
   INTEGER(8)::IT(IM/2), JC(MM*(2*NM-MM-1)/16+MM),IPOW,IFLAG
   REAL(8)::P(JM/2,2*MM+5)
   REAL(8)::C(MM+1+MM*(2*MM-MM+1))
   REAL(8)::D((MM+1+MM*(2*MM-MM+1)),2)
   REAL(8)::XG(0:IM-1),YG(JM)
   REAL(8)::RN(LM,2),DRN(LM),DOM(LM)
   REAL(8)::CD(LM),SD(LM)
   REAL(8)::CP(LM),SP(LM)
   REAL(8)::AVTT(LM),AVTC(LM),DIV(LM),TMP(LM)
   REAL(8),DIMENSION(MM*MM+4*MM+2)::STFC,STFT,STFTX,STFCX,VPO,VPOX,TMPX,TMPY,SXTMP
   REAL(8),DIMENSION(MM*MM+4*MM+2)::VPOX2,STFTX2,STFCX2,STFTY,STFCY,VPOY
   REAL(8),DIMENSION(MM*MM+4*MM+2)::AVTINIS,GEOINIS,S11,S13,S17,S18,S16
   REAL(8),DIMENSION(MM*MM+4*MM+2)::S1,S2,S3,S4,S5,S6,S7,S9,SX9,SX11,SX13,SLTMP,SL152
   REAL(8),DIMENSION(MM*MM+4*MM+2)::S8,SL92,SX92,SX112,SX132,S152,S15,SL15
   REAL(8),DIMENSION(MM*MM+4*MM+2)::S10,S12,S14,SY10,SY12,SY14,SY102,SY122,SY142
   REAL(8),DIMENSION(LM)::DAVTT,DAVTC,DDIV,DTMP
   REAL(8),DIMENSION(0:IM-1,JM)::G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G17,G18
   REAL(8),DIMENSION(0:IM-1,JM)::GY6,GY5,GY8,G,AVTINIG,EYE,G14,G15,G16,G152
   REAL(8),DIMENSION(0:IM-1,JM)::GTMPX,GTMPY
   REAL(8)::TIMESERIES,Lamda0,Period,IniShift,Adjust

   INTEGER(8)::IU = 1
   !---- OPEN SUBROUTINE PACKAGE AND INITIALIZE VARIABLES -------------
   
   INTEGER(8),PARAMETER::LEV=4                  !高階粘性項のラプラシアンの階数
   INTEGER(8)::ITM=1              !時間発展するステップ数
   INTEGER(8)::NDV=1                 !Runge-Kuttaでのステップ分割数
   !REAL(8)::DT=1/8D0                 !ファイル出力の時間間隔
   REAL(8)::SIGMA = 1D0
   
   REAL(8)::DNU = 10D0/(1D0*MM*(MM+1))**LEV !高階粘性の係数  !10D0/(1D0*MM*(MM+1)-2)**LEV
   
   !----/ガウシアン型の擾乱のパラメター(ジオポテンシャルに与える) /
   
   REAL(8)::OMEGA = 2*PI
   REAL(8)::BARPHI = 30D0
   REAL(8)::DA=1D-1                   !擾乱の振幅
   REAL(8)::DB=10D0                   !擾乱の幅の逆数のルートに相当するパラメター
   REAL(8)::X0=PI                     !擾乱の中心位置のλ座標
   REAL(8)::Y0=PI/6D0                      !擾乱の中心位置のφ座標

   REAL(8)::KAP = 2/7D0

   !------/Held&Suarez
   REAL(8),parameter::kf = 1D0 !(/day)
   REAL(8),parameter::ks = 1/4D0 !(/day)
   REAL(8),parameter::Rd = 287.05 !(J/K/kg) !(K*m^2/s^2)
   REAL(8),parameter::a = 6.371D6 !(m)
   REAL(8),parameter::day = 86400 !(s)
   REAL(8),parameter::TBAR = 295*Rd*(day*day)/(a*a) !()
   REAL(8),parameter:: S = 5/7D0*TBAR
   REAL(8),parameter::That_eq00 = 0D0                     
   !That_eq02 = (-32*SQRT(5D0))*Rd*(day*day)/(a*a) !() about -76*      -80*sqrt(2)   *sqrt(2)/sqrt(5)
   !That_eq00 = (-10*(SQRT(2D0) - 1D0)/SQRT(2D0))*Rd*(day*day)/(a*a) !()
   REAL(8)::That_eq02 ! 
   INTEGER(8)::Md 
   INTEGER(8)::Nd 
   REAL(8)::new_yF(n)

   contains
   
   SUBROUTINE SXOPEN(LEV_,DNU_,NDV_)
      INTEGER(8)::LEV_,NDV_
      REAL(8)::DNU_,DTDNU,DTINT
      INTEGER(8)::J,L_,I_

      CALL SXINI1(MM,MM+1,IM,IT,T,R)
      CALL SXINI2(MM,MM+1,JM,1_8,P,R,JC)    
      CALL SXINIC(MM,MM,C)
      CALL SXINID(MM,MM,D)
  
      DO J=1,JM/2
          YG(JM/2+J)= ASIN(P(J,1))
          YG(JM/2-J+1)=-ASIN(P(J,1))
      END DO

      DO I_=0,IM-1
          XG(I_)=2*PI*I_/IM
      END DO

      RETURN
   
   END SUBROUTINE SXOPEN

end module definition

!*************************************************************************
 include 'NewtonHook.f90'
 include 'GMRESm.f90'
  PROGRAM MAIN
 !*************************************************************************
    use orbit
    use definition
    use newton
    implicit none
    external :: getrhs, multJ, multJp, saveorbit 
    double precision, external :: dotprod
    double precision :: d2
    intEGER(8) :: L_,I_,N_,M_
    real(8)::VAR2(4*LM),Shift2,norm_x
    
    allocate(new_x(n))
    allocate(new_fx(n))  
 
 ! 		 Initial guesses for PO, initial T,X,Y,Z, num timesteps
 !   new_x = (/1.55, -13.7, -19.5, 27./)
    !OPEN(1_8,FILE='Savedata8_124_400.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata8_120_400.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata12_108.5_400.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Fin3.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata11_112.6_1500.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Fin4.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata8_184_800_L3.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Fin6.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata21_400_93.dat',FORM='UNFORMATTED')
    !OPEN(1_8,FILE='Savedata21_750_93.dat',FORM='UNFORMATTED')
    OPEN(1_8,FILE='Savedata21_900_93.dat',FORM='UNFORMATTED')
    READ(1_8) VAR
    CLOSE(1_8)

    !do L_ = 1,LM
      !CALL SXL2NM(MM,L_,N_,M_)
      !print*,M_,N_,VAR(L_),VAR(LM+L_),VAR(2*LM+L_),VAR(3*LM+L_)
    !END do
    !stop
    !OPEN(400_8,FILE='PeriLam6.dat',FORM='UNFORMATTED')
    !READ(400_8) Period,Lamda0
    !CLOSE(400_8)
    OPEN(1_8,FILE='PeriLam6.dat',FORM='UNFORMATTED')
    READ(1_8) Period,Lamda0
    CLOSE(1_8)

    !Period = 133d0 
    !Period = 6d0/300d0*566d0
    !Period = 36.5d0
    !Period = 37.0723305
    !Period = 37.0723305*2  !Semi
    !Period = 37.5
    !Period = 17.5
    !Period = 36.5
    !Period = Period/18d0
    Period = 21
    !Period = 5.5

    Md = 1
    !Nd = 3
    !Nd = 4!
    !Nd = 3
    Nd = 2
    !That_eq02 = (-120*SQRT(5D0))*Rd*(day*day)/(a*a)
    !That_eq02 = (-124*SQRT(5D0))*Rd*(day*day)/(a*a)
    !That_eq02 = (-108.5*SQRT(5D0))*Rd*(day*day)/(a*a)
    !That_eq02 = (-112.6*SQRT(5D0))*Rd*(day*day)/(a*a)
    !That_eq02 = (-184*SQRT(5D0))*Rd*(day*day)/(a*a)
    That_eq02 = (-93*SQRT(5D0))*Rd*(day*day)/(a*a)

    Adjust = 0!-0.5*PI

    CALL SXOPEN(LEV,DNU,NDV)

    CALL SXNM2L(MM,Nd,Md,L_)
    IniShift = atan(VAR(3*LM + L_+1)/VAR(3*LM + L_))/dble(Md)   

    DO I_= 1,4
        CALL SHIFT(VAR((I_-1)*LM+1:I_*LM),IniShift)
    END DO

    !print*,"here"
    print*,VAR(3*LM + L_+1)

    OPEN(1_8,FILE='IniShift.dat',FORM='UNFORMATTED')
    Write(1_8) VAR
    CLOSE(1_8)

    print*,"Ini12",VAR(3*LM + L_+1),VAR(3*LM + L_)

    N_ = 4
    M_ = 2
    CALL SXNM2L(MM,N_,M_,L_)

    !print*,"Ini24",VAR(3*LM + L_+1),VAR(3*LM + L_)

    CALL new_xSet(new_x,Period,Var,MM,Md,Nd)
   
    ndts = nint(48*Period)
    !print*,ndts
    fixT = .false. 
    !fixT = .True. 
    !fixLamda = .True.   
                   ! scale params by |x|
    d2 = dotprod(-1,new_x,new_x)
    tol  = rel_err * dsqrt(d2)
    del  = del     * dsqrt(d2)
    mndl = mndl    * dsqrt(d2)
    mxdl = mxdl    * dsqrt(d2)

    info = 1
    call newtonhook(getrhs, multJ, multJp, saveorbit, dotprod, &
                   mgmres, n, gtol, tol, del, mndl, mxdl, nits, info)

    print*, 'new_x = '
    print*, real(new_x)

    CALL new_xSetInv(new_x,Period,Var,Md,Nd)

    OPEN(2_8,FILE='Fin11.dat',FORM='UNFORMATTED')
    Write(2_8) VAR
    CLOSE(2_8)

    CALL steporbit(ndts,new_x, new_yF) 
    CALL new_xSetInvOrbit(new_yF,Var2,Md,Nd)
    OPEN(20_8,FILE='FinOrbit11.dat',FORM='UNFORMATTED')
    Write(20_8) VAR2
    CLOSE(20_8)

    OPEN(30_8,FILE='Peri11.dat',FORM='UNFORMATTED')
    Write(30_8) new_x(1)
    CLOSE(30_8)

                   ! check solution
    if(fixT) then
       d2 = sqrt((28d0-1d0)*(8d0/3d0)) 
       print*, 'Eqm = '
       print*, real(d2), real(d2), 28.-1.
    else    
       call steporbit(ndts,new_x, new_fx) ;
       print*, 'f(new_x) = '
       print*, real(new_fx(2:))
    end if

    print*,"result"

    norm_x = dsqrt(dotprod(-1,new_x,new_x))

    CALL new_xSetInv(new_x,Period,Var,Md,Nd)

    CALL getrhs(n,new_x,new_fx)

    CALL new_xSetInvOrbit(new_fx,Var2,Md,Nd)
    
    print*,Period,Lamda0

    do n_ = 1,n
      print*,n_,(new_fx(n_)**2)!/norm_x
    end do

    call saveorbit()

    stop
 
  contains
 
 !*************************************************************************
  END PROGRAM MAIN
 !*************************************************************************

  SUBROUTINE new_xSet(new_x,Period,Var,MM,Md,Nd)
   implicit none
   integer(8),intent(in)::MM,Nd,Md
   real(8),intent(in)::Period,Var(4*(MM+1)*(MM+1))
   real(8),intent(out)::new_x(4*(MM+1)*(MM+1) -5 + 1)
   integer(8)::LM,L

   LM = (MM+1)*(MM+1)
   CALL SXNM2L(MM,Nd,-Md,L)

   new_x(1) = Period
   new_x(2:LM) = Var(2:LM)
   new_x(LM+1:2*LM-1) = Var(LM+2:2*LM)
   new_x(2*LM:3*LM-2) = Var(2*LM+2:3*LM)
   new_x(3*LM-1:3*LM+L-4) = Var(3*LM+2:3*LM+L-1)
   new_x(3*LM+L-3:4*LM-4) = Var(3*LM+L+1:4*LM)

  END SUBROUTINE new_xSet

  SUBROUTINE new_xSetInv(new_x,Period,Var_,Md,Nd) 
   use orbit
   implicit none
   integer(8),intent(in)::Md,Nd
   real(8),intent(in)::new_x(4*(MM+1)*(MM+1) -5 + 1)
   real(8),intent(out)::Period,Var_(4*(MM+1)*(MM+1))
   integer(8)::LM,L_

   LM = (MM+1)*(MM+1)

   Period = new_x(1)
   Var_(1) = Var(1)
   Var_(2:LM) = new_x(2:LM)
   Var_(LM+1) = Var(LM+1)     
   Var_(LM+2:2*LM) = new_x(LM+1:2*LM-1)
   Var_(2*LM+1) = Var(2*LM+1)
   Var_(2*LM+2:3*LM) = new_x(2*LM:3*LM-2)
   CALL SXNM2L(MM,Nd,-Md,L_)
   Var_(3*LM+1) = Var(3*LM+1)
   Var_(3*LM+2:3*LM+L_-1) = new_x(3*LM-1:3*LM+L_-4)
   Var_(3*LM+L_) = Var(3*LM+L_)
   Var_(3*LM+L_+1:4*LM) = new_x(3*LM+L_-3:4*LM-4)

  END SUBROUTINE new_xSetInv

  SUBROUTINE new_xSetInvOrbit(new_x,Var_,Md,Nd) 
   use orbit
   implicit none
   integer(8),intent(in)::Md,Nd
   real(8),intent(in)::new_x(4*(MM+1)*(MM+1) -5 + 1)
   real(8),intent(out)::Var_(4*(MM+1)*(MM+1))
   real(8)::Period_TmnI
   integer(8)::LM,L_

   LM = (MM+1)*(MM+1)

   Period_TmnI = new_x(1)

   Var_(1) = Var(1)
   Var_(2:LM) = new_x(2:LM)
   Var_(LM+1) = Var(LM+1)     
   Var_(LM+2:2*LM) = new_x(LM+1:2*LM-1)
   Var_(2*LM+1) = Var(2*LM+1)
   Var_(2*LM+2:3*LM) = new_x(2*LM:3*LM-2)
   CALL SXNM2L(MM,Nd,-Md,L_)
   Var_(3*LM+1) = Var(3*LM+1)
   Var_(3*LM+2:3*LM+L_-1) = new_x(3*LM-1:3*LM+L_-4)
   Var_(3*LM+L_) = Period_TmnI
   Var_(3*LM+L_+1:4*LM) = new_x(3*LM+L_-3:4*LM-4)      

  END SUBROUTINE new_xSetInvOrbit

  !-------------------------------------------------------------------------
!  function to be minimised   
!-------------------------------------------------------------------------
subroutine getrhs(n_,x_, y_)
use orbit
use definition
implicit none
integer,          intent(in)  :: n_
double precision, intent(in)  :: x_(n)
double precision, intent(out) :: y_(n)
double precision :: y2_(n),Var_(4*LM),y3_(n)
integer(8)::I_,L_,Nd_,Md_
real(8)::OrbitShift


call steporbit(ndts,x_, y2_) !Period in x_,TmnI in y2_(1)

call new_xSetInvOrbit(y2_,Var_,Md,Nd)

!CALL SXNM2L(MM,Nd,Md,L_)

CALL SXNM2L(MM,Nd,Md,L_)

OrbitShift = atan2(VAR_(3*LM + L_+1),VAR_(3*LM + L_))/dble(Md) - PI!+ Adjust

DO I_= 1,4
   CALL SHIFT(VAR_((I_-1)*LM+1:I_*LM),OrbitShift)
END DO

!print*,VAR_(3*LM + L_)

OPEN(3_8,FILE='Orbit.dat',FORM='UNFORMATTED')
Write(3_8) VAR_
CLOSE(3_8)

!print*,"Orbit12",VAR_(3*LM + L_+1),VAR_(3*LM + L_)

!Nd_ = 4
!Md_ = 2
!CALL SXNM2L(MM,Nd_,Md_,L_)
!print*,"Orbit24",VAR(3*LM + L_+1),VAR(3*LM + L_)

!stop
!print*,VAR_(3*LM + L_+1)

call new_xSet(y3_,Period,Var_,MM,Md,Nd) !Period in y3_(1)

y_ = y3_ - x_      ! diff
y_(1) = 0d0			 ! constraints, rhs=0   

end subroutine getrhs

!-------------------------------------------------------------------------
!  Action of Jacobian and constraints on update.
!  Use approximation  dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
!-------------------------------------------------------------------------

subroutine multJ(n_,x, y)
use newton
use orbit,    only : n, epsJ, ndts, fixT
implicit none
integer,          intent(in)  :: n_
double precision, intent(in)  :: x(n_)
double precision, intent(out) :: y(n_)   
double precision, external :: dotprod
double precision :: eps, s(n_), s2(n_),dt,dlamda
             ! (F(x0+eps.x)-F(x0))/eps
!print*,x
eps = dsqrt(dotprod(1,x,x))
if(eps==0d0)  stop 'multJ: eps=0 (1)'
eps = epsJ * dsqrt(dotprod(1,new_x,new_x)) / eps  
if(eps==0d0)  stop 'multJ: eps=0 (2)' 
y = new_x + eps*x 
call getrhs(n_,y, s) !
y = (s - new_fx) / eps

if(fixT) then		! no extra constraint if T fixed  !fixT=.false.
   y(1) = 0d0
else				! contstraint, 
         ! no update in trajectory direction
   call steporbit(1,new_x, s)  !proceed 1 step
   dt = new_x(1)/ndts 
   s2 = (s - new_x) / dt !s is Derivative function
   y(1) = dotprod(-1,s2,x) !product !x is direction


!y(2) = 0d0
!----------
!   dlamda = new_x(2)/ndts
!   s = (s - new_x) / dlamda
!   y(2) = dotprod(-1,s,x)
!----------

end if

end subroutine multJ
!-------------------------------------------------------------------------
!  preconditioner for multJ.  Empty - no preconditioner required
!-------------------------------------------------------------------------
subroutine multJp(n, x)
implicit none
integer,          intent(in)    :: n
double precision, intent(inout) :: x(n)
end subroutine multJp

!-------------------------------------------------------------------------
!  called at each newton iteration   
!-------------------------------------------------------------------------
subroutine saveorbit()
use newton
use orbit
implicit none
double precision :: norm_x
double precision, external :: dotprod

print*, 'newton: iteration ', new_nits

norm_x = dsqrt(dotprod(-1,new_x,new_x))
print*, 'relative error = ', new_tol/norm_x
print*,new_x(1)!,new_x(2)

!  SAVE current solution, new_x

end subroutine saveorbit

!-------------------------------------------------------------------------
! dot product.  can flag to exclude parameter T.  Could include weights
!-------------------------------------------------------------------------
double precision function dotprod(n_,a,b)
use orbit
implicit none
integer,          intent(in) :: n_
double precision, intent(in) :: a(n), b(n)
integer :: n1
n1 = 1
if(n_==-1) n1 = 2
dotprod = dot_product(a(n1:n),b(n1:n))
end function dotprod

subroutine steporbit(ndts_,x_, y_)
use orbit
implicit none
integer,          intent(in)  :: ndts_
double precision, intent(in)  :: x_(n)
double precision, intent(out) :: y_(n)
double precision, save :: dt_

!print*,"Period=",x_(1)

if(ndts_/=1) then		! Set timestep size dt=T/ndts_
   dt_ = x_(1) / dble(ndts_)	! If only doing one step to calc \dot{x},
end if			! then use previously set dt.

!print*,x_(1),ndts_,dt_

CALL SAMAIN(ndts_,x_,y_,dt_)

end subroutine steporbit

SUBROUTINE SAMAIN(ndts_,x_,y_,dt_)
   
!     IU:The unit number of output file
!     CF:The name of output file
      use definition
      IMPLICIT NONE
      integer,          intent(in)  :: ndts_
      double precision, intent(in)  :: x_(n),dt_
      real(8)::Var_(4*LM)
      double precision, intent(out) :: y_(n)
      REAL(8) :: DTDNU
      integer(8)::I_,L_
      EXTERNAL TIME_DEVL,TIME_DERIVN
      
      CALL new_xSetInv(x_,Period,Var_,Md,Nd) 
      
      DTDNU=0.5D0*dt_/NDV*DNU
      DRN(1)=1
      DO L_=1,LM
         DRN(L_)=EXP(-DTDNU*(ABS(D(L_,1))-2)**LEV) ! 高階粘性のための設定
      END DO

      !OPEN(400_8,FILE='Orbit.dat',FORM='UNFORMATTED')
      !WRITE(400_8) VAR_

      DO II=1,ndts_
         CALL TXRKNU(LM*4,NDV,dt_,TIM,Var_,WORK,TIME_DEVL,TIME_DERIVN)
      END DO

      !WRITE(400_8) VAR_
      !CLOSE(400_8)
      !stop
      CALL SXNM2L(MM,Nd,-Md,L_)
      CALL new_xSet(y_,Var_(3*LM + L_),Var_,MM,Md,Nd) !Var -> y_

      
      END

!******************************************************************
!*     OPEN SUBROUTINE PACKAGE
!****************************************************************

!-----------------------------------------------------------------
!     CALCULATION OF d(VAR)/dt
!-----------------------------------------------------------------
      SUBROUTINE TIME_DERIVN(TIME,VAR_,DVAR_)
      use definition
      IMPLICIT NONE
      INTEGER(8)::L_,I_,J_
      REAL(8) :: TIME
      REAL(8) :: VAR_(LM*4), DVAR_(LM*4)

      DO L_ = 1,LM
          AVTT(L_) = VAR_(L_)
          AVTC(L_) = VAR_(LM+L_)
          DIV(L_) = VAR_(2*LM+L_)
          TMP(L_) = VAR_(3*LM+L_)
      END DO
      !print*,AVT(2)

      CALL SXCLAP(MM,MM,AVTT,STFT,D,2_8) !qnm → psinm
      CALL SXCLAP(MM,MM,AVTC,STFC,D,2_8)
      CALL SXCLAP(MM,MM,DIV,VPO,D,2_8) !Dnm → χnm

!     S1 Unm,Vnm
      CALL SXCS2X(MM,MM,VPO,VPOX)
      CALL SXCS2X(MM,MM,STFT,STFTX)
      CALL SXCS2X(MM,MM,STFC,STFCX)
      CALL SXCS2Y(MM,MM,STFT,STFTY,C) 
      CALL SXCS2Y(MM,MM,STFC,STFCY,C) 
      CALL SXCS2Y(MM,MM,VPO,VPOY,C) 
      CALL SXCRPK(MM,MM,MM+1,VPOX,VPOX2)
      CALL SXCRPK(MM,MM,MM+1,STFTX,STFTX2)
      CALL SXCRPK(MM,MM,MM+1,STFCX,STFCX2)

      DO L_=1,MM*MM+4*MM+2
            S1(L_) = -STFTY(L_) !Ut
            S2(L_) = VPOX2(L_) - STFCY(L_) !Uc
            S3(L_) = STFTX2(L_)  !Vt
            S4(L_) = STFCX2(L_) + VPOY(L_) !Vc
      END DO

      AVTT(2) = AVTT(2) + 2*OMEGA/SQRT(3D0)

      IPOW=0_8
      CALL SXTS2V(MM,NM,MM+1,IM,JM,S1,S2,G1,G2,IT,T,P,R,JC,W2,IPOW) !Ut,Uc
      CALL SXTS2V(MM,NM,MM+1,IM,JM,S3,S4,G3,G4,IT,T,P,R,JC,W2,IPOW) !Vt,Vc

      CALL SXTS2V(MM,NM,MM,IM,JM,AVTT,AVTC,G5,G6,IT,T,P,R,JC,W2,IPOW) !Xit,Zetac
      CALL SXTS2V(MM,NM,MM,IM,JM,DIV,TMP,G7,G8,IT,T,P,R,JC,W2,IPOW) !Deltac,T
      
!/ 非線形項の計算
      DO J_ = 1,JM
         DO I_ = 0,IM-1
            G9(I_,J_) = G1(I_,J_)*G5(I_,J_) + G2(I_,J_)*G6(I_,J_) + G7(I_,J_)*G4(I_,J_)
            G10(I_,J_) = G3(I_,J_)*G5(I_,J_) + G4(I_,J_)*G6(I_,J_) - G7(I_,J_)*G2(I_,J_)
            G11(I_,J_) = G2(I_,J_)*G5(I_,J_) + G1(I_,J_)*G6(I_,J_)
            G12(I_,J_) = G4(I_,J_)*G5(I_,J_) + G3(I_,J_)*G6(I_,J_)
            G13(I_,J_) = G4(I_,J_)*G5(I_,J_) + G3(I_,J_)*G6(I_,J_)  !!=G12
            G14(I_,J_) = G2(I_,J_)*G5(I_,J_) + G1(I_,J_)*G6(I_,J_)  !!=G11  
            !G15(I,J) = SQRT(3D0)/12D0*G8(I,J)  
            G152(I_,J_) = G1(I_,J_)*G2(I_,J_) + G3(I_,J_)*G4(I_,J_)
            G16(I_,J_) = SQRT(3D0)*G7(I_,J_)*S !+ KAP/2D0*G8(I,J))
         END DO
      END DO
      
      IPOW=2_8           
      CALL SXTV2S(MM,NM,MM,IM,JM,S9,S11,G9,G11,IT,T,P,R,JC,W2,IPOW) 
      CALL SXTV2S(MM,NM,MM,IM,JM,S13,S152,G13,G152,IT,T,P,R,JC,W2,IPOW)
    
      CALL SXCS2X(MM,MM,S9,SX9) !d(ut*xit ...)/dlam/cos
      CALL SXCS2X(MM,MM,S11,SX11) !
      CALL SXCS2X(MM,MM,S13,SX13)  

      IFLAG = 1_8
      CALL SXCLAP(MM,MM,S152,SL152,D,IFLAG) !(Lap_E)nm SL152
      CALL SXCLAP(MM,MM,TMP,SLTMP,D,IFLAG) !(Lap_TMP)nm SLTMP

!緯度微分に関わる項
      IPOW=2_8
      CALL SXTV2S(MM,NM,MM+1,IM,JM,S10,S12,G10,G12,IT,T,P,R,JC,W2,2_8)
      CALL SXTG2S(MM,NM,MM+1,IM,JM,S14,G14,IT,T,P,R,JC,W1,2_8)
      CALL SXCY2S(MM,MM,S10,SY10,C)
      CALL SXCY2S(MM,MM,S12,SY12,C)
      CALL SXCY2S(MM,MM,S14,SY14,C)
 
      !配列の大きさを揃える。
      !CALL SXCRPK(MM,MM+1,MM,SY10,SY102) 
      !CALL SXCRPK(MM,MM+1,MM,SY12,SY122) 
      !CALL SXCRPK(MM,MM+1,MM,SY14,SY142)
!T 
      CALL SXCS2X(MM,MM,TMP,TMPX) !dT/dlam

      CALL SXCS2Y(MM,MM,TMP,TMPY,C) !cos(phi)*dT/dphi

      IPOW = 0_8
      !CALL SXTS2V(MM,NM,MM,IM,JM,TMPX,TMPY,GTMPX,GTMPY,IT,T,P,R,JC,W2,IPOW)
      CALL SXTS2G(MM,NM,MM,IM,JM,TMPX,GTMPX,IT,T,P,R,JC,W1,0_8)
      CALL SXTS2G(MM,NM,MM+1,IM,JM,TMPY,GTMPY,IT,T,P,R,JC,W1,0_8)

      DO J_ = 1,JM
         DO I_ = 0,IM-1
            G17(I_,J_) = G1(I_,J_)*GTMPX(I_,J_) !ut*cos *dT/dlam
            G18(I_,J_) = G3(I_,J_)*GTMPY(I_,J_) !vt*cos *cos(phi)*dT/dphi
         END DO
      END DO

      CALL SXTG2S(MM,NM,MM,IM,JM,S16,G16,IT,T,P,R,JC,W1,0_8) !SQRT(3D0)*G7(I,J)*(S(I,J) + KAP/2D0*G8(I,J))
      CALL SXTG2S(MM,NM,MM,IM,JM,S17,G17,IT,T,P,R,JC,W1,2_8) !ut/cos *dT/dlam
      CALL SXTG2S(MM,NM,MM,IM,JM,S18,G18,IT,T,P,R,JC,W1,2_8) !vt*dT/dphi

      DO L_=1,LM
            DAVTT(L_) = -SX9(L_) - SY10(L_)
            DAVTC(L_) = -SX11(L_) - SY12(L_)
            DDIV(L_) = SX13(L_) - SY14(L_) - SQRT(3D0)/12D0*SLTMP(L_) - SL152(L_)
            DTMP(L_) = -S17(L_) - S18(L_) - S16(L_)
      END DO

      DO L_ = 1,LM
            DVAR_(L_) = DAVTT(L_) - kf*0.15_8*VAR_(L_)  !up
            DVAR_(LM+L_) = DAVTC(L_) - kf*0.15_8*AVTC(L_)
            DVAR_(2*LM+L_) = DDIV(L_) - kf*0.15_8*DIV(L_)
            DVAR_(3*LM+L_) = DTMP(L_) -ks*(TMP(L_))
            !print*, L,DAVTT(L),DAVTC(L),DDIV(L),DTMP(L)
      END DO

      !print*,DDIV(3),DIV(3),DDIV(5),DIV(5)

      DVAR_(3*LM+1) = DVAR_(3*LM+1) + ks*That_eq00
      DVAR_(3*LM+3) = DVAR_(3*LM+3) + ks*That_eq02  !up

      !stop
      END SUBROUTINE TIME_DERIVN

      SUBROUTINE TIME_DEVL(TIME_,DTIME_,VAR_)
            use definition
            IMPLICIT NONE        
            INTEGER(8) :: L_
            REAL(8) :: VAR_(4*LM),TIME_,DTIME_
        
            DO L_=1,LM
               VAR_(L_)=VAR_(L_)*DRN(L_) ! 高階粘性による減衰
               VAR_(LM+L_)=VAR_(LM+L_)*DRN(L_)
               VAR_(2*LM+L_)=VAR_(2*LM+L_)*DRN(L_)
               VAR_(3*LM+L_)=VAR_(3*LM+L_)*DRN(L_)
            END DO
         
      END SUBROUTINE TIME_DEVL
   
      SUBROUTINE ENEENS(I_,VAR_)
            use definition
            IMPLICIT NONE
            REAL(8) :: ENE,ET
            REAL(8) :: VAR_(4*LM)
            INTEGER(8) :: I_, L_
      
            ENE=0
            ET = 0
            DO L_=2,MM+1
               ENE=ENE-D(L_,2)*VAR_(L_)*VAR_(L_)-D(L_,2)*VAR_(LM+L_)*VAR_(LM+L_)-D(L_,2)*VAR_(2*LM+L_)*VAR_(2*LM+L_) ! ZetaT , ZetaC, Delta
               ET = ET + VAR_(3*LM+L_)*VAR_(3*LM+L_) !T
            END DO
            DO L_=MM+2,LM
               ENE=ENE-2*D(L_,2)*VAR_(L_)*VAR_(L_)-2*D(L_,2)*VAR_(LM+L_)*VAR_(LM+L_)-2*D(L_,2)*VAR_(2*LM+L_)*VAR_(2*LM+L_)
               ET = ET + 2*VAR_(3*LM+L_)*VAR_(3*LM+L_)
            END DO
            ENE=ENE/2 + ET/2D0/(12D0*S)
            WRITE(6,'(A,I3,A,F21.15)') 'step=',I_,' energy=',ENE
        
      END SUBROUTINE ENEENS

      SUBROUTINE SHIFT(VAR_,Lamda_)
            use definition
            IMPLICIT NONE
            REAL(8),dimension(LM)::VAR_,VAR2
            REAL(8)::Lamda_
            INTEGER(8)::M_,N_,L_

            VAR2 = VAR_
            DO M_ = 1,MM
               DO N_ = M_,MM
                  CALL SXNM2L(MM,N_,M_,L_)
                  VAR2(L_) = VAR_(L_)*COS(dble(M_)*Lamda_) + VAR_(L_+1)*SIN(dble(M_)*Lamda_)
                  VAR2(L_+1) = - VAR_(L_)*SIN(dble(M_)*Lamda_) + VAR_(L_+1)*COS(dble(M_)*Lamda_)
               END DO
            END DO

            VAR_ = VAR2

      END SUBROUTINE SHIFT



      !---------------------------------------------------------------------
! 非線形項と線形項を分離して Runge-Kutta法を使う
! (ISPACK1の TDRKNU の整数型変数のところを INTEGER(8)に変えたもの)
!---------------------------------------------------------------------
SUBROUTINE TXRKNU(N,M,H,X,Y,W,SUBL,SUBN)

      IMPLICIT NONE
      INTEGER(8) :: N,M,I
      REAL(8) :: H,X,Y(N),W(N,3),DX
      EXTERNAL SUBL,SUBN
    
      DX=H/M
      DO I=1,M
         CALL TXRKNL(N,DX,X,Y,W,SUBL,SUBN)
      END DO
    
    END SUBROUTINE TXRKNU
    !--------------------------------------        
    SUBROUTINE TXRKNL(N,DX,X,Y,W,SUBL,SUBN)
    
      IMPLICIT NONE
      INTEGER(8) :: N, I
      REAL(8) :: DX,X,Y(N),W(N,3)
      EXTERNAL SUBL,SUBN
    
      CALL SUBN(X,Y,W(1,1))
    
      DO I=1,N
         W(I,2)=Y(I)+(DX/2)*W(I,1)
         W(I,3)=Y(I)+(DX/6)*W(I,1)
      END DO
    
      CALL SUBL(X,DX/2,Y)
      CALL SUBL(X,DX/2,W(1,2))
      CALL SUBL(X,DX/2,W(1,3))
      CALL SUBN(X+DX/2,W(1,2),W(1,1))
    
      DO I=1,N
         W(I,2)=Y(I)+(DX/2)*W(I,1)
         W(I,3)=W(I,3)+(DX/3)*W(I,1)
      END DO
    
      CALL SUBN(X+DX/2,W(1,2),W(1,1))
    
      DO I=1,N
         W(I,2)=Y(I)+DX*W(I,1)
         W(I,3)=W(I,3)+(DX/3)*W(I,1)
      END DO
    
      CALL SUBL(X+DX/2,DX/2,W(1,2))
      CALL SUBL(X+DX/2,DX/2,W(1,3))
      CALL SUBN(X+DX,W(1,2),W(1,1))
    
      X=X+DX
      DO I=1,N
         Y(I)=W(I,3)+(DX/6)*W(I,1)
      END DO
    
    END SUBROUTINE TXRKNL



 