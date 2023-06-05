Program Main
  Implicit None
  real, parameter :: PI = 4.0*atan(1.0)
  character(*), parameter:: InputFile='input.txt', InputParticleFile ='inputParticle.txt' ! names of input and output files
  character(*), parameter:: OutputParticleFile='OutputParticle.dat', TimeForce = 'TimeForce.txt'
  character MeshFile*30, SolutionFile*30, ResultFile*30, str       ! name of file with computational mesh
  integer, parameter:: IO = 12, I1 = 13 ! input-output unit
  integer :: GGI, InterpolateOrder, I, J
  Real :: rtmp
  real :: Magnitude2D

  !MESH VAR
  integer ::  NI, NJ
  real,allocatable,dimension(:,:):: X,Y,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays
  real :: LRef

  !PARTICLE VAR
  integer :: BoundStatus
  integer :: IP, JP, IPNext, JPNext
  integer :: ITime, NTime
  real :: diameter, densityFlow, densityParticle, viscosity, StokesNumber, Square
  real :: deltaTime
  real, Dimension(2) :: r_m, r_mNext, V_m, V_mNext, tempV
  real, Dimension(2) :: Force, FDrag, FMag
  real :: w_m, w_mNext, Moment
  real :: ReV, ReW
  real :: addMass, mass, inertia
  real :: f, DragCoef, MomentCoef,  magnV, magnW

  !FLOW VAR
  real, Dimension(2) :: Vflow, CellCenterParticle
  real,allocatable,dimension(:,:,:)::V
  real,allocatable,dimension(:,:)::P, T, Vmagn
  real,allocatable,dimension(:,:) :: RotV
  real,allocatable,dimension(:,:,:)::GradU, GradV
  real :: w
  real :: ReFlow

!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) SolutionFile  ! read name of file with computational mesh
  READ(IO,*) ResultFile
  READ(IO,*) GGI
  CLOSE(IO)

  !===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputParticleFile
  OPEN(IO,FILE=InputParticleFile)
  READ(IO,*) Lref
  READ(IO,*) ReFlow
  READ(IO,*) densityFlow, densityParticle  ! read name of file with computational mesh
  READ(IO,*) viscosity  ! read name of file with computational mesh
  READ(IO,*) diameter
  READ(IO,*) r_m(1), r_m(2)
  READ(IO,*) V_m(1), V_m(2)
  READ(IO,*) w_m
  READ(IO,*) deltaTime
  READ(IO,*) NTime
  READ(IO,*) InterpolateOrder
  CLOSE(IO)

  diameter = diameter * Lref

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ
  WRITE(*,*) 'GGI = ', GGI

  WRITE(*,*)"LRef = ", LRef, " ReFlow = ", ReFlow
  WRITE(*,*)"densityFlow = ", densityFlow,"densityParticle = ", densityParticle
  WRITE(*,*)"viscosity = ", viscosity, "diameter = ",  diameter
  WRITE(*,*)"Vx_0 = ", V_m(1), "Vy_0 = ", V_m(2)

  StokesNumber = ReFLow*densityParticle/densityFlow*(diameter / LRef)**2.0 / 18.0

  WRITE(*,*) 'StokesNumber = ', StokesNumber

!=== ALLOCATE ALL ARRAYS ===
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(T(0:NI,0:NJ))
  allocate(V(0:NI,0:NJ,2), Vmagn(0:NI,0:NJ))   ! Velocity
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces

  allocate(RotV( 0:NI,0:NJ))
  allocate(GradU( 0:NI,0:NJ,2), GradV( 0:NI,0:NJ,2))

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
    READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  CLOSE(IO)

!=== READ SOLUTION ===
  WRITE(*,*) 'Read solution'
  OPEN (IO, FILE = SolutionFile)
  READ(IO,*) str
  READ(IO,*) str
  READ(IO,*) ((rtmp,rtmp, V(I,J,1),V(I,J,2),rtmp,P(I,J),rtmp,rtmp, I=0, NI), J=0, NJ)
  CLOSE(IO)


  DO J=1,NJ
    DO I=1, NI
        X(I,J) = X(I,J) * LRef
        Y(I,J) = Y(I,J) * Lref
        V(I,J,:) = V(I,J,:) * ReFlow * viscosity/LRef/densityFlow
    END DO
  END DO

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
  Vmagn(:,:) = sqrt(V(:,:,1)**2.0 + V(:,:,2)**2.0)

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate gradient U and V'

  If (GGI .EQ. 0) THEN
    Call B_CalcGradientGG(NI,NJ,V(0:NI,0:NJ,1),CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradU)
    Call B_CalcGradientGG(NI,NJ,V(0:NI,0:NJ,2),CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradV)
  ELSE
    Call B_CalcGradientGGI(NI,NJ,V(0:NI,0:NJ,1),CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradU)
    Call B_CalcGradientGGI(NI,NJ,V(0:NI,0:NJ,2),CellVolume,IFaceVector,JFaceVector,IFaceCenter,JFaceCenter,CellCenter, GradV)
  END IF

!=== CALCULATE ROTOR ===
  WRITE(*,*) 'Calculate Rotor V'
    Call B_CalcRotorV(NI,NJ,V,CellVolume,IFaceVector,JFaceVector, IFaceCenter,JFaceCenter,CellCenter,RotV)

  WRITE(*,*) 'Output fields to file: ', ResultFile
  Open(IO,FILE=ResultFile)
  Call B_OutputResult(IO,NI,NJ,X,Y,P,V,Vmagn,GradU,GradV,RotV)
  Close(IO)


Square = 0.25*PI*(diameter)**2.0
mass = PI*(diameter)**3.0*densityParticle/6.0
addMass = 0.5*PI*(diameter)**3.0*densityFlow/6.0
inertia = PI*diameter**5.0*densityParticle/60.0

OPEN(IO,FILE =  OutputParticleFile)
OPEN(I1,FILE =  TimeForce)
    WRITE(IO,*)'ITime, x, y, u, v, w, Fx, Fy, Mz'
  DO ITime = 1, NTime
    IP = -1
    JP = -1
    IPNext = -1
    JPNext = -1

    CALL ParticleInCell(NI,NJ, IFaceVector,JFaceVector,IFaceCenter, JFaceCenter, r_m, IP, JP)
    !WRITE(*,*) "ITime = ", ITime, "IP = ", IP, "JP = ", JP

    IF(InterpolateOrder .EQ. 1) THEN
        Vflow(:) = V(IP,JP,:)
    ELSE
        CellCenterParticle(:) = CellCenter(IP,JP,:) - r_m(:)
        Vflow(1) = V(IP,JP,1)  + DOT_PRODUCT(CellCenterParticle(:), GradU(IP,JP,:))
        Vflow(2) = V(IP,JP,2)  + DOT_PRODUCT(CellCenterParticle(:), GradV(IP,JP,:))
    ENDIF

    magnV = Magnitude2D(Vflow(:)-V_m(:))
    ReV = diameter * densityFlow*magnV/viscosity
    f =1.0+0.179*sqrt(ReV)+0.013*ReV
    DragCoef = 24.0*f/ReV
    FDrag(:) = 0.5*DragCoef*densityFlow*Square*magnV*(Vflow(:)-V_m(:))

    w =0.5*RotV(IP,JP)
    tempV(1) = (V_m(2)-Vflow(2))
    tempV(2) = (V_m(1)-Vflow(1))
    FMag(:) = PI*diameter**3.0*densityFlow*(w_m-w)*tempV(:)/8.0

    Force(:) = FDrag(:) + FMag(:)

    magnW=abs(w_m-w)
    ReW=diameter**2.0*magnW*densityFlow/viscosity
    MomentCoef = 64.0/ReW
    Moment = MomentCoef*PI*diameter**5.0*magnW*(w-w_m)/64.0


    r_mNext(:) = r_m(:) + deltaTime*V_m(:)
    V_mNext(:) = V_m(:) + deltaTime*Force(:)/(mass+addMass)
    w_mNext = w_m + deltaTime * Moment/inertia

    WRITE(IO,*)ITime, r_m(:), V_m(:), w_m, Force(:), Moment
    WRITE(I1,*)ITime*deltaTime, Magnitude2D(FDrag(:))/mass, Magnitude2D(FMag(:))/mass, Magnitude2D(addMass*(V_mNext(:)-V_m(:))/deltaTime)/mass

    CALL ParticleInCell(NI,NJ, IFaceVector,JFaceVector,IFaceCenter, JFaceCenter, r_mNext, IPNext, JPNext)

    IF (IPNext .EQ. -1 .AND.  JPNext .EQ. -1)THEN
        CALL C_Boundary(X,Y,IFaceVector,JFaceVector,NI,NJ,r_m,V_m,w_m,r_mNext,V_mNext,w_mNext, diameter, BoundStatus)
        IF(BoundStatus .EQ. 3 .OR. BoundStatus .EQ. 4) THEN
            EXIT
        END IF
    END IF
    r_m(:) = r_mNext(:)
    V_m(:) = V_mNext(:)
    w_m = w_mNext
  END DO
CLOSE(IO)
CLOSE(I1)


  Deallocate (X,Y,CellVolume,IFaceVector,JFaceVector, IFaceCenter,JFaceCenter,CellCenter)
  Deallocate (P, V, T)
  Deallocate (GradU, GradV)


END PROGRAM Main
