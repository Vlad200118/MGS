Subroutine C_Boundary(X,Y,IFaceVector,JFaceVector,NI,NJ,r_m,V_m,w_m,r_mNext,V_mNext,w_mNext, diameter, BoundStatus)
    implicit none
    real, parameter :: recoveryFactor = 0.7
    integer :: I, J, NI, NJ, BoundStatus
    logical :: isCrossed
    Real :: w_m, w_mNext
    Real :: diameter
    Real,Dimension(NI,NJ):: X,Y
    Real,Dimension(NI,NJ-1,2) :: IFaceVector
    Real,Dimension(NI-1,NJ,2) :: JFaceVector
    Real,Dimension(2) :: r_m, V_m, r_mNext, r_mNext_new, V_mNext, rCross
    Real,Dimension(2) :: r_Bound, r_BoundNext

    BoundStatus = 0

    DO I = 1, NI - 1

        r_Bound(1) = X(I,1)
        r_Bound(2) = Y(I,1)
        r_BoundNext(1) = X(I+1,1)
        r_BoundNext(2) = Y(I+1,1)
        CALL Cross_Edges (r_m, r_mNext, r_Bound, r_BoundNext, rCross, isCrossed)
        IF(isCrossed) THEN
            BoundStatus = 1

            CALL Reflect (r_mNext, rCross, JFaceVector(I,1,:), r_mNext_new)
            WRITE(*,*) 'Particle cross bounary ', BoundStatus, ' in coordinates ', rCross(1), rCross(2)

            r_mNext(:) = r_mNext_new(:)

            CALL ReflectParameters (V_m, w_m, JFaceVector(I,1,:), diameter, V_mNext, w_mNext)

            return
        END IF

        r_Bound(1) = X(I,NJ)
        r_Bound(2) = Y(I,NJ)
        r_BoundNext(1) = X(I+1,NJ)
        r_BoundNext(2) = Y(I+1,NJ)
        CALL Cross_Edges (r_m, r_mNext, r_Bound, r_BoundNext, rCross, isCrossed)

        IF(isCrossed) THEN
            BoundStatus = 2
            CALL Reflect (r_mNext, rCross, JFaceVector(I,NJ,:), r_mNext_new)
            WRITE(*,*) 'Particle cross bounary ', BoundStatus, ' in coordinates ', rCross(1), rCross(2)

            r_mNext(:) = r_mNext_new(:)

            CALL ReflectParameters (V_m, w_m, JFaceVector(I,NJ,:), diameter, V_mNext, w_mNext)

            return
        END IF
    END DO

    DO J = 1, NJ - 1

        r_Bound(1) = X(1,J)
        r_Bound(2) = Y(1,J)
        r_BoundNext(1) = X(1,J+1)
        r_BoundNext(2) = Y(1,J+1)
        CALL Cross_Edges (r_m, r_mNext, r_Bound, r_BoundNext, rCross, isCrossed)
        IF(isCrossed) THEN
            BoundStatus = 3
            WRITE(*,*) 'Particle cross bounary ', BoundStatus, ' in coordinates ', rCross(1), rCross(2)
            return
        END IF

        r_Bound(1) = X(NI,J)
        r_Bound(2) = Y(NI,J)
        r_BoundNext(1) = X(NI,J+1)
        r_BoundNext(2) = Y(NI,J+1)
        CALL Cross_Edges (r_m, r_mNext, r_Bound, r_BoundNext, rCross, isCrossed)
        IF(isCrossed) THEN
            BoundStatus = 4
            WRITE(*,*) 'Particle cross bounary ', BoundStatus, ' in coordinates ', rCross(1), rCross(2)
            return
        END IF
    END DO
End Subroutine

Subroutine Cross_Edges (r_m, r_mNext, r_Bound, r_BoundNext, rCross, isCrossed)
implicit none
    logical :: isCrossed
    Real,Dimension(2) :: r_m, r_mNext, r_Bound, r_BoundNext, rCross
    Real :: det, detS, detT
    Real :: a11, a12, a21, a22
    Real :: S, T

    isCrossed = .false.
    a11 = r_mNext(1) - r_m(1)
    a12 = r_Bound(1) - r_BoundNext(1)
    a21 = r_mNext(2) - r_m(2)
    a22 = r_Bound(2) - r_BoundNext(2)

    det = a11*a22-a12*a21
    detS = a22 * (r_Bound(1) - r_m(1)) - a12*(r_Bound(2) - r_m(2))
    detT = a11*(r_Bound(2) - r_m(2)) - a21*(r_Bound(1) - r_m(1))

    IF(det .EQ. 0) THEN
        return
    ENDIF
    S = detS / det
    T = detT / det

    IF(S .GE. 0.0  .AND. S .LE. 1.0 .AND. T .GE. 0.0  .AND. T .LE. 1.0 ) THEN
        !WRITE(*,*)S, T
        rCross(1) = r_m(1) + S * a11
        rCross(2) = r_m(2) + S * a21
        isCrossed = .true.
        return
    END IF

End Subroutine

Subroutine ReflectParameters (V_m, w_m, FaceVector, diameter, V_mNext, w_mNext)
implicit none
    real, parameter :: recoveryFactor = 0.7
    REAL :: Magnitude2D, UctMagn
    REAL :: diameter, w_m, w_mNext
    Real,Dimension(2) :: V_m, V_mNext, FaceVector
    Real,Dimension(2) :: normal, tau, U, Uc, Uct

    !¬нутрен€€ нормаль
    normal(:) = FaceVector(:)/Magnitude2D(FaceVector(:))

    tau(1) =  normal(2)
    tau(2) = - normal(1)

    U(:) =  V_m(:)

    Uc(1) = U(1) - 0.5*diameter*w_m*normal(2)
    Uc(2) = U(2) + 0.5*diameter*w_m*normal(1)

    Uct(:) = Uc(:) - DOT_PRODUCT(Uc(:),normal(:))*normal(:)
    UctMagn = Magnitude2D(Uct(:))


    V_mNext(:) = V_m(:) - (1+recoveryFactor)*DOT_PRODUCT(U(:),normal(:))*normal(:) - 2.0*UctMagn*tau(:)/7.0
    w_mNext = w_m - 10.0 * UctMagn / 7.0 / diameter

End Subroutine

Subroutine Reflect (r_mNext, r_bound, FaceVector, r_mNext_new)
implicit none
    real :: Magnitude2D
    Real,Dimension(2) :: r_mNext, r_bound, FaceVector, r_mNext_new
    Real,Dimension(2) :: reflectVector, reflectVectorTN, difVector, tau, normal
!!Ќужна внутренн€€ нормаль
    difVector(:) = r_mNext(:) - r_bound(:)
    normal(:) = FaceVector(:) / Magnitude2D(FaceVector(:))
    tau(1) = normal(2)
    tau(2) = -normal(1)
    reflectVectorTN(1) = DOT_PRODUCT(difVector(:), tau(:))
    reflectVectorTN(2) = -DOT_PRODUCT(difVector(:), normal(:))
    reflectVector(1) = reflectVectorTN(1)*tau(1) + reflectVectorTN(2)*normal(1)
    reflectVector(2) = reflectVectorTN(1)*tau(2) + reflectVectorTN(2)*normal(2)

    r_mNext_new(:) = r_bound(:) + reflectVector(:)

End Subroutine
