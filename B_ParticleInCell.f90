Subroutine ParticleInCell(NI,NJ, IFaceVector,JFaceVector,IFaceCenter, JFaceCenter, coordParticle, IP, JP)
    implicit none
    integer:: I, J, NI, NJ, IP, JP,IFace
    real :: ScalarProduct
    logical :: IsInCell
    Real,Dimension(NI,NJ-1,2) :: IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) :: JFaceCenter, JFaceVector
    Real,Dimension(2) :: PB, FaceCenter, FaceVector
    Real,Dimension(2) :: coordParticle


    DO J = 1, NJ-1
        DO I = 1, NI-1
        IsInCell = .true.
            DO IFace=1, 4 !Цикл по соседям
                select case (IFace)
                case(1)!нижняя грань
                    FaceCenter(:) =  JFaceCenter(I,J,:)
                    FaceVector(:) = -JFaceVector(I,J,:)!минус из-за внешней нормали

                case(2)!верхняя грань
                    FaceCenter(:) = JFaceCenter(I,J+1,:)
                    FaceVector(:) = JFaceVector(I,J+1,:)

                case(3)!грань слева
                    FaceCenter(:) =  IFaceCenter(I,J,:)
                    FaceVector(:) = -IFaceVector(I,J,:)!минус из-за внешней нормали

                case(4)!грань справа
                    FaceCenter(:) = IFaceCenter(I+1,J,:)
                    FaceVector(:) = IFaceVector(I+1,J,:)
                end select
                PB(1) = FaceCenter(1) - coordParticle(1)
                PB(2) = FaceCenter(2) - coordParticle(2)
                ScalarProduct = DOT_PRODUCT(PB(:), FaceVector(:))

                IsInCell = IsInCell .and. ScalarProduct .GE. 0
            END DO
            IF (IsInCell) THEN
                IP = I
                JP = J
                EXIT
            END IF
        END DO
        IF (IsInCell) THEN
            EXIT
        END IF
    END DO

End Subroutine
