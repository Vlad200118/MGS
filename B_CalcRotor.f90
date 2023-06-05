Subroutine B_CalcRotorV (NI, NJ, V, CellVolume,IFaceVector,JFaceVector, IFaceCenter, &
                            JFaceCenter,CellCenter, RotV)
    Implicit none
    Integer :: NI, NJ ,I, J, IFace, I_neighbor, J_neighbor
    Real :: LinearInterpolation
    Real,Dimension(NI-1,NJ-1) :: CellVolume
    Real,Dimension(0:NI,0:NJ,2) :: V
    Real, Dimension (2) :: VLinearInterpolation
    Real,Dimension(0:NI,0:NJ) :: RotV
    Real,Dimension(NI,NJ-1,2) :: IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) :: JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) :: CellCenter
    Real,Dimension(2) :: FaceCenter, FaceVector

    DO J =1, NJ-1
        DO I=1, NI-1
            DO IFace=1, 4 !÷икл по сосед€м
                select case (IFace)
                case(1)!сосед снизу
                    I_neighbor = I
                    J_neighbor = J-1
                    FaceCenter(:) =  JFaceCenter(I,J,:)
                    FaceVector(:) = -JFaceVector(I,J,:)!минус из-за внешней нормали

                case(2)!сосед сверху
                    I_neighbor = I
                    J_neighbor = J+1
                    FaceCenter(:) = JFaceCenter(I,J+1,:)
                    FaceVector(:) = JFaceVector(I,J+1,:)

                case(3)!сосед слева
                    I_neighbor = I-1
                    J_neighbor = J
                    FaceCenter(:) =  IFaceCenter(I,J,:)
                    FaceVector(:) = -IFaceVector(I,J,:)!минус из-за внешней нормали

                case(4)!сосед справа
                    I_neighbor = I+1
                    J_neighbor = J
                    FaceCenter(:) = IFaceCenter(I+1,J,:)
                    FaceVector(:) = IFaceVector(I+1,J,:)

                end select

                    VLinearInterpolation(1) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),V(I,J,1), V(I_neighbor,J_neighbor,1))!V_linear_x

                    VLinearInterpolation(2) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),V(I,J,2), V(I_neighbor,J_neighbor,2))!V_linear_y


                    !! Ќужно обработать границу!!!
                    RotV(I,J) = RotV(I,J) + VLinearInterpolation(2)*FaceVector(1)-VLinearInterpolation(1)*FaceVector(2)

            END DO
            RotV(I,J)=RotV(I,J)/CellVolume(I,J)

        END DO
    END DO

End subroutine

