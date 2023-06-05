Subroutine B_CalcGradientGGI(NI,NJ,P,CellVolume,IFaceVector,JFaceVector,IFaceCenter, &
                            JFaceCenter,CellCenter, GradP)
    Implicit none
    Integer:: NI, NJ ,I, J, IFace, I_neighbor, J_neighbor, ITER, NITER
    Real:: PLinearInterpolation, LinearInterpolation
    Real,Dimension(NI-1,NJ-1):: CellVolume
    Real,Dimension(0:NI,0:NJ)::P
    Real,Dimension(0:NI,0:NJ,2)::GradP, GradPLast
    Real,Dimension(NI,NJ-1,2) ::IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) ::JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) ::CellCenter
    Real,Dimension(2)::FaceCenter, FaceVector, EMVector, PointE, GradLinearInterpolationE

    NITER = 10
    GradP(:,:,:) = 0.0
    DO ITER = 1, NITER
    GradPLast(:,:,:) = GradP(:,:,:)
    GradP(:,:,:) = 0.0
    DO J =1, NJ-1
        DO I=1, NI-1
            DO IFace=1, 4 !Цикл по соседям
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
                !!Находим координаты точки Е
                    PointE(1) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                           CellCenter(I_neighbor,J_neighbor,:),CellCenter(I,J,1), CellCenter(I_neighbor,J_neighbor,1))!!!RE
                    PointE(2) = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                           CellCenter(I_neighbor,J_neighbor,:),CellCenter(I,J,2), CellCenter(I_neighbor,J_neighbor,2))!!!RE
                !!Находим вектор EM
                    EMVector(:) = FaceCenter(:) - PointE(:)

                    PLinearInterpolation = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),P(I,J), P(I_neighbor,J_neighbor))!PE

                    GradLinearInterpolationE(1) = LinearInterpolation(EMVector(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:), GradPLast(I,J,1), GradPLast(I_neighbor,J_neighbor,1))!GradPX в точке E

                    GradLinearInterpolationE(2) = LinearInterpolation(EMVector(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:), GradPLast(I,J,2), GradPLast(I_neighbor,J_neighbor,2))!GradPY в точке E

                    GradP(I,J,:) = GradP(I,J,:) + (PLinearInterpolation + &
                                    DOT_PRODUCT(GradLinearInterpolationE(:), EMVector(:)))*FaceVector(:)

            END DO
            GradP(I,J,:)=GradP(I,J,:)/CellVolume(I,J)

        END DO
    END DO
    END DO
End Subroutine 
