Subroutine B_CalcGradientGG(NI,NJ,P,CellVolume,IFaceVector,JFaceVector,IFaceCenter, &
                            JFaceCenter,CellCenter, GradP)
    Implicit none
    Integer:: NI, NJ ,I, J, IFace, I_neighbor, J_neighbor
    Real:: PLinearInterpolation, LinearInterpolation
    Real,Dimension(NI-1,NJ-1):: CellVolume
    Real,Dimension(0:NI,0:NJ)::P
    Real,Dimension(0:NI,0:NJ,2)::GradP
    Real,Dimension(NI,NJ-1,2) ::IFaceCenter, IFaceVector
    Real,Dimension(NI-1,NJ,2) ::JFaceCenter, JFaceVector
    Real,Dimension(0:NI,0:NJ,2) ::CellCenter
    Real,Dimension(2)::FaceCenter, FaceVector

    GradP(:,:,:) = 0.0
    DO J =1, NJ-1
        DO I=1, NI-1
            DO IFace=1, 4 !���� �� �������
                select case (IFace)
                case(1)!����� �����
                    I_neighbor = I
                    J_neighbor = J-1
                    FaceCenter(:) =  JFaceCenter(I,J,:)
                    FaceVector(:) = -JFaceVector(I,J,:)!����� ��-�� ������� �������

                case(2)!����� ������
                    I_neighbor = I
                    J_neighbor = J+1
                    FaceCenter(:) = JFaceCenter(I,J+1,:)
                    FaceVector(:) = JFaceVector(I,J+1,:)

                case(3)!����� �����
                    I_neighbor = I-1
                    J_neighbor = J
                    FaceCenter(:) =  IFaceCenter(I,J,:)
                    FaceVector(:) = -IFaceVector(I,J,:)!����� ��-�� ������� �������


                case(4)!����� ������
                    I_neighbor = I+1
                    J_neighbor = J
                    FaceCenter(:) = IFaceCenter(I+1,J,:)
                    FaceVector(:) = IFaceVector(I+1,J,:)

                end select

                    PLinearInterpolation = LinearInterpolation(FaceCenter(:), CellCenter(I,J,:), &
                                   CellCenter(I_neighbor,J_neighbor,:),P(I,J), P(I_neighbor,J_neighbor))!PE

                    GradP(I,J,:) = GradP(I,J,:) + PLinearInterpolation*FaceVector(:)

            END DO
            GradP(I,J,:)=GradP(I,J,:)/CellVolume(I,J)

        END DO
    END DO
End Subroutine
