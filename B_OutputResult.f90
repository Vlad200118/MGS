Subroutine B_OutputResult(IO,NI,NJ,X,Y,P,V,Vmagn,GradU,GradV,RotV)
  integer :: NI, NJ
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P, Vmagn
  Real,Dimension(0:NI,0:NJ):: RotV
  Real,Dimension(0:NI,0:NJ,2):: V, GradU, GradV

  Write(IO,*) 'VARIABLES = "X", "Y","V_x","V_y","Vmagn","P","GradU_x","GradU_y","GradV_x","GradV_y","RotV"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ)
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)

  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') Vmagn(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)

  Write(IO,'(100F14.7)') GradU(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradU(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') GradV(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradV(1:NI-1,1:NJ-1,2)

  Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)

End Subroutine

