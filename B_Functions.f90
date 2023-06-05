Function LinearInterpolation(TargetPoint, PointIn, PointOut, funIn, funOut)
    implicit none
    real, dimension(2) :: TargetPoint, PointIn, PointOut
    real::LenIn, LenOut, funIn, funOut, LinearInterpolation

    LenIn  = Norm2(PointIn(:) - TargetPoint(:))!¬нутри €чейки
    LenOut = Norm2(PointOut(:) - TargetPoint(:))!—наружи €чейки €чейки

    LinearInterpolation = (funIn*LenOut + funOut*LenIn)/(LenOut + LenIn)
End Function

pure real function Magnitude2D(Vector)
  real,Dimension(2), intent(in) :: Vector
  Magnitude2D = sqrt(Vector(1)*Vector(1) + Vector(2)*Vector(2))
end function



