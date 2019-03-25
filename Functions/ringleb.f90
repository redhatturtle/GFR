pure function ringleb_solution(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref, aref, uref, machref
  use ovar, only : gam, rhoref, tref, pref, itestcase
  use ovar, only : ringleb_qmin, ringleb_kmax
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc, mb, me, ne
  real(wp) :: gm1, rhoa, rhoa2
  real(wp) :: Vmin, Vmax, q, k
  real(wp) :: a, J, rho, p, Eint, Ek, theta
  !
  !.. Local Arrays ..
  real(wp), dimension(1:3) :: vel
  !
continue
  !
  nc = 1
  mb = nc + 1
  me = nc + size(xyz)
  ne = me + 1
  !
  ! Some needed constants
  !
  gm1   = gam - one
  rhoa  = rhoref * aref
  rhoa2 = rhoa * aref
  !
  Vmin = ringleb_qmin
  Vmax = ringleb_kmax
  q    = v_from_xy_bisec( xyz(1), xyz(2), Vmin, Vmax )
  !
  a    = sqrt( one - (gm1/two)*q**two )
  rho  = a**(two/gm1)
  J    = (one/a) + ( one/(three*a**three) ) +( one/(five*a**five) ) - half*log( (one+a)/(one-a) )
  p    = a**(two*gam/gm1)/gam
  Eint = p/(gm1*rho)
  Ek   = half*q**two
  !
  k    = sqrt(   one/( half/q**two + rho*(xyz(1)+half*J) )   )
  !
  ! Compute the analytical value for velocity components
  !
  theta  = asin(q/k)
  vel(1) = q*cos(theta)
  vel(2) = q*sin(theta)
  vel(3) = zero
  if ( xyz(2) > zero ) vel(1) = -vel(1)
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(nc)    = rho
  return_value(mb:me) = rho * vel(:)
  return_value(ne)    = rho * (Eint + Ek)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function ringleb_solution
!
!###############################################################################
!
pure function V_from_xy_newt(x, y, Vini) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x, y, Vini
  !
  !.. Function Return Value ..
  real(wp)             :: return_value
  !
  !.. Local Scalars ..
  integer  :: i, iterMax
  real(wp) :: tol, gm1
  real(wp) :: v, a, rho, J, F, dF
  !
  iterMax = 100
  tol     = ten*epsilon(Vini)
  gm1 = gam - one
  !
  ! Initial guess
  !
  v = Vini
  !
  ! Start Newton iterations
  !
  do i=1, iterMax
    a   = sqrt( one - (gm1/two)*V**two )
    rho = a**(two/gm1)
    J   = (one/a) + ( one/(three*a**three) ) +( one/(five*a**five) ) - half*log( (one+a)/(one-a) )
    !
    F  = ( one/(two*rho*v**two) )**two - (x+half*J)**two - y**two
    dF = -one/( rho**two * v**five )
    !
    ! Iterate Newton method
    !
    v = v - ( F/dF )
    !
    if ( abs( F/dF ) < tol ) exit
  enddo
  !
  return_value = v
  !
end function V_from_xy_newt
!
!###############################################################################
!
pure function v_from_xy_bisec(x, y, Vmin, Vmax) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x, y, Vmin, Vmax
  !
  !.. Function Return Value ..
  real(wp)             :: return_value
  !
  !.. Local variables ..
  integer  :: iter, iterMax
  real(wp) :: tol, gm1
  real(wp) :: Vn, Vl, Vh, a, rho, J
  real(wp) :: yn, Ymin, Ymax
  !
  iterMax = 100
  tol  = ten*epsilon(Vmin)
  gm1 = gam - one
  !
  ! Begin calculating the function for V=Vmin and V=Vmax
  !
  Vl   = Vmin
  a    = sqrt( one - (gm1/two)*Vl**two )
  rho  = a**(two/gm1)
  J    = (one/a) + ( one/(three*a**three) ) +( one/(five*a**five) ) - log( (one+a)/(one-a) )/two
  Ymin = (x+half*J)**two + y**two - ( one/(two*rho*Vl**two) )**two
  !
  Vh   = Vmax
  a    = sqrt( one - (gm1/two)*Vh**two )
  rho  = a**(two/gm1)
  J    = (one/a) + ( one/(three*a**three) ) +( one/(five*a**five) ) - log( (one+a)/(one-a) )/two
  Ymax = (x+half*J)**two + y**two - ( one/(two*rho*Vh**two) )**two
  !
  ! Start bisection method iterations
  !
  do iter=1, iterMax
    Vn  = half*( Vl + Vh )
    a   = sqrt( one - (gm1/two)*Vn**two )
    rho = a**(two/gm1)
    J   = (one/a) + ( one/(three*a**three) ) +( one/(five*a**five) ) - log( (one+a)/(one-a) )/two
    !
    yn = (x+half*J)**two + y**two - ( one/(two*rho*Vn**two) )**two
    !
    if ( yMin*yn > 0 ) Vl = Vn
    if ( yMax*yn > 0 ) Vh = Vn
    !
    ! Check convergence
    !
    if ( abs( Vh-Vl ) < tol ) exit
  enddo
  !
  return_value = Vn
  !
end function V_from_xy_bisec
