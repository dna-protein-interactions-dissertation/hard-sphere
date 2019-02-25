program hard_potential

use randgen
implicit none

integer, parameter :: long = selected_real_kind(15, 50)
double precision, parameter :: PI=3.141592653589793238462

! Defining our parameters

real(kind=long), parameter :: D = 25e-12
real(kind=long), parameter :: b = 1000
real(kind=long), parameter :: u = 1000
real(kind=long), parameter :: r = 10e-9
real(kind=long), parameter :: L = 400e-9
real(kind=long), parameter :: lambda = (3*sqrt((u+1)*b))/(8*r**2*sqrt(2*PI*D))
real(kind=long), parameter :: sigma = sqrt((u+1)*D*b)
real(kind=long), parameter :: dt = 1e-9
real(kind=long), parameter :: p_new = (3*b*(u+1)*dt*L**2)/(16*PI*r**2)

! Defining the variables we need to use

integer :: nsolv, nloc, length
real(kind=long), dimension(3) :: seed, dp, dv, Xcol, Xpar, Xper, Ycol, Ypar, Yper, cax
real(kind=long), allocatable, dimension(:, :) :: Y, nY
real(kind=long), dimension(6) :: X, nX, temp
real(kind=long) :: square_distance, t_tot = 0, ac, bc, cc, tau, mcax, trunc_normal, rand_cerf, xi, ike,ake, &
xmom, ymom, zmom

! Defining the dummy variables we use for indexing

integer :: ii, jj, kk, mm, nadd=0, nlost=0, qq, t
character(len=1024) :: filename

! Initializing our system

X = [real(kind=long) :: L/2, L/2, L/2, 0, 0, 0]

nsolv = random_poisson(real(nint(L**3*lambda)), .TRUE.)
length = nint(nsolv + 100*sqrt(real(nsolv)))
allocate(Y(length, 7))
allocate(nY(length, 7))

print*, sigma

do ii=1, nsolv

	call random_number(seed)
	Y(ii, 1) = 1
	Y(ii, 2:4) = L*seed
	Y(ii, 5) = sigma*random_normal()
	Y(ii, 6) = sigma*random_normal()
	Y(ii, 7) = sigma*random_normal()
	
	if (square_distance(X(1:3), Y(ii, 2:4)) < r**2) then
		Y(ii, 1) = 0
	end if
	
end do

do jj=nsolv+1, length
	Y(jj, :) = [0,0,0,0,0,0,0]
end do


! Looping until either a predetermined time or some logical condition is met

do while(t_tot < 0.001)

nX(1:3) = X(1:3) + X(4:6)*dt
nx(4:6) = X(4:6)

nY(:, 1) = Y(:, 1)
nY(:, 2:4) = Y(:, 2:4) + Y(:, 5:7)*dt
nY(:, 5:7) = Y(:, 5:7)

do kk=1, length
	if (Y(kk, 1)>0 .AND. square_distance(nX(1:3), nY(kk, 2:4)) < r**2) then
	
		! a collision has occured, find out when
		ike = Y(kk, 5)**2 + Y(kk, 6)**2 + Y(kk, 7)**2 + 1000*(X(4)**2 + X(5)**2 + X(6)**2)
		xmom = Y(kk, 5) + 1000*X(4)
		ymom = Y(kk, 6) + 1000*X(5)
		zmom = Y(kk, 7) + 1000*X(6)
		!print*, xmom, ymom, zmom
		dp = Y(kk, 2:4) - X(1:3) ! difference in prior positions
		dv = Y(kk, 5:7) - X(4:6) ! difference in prior velocities
		ac = dot_product(dv, dv)
		bc = 2*dot_product(dv, dp)
		cc = dot_product(dp, dp) - r**2
		tau = (-bc - (bc**2 - 4*ac*cc)**0.5)/(2*ac)
		print*, tau
		print*, nadd-nlost
		! now calculate quantities needed for update rule
		Xcol = X(1:3) + tau*X(4:6)
		Ycol = Y(kk, 2:4) + tau*Y(kk, 5:7)
		cax = Ycol - Xcol
		mcax = sqrt(dot_product(cax, cax))
		cax = cax/mcax
		Xper = (dot_product(cax, X(4:6)))*cax
		Xpar = X(4:6) - Xper
		Yper = (dot_product(cax, Y(kk, 5:7)))*cax
		Ypar = Y(kk, 5:7) - Yper
		
		! now update nY, nX according to formulas in notes
		
		nY(kk, 5:7) = Ypar + ((1-u)/(1+u))*Yper + ((2*u)/(1+u))*Xper
		nX(4:6) = Xpar + ((u-1)/(1+u))*Xper + (2/(1+u))*Yper
		
		nY(kk, 2:4) = Ycol + nY(kk, 5:7)*(dt-tau)
		nX(1:3) = Xcol + nX(4:6)*(dt-tau)	
		ake = nY(kk, 5)**2 + nY(kk, 6)**2 + nY(kk, 7)**2 + 1000*(nX(4)**2 + &
		nX(5)**2 + nX(6)**2)
		xmom = nY(kk, 5) + 1000*nX(4)
		ymom = nY(kk, 6) + 1000*nX(5)
		zmom = nY(kk, 7) + 1000*nX(6)
		!print*, xmom, ymom, zmom
	end if
end do

X = nX
Y = nY

! Now we just need to apply boundary conditions

! first we remove particles which have left the domain
!nlost = 0
do mm=1, length
	
	if (maxval(Y(mm, 2:4)) > L .OR. minval(Y(mm, 2:4)) < 0) then
	
		if (Y(mm, 1)==1) then
			nlost=nlost+1
			!print*, mm
		end if
		Y(mm, 1) = 0
	end if

end do

! now we introduce particles into the domain

call random_number(temp)
!nadd=0

if (temp(1) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 1
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, dt*sigma*xi*2**0.5, L*seed(1), L*seed(2), sigma*trunc_normal&
		(sigma*xi*2**0.5), sigma*random_normal(), sigma*random_normal()]
	!print*, int(sum(Y(:, 1)))

end if

if (temp(2) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 2
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, L*seed(1), dt*sigma*xi*2**0.5, L*seed(2), sigma*random_normal(),&
		sigma*trunc_normal(sigma*xi*2**0.5), sigma*random_normal()]
	!print*, int(sum(Y(:, 1)))
end if

if (temp(3) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 3
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, L*seed(1), L*seed(2), dt*sigma*xi*2**0.5, sigma*random_normal(),&
	sigma*random_normal(), sigma*trunc_normal(sigma*xi*2**0.5)]
	!print*, int(sum(Y(:, 1)))
end if

if (temp(4) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 4
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, L - dt*sigma*xi*2**0.5, L*seed(1), L*seed(2), -sigma*trunc_normal&
		(sigma*xi*2**0.5), sigma*random_normal(), sigma*random_normal()]
	!print*, int(sum(Y(:, 1)))
end if

if (temp(5) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 5
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, L*seed(1), L-dt*sigma*xi*2**0.5, L*seed(2), sigma*random_normal(),&
		-sigma*trunc_normal(sigma*xi*2**0.5), sigma*random_normal()]
	!print*, int(sum(Y(:, 1)))
end if

if (temp(6) < p_new) then
	nadd=nadd+1
	xi = rand_cerf()
	call random_number(seed)
	nloc = minloc(Y(:, 1), DIM=1)
	!print*, nloc, 6
	!print*, int(sum(Y(:, 1)))
	Y(nloc, :) = [real(kind=long):: 1, L*seed(1), L*seed(2), L-dt*sigma*xi*2**0.5, sigma*random_normal(),&
	sigma*random_normal(), -sigma*trunc_normal(sigma*xi*2**0.5)]
	!print*, int(sum(Y(:, 1)))
end if

t_tot = t_tot + dt

!print*, int(sum(Y(:, 1)))
!print*, nsolv -nlost + nadd - int(sum(Y(:, 1)))

! if (mod(nint(t_tot/dt), 100)==0) then
! write(filename,"(A14,I0.3)") "solvent",nint(t_tot/dt)
! open(unit=1, file=filename, status='replace')
! do qq=1, 10000
	! write(1, *) (Y(qq, t), ",", t=1, 7)
! end do
! end if

end do

end program

function square_distance(X, Y)

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	real(kind=long), dimension(3), intent(in) :: X, Y
	real(kind=long) :: square_distance
	
	square_distance = dot_product(Y-X, Y-X)
	
end function

function rand_cerf()

	implicit none
	integer, parameter :: long = selected_real_kind(15, 50)
	real(kind=long) :: rand_cerf, temp(2), rr
	real(kind=long), parameter:: a1=0.532
	real(kind=long), parameter:: a2=0.814
	
20  call random_number(temp)
    rr = -a1*log(temp(1))
    if(temp(2)*temp(1) < a2*ERFC(rr)) then
        rand_cerf=rr
    else
        go to 20
    end if
	
end function

function trunc_normal(A)
    implicit none
    integer, parameter :: long = selected_real_kind(15, 50)
    real(kind=long), intent(in) :: A
    real(kind=long) :: z1, z2, z3, cc, aa, trunc_normal
    double precision, dimension(2) :: temp
    aa=0.5*(A+(A**2 + 4)**0.5)
10  call random_number(temp)
    z1=temp(1)
    z2=A-log(z1)/a
    z3=temp(2)
    if (z3<exp(-0.5*(z2-a)**2)) then
        cc=z2
    else
        go to 10
    end if
    trunc_normal=cc
end function
