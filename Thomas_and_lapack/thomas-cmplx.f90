program thomas
implicit none

integer, parameter :: ikind=selected_int_kind(8)
integer (kind=ikind) :: i,ndim,nmid
complex*16 :: denom
complex*16, allocatable :: a(:),b(:),c(:),d(:),cp(:),dp(:),x(:),r(:)
integer*4 itime_start, itime_end, itime

read(5,*) ndim
allocate (a(ndim),b(ndim),c(ndim),d(ndim),cp(ndim),dp(ndim),x(ndim),r(ndim))

do i=1,ndim
  b(i) = cmplx(sqrt(dble(i)),1.0d0/dble(i))
  a(i) = cmplx(sin(dble(i)),-cos(dble(i)*0.5d0))
  c(i) = cmplx(cos(dble(i)),sin(dble(i)*1.2d0))
  d(i) = cmplx(log(dble(i)+1.0d0),1.0/sqrt(dble(i)))
end do

cp(1) = c(1)/b(1)
dp(1) = d(1)/b(1)

print *,'starting now'
call flush(6)

call system_clock(itime_start,itime)
floop: do i=2,ndim
denom = 1.0d0/(b(i)-a(i)*cp(i-1))
cp(i) = c(i)*denom
dp(i) = (d(i)-a(i)*dp(i-1))*denom
end do floop

x(ndim) = dp(ndim)
bloop: do i=ndim-1,1,-1
x(i) = dp(i)-cp(i)*x(i+1)
end do bloop
call system_clock(itime_end)
write(*,'(1x,a,i10,a,f11.3,a)') 'Time elapsed for ndim = ',ndim,&
      ': ',real(itime_end-itime_start)/real(itime), ' seconds.'
call flush(6)


printloop: do i=1,min(10,ndim)
write(*,1000) i,a(i),b(i),c(i),d(i),x(i)
1000 format(i10,1p10e16.8)
end do printloop

!check for n=5
if(ndim==5) then
 r(1) =           b(1)*x(1)+c(1)*x(2)-d(1)
 r(2) = a(2)*x(1)+b(2)*x(2)+c(2)*x(3)-d(2)
 r(3) = a(3)*x(2)+b(3)*x(3)+c(3)*x(4)-d(3)
 r(4) = a(4)*x(3)+b(4)*x(4)+c(4)*x(5)-d(4)
 r(5) = a(5)*x(4)+b(5)*x(5)          -d(5)
 print *,'check the residuals'
 do i=1,ndim
  write(*,1000) i,r(i)
 end do 
endif

if(ndim.gt.5) then
 nmid = ndim/2
 r(nmid) =        a(nmid)*x(nmid-1)+b(nmid)*x(nmid)+c(nmid)*x(nmid+1)-d(nmid)
 r(ndim) =        a(ndim)*x(ndim-1)+b(ndim)*x(ndim)                  -d(ndim)
 write(*,1000)    nmid,r(nmid)
 write(*,1000)    ndim,r(ndim)
endif

end program thomas
