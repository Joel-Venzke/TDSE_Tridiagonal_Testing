      program trilapackcmplx
      implicit double precision (a-h,o-z)
      complex*16, allocatable :: a(:),b(:),c(:),d(:,:)
      integer*4 itime_start, itime_end, itime
    
      read(5,*) ndim
      allocate (a(ndim),b(ndim),c(ndim),d(ndim,1))

      do i=1,ndim
       b(i) = cmplx(sqrt(dble(i)),1.0d0/dble(i))
       if (i.gt.1) a(i-1) = cmplx(sin(dble(i)),-cos(dble(i)*0.5d0))
       c(i) = cmplx(cos(dble(i)),sin(dble(i)*1.2d0))
       d(i,1) = cmplx(log(dble(i)+1.0d0),1.0/sqrt(dble(i)))
      end do
    
      call system_clock(itime_start,itime)
      call zgtsv(ndim,1,a,b,c,d,ndim,ifail)
      call system_clock(itime_end)
      write(*,'(1x,a,i10,a,f11.3,a)') 'Time elapsed for ndim = ',ndim,
     > ': ',real(itime_end-itime_start)/real(itime), ' seconds.'
      call flush(6)
      do i=1,min(10,ndim)
       write(6,1000) i,d(i,1)
      end do
      do i=max(1,ndim-10),ndim
       write(6,1000) i,d(i,1)
      end do
1000  format(i10,1p10e16.8)
      stop
      end
