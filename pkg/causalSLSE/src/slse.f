      subroutine mypnorm(x, n, mu, sig, p)

      integer i, n
      double precision x(n), mu, sig, p(n)

      do i=1,n
         call fpnorm(p(i), x(i), mu, sig, 0, 0)
      end do

      end

      subroutine myhat(x, n, k, qraux, tol, hat)
      
      integer n, k, info, i, j
      double precision x(n,k), qraux(k), hat(n), sigma(n)
      double precision dummy, tol
      
      do i = 1,n
        hat(i) = 0.0d0
      end do
      
      do j = 1,k
         do i = 1,n
            sigma(i) = 0.0d0
         end do
         sigma(j) = 1.0d0
         call dqrsl(x, n, n, k, qraux, sigma, sigma, dummy,
     c        dummy, dummy, dummy, 10000, info)
         do i = 1, n
            hat(i) = hat(i)+sigma(i)*sigma(i)
         end do
      end do
      do i = 1, n
        if(hat(i) .ge. 1.0d0 - tol) hat(i) = 1.0d0
      end do
      end  

      subroutine myvcov(x, xqr, qraux, e, tol, n, k, rk, type, v)

      integer n, k, rk, i, j, type, info
      double precision x(n, k), xqr(n,k), v(rk, rk)
      double precision sig, e(n), qraux(k), tol, e2(n), bread(rk, rk)
      double precision xxi(rk*(rk+1)/2), hat(n), meat(rk, rk)

      do i=1,rk
         do j=i,rk
            xxi(i+(j-1)*j/2) = xqr(i,j)
         end do
      end do         
      call dpptri('u', rk, xxi, info)
      if (type .gt. 1) then
         call myhat(xqr, n, k, qraux, tol, hat)      
      end if
      sig = 1.0d0
      
      if (type .eq. -1) then
         sig = sum(e*e)/real(n-rk)
      end if

      do i=1,rk
         do j=i,rk
            bread(i,j) = xxi(i+(j-1)*j/2)*sig
            bread(j,i) = bread(i,j)
         end do
      end do

      if (type .eq. -1) then
         v = bread
      else
         e2 = abs(e)
         if (type .eq. 2) then
            e2 = e2/sqrt((1-hat))
         end if
         if (type .eq. 3) then
            e2 = e2/(1-hat)
         end if
         do i=1,k
            x(:,i) = x(:,i)*e2
         end do
         do i=1,rk
            do j=i,rk
               meat(i,j) = sum(x(:,j)*x(:,i))
               meat(j,i) = meat(i,j)
            end do
         end do
         v = matmul(bread, meat)
         v = matmul(v, bread)
         if (type .eq. 1) then
            v = v*n/(n-k)
         end if
      end if
      end

      subroutine myls(y, x, n, k, tol, type, rk, pv, e, b, vcov)
      integer n, k, rk, pv(k), info, i, j, type
      double precision y(n), x(n,k), b(k,1), qraux(k), sig
      double precision tol,  work(2*k), e(n,1), qty(n,1), xxi(k*(k+1)/2)
      double precision hat(n), b2(n, k), sigma(n, k) 
      double precision vcov(k,k), xqr(n,k), dd

      se = 0.0d0
      xqr = x
      
      call dqrls(xqr,n,k,y,1,tol,b,e,qty,rk,pv,qraux,work)
      call myvcov(x, xqr, qraux, e, tol, n, k, rk, type, vcov)

      end   

      subroutine splinei(x, n, k, nk, mnk, sx)

      integer n, nk, i
      double precision x(n), k(mnk), sx(n, nk+1)
      sx = 0.0d0
      if (nk .gt. 0) then
         sx = 0.0d0
         do i=1,(nk+1)
            if (i .eq. 1) then
               where (x .le. k(1))
                  sx(:,i) = x
               elsewhere    
                  sx(:,i) = k(i)
               end where
            else if (i .eq. (nk+1)) then
               where (x .gt. k(nk))
                  sx(:,i) = x-k(nk)
               end where
            else
               where ((x .ge. k(i-1)) .and. (x .le. k(i)))
                  sx(:,i) = x-k(i-1)
               end where
               where (x .gt. k(i))
                  sx(:, i) = k(i)-k(i-1)
               end where
            end if
         end do
      else
         sx(:,1) = x
      end if
      end

      subroutine spline(x, n, p, k, nk, mnk, tnk, sx)

      integer n, p, mnk, tnk, nk(p), i, ind1, ind2
      double precision x(n,p), k(mnk, p), sx(n, tnk+p)

      ind2 = 0
      do i=1,p
         ind1 = ind2+1
         ind2 = ind1+nk(i)
         call splinei(x(:,i), n, k(:,i), nk(i), mnk, sx(:,ind1:ind2))
      end do      
      end

      subroutine splinefit(y, x, n, p, k, nk, mnk, tnk, tol,
     *                     rk, b, pv, rsd, eff)
      integer n, p, mnk, tnk, nk(p), pv(tnk+p+1), rk
      double precision x(n,p), k(mnk, p), sx(n, tnk+p+1), eff(n), tol
      double precision b(tnk+p+1,1), rsd(n), work(2*(tnk+p+1))
      double precision y(n), qr(tnk+p+1) 

      sx(:,1) = 1.0d0
      call spline(x, n, p, k, nk, mnk, tnk, sx(:,2:(tnk+p+1)))
      call dqrls(sx,n,tnk+p+1,y,1,tol,b,rsd,eff,rk,pv,qr,work)      
      end
 

      subroutine modelfit(y0, y1, x0, x1, n0, n1, p, tol, 
     *                    k0, nk0, mnk0, tnk0, rk0,
     *                    k1, nk1, mnk1, tnk1, rk1,
     *                    w0, w1, b0, pv0, b1, pv1, 
     *                    e0, e1, eff0, eff1, bic, aic)
      integer n0, n1, p, mnk0, tnk0, nk0(p), pv0(tnk0+p+1), rk0
      integer mnk1, tnk1, nk1(p), pv1(tnk1+p+1), rk1
      integer w0(mnk0,p), w1(mnk1,p), i, tnk0s, tnk1s
      integer nk0s(p), nk1s(p)
      double precision x0(n0,p), x1(n1,p), y0(n0), y1(n1)
      double precision tol, e0(n0), e1(n1), eff0(n0), eff1(n1)
      double precision k0(mnk0, p), b0(tnk0+p+1,1), ssr
      double precision k1(mnk1, p), b1(tnk1+p+1,1), aic, bic, ll, pi
      double precision k1s(mnk1,p), k0s(mnk0,p)

      k1s = k1
      k0s=k0
      nk0s=nk0
      nk1s=nk1
      tnk0s=tnk0
      tnk1s=tnk1

      do i=1,p
         if (w0(1,i) .eq. 0) then
            nk0s(i) = 0
         else
            nk0s(i) = count(w0(:,i) .gt. 0)
            k0s(1:nk0s(i),i) = k0s(w0(1:nk0s(i),i),i)
         end if
         if (w1(1,i) .eq. 0) then
            nk1s(i) = 0
         else
            nk1s(i) = count(w1(:,i) .gt. 0)
            k1s(1:nk1s(i),i) = k1s(w1(1:nk1s(i),i),i)
         end if
      end do
      tnk0s = sum(nk0s)
      tnk1s = sum(nk1s)
      
      call splinefit(y0, x0, n0, p, k0s, nk0s, mnk0, tnk0s,
     *     tol, rk0, b0, pv0, e0, eff0)
      call splinefit(y1, x1, n1, p, k1s, nk1s, mnk1, tnk1s,
     *     tol, rk1, b1, pv1, e1, eff1)

      pi = 4.d0*datan(1.d0)
      ssr = sum(e1**2) + sum(e0**2)
      ll = real(n1+n0)*(log(2.0d0 * pi) + 1 - log(real(n1+n0)) +
     *     log(ssr))
      aic = ll+2.0d0*real(rk0+rk1+1)
      bic = ll+log(real(n1+n0))*real(rk0+rk1+1)     
      end
 

      subroutine selic(y0, y1, x0, x1, n0, n1, p, tol, 
     *                 k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1, pval0,
     *     pval1, spval, npval, bic, aic,
     *     w0bic, w1bic, w0aic, w1aic)
      integer n0, n1, p, mnk0, tnk0, nk0(p), pv0(tnk0+p+1)
      integer mnk1, tnk1, nk1(p), pv1(tnk1+p+1), npval
      integer w0(mnk0, p), w1(mnk1, p), j, i, l, s
      integer w0bic(mnk0, p), w1bic(mnk1, p)
      integer w0aic(mnk0, p), w1aic(mnk1, p)
      double precision x0(n0,p), x1(n1,p), y0(n0), y1(n1), tol
      double precision k0(mnk0, p), spval(npval)
      double precision pval0(mnk0, p), pval1(mnk1, p)
      double precision k1(mnk1, p), aic(npval+1), bic(npval+1)
      double precision aicsel, bicsel

      w0(:,:) = 0
      w1(:,:) = 0
      call modelfit0(y0, y1, x0, x1, n0, n1, p, tol, k0, nk0, mnk0,
     *     tnk0, k1, nk1, mnk1, tnk1, w0, w1, bic(1), aic(1))
      bicsel = bic(1)
      aicsel = aic(1)
      w0bic = w0
      w1bic = w1
      w0aic = w0
      w1aic = w1
      do i=1,npval
         do j=1,p
            if (nk0(j) .gt. 0) then
               l = 1
               do s=1,nk0(j)
                  if (pval0(s,j) .le. spval(i)) then
                     w0(l,j) = s
                     l = l+1
                  end if
               end do
            end if
            if (nk1(j) .gt. 0) then
               l = 1
               do s=1,nk1(j)
                  if (pval1(s,j) .le. spval(i)) then
                     w1(l,j) = s
                     l = l+1
                  end if
               end do
            end if            
         end do
         call modelfit0(y0, y1, x0, x1, n0, n1, p, tol, k0, nk0, mnk0,
     *        tnk0, k1, nk1, mnk1, tnk1, w0, w1, bic(i+1), aic(i+1))
         if (aic(i+1) .le. aicsel) then
            aicsel = aic(i+1)
            w0aic = w0
            w1aic  = w1
         end if
         if (bic(i+1) .le. bicsel) then
            bicsel = bic(i+1)
            w0bic = w0
            w1bic  = w1
         end if         
      end do
      end

      subroutine modelfit0(y0, y1, x0, x1, n0, n1, p, tol, 
     *     k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1, 
     *     w0, w1, bic, aic)
      integer n0, n1, p, mnk0, tnk0, nk0(p), pv0(tnk0+p+1)
      integer mnk1, tnk1, nk1(p), pv1(tnk1+p+1)
      integer w0(mnk0,p), w1(mnk1,p), i, tnk0s, tnk1s
      integer nk0s(p), nk1s(p), rk0, rk1
      double precision x0(n0,p), x1(n1,p), y0(n0), y1(n1)
      double precision tol, e0(n0), e1(n1), eff0(n0), eff1(n1)
      double precision k0(mnk0, p), b0(tnk0+p+1,1), ssr
      double precision k1(mnk1, p), b1(tnk1+p+1,1), aic, bic, ll, pi
      double precision y0f(n0), x0f(n0,p), y1f(n1), x1f(n1,p)
      double precision k1s(mnk1,p), k0s(mnk0,p)
      
      k1s = k1
      k0s=k0
      nk0s=nk0
      nk1s=nk1
      tnk0s=tnk0
      tnk1s=tnk1

      do i=1,p
         if (w0(1,i) .eq. 0) then
            nk0s(i) = 0
         else
            nk0s(i) = count(w0(:,i) .gt. 0)
            k0s(1:nk0s(i),i) = k0s(w0(1:nk0s(i),i),i)
         end if
         if (w1(1,i) .eq. 0) then
            nk1s(i) = 0
         else
            nk1s(i) = count(w1(:,i) .gt. 0)
            k1s(1:nk1s(i),i) = k1s(w1(1:nk1s(i),i),i)
         end if
      end do
      tnk0s = sum(nk0s)
      tnk1s = sum(nk1s)
      x0f=x0
      x1f=x1
      y0f=y0
      y1f=y1
      do i=1,(tnk0+p+1)
         pv0(i) = i
      end do
      do i=1,(tnk1+p+1)
         pv1(i) = i
      end do
      
      call splinefit(y0f, x0f, n0, p, k0s, nk0s, mnk0, tnk0s,
     *     tol, rk0, b0, pv0, e0, eff0)
      call splinefit(y1f, x1f, n1, p, k1s, nk1s, mnk1, tnk1s,
     *     tol, rk1, b1, pv1, e1, eff1)

      pi = 4.d0*datan(1.d0)
      ssr = sum(e1**2) + sum(e0**2)
      ll = real(n1+n0)*(log(2.0d0 * pi) + 1 - log(real(n1+n0)) +
     *     log(ssr))
      aic = ll+2.0d0*real(rk0+rk1+1)
      bic = ll+log(real(n1+n0))*real(rk0+rk1+1)     
      end
 

