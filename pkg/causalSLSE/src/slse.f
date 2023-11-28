c     This generates the hat values h_t of the LSE.
c     x is the k first column of the QR decomposition
c     returned by dqrls, where k is the rank of the 
c     regression fit.
      
      subroutine hatvals(x, n, k, qraux, tol, hat)
      
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
     *        dummy, dummy, dummy, 10000, info)
         do i = 1, n
            hat(i) = hat(i)+sigma(i)*sigma(i)
         end do
      end do
      do i = 1, n
        if(hat(i) .ge. 1.0d0 - tol) hat(i) = 1.0d0
      end do
      end  

c     LSE with vcov
c     vcov(1:rk, 1:rk) is the covariance matrix of the feasible estimator
c     The type is classical when tv = -1 and HCi for tv=i, i=0,1,2,3.
c     If tv .eq. -10, no v is computed
      
      subroutine lse(y, x, n, k, tol, tv, rk, piv, e, b, v)
      integer n, k, rk, piv(k), info, i, j, tv
      double precision y(n), x(n,k), b(k,1), qraux(k)
      double precision tol,  work(2*k), e(n,1), qty(n,1)
      double precision v(k,k), xqr(n,k)

      se = 0.0d0
      xqr = x
      do i=1,k
         piv(i) = i
      end do

      call dqrls(xqr,n,k,y,1,tol,b,e,qty,rk,piv,qraux,work)
      if (tv .ne. -10) then
         call vcovhc(x, xqr, qraux, e, tol, n, k, rk, tv, piv,
     *        v(1:rk, 1:rk))
      end if
      end   

c     The function computes the LSE covariance matrix using HCCM
c     It requires the residuals e, the matrix of covariates x
c     and the following output from dqrls:
c     - xqr: the QR decomposition
c     - piv: the pivot indicator of 1:k
c     - rk: the rank of the regression.
c     - qraux: info from the QR decomposition
c     v(1:rk, 1:rk) is the covariance matrix of the feasible estimates b(piv(1:rk))
c     the bread (X'X)^{-1} is computed by inverting the upper triangular matrix R
c     using the Lapack subroutine dpptri.
      
      subroutine vcovhc(x, xqr, qraux, e, tol, n, k, rk, tv, piv, v)
      integer n, k, rk, i, j, tv, info, piv(k)
      double precision x(n, k), xqr(n,k), v(rk, rk)
      double precision sig, e(n), qraux(k), tol, e2(n), bread(rk, rk)
      double precision xxi(rk*(rk+1)/2), hat(n), meat(rk, rk)

      do i=1,rk
         do j=i,rk
            xxi(i+(j-1)*j/2) = xqr(i,j)
         end do
      end do         
      call dpptri('u', rk, xxi, info)
      if (tv .gt. 1) then
         call hatvals(xqr(:,1:rk), n, rk, qraux, tol, hat)      
      end if
      sig = 1.0d0
      
      if (tv .eq. -1) then
         sig = sum(e*e)/dble(n-rk)
      end if

      do i=1,rk
         do j=i,rk
            bread(i,j) = xxi(i+(j-1)*j/2)*sig
            bread(j,i) = bread(i,j)
         end do
      end do

      if (tv .eq. -1) then
         v = bread
      else
         e2 = abs(e)
         if (tv .eq. 2) then
            e2 = e2/sqrt((1-hat))
         end if
         if (tv .eq. 3) then
            e2 = e2/(1-hat)
         end if
         do i=1,rk
            x(:,piv(i)) = x(:,piv(i))*e2
         end do
         do i=1,rk
            do j=i,rk
               meat(i,j) = sum(x(:,piv(j))*x(:,piv(i)))
               meat(j,i) = meat(i,j)
            end do
         end do
         v = matmul(bread, meat)
         v = matmul(v, bread)
         if (tv .eq. 1) then
            v = v*n/(n-k)
         end if
      end if
      end

c     This function generate the spline functions for one x
c     It is called by the spline subroutine.
      
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

c     This function generate the basis functions from the knots
c     It is called separately for the treated and nontreated
c     On return, sx is the n x (tnk+p) matrix of basis function,
c     where tnk is the total number of knots and p is the number of covariates
      
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

c     This function generate the basis function from the knots
c     and estimate the model by least squares. The matrix x does not
c     contain an intercept but it is added to the spline matrix.
c     The function returns the AIC, BIC and vcov of the LS coefficients
c     using the classical vcov if vt=-1 and HCi (i=0,1,2,3) for vt=i.
c     If vt .eq. -10, no v is computed and it isa set to 0.      
c     It is called separately for the treated and nontreated
      
      subroutine splinefit(y, x, n, p, k, nk, mnk, tnk, tol,
     *                     rk, b, piv, rsd, vt, v)
      integer n, p, mnk, tnk, nk(p), piv(tnk+p+1), rk, vt
      double precision x(n,p), k(mnk, p), sx(n, tnk+p+1), tol
      double precision b(tnk+p+1,1), rsd(n)
      double precision y(n), v(tnk+p+1,tnk+p+1)

      sx(:,1) = 1.0d0
     
      call spline(x, n, p, k, nk, mnk, tnk, sx(:,2:(tnk+p+1)))
      call lse(y, sx, n, tnk+p+1, tol, vt, rk, piv, rsd, b, v)

      end
     
c     This subroutine updates the list of knots using a selection
c     The selection index is w, which is mnk x p
c     where mnk is the maximum number of knots and p is the number of cavariates.      
c     For each covariate:
c     If the number of knots was originally 0, nothing change.
c     If w(1,i) >= mnk, nothing change.
c     If w(1,i) == 0, all knots are removed.
c     It is assumed that the knots selection index w(:,i) is sorted and at the beginning.  
      
      subroutine updatek(k, p, nk, mnk, w, kf, nkf)
      
      integer p, nk(p), mnk, w(mnk, p), i, nkf(p)
      double precision k(mnk, p), kf(mnk, p)

      nkf = nk
      kf = k
      
      do i=1,p
         if ((nk(i) .gt. 0)) then
            if (w(1,i) .le. mnk) then
               where (w(:,i) .gt. 0)
                  kf(:,i) = kf(w(:,i),i)
               end where
               nkf(i) = count(w(:,i).gt.0)
            end if
         end if
      end do
      end

c     This function estimate the coefficient of a model
c     It creates the splines for each group and estimate the model for each group
c     as well. The estimation is done separately for each group.
c     It returns the LSE, the AIC and BIC, and the covariance matrix of the two
c     vectors of LSE estimates based on the type vt. The latter is classical one
c     if vt=-1 and HCi (i=0,1,2,3) for vt=i.
c     If vt .eq. -10, no v is computed and it isa set to 0.      
      
      subroutine modelfit(y0, y1, x0, x1, p, n0, n1, tol, 
     *                    k0, nk0, mnk0, tnk0, rk0, piv0,
     *                    k1, nk1, mnk1, tnk1, rk1, piv1,
     *                    w0, w1, vt, tnk0s, tnk1s, nk0s, nk1s,
     *                    b0, b1, v0, v1, bic, aic)
      integer p, n0, nk0(p), mnk0, tnk0, n1, nk1(p), mnk1, tnk1
      integer w0(mnk0,p), w1(mnk1,p)
      integer piv0(tnk0+p+1), rk0, piv1(tnk1+p+1), rk1
      integer nk0s(p), nk1s(p), tnk0s, tnk1s, i, vt
      double precision k1s(mnk1,p), k0s(mnk0,p), k1(mnk1,p), k0(mnk0,p)
      double precision y0(n0), y1(n1), x0(n0,p), x1(n1,p)
      double precision b0(tnk0+p+1,1), b1(tnk1+p+1,1)
      double precision tol, aic, bic, ll, pi, v1(tnk1+p+1,tnk1+p+1)
      double precision e0(n0), e1(n1), v0(tnk0+p+1,tnk0+p+1)

      call updatek(k0, p, nk0, mnk0, w0, k0s, nk0s)
      call updatek(k1, p, nk1, mnk1, w1, k1s, nk1s)
      tnk0s = sum(nk0s)
      tnk1s = sum(nk1s)
      
      call splinefit(y0, x0, n0, p, k0s, nk0s, mnk0, tnk0s,
     *     tol, rk0, b0(1:(tnk0s+p+1),1), piv0(1:(tnk0s+p+1)),
     *     e0, vt, v0(1:(tnk0s+p+1),1:(tnk0s+p+1)))
      call splinefit(y1, x1, n1, p, k1s, nk1s, mnk1, tnk1s,
     *     tol, rk1, b1(1:(tnk1s+p+1),1), piv1(1:(tnk1s+p+1)),
     *     e1, vt, v1(1:(tnk1s+p+1),1:(tnk1s+p+1)))
     
      pi = 4.d0*datan(1.d0)
      ssr = sum(e1**2) + sum(e0**2)
      ll = dble(n1+n0)*(log(2.0d0 * pi) + 1 - log(dble(n1+n0)) +
     *     log(ssr))
      aic = ll+2.0d0*dble(rk0+rk1+1)
      bic = ll+log(dble(n1+n0))*dble(rk0+rk1+1)
      end

c     Return the pvalues of the test H0: 'slopes adjacent to a knot is the same'
c     for all knots of a given covariate (pi=1,...,p)
c     All p-values are set to 0 when the number of knots is equal to 0, and it is set to
c     2 when it cannot be computed (for example when b is missing). This number means that
c     the p-value if missing.
      
      subroutine testknoti(b, v, n, nk, mnk, tnk, p, rk, piv, pi, pval)

      integer mnk, tnk, p, n, nk(p), rk, piv(tnk+p+1), pi, i, j1, j2
      double precision b(tnk+p+1), v(tnk+p+1, tnk+p+1)
      double precision pval(mnk), b1, b2, vi, test, sortb(tnk+p+1)
      double precision numna, sortv(tnk+p+1, tnk+p+1), pvali

      pval = 0.0d0     
      numna = maxval(b)+1
      sortb = numna+1
      sortb(piv(1:rk)) = b(1:rk)
      sortv(piv(1:rk), piv(1:rk)) = v(1:rk, 1:rk)
          
      if (nk(pi) .gt. 0) then
         if (pi .eq. 1) then
            j1 = 2
         else
            j1 = pi+1+sum(nk(1:(pi-1)))
         end if
         do i=1,nk(pi)
            j2 = j1+1
            if ((sortb(j1).gt.numna) .or. (sortb(j2).gt.numna)) then
               pval(i) = 2.0d0
            else
               b1 = b(j1)
               b2 = b(j2)
               vi = sortv(j1,j1)+sortv(j2,j2)-2*sortv(j1,j2)
               test = (b2-b1)*(b2-b1)/vi
               call fpf(pvali, test, 1.0d0, dble(n-rk))
               pval(i) = 1.0d0-pvali
            end if
            j1 = j1+1
         end do 
      end if
      end

c     The following 2 functions compute the p-values using the Backward and Forward methods.
c     The p-values are stored in an mnk x p array. The number of valid pvalues per column
c     is determined by nk, which represents the number of knots per covariate.
c     A p-value greater
c     than 1 is considered as an NA. 
      
      subroutine pvalb(y, x, k, tol, n, p, nk, mnk, tnk, vt, pval)

      integer n, p, nk(p), mnk, tnk, rk, vt, piv(tnk+p+1), i
      double precision y(n), x(n,p), pval(mnk, p), tol, e(n) 
      double precision k(mnk, p), b(tnk+p+1,1), v(tnk+p+1, tnk+p+1)
     
      call splinefit(y, x, n, p, k, nk, mnk, tnk, tol, rk, b(:,1),
     *     piv, e, vt, v) 
      
      do i=1,p
         call testknoti(b(:,1), v, n, nk, mnk, tnk, p, rk, piv, i,
     *        pval(:,i))
      end do
      
      end

      subroutine pvalf(y, x, k, tol, n, p, nk, mnk, tnk, vt, pval)

      integer n, p, nk(p), mnk, tnk, rk, vt, piv(tnk+p+1), i, w(mnk,p),
     *     j, nks(p), tnks
      double precision y(n), x(n,p), pval(mnk, p), tol, e(n)
      double precision k(mnk, p), b(tnk+p+1,1), v(tnk+p+1, tnk+p+1)
      double precision ks(mnk, p), pval0(mnk)

      w = 0      
      do i=1,p
         if (nk(i) .gt. 0) then
            if (nk(i) .le. 2) then
               w(1,i) = 1
               if (nk(i) .eq. 2) then
                  w(2,i) = 2
               end if
               call updatek(k, p, nk, mnk, w, ks, nks)
               tnks = sum(nks)
               call splinefit(y, x, n, p, ks, nks, mnk, tnks, tol,
     *              rk, b(1:(tnks+p+1),1), piv(1:(tnks+p+1)), e, vt,
     *              v(1:(tnks+p+1),1:(tnks+p+1)))
               call testknoti(b(1:(tnks+p+1),1),
     *              v(1:(tnks+p+1),1:(tnks+p+1)), n, nks, mnk, tnks, p,
     *              rk, piv(1:(tnks+p+1)), i, pval(:,i))
            else
               do j=1,nk(i)
                  w(:,i) = 0
                  if (j .eq. 1) then
                     w(1,i) = 1
                     w(2,i) = 2
                  else if (j .eq. nk(i)) then
                     w(1,i) = nk(i)-1
                     w(2,i) = nk(i)
                  else
                     w(1,i) = j-1
                     w(2,i) = j
                     w(3,i) = j+1
                  end if
                  call updatek(k, p, nk, mnk, w, ks, nks)
                  tnks = sum(nks)
                  call splinefit(y, x, n, p, ks, nks, mnk, tnks, tol,
     *                 rk, b(1:(tnks+p+1),1), piv(1:(tnks+p+1)), e, vt,
     *                 v(1:(tnks+p+1),1:(tnks+p+1)))
                  call testknoti(b(1:(tnks+p+1),1),
     *                 v(1:(tnks+p+1),1:(tnks+p+1)), n, nks, mnk, tnks,
     *                 p, rk, piv(1:(tnks+p+1)), i, pval0)
                  if (j .eq. 1) then
                     pval(j,i) = pval0(1)
                  else 
                     pval(j,i) = pval0(2)
                  end if
               end do
            end if
         end if
         w(:,i) = 0
      end do
      end

c     This is the selection method function for causalSLSE models.
c     To avoid having to recompute the p-values for
c     different selection method, it is possible to have more than one selection methods
c     - pvm is the p-value method: 1 for Backward and 2 for Forward.
c     - selm is the selection method. 1 = PVT only and 2 = PVT, AIC and BIC
c     - t1 and t0 are the PV threshold.
c     mnk0 and mnk1 cannot be 0. 
      

      subroutine selcmodel(y0,y1,x0,x1,n0,n1,p,tol,t0,t1,pvm,vt,selm, 
     *     k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1, pval0, 
     *     pval1, bic, aic, w0bic, w0aic, w0pvt, w1bic, w1aic, w1pvt,
     *     npval)

      integer n0, n1, p, mnk0, tnk0, nk0(p), mnk1, tnk1, nk1(p), pvm, vt
      integer selm, w0pvt(mnk0, p), w1pvt(mnk1, p), npval
      integer w0bic(mnk0, p), w1bic(mnk1, p)
      integer w0aic(mnk0, p), w1aic(mnk1, p)
      double precision x0(n0,p), x1(n1,p), y0(n0), y1(n1), tol, t0, t1
      double precision k0(mnk0, p), pval0(mnk0, p), pval1(mnk1, p)
      double precision k1(mnk1, p), aic(tnk0+tnk1+1), bic(tnk0+tnk1+1)

      if (pvm .eq. 1) then
         call pvalb(y0, x0, k0, tol, n0, p, nk0, mnk0, tnk0, vt, pval0)
         call pvalb(y1, x1, k1, tol, n1, p, nk1, mnk1, tnk1, vt, pval1)
      else 
         call pvalf(y0, x0, k0, tol, n0, p, nk0, mnk0, tnk0, vt, pval0)
         call pvalf(y1, x1, k1, tol, n1, p, nk1, mnk1, tnk1, vt, pval1)
      end if

      call selpvt(p, nk0, mnk0, t0, pval0, w0pvt)
      call selpvt(p, nk1, mnk1, t1, pval1, w1pvt)

      if (selm .eq. 2) then
         call selicc(y0, y1, x0, x1, n0, n1, p, tol, 
     *        k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1, pval0,
     *        pval1, bic, aic, w0bic, w1bic, w0aic, w1aic, npval)
      end if
      end


c     This is the selection method function for slse models.
c     To avoid having to recompute the p-values for
c     different selection method, it is possible to have more than one selection methods
c     - pvm is the p-value method: 1 for Backward and 2 for Forward.
c     - selm is the selection method. 1 = PVT only and 2 = PVT, AIC and BIC
c     - t1 and t0 are the PV threshold.
c     Note that mnk and tnk cannot be 0. If they are 0, we set them to 1.
      

      subroutine selmodel(y, x, n, p, tol, t, pvm, vt, selm, k, nk, mnk,
     *     tnk, pval, bic, aic, wbic, waic, wpvt, npval)

      integer n, p, mnk, tnk, nk(p), pvm, vt, selm, wpvt(mnk, p),  npval
      integer wbic(mnk, p), waic(mnk, p)
      double precision x(n,p), y(n), tol, t, k(mnk, p), pval(mnk, p)
      double precision aic(tnk+1), bic(tnk+1), spval(tnk)

      if (pvm .eq. 1) then
         call pvalb(y, x, k, tol, n, p, nk, mnk, tnk, vt, pval)
      else 
         call pvalf(y, x, k, tol, n, p, nk, mnk, tnk, vt, pval)
      end if
            
      call selpvt(p, nk, mnk, t, pval, wpvt)

      if (selm .eq. 2) then
         call vecpval(pval, nk, mnk, tnk, p, spval, npval)         
         call selic(y, x, n, p, tol, k, nk, mnk, tnk,
     *        pval, bic, aic, wbic, waic, spval(1:npval), npval)
      end if
      end

      
c     This is the PVT selection method. Once the p-values are computed, it is just a knot
c     selection based on a threshold. No additional estimation is needed
c     The threshold is t. Knots with p-values less than t are kept. Since a p-value greater
c     than 1 is a missing value, knots with such p-values are removed.
c     Note that there is no need to have one function for both groups.
c     We can apply it to each group
c     separately. The selected model is returned through the knot selection matrix w.
      
      subroutine selpvt(p, nk, mnk, t, pval, w)

      integer p, nk(p), mnk, w(mnk, p), i, j, l
      double precision pval(mnk,p), t
      
      do i=1,p
         w(:,i) = 0
         if (nk(i) .gt. 0) then
            l = 1
            do j=1,nk(i)
               if (pval(j,i) .le. t) then
                  w(l,i) = j
                  l = l+1
               end if
            end do
         end if
      end do
      end

c     The function select the best model based on BIC and AIC.
c     The selected model is returned  through the knot selection matrices
c     waic for the AIC criterion and wbic for the BIC.
      
      subroutine selic(y, x, n, p, tol, k, nk, mnk, tnk,
     *     pval, bic, aic, wbic, waic, spval, npval)
      
      integer n, p, mnk, tnk, nk(p), npval, j, i, l, s
      integer w(mnk, p), wbic(mnk, p), waic(mnk, p)
      double precision x(n,p), y(n), tol, minpv, k(mnk, p), spval(npval)
      double precision pval(mnk, p), aic(tnk+1), bic(tnk+1)
      double precision aicsel, bicsel
     
      w(:,:) = 0

      call modfitsel(y, x, p, n, tol, k, nk, mnk, tnk,
     *     w, bic(1), aic(1))

      bicsel = bic(1)
      aicsel = aic(1)
      wbic = w
      waic = w

      do i=1,npval
         minpv = spval(i)
         do j=1,p
            if (nk(j) .gt. 0) then
               l = 1
               do s=1,nk(j)
                  if (pval(s,j) .le. minpv) then
                     w(l,j) = s
                     l = l+1
                  end if
               end do
            end if
         end do
         call modfitsel(y, x, p, n, tol, k, nk, mnk, tnk,
     *        w, bic(i+1), aic(i+1))        
         if (aic(i+1) .lt. aicsel) then
            aicsel = aic(i+1)
            waic = w
         end if
         if (bic(i+1) .lt. bicsel) then
            bicsel = bic(i+1)
            wbic = w
         end if         
      end do
      end
      
c     The function select the best causal model based on BIC and AIC.
c     The selected model is returned  through the knot selection matrices
c     w0aic and w1aic for the AIC criterion and w0bic and w1bic for the BIC.
      
      subroutine selicc(y0, y1, x0, x1, n0, n1, p, tol, 
     *     k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1, pval0,
     *     pval1, bic, aic, w0bic, w1bic, w0aic, w1aic, npval)
      
      integer n0, n1, p, mnk0, tnk0, nk0(p), mnk1, tnk1, nk1(p), npval
      integer w0(mnk0, p), w1(mnk1, p), j, i, l, s
      integer w0bic(mnk0, p), w1bic(mnk1, p)
      integer w0aic(mnk0, p), w1aic(mnk1, p)
      double precision x0(n0,p), x1(n1,p), y0(n0), y1(n1), tol, minpv
      double precision k0(mnk0, p), spval(tnk0+tnk1)
      double precision pval0(mnk0, p), pval1(mnk1, p)
      double precision k1(mnk1, p), aic(tnk0+tnk1+1), bic(tnk0+tnk1+1)
      double precision aicsel, bicsel
     
      call vecpvalc(pval0, nk0, mnk0, tnk0, pval1, nk1, mnk1, tnk1, p,
     *     spval, npval)

      w0(:,:) = 0
      w1(:,:) = 0
      
      call cmodfitsel(y0, y1, x0, x1, p, n0, n1, tol,
     *     k0, nk0, mnk0, tnk0,  k1, nk1, mnk1, tnk1,
     *     w0, w1, bic(1), aic(1)) 
      bicsel = bic(1)
      aicsel = aic(1)
      w0bic = w0
      w1bic = w1
      w0aic = w0
      w1aic = w1
      do i=1,npval
         minpv = spval(i)
         do j=1,p
            if (nk0(j) .gt. 0) then
               l = 1
               do s=1,nk0(j)
                  if (pval0(s,j) .le. minpv) then
                     w0(l,j) = s
                     l = l+1
                  end if
               end do
            end if
            if (nk1(j) .gt. 0) then
               l = 1
               do s=1,nk1(j)
                  if (pval1(s,j) .le. minpv) then
                     w1(l,j) = s
                     l = l+1
                  end if
               end do
            end if            
         end do
         call cmodfitsel(y0, y1, x0, x1, p, n0, n1, tol,
     *        k0, nk0, mnk0, tnk0,  k1, nk1, mnk1, tnk1,
     *        w0, w1, bic(i+1), aic(i+1)) 
         if (aic(i+1) .lt. aicsel) then
            aicsel = aic(i+1)
            w0aic = w0
            w1aic  = w1
         end if
         if (bic(i+1) .lt. bicsel) then
            bicsel = bic(i+1)
            w0bic = w0
            w1bic  = w1
         end if         
      end do
      end

c     This function vectorizes the pvalues of causal models and return the sorted ones
c     The output npval indicates how many p-values are valid (not NA's)
      
      subroutine vecpvalc(pval0, nk0, mnk0, tnk0, 
     *     pval1, nk1, mnk1, tnk1, p, spval, npval)

      integer p, nk0(p), nk1(p), mnk0, mnk1, tnk0, tnk1
      integer npval0, npval1, npval, nna, i, j, l
      double precision pval0(mnk0,p), pval1(mnk1,p)
      double precision spval0(tnk0), spval1(tnk1), spval(tnk1+tnk0)
      
      call vecpval(pval0, nk0, mnk0, tnk0, p, spval0, npval0)
      call vecpval(pval1, nk1, mnk1, tnk1, p, spval1, npval1)
      npval = npval0+npval1
      spval(1:npval1) = spval1(1:npval1)
      spval((npval1+1):(npval1+npval0)) = spval0(1:npval0)
      call qsort3(spval, 1, npval)
      end

c     This function vectorizes the pvalues of models and return the sorted ones
c     The output npval indicates how many p-values are valid (not NA's)
      
      subroutine vecpval(pval, nk, mnk, tnk, p, spval, npval)

      integer p, nk(p), mnk, tnk, npval, nna, i, j, l
      double precision pval(mnk,p), spval(tnk)

      nna = count(pval .gt. 1.0d0)
      npval = tnk-nna

      l = 1
      do i=1,p
         if (nk(i) .gt. 0) then
            do j=1,nk(i)
               if (pval(j,i) .le. 1.0d0) then
                  spval(l) = pval(j,i)
                  l = l+1
               end if
            end do
         end if
      end do
      call qsort3(spval, 1, npval)
      end
      
c     This is a simplified version of the modelfit function for causal model
c     It only returns what is needed for the selection:
c     - new knots: k0s and k1s
c     - new number of knots: nk0s, nk1s
c     - AIC and BIC
c     It does not compute the variance of LSE, so it is faster
      
      subroutine cmodfitsel(y0, y1, x0, x1, p, n0, n1, tol, 
     *     k0, nk0, mnk0, tnk0, k1, nk1, mnk1, tnk1,
     *     w0, w1, bic, aic)
      integer p, n0, nk0(p), mnk0, tnk0, n1, nk1(p), mnk1, tnk1
      integer w0(mnk0,p), w1(mnk1,p), rk0, rk1
      integer nk0s(p), nk1s(p), tnk0s, tnk1s, i
      integer piv0(tnk0+p+1), piv1(tnk1+p+1)
      double precision k1s(mnk1,p), k0s(mnk0,p), k1(mnk1,p), k0(mnk0,p)
      double precision y0(n0), y1(n1), x0(n0,p), x1(n1,p)
      double precision tol, aic, bic, ll, pi
      double precision e0(n0), e1(n1), b0(tnk0+p+1), b1(tnk1+p+1)
      double precision v0(tnk0+p+1,tnk0+p+1), v1(tnk1+p+1,tnk1+p+1)

      call updatek(k0, p, nk0, mnk0, w0, k0s, nk0s)
      call updatek(k1, p, nk1, mnk1, w1, k1s, nk1s)
      tnk0s = sum(nk0s)
      tnk1s = sum(nk1s)
      
      call splinefit(y0, x0, n0, p, k0s, nk0s, mnk0, tnk0s,
     *     tol, rk0, b0(1:(tnk0s+p+1)), piv0(1:(tnk0s+p+1)), e0, -10,
     *     v0(1:(tnk0s+p+1), 1:(tnk0s+p+1)))
      call splinefit(y1, x1, n1, p, k1s, nk1s, mnk1, tnk1s,
     *     tol, rk1, b1(1:(tnk1s+p+1)), piv1(1:(tnk1s+p+1)), e1, -10,
     *     v1(1:(tnk1s+p+1), 1:(tnk1s+p+1)))
      
      pi = 4.d0*datan(1.d0)
      ssr = sum(e1**2) + sum(e0**2)
      ll = dble(n1+n0)*(log(2.0d0 * pi) + 1 - log(dble(n1+n0)) +
     *     log(ssr))
      aic = ll+2.0d0*dble(rk0+rk1+1)
      bic = ll+log(dble(n1+n0))*dble(rk0+rk1+1)
      end


c     This is a simplified version of the modelfit function for all Spline models
c     It only returns what is needed for the selection:
c     - new knots: k0s and k1s
c     - new number of knots: nk0s, nk1s
c     - AIC and BIC
c     It does not compute the variance of LSE, so it is faster
      
      subroutine modfitsel(y, x, p, n, tol, k, nk, mnk, tnk,
     *     w, bic, aic)
      integer p, n, nk(p), mnk, tnk, w(mnk0,p), rk
      integer nks(p), tnks, piv(tnk+p+1), i 
      double precision ks(mnk,p), k(mnk,p), y(n), x(n,p)
      double precision tol, aic, bic, ll, pi
      double precision e(n), b(tnk+p+1), v(tnk+p+1,tnk+p+1)

      call updatek(k, p, nk, mnk, w, ks, nks)
      tnks = sum(nks)
      
      call splinefit(y, x, n, p, ks, nks, mnk, tnks,
     *     tol, rk, b(1:(tnks+p+1)), piv(1:(tnks+p+1)), e, -10,
     *     v(1:(tnks+p+1), 1:(tnks+p+1)))
      
      pi = 4.d0*datan(1.d0)
      ssr = sum(e**2)
      ll = dble(n)*(log(2.0d0 * pi) + 1 - log(dble(n)) +
     *     log(ssr))
      aic = ll+2.0d0*dble(rk+1)
      bic = ll+log(dble(n))*dble(rk+1)
      end
      

