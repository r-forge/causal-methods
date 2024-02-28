      subroutine edist(x1,x2,n,k,dist)
      integer n, i
      double precision x1(k), x2(n,k), dist(n)
      do i=1,n
         call exprsum((x1-x2(i,:))**2, k, dist(i))
      end do
      end

      subroutine allmin(x,n,nm,loc)
      integer n, nm, loc(n), indi
      double precision x(n), tmp1, tmp2
      logical sel(n)

      nm = 0
      sel = .true.      
      do
         indi = minloc(x,dim=1,mask=sel)
         sel(indi) = .false.                     
         if (nm == 0) then
            tmp1 = x(indi)
            nm = 1
            loc(1) = indi
            cycle
         end if
         tmp2 = x(indi)
         if (tmp1 == tmp2) then
            nm = nm+1
            loc(nm) = indi
            tmp1 = tmp2
         else
            exit
         end if
      end do      
      end

c     x1 is for group 1 and x2 fro group 2 (either could be treated)
c     Here units in x2 are matched to THE unit in x1      
      
      subroutine findnni(x1, x2, tol, n, k, minn, nn, ind)
      integer n, k, minn, nn, ind(n), indi
      double precision x1(k), x2(n,k), dist(n), tmp1, tmp2, tol
      logical sel(n)
      
      nn = 0
      sel = .true.
      ind = 0
      
      call edist(x1,x2,n,k,dist)
      do
         indi = minloc(dist,dim=1,mask=sel)
         if (indi == 0) then
            exit
         end if
         sel(indi) = .false.                     
         if (nn==0) then
            tmp1 = sqrt(dist(indi))
            ind(1) = indi
            nn = 1
            cycle
         end if         
         tmp2 = sqrt(dist(indi))
         if (nn < minn) then
            nn = nn+1 
            ind(nn) = indi
         else
            if (abs(tmp2-tmp1)<tol) then
               nn = nn+1
               ind(nn) = indi
            else
               exit
            end if            
         end if
         tmp1 = tmp2
      end do
      end

c     Here units in x2 are matched to every unit in x1
c     Here we also estimate the missing Y(2) [ Y(2)|Z=1 ] by averaging the Y's
c     from group 2 matched to group 1.
c     Y(2) is therefore n1 x 1 using Yi from group 2.
c     Also, rk is the weights attached to y2. it sums to n1
      
      subroutine findnn(x1, x2, y2, tol, n1, n2, k, minn, rk, nn, ind,
     *     y2h)
      integer n1, n2, k, minn, nn(n1), ind(n1,n2), i
      double precision x1(n1,k), x2(n2,k), tol, rk(n2), y2(n2), y2h(n1)

      rk = 0.0d0
      do i=1,n1
         call findnni(x1(i,:), x2, tol, n2, k, minn, nn(i), ind(i,:))
         rk(ind(i,1:nn(i))) = rk(ind(i,1:nn(i))) + 1.0d0/nn(i)
         call exprsum(y2(ind(i,1:nn(i))), nn(i), y2h(i))
         y2h(i) = y2h(i)/nn(i)
      end do
      end

c     kern=1 means Gaussian, else it is Epanechnikok
      
      subroutine kernel(x, n, h, kern, valk)
      integer n, kern
      double precision x(n), h, valk(n), cst
         
      if (kern == 1) then
         cst = 4.0d0*atan(1.0d0)
         cst = 1.0d0/sqrt(2.0d0*cst)
         valk = cst*exp(-0.5*(x/h)**2)
      else
         where (abs(x)>h)
            valk = 0.0d0
         else where
            valk = 0.75*(1.0d0-(x/h)**2)
         end where
      end if
      end

c     It computes the missing Y(2) for a single individual in group 1
c     using local linear regression matching
c     p1 is the propensity score for the individual in group 1
c     p2 are the n propensity scores in group 2
c     y2 is the associated observed Y(2) in group 2
c     kern is 1 for Gaussian and for any other value it is the Epanechnikov Kernel
c     the output y2h is the estimated Y(2) for the group 1 individual
      
      subroutine llri(p1, p2, y2, n, h, kern, info, y2h)
      integer n, info, kern
      double precision p1, p2(n), y2(n), y2h, h, x1, r, c
      double precision valk(n), dist(n), s0, t0, s2, t1

      info = 0
      y2h = 0.0d0
      if (h >= 10000) then
         dist = 0.0d0
      else
         dist = p1-p2
      end if
      call kernel(dist, n, h, kern, valk)
      if (kern == 1) then
         c = 0.35
      else
         c = 0.3125
      end if      
      call exprsum(valk, n, s0)
      call exprsum(valk*y2, n, t0)
      if (abs(s0) <= epsilon(p1)) then
         info = 1
      else
         call exprsum(valk*p2, n, x1)
         x1 = x1/s0
         dist = p2-x1
         call exprsum(valk*dist**2, n, s2)
         call exprsum(valk*dist*y2, n, t1)
         if ((h>0) .and. (h<10000)) then
            r = h*abs(p1-x1)*c
         else
            r = 999999.0d0
         end if
         y2h = t0/s0 + t1*(p1-x1)/(s2+r)
      end if
      end

      subroutine ksum(x,n,s)
      integer n, i
      double precision x(n), s, c, y, t
      s = 0.0d0
      c = 0.0d0
      do i=1,n
         y = x(i) - c
         t = s + y
         c = (t-s) - y
         s = t
      end do
      end
      
      subroutine llr(p1, p2, y2, n1, n2, h, kern, info, y2h)
      integer n1, n2, info(n1), i, kern      
      double precision p1(n1), p2(n2), y2(n2), y2h(n1), h

      do i=1,n1
         call llri(p1(i), p2, y2, n2, h, kern, info(i), y2h(i))
      end do

      end

      subroutine cvfct(p,y,n,h,kern,nh, cv)
      integer n, i, j, ind(n), info, nh, ind2(n-1), kern
      double precision p(n), y(n), h(nh), cv(nh), yi, yh
      
      ind = (/(i, i=1,n, 1)/)
      cv = 0.0d0
      do j=1,nh
         do i=1,n
            ind2 = pack(ind, ind/=i)
            call llri(p(i), p(ind2), y(ind2), n-1, h(j), kern, info, yh)
            cv(j) = cv(j) + (y(i)-yh)**2
         end do
         cv(j) = cv(j)/n
      end do
      end  
      
      subroutine gridcv(p,y,n,kern,lh,uh,nh,info,h,cv)
      integer n, ind(n), nh, wmin, info, i, kern
      double precision p(n), y(n), h(3), cv(3), uh, lh, by
      double precision h2(nh), cv2(nh)
      info = 0
      by = (uh-lh)/(nh-1)      
      do i=1,nh
         h2(i) = lh + by*(dble(i)-1.0d0)
      end do
      call cvfct(p,y,n,h2,kern,nh,cv2)
      wmin = minloc(cv2,1)
      if (wmin==1 .or. wmin==nh) then
         info = 1
         h = h2(wmin)
         cv = cv2(wmin)
      else
         cv = cv2((/wmin-1, wmin, wmin+1/))
         h = h2((/wmin-1, wmin, wmin+1/))
      end if
      end

c     it is assumed that cvg(2) is smaller than cvg(1) and cvg(3)
c     it is assumed that hg(1)<hg(2)<hg(3)
c     the above is satisfied if hg ang cvg comes from gridcv subroutine
      
      subroutine brentcv(p1,y1,n,kern,tol,maxit,h,cv,
     *     info, minh,mincv)
      integer i, n, maxit, info, kern
      double precision p1(n), y1(n), tol, h(3), cv(3), minh, mincv
      double precision cgold, zeps, a, b, c, x, w, v, fw, fv, fx, xm
      double precision tol1, tol2, e, r, q, p, etemp, d, u, fu, u2(1)
      double precision fu2(1)
      
      cgold = 0.3819660; e = 0.0d0; d = 0.0d0;  mincv = cv(2)
      minh = h(2); info = 0
      zeps = epsilon(cgold)*(1.0e-3)
      a = h(1); b = h(3)
      x = h(2); w = h(2); v = h(2)
      fw = cv(2); fv = cv(2); fx = cv(2)
      do i=1,maxit
         xm = 0.5*(a+b)
         tol1 = abs(x)*tol+zeps
         tol2 = 2.0d0*tol1
         if (abs(x-xm) <= (tol2-0.5*(b-a))) then
            mincv=fx
            minh=x
            exit
         end if
         if (abs(e) > tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q - (x-w)*r
            q = 2.0d0*(q-r)
            if (q > 0.0d0) then
               p = -p
            end if
            q = abs(q); etemp = e; e = ed
            if (abs(p)>=abs(0.5*q*etemp) .or. p<=q*(a-x) .or.
     *           p>=q*(b-x)) then
               e = merge(a-x, b-x, x>=xm)
               d = cgold*e
            else
               d = p/q
               u = x+d
               if ( (u-a)<tol2 .or. (b-u)<tol2 )then
                  d = sign(tol1, xm-x)
               end if
            end if
         else
            e = merge(a-x, b-x, x>=xm)
            d = cgold*e
         end if
         u = merge(x+d, x+sign(tol1,d), abs(d) >= tol1)
         u2 = u;          
         call cvfct(p1,y1,n,u2,kern,1,fu2)
         fu = fu2(1)
         if (fu <= fx) then
            if (u >= x) then
               a = x
            else
               b = x
            end if
            call shift3(v,w,x,u)
            call shift3(fv,fw,fx,fu)
         else
            if (u < x) then
               a=u
            else
               b=u
            end if
            if (fu <= fw .or. w == x) then
               v=w; w=u; fv=fw; fw=fu
            else if (fu <=fv .or. v == x .or. v == w) then
               v=u; fv=fu
            end if
         end if
         if (i >= maxit) then
            info = 1
         end if
      end do
      end

      

      subroutine shift3(a,b,c,d)
      double precision a, b, c, d
      a=b; b=c; c=d
      end
   

      subroutine exprsum(x, n, value)
      integer n, i
      double precision x(n), value
      real*16 s

      s = 0.0
      do i =1,n
         s = s+x(i)
      end do
      value = dble(s)

      end

