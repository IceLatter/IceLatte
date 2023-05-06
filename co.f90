MODULE COMPU
CONTAINS
SUBROUTINE LATTICECONSTANTink(lr,lk)
IMPLICIT NONE
INTEGER, PARAMETER              ::             dp = SELECTED_REAL_KIND(15,14)      
REAL (KIND=dp)                  ::             lr(3,3)                             
REAL (KIND=dp)                  ::             lk(3,3)                             
REAL (KIND=dp)                  ::             tr(4)                               
REAL (KIND=dp)                  ::             pi                                  

tr(1) = lr(2,2)*lr(3,3)-lr(3,2)*lr(2,3)
tr(2) = lr(3,1)*lr(2,3)-lr(2,1)*lr(3,3)
tr(3) = lr(2,1)*lr(3,2)-lr(3,1)*lr(2,2)
tr(4) = lr(1,1)*tr(1)+lr(1,2)*tr(2)+lr(1,3)*tr(3)
lk(1,1) = tr(1)/tr(4)
lk(1,2) = tr(2)/tr(4)
lk(1,3) = tr(3)/tr(4)
tr(1) = lr(3,2)*lr(1,3)-lr(1,2)*lr(3,3)
tr(2) = lr(1,1)*lr(3,3)-lr(3,1)*lr(1,3)
tr(3) = lr(3,1)*lr(1,2)-lr(1,1)*lr(3,2)
tr(4) = lr(2,1)*tr(1)+lr(2,2)*tr(2)+lr(2,3)*tr(3)
lk(2,1) = tr(1)/tr(4)
lk(2,2) = tr(2)/tr(4)
lk(2,3) = tr(3)/tr(4)
tr(1) = lr(1,2)*lr(2,3)-lr(2,2)*lr(1,3)
tr(2) = lr(2,1)*lr(1,3)-lr(1,1)*lr(2,3)
tr(3) = lr(1,1)*lr(2,2)-lr(2,1)*lr(1,2)
tr(4) = lr(3,1)*tr(1)+lr(3,2)*tr(2)+lr(3,3)*tr(3)
lk(3,1) = tr(1)/tr(4)
lk(3,2) = tr(2)/tr(4)
lk(3,3) = tr(3)/tr(4)
pi = DACOS(-1.0d0)
lk = lk * 2.0d0 * pi
RETURN
END SUBROUTINE LATTICECONSTANTink
SUBROUTINE MULTIPLICATION_MATRIX(di,ou,i1,i2)
IMPLICIT NONE
INTEGER, PARAMETER              ::             dp = SELECTED_REAL_KIND(15,14)      
INTEGER                         ::             di                                  
COMPLEX (KIND=dp), ALLOCATABLE  ::             ou(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             i1(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             i2(:,:)                             
INTEGER                         ::             i, j, k
ou = DCMPLX(0.0d0, 0.0d0)
DO i = 1, di, 1
   DO j = 1, di, 1
      DO k = 1, di, 1
         ou(i,j) = ou(i,j)+i1(j,k)*i2(k,j)
      END DO
   END DO
END DO
RETURN
END SUBROUTINE MULTIPLICATION_MATRIX
END MODULE COMPU
MODULE TRANS
CONTAINS
SUBROUTINE TRANSPORT(n1,n2,n3,dk,ra,ch,it,nm,te,nb,po,tv,lr,ei,ks,ol,t2,t3)
USE COMPU
IMPLICIT NONE
INTEGER, PARAMETER              ::             dp = SELECTED_REAL_KIND(15,14)      
INTEGER                         ::             nk                                  
INTEGER                         ::             n1                                  
INTEGER                         ::             n2                                  
INTEGER                         ::             n3                                  
REAL (KIND=dp)                  ::             dk(3)                               
INTEGER                         ::             ra(2)                               
REAL (KIND=dp)                  ::             ch(3)                               
REAL (KIND=dp)                  ::             it                                  
INTEGER                         ::             nm                                  
REAL (KIND=dp)                  ::             te                                  
INTEGER                         ::             nb                                  
INTEGER                         ::             po(3)                               
INTEGER, ALLOCATABLE            ::             tv(:,:)                             
REAL (KIND=dp)                  ::             lr(3,3)                             
REAL (KIND=dp), ALLOCATABLE     ::             ei(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             ks(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             ol(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             t2(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             t3(:)                               
INTEGER                         ::             i, j, k, l, m
REAL (KIND=dp)                  ::             kc(3)                               
COMPLEX (KIND=dp)               ::             ic                                  
COMPLEX (KIND=dp), ALLOCATABLE  ::             sr(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             ll(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             ls(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             rl(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             rs(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             th(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             rg(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             ag(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             m1(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             m2(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             m3(:,:)                             
COMPLEX (KIND=dp)               ::             pf                                  
COMPLEX (KIND=dp), ALLOCATABLE  ::             l1(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             l2(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             r1(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             r2(:,:)                             
REAL (KIND=dp)                  ::             fd(3)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             tr(:,:)                             
INTEGER                         ::             lwork, info
INTEGER, ALLOCATABLE            ::             ipiv(:)
COMPLEX (KIND=dp), ALLOCATABLE  ::             work(:)
ic = DCMPLX(0.0d0, 1.0d0)
kc(1) = (n1-1)*dk(1)
kc(2) = (n2-1)*dk(2)
kc(3) = (n3-1)*dk(3)
k = ra(2) - ra(1) + 1
ALLOCATE (sr(k,k))
ALLOCATE (ls(k,k))
ALLOCATE (rs(k,k))
ALLOCATE (rg(k,k))
ALLOCATE (ag(k,k))
ALLOCATE (th(nb,nb))
ALLOCATE (ll(3,k,k))
ALLOCATE (rl(3,k,k))
ALLOCATE (m1(k,k))
ALLOCATE (m2(k,k))
ALLOCATE (m3(k,k))
ALLOCATE (ipiv(k))
lwork = k
ALLOCATE (work(lwork))
ALLOCATE (l1(k,k))
ALLOCATE (l2(k,k))
ALLOCATE (r1(k,k))
ALLOCATE (r2(k,k))
ALLOCATE (tr(k,k))
t2 = DCMPLX(0.0d0, 0.0d0)
t3 = DCMPLX(0.0d0, 0.0d0)
write(unit=10,fmt=*)'before openmp'
!$OMP PARALLEL DO PRIVATE(i,j,l,pf,th,ll,sr,rl,ls,rs,rg,ag,l1,l2,r1,r2,ipiv,work,lwork,info,fd,tr,m1,m2,m3)&
!$OMP             SHARED(nm,nb,po,ic,kc,tv,lr,ol,ks,k,ei,ra,ch,te) REDUCTION(+:t2,t3)
DO i = 1, nm, 1
   th = DCMPLX(0.0d0, 0.0d0)
   DO j = 1, nb, 1
      DO l = 1, (po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1), 1
         pf = EXP(ic*(kc(1)*(tv(1,l)*lr(1,1)+tv(2,l)*lr(2,1)+tv(3,l)*lr(3,1))&
              +kc(2)*(tv(1,l)*lr(1,2)+tv(2,l)*lr(2,2)+tv(3,l)*lr(3,2))&
              +kc(3)*(tv(1,l)*lr(1,3)+tv(2,l)*lr(2,3)+tv(3,l)*lr(3,3))))
         th(:,j) = th(:,j)+pf*((ei(i)+ic)*ol(:,j,l)-ks(:,j,l))
      END DO
   END DO
   ll(1,:,:) = th(ra(1):ra(2),ra(1)-k:ra(2)-k)                                     
   m1 = th(ra(1)-k:ra(2)-k,ra(1)-k:ra(2)-k)                                        
   CALL ZGETRF(ra(2)-ra(1)+1,ra(2)-ra(1)+1,m1,ra(2)-ra(1)+1,ipiv,info)
   CALL ZGETRI(ra(2)-ra(1)+1,m1,ra(2)-ra(1)+1,ipiv,work,lwork,info)
   ll(2,:,:) = m1
   ll(3,:,:) = th(ra(1)-k:ra(2)-k,ra(1):ra(2))                                     
   sr = th(ra(1):ra(2),ra(1):ra(2))                                                
   rl(1,:,:) = th(ra(1):ra(2),ra(1)+k:ra(2)+k)                                     
   m1 = th(ra(1)+k:ra(2)+k,ra(1)+k:ra(2)+k)                                        
   CALL ZGETRF(ra(2)-ra(1)+1,ra(2)-ra(1)+1,m1,ra(2)-ra(1)+1,ipiv,info)
   CALL ZGETRI(ra(2)-ra(1)+1,m1,ra(2)-ra(1)+1,ipiv,work,lwork,info)
   rl(2,:,:) = m1
   rl(3,:,:) = th(ra(1)+k:ra(2)+k,ra(1):ra(2))                                     
   m2 = ll(1,:,:)
   m3 = ll(2,:,:)
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   m2 = m1
   m3 = ll(3,:,:)
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   ls = m1                                                                         
   m2 = rl(1,:,:)
   m3 = rl(2,:,:)
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   m2 = m1
   m3 = rl(3,:,:)
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   rs = m1                                                                         
   rg = sr-ls-rs
   CALL ZGETRF(ra(2)-ra(1)+1,ra(2)-ra(1)+1,rg,ra(2)-ra(1)+1,ipiv,info)
   CALL ZGETRI(ra(2)-ra(1)+1,rg,ra(2)-ra(1)+1,ipiv,work,lwork,info)
   DO j = 1, ra(2)-ra(1)+1, 1
      ag(j,:) = DCONJG(rg(:,j))
   END DO
   l1 = ls
   DO j = 1, ra(2)-ra(1)+1, 1
      l2(j,:) = DCONJG(l1(:,j))
   END DO
   m2 = DCMPLX(0.0d0, 0.0d0)
   DO j = 1, ra(2)-ra(1)+1
      m2(j,j) = DCMPLX(0.0d0, 1.0d0)
   END DO
   m3 = l1 - l2
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   ls = m1
   r1 = rs
   DO j = 1, ra(2)-ra(1)+1, 1
      r2(j,:) = DCONJG(r1(:,j))
   END DO
   m2 = DCMPLX(0.0d0, 0.0d0)
   DO j = 1, ra(2)-ra(1)+1, 1
      m2(j,j) = DCMPLX(0.0d0, 1.0d0)
   END DO
   m3 = r1 - r2
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   rs = m1
   fd(1) = 1.0d0/(EXP((ei(i)-ch(1))/te)+1.0d0)
   fd(2) = 1.0d0/(EXP((ei(i)-ch(3))/te)+1.0d0)
   IF (fd(1) > fd(2)) THEN
      fd(3) = fd(1)
      fd(1) = fd(2)
      fd(2) = fd(3)
      fd(3) = fd(2)-fd(1)
   ELSE IF (fd(1) == fd(2)) THEN
           fd(3) = 1.0d0
   ELSE
      fd(3) = fd(2)-fd(1)
   END IF
   m2 = ls
   m3 = ag
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   m2 = m1
   m3 = rs
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   m2 = m1
   m3 = rg
   CALL MULTIPLICATION_MATRIX(ra(2)-ra(1)+1,m1,m2,m3)
   tr = m1
   DO j = 1, ra(2)-ra(1)+1, 1
      t2(i) = t2(i)+tr(j,j)
   END DO
   t2(i) = t2(i)*fd(3)
   tr = rg - ag
   DO j = 1, ra(2)-ra(1)+1, 1
      t3(i) = t3(i)+tr(j,j)
   END DO
END DO
!$OMP END PARALLEL DO
write(unit=10,fmt=*)'after openmp'
DEALLOCATE (sr)
DEALLOCATE (ll)
DEALLOCATE (ls)
DEALLOCATE (rg)
DEALLOCATE (ag)
DEALLOCATE (rl)
DEALLOCATE (rs)
DEALLOCATE (th)
DEALLOCATE (ipiv)
DEALLOCATE (work)
DEALLOCATE (l1)
DEALLOCATE (l2)
DEALLOCATE (r1)
DEALLOCATE (r2)
DEALLOCATE (tr)
DEALLOCATE (m1)
DEALLOCATE (m2)
DEALLOCATE (m3)
RETURN
END SUBROUTINE TRANSPORT
END MODULE TRANS
PROGRAM QUANTUMTRANSPORT_NEGF
USE MPI
USE COMPU
USE TRANS
IMPLICIT NONE
INTEGER, PARAMETER              ::             dp = SELECTED_REAL_KIND(15,14)      
INTEGER                         ::             ra(2)                               
REAL (KIND=dp)                  ::             ch(3)                               
REAL (KIND=dp)                  ::             it                                  
REAL (KIND=dp)                  ::             mi                                  
INTEGER                         ::             nm                                  
REAL (KIND=dp)                  ::             te                                  
INTEGER                         ::             km(3)                               
INTEGER                         ::             ps(3)                               
REAL (KIND=dp)                  ::             lr(3,3)                             
REAL (KIND=dp)                  ::             lk(3,3)                             
INTEGER                         ::             na                                  
REAL (KIND=dp), ALLOCATABLE     ::             ap(:,:)                             
INTEGER, ALLOCATABLE            ::             ns(:,:)                             
INTEGER                         ::             po(3)                               
INTEGER                         ::             nb                                  
REAL (KIND=dp)                  ::             p1(3)                               
INTEGER, ALLOCATABLE            ::             p2(:,:)                             
INTEGER, ALLOCATABLE            ::             tv(:,:)                             
COMPLEX (KIND=dp), ALLOCATABLE  ::             ok(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             oo(:,:,:)                           
INTEGER                         ::             i, j, k, l, m, n                    
INTEGER                         ::             tc(4)                               
REAL (KIND=dp)                  ::             ti                                  
REAL (KIND=dp)                  ::             dk(3)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             ks(:,:,:)                           
COMPLEX (KIND=dp), ALLOCATABLE  ::             ol(:,:,:)                           
REAL (KIND=dp), ALLOCATABLE     ::             ei(:)                               
REAL (KIND=dp)                  ::             pi                                  
COMPLEX (KIND=dp), ALLOCATABLE  ::             t2(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             c2(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             c3(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             t3(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             c4(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             c5(:)                               
COMPLEX (KIND=dp), ALLOCATABLE  ::             cc(:)                               
INTEGER                         ::             world_size, world_rank, ierr        
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,world_size,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,ierr)
IF (world_rank == 0) THEN
   CALL SYSTEM_CLOCK(count_rate=tc(1))
   ti = REAL(tc(1))
END IF
IF (world_rank == 0) THEN
   CALL SYSTEM_CLOCK(tc(2))
   CALL SYSTEM_CLOCK(tc(3))
END IF
IF (world_rank == 0) THEN
   OPEN (UNIT=3, FILE='INPUT.dat', STATUS='OLD')
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) ra
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) ch(1)
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) ch(2)
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) ch(3)
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) it
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) mi
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) nm
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) te
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) km
   READ (UNIT=3, FMT=*)
   READ (UNIT=3, FMT=*) ps
   CLOSE (UNIT=3)
   ch = ch * 3.67493036006967667962692106984d-2
   it = it * 3.67493036006967667962692106984d-2
   mi = mi * 3.67493036006967667962692106984d-2
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ra,2,MPI_INT,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ch,3,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(it,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(mi,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nm,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(te,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(km,3,MPI_INT,0,MPI_COMM_WORLD,ierr)
IF (world_rank == 0) THEN
   OPEN (UNIT=11, FILE='final.config', STATUS='OLD')
   READ (UNIT=11, FMT=*) na
   READ (UNIT=11, FMT=*)
   DO i = 1, 3, 1
      READ (UNIT=11, FMT=*) p1
      lr(i,:) = p1
   END DO
   READ (UNIT=11, FMT=*)
   ALLOCATE (ap(na,6))
   ALLOCATE (ns(3,na))
   DO i = 1, na, 1
      READ (UNIT=11, FMT=*) j, p1, k, l, m
      ap(i,4:6) = p1
      ns(1,i) = i
   END DO
   CLOSE (UNIT=11)
   DO i = 1, na, 1
      ap(i,1) = ap(i,4)*lr(1,1)+ap(i,5)*lr(2,1)+ap(i,6)*lr(3,1)
      ap(i,2) = ap(i,4)*lr(1,2)+ap(i,5)*lr(2,2)+ap(i,6)*lr(3,2)
      ap(i,3) = ap(i,4)*lr(1,3)+ap(i,5)*lr(2,3)+ap(i,6)*lr(3,3)
   END DO
   DO i = 1, na, 1
      DO j = i+1, na, 1
         IF (ap(i,ps(1)) > ap(j,ps(1))) THEN
             k = ns(1,i)
             ns(1,i) = ns(1,j)
             ns(1,j) = k
         ELSE IF (ap(i,ps(1)) == ap(j,ps(1))) THEN
                 IF (ap(i,ps(2)) > ap(j,ps(2))) THEN
                    k = ns(1,i)
                    ns(1,i) = ns(1,j)
                    ns(1,j) = k
                 ELSE IF (ap(i,ps(2)) == ap(i,ps(2))) THEN
                    IF (ap(i,ps(3)) > ap(j,ps(3))) THEN
                       k = ns(1,i)
                       ns(1,i) = ns(1,j)
                       ns(1,j) = k
                    END IF
                 END IF
         END IF
      END DO
   END DO
   DEALLOCATE (ap)
   ALLOCATE (p2(2,na))
   OPEN (UNIT=12, FILE='OUT.GAUSSIAN_BASIS_INDEX', FORM='UNFORMATTED', STATUS='OLD')
   READ (UNIT=12) i, j
   READ (UNIT=12) p2(:,1:na)
   DO i = 1, na, 1
      ns(2:3,i) = p2(:,i)
   END DO
   CLOSE (UNIT=12)
   DEALLOCATE (p2)
END IF
IF (world_rank == 0) THEN
   OPEN (UNIT=4, FILE='OUT.GAUSSIAN_H_T', FORM='UNFORMATTED', STATUS='OLD')
   OPEN (UNIT=6, FILE='OUT.GAUSSIAN_S_T', FORM='UNFORMATTED', STATUS='OLD')
   READ (UNIT=4) po, nb, lr
   READ (UNIT=6) po, nb, lr
   ALLOCATE (ok(nb,nb,(po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1)))
   ALLOCATE (oo(nb,nb,(po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1)))
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(po,3,MPI_INT,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nb,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(lr,3*3,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
ALLOCATE (tv(3,(po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1)))
ALLOCATE (ks(nb,nb,(po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1)))
ALLOCATE (ol(nb,nb,(po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1)))
IF (world_rank == 0) THEN
   DO i = 1, po(1)*2+1, 1
      DO j = 1, po(2)*2+1, 1
         DO k = 1, po(3)*2+1, 1
            l = (i-1)*(po(2)*2+1)*(po(3)*2+1)+(j-1)*(po(3)*2+1)+k
            READ (UNIT=4) tv(:,l)
            READ (UNIT=6) tv(:,l)
            DO m = 1, nb, 1
               READ (UNIT=4) ok(:,m,l)
               READ (UNIT=6) oo(:,m,l)
            END DO
         END DO
      END DO
   END DO
   DO i = 1, (po(1)*2+1)*(po(2)*2+1)*(po(3)*2+1), 1
      l = 0
      DO j = 1, na, 1
         DO k = 1, na, 1
            IF (ns(1,k) == j) THEN
               ks(:,l+1:l+ns(3,k)-ns(2,k)+1,i) = ok(:,ns(2,k):ns(3,k),i)
               ol(:,l+1:l+ns(3,k)-ns(2,k)+1,i) = oo(:,ns(2,k):ns(3,k),i)
               ks(l+1:l+ns(3,k)-ns(2,k)+1,:,i) = ok(ns(2,k):ns(3,k),:,i)
               ol(l+1:l+ns(3,k)-ns(2,k)+1,:,i) = oo(ns(2,k):ns(3,k),:,i)
               l = l + ns(3,k) - ns(2,k) + 1
            END IF
         END DO
      END DO
   END DO
   CLOSE (UNIT=4)
   CLOSE (UNIT=6)
   DEALLOCATE (ns)
   DEALLOCATE (ok)
   DEALLOCATE (oo)
END IF
print*,'l4'
IF (world_rank == 0) THEN
   CALL LATTICECONSTANTink(lr,lk)
   IF (km(1) == 1) THEN
      dk(1) = 0.0d0
   ELSE
      dk(1) = DSQRT(lk(1,1)**2+lk(1,2)**2+lk(1,3)**2)/DBLE(km(1)-1)
   END IF
   IF (km(2) == 1) THEN
      dk(2) = 0.0d0
   ELSE
      dk(2) = DSQRT(lk(2,1)**2+lk(2,2)**2+lk(2,3)**2)/DBLE(km(2)-1)
   END IF
   IF (km(3) == 1) THEN
      dk(3) = 0.0d0
   ELSE
      dk(3) = DSQRT(lk(3,1)**2+lk(3,2)**2+lk(3,3)**2)/DBLE(km(3)-1)
   END IF
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(dk,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
ALLOCATE (ei(nm))
IF (world_rank == 0) THEN
   DO i = 1, nm, 1
      ei(i) = mi+(i-1)*(ch(2)-mi)/DBLE(nm-1)
   END DO
END IF
print*,'l2'
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ei,nm,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
IF (world_rank == 0) THEN
   CALL SYSTEM_CLOCK(tc(4))
   OPEN (UNIT=10, FILE='output.dat', STATUS='UNKNOWN')
   WRITE (UNIT=10, FMT='(A81)') 'Time for reading input parameters, atom positions, re-organising the Hamiltonian'
   WRITE (UNIT=10, FMT='(A62,F20.4,A7)') 'matrice, generating k points and broadcasting all parameters:', (tc(4)-tc(3))/ti, 'seconds'
END IF
IF (world_rank == 0) THEN
   CALL SYSTEM_CLOCK(tc(2))
   CALL SYSTEM_CLOCK(tc(3))
END IF
ALLOCATE (t2(nm))
ALLOCATE (t3(nm))
write(unit=10,fmt=*)'before mpi'
IF (world_rank == 0) THEN
   ALLOCATE (c2(nm))
   c2 = DCMPLX(0.0d0, 0.0d0)
   ALLOCATE (c3(nm))
   ALLOCATE (c4(nm))
   c4 = DCMPLX(0.0d0, 0.0d0)
   ALLOCATE (c5(nm))
END IF
write(unit=10,fmt=*)'before mpi'
DO i = 1, km(1), 1
   DO j = 1, km(2), 1
      DO k = 1, km(3), 1
         IF (MOD((i-1)*km(2)*km(3)+(j-1)*km(3)+k-1, world_size) /= world_rank) CYCLE
         CALL TRANSPORT(i,j,k,dk,ra,ch,it,nm,te,nb,po,tv,lr,ei,ks,ol,t2,t3)
         IF (world_rank == 0) THEN
            c2 = c2 + t2
            c4 = c4 + t3
            DO l = 1, world_size-1, 1
               IF ((i-1)*km(2)*km(3)+(j-1)*km(3)+k-1+l == km(1)*km(2)*km(2)) EXIT
               m = (i-1)*km(2)*km(3)+(j-1)*km(3)+k + l + 100000
               CALL MPI_RECV(c3,nm,MPI_DOUBLE_COMPLEX,l,m,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
               c2 = c2 + c3
               n = (i-1)*km(2)*km(3)+(j-1)*km(3)+k + l + 1000000
               CALL MPI_RECV(c5,nm,MPI_DOUBLE_COMPLEX,l,n,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
               c4 = c4 + c5
            END DO
         ELSE
            m = (i-1)*km(2)*km(3)+(j-1)*km(3)+k + 100000
            CALL MPI_SEND(t2,nm,MPI_DOUBLE_COMPLEX,0,m,MPI_COMM_WORLD,ierr)
            n = (i-1)*km(2)*km(3)+(j-1)*km(3)+k + 1000000
            CALL MPI_SEND(t3,nm,MPI_DOUBLE_COMPLEX,0,n,MPI_COMM_WORLD,ierr)
         END IF
      END DO
   END DO
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
write(unit=10,fmt=*)'after mpi'
IF (world_rank == 0) THEN
   IF ((dk(1) /= 0) .AND. (dk(2) /= 0) .AND. (dk(3) /= 0)) THEN
      c2 = c2*dk(1)*dk(2)*dk(3)
      c4 = c4*dk(1)*dk(2)*dk(3)
   ELSE IF ((dk(1) /= 0) .AND. (dk(2) == 0) .AND. (dk(3) /= 0)) THEN
           c2 = c2*dk(1)*dk(3)
           c4 = c4*dk(1)*dk(3)
   ELSE IF ((dk(1) /= 0) .AND. (dk(2) /= 0) .AND. (dk(3) == 0)) THEN
           c2 = c2*dk(1)*dk(2)
           c4 = c4*dk(1)*dk(2)
   ELSE IF ((dk(1) == 0) .AND. (dk(2) /= 0) .AND. (dk(3) /= 0)) THEN
           c2 = c2*dk(2)*dk(3)
           c4 = c4*dk(2)*dk(3)
   ELSE IF ((dk(1) /= 0) .AND. (dk(2) == 0) .AND. (dk(3) == 0)) THEN
           c2 = c2*dk(1)
           c4 = c4*dk(1)
   ELSE IF ((dk(1) == 0) .AND. (dk(2) /= 0) .AND. (dk(3) == 0)) THEN
           c2 = c2*dk(2)
           c4 = c4*dk(2)
   ELSE IF ((dk(1) == 0) .AND. (dk(2) == 0) .AND. (dk(3) /= 0)) THEN
           c2 = c2*dk(3)
           c4 = c4*dk(3)
   END IF
END IF
IF (world_rank == 0) THEN
   CALL SYSTEM_CLOCK(tc(4))
   WRITE (UNIT=10, FMT='(A45,F8.4,A7)') 'Time for computing the transport properties:', (tc(4)-tc(3))/ti, 'seconds'
END IF
IF (world_rank == 0) THEN
   pi = DACOS(-1.0d0)
   ALLOCATE (cc(nm))
   cc = c2
   c2 = c2*(ch(2)-mi)/DBLE(nm-1)
   cc(1) = cc(1)/2.0d0
   DO i = nm, 1, -1
      cc(i) = cc(i)/2.0d0
      DO j = 1, i, 1
         cc(i) = cc(i) + cc(j)
      END DO
   END DO
   c4 = c4/2.0d0/pi*(-1.0d0)
END IF
IF (world_rank == 0) THEN
   OPEN (UNIT=8, FILE='charge_current.dat', STATUS='UNKNOWN')
   OPEN (UNIT=9, FILE='dos.dat', STATUS='UNKNOWN')
   WRITE (UNIT=8, FMT=*) 'Charge current with the unit of 2e^2/h'
   WRITE (UNIT=8, FMT='(A60)') 'The integrate of the charge current v.s. Fermi energy level'
   DO i = 1, nm, 1
      WRITE (UNIT=8, FMT=*) ei(i), c2(i)
   END DO
   WRITE (UNIT=8, FMT=*)
   WRITE (UNIT=8, FMT='(A59)') 'The integral of the charge current v.s. Fermi energy level'
   DO i = 1, nm, 1
      WRITE (UNIT=8, FMT=*) ei(i), cc(i)
   END DO
   DO i = 1, nm, 1
      WRITE (UNIT=9, FMT=*) 'DOS with the unit of States/eV'
      WRITE (UNIT=9, FMT=*) ei(i), IMAG(c4(i))
   END DO
   DEALLOCATE (c2)
   DEALLOCATE (c3)
   DEALLOCATE (c4)
   DEALLOCATE (c5)
   DEALLOCATE (cc)
   CLOSE (UNIT=8)
   CLOSE (UNIT=9)
   CLOSE (UNIT=10)
END IF
CALL MPI_FINALIZE(ierr)
DEALLOCATE (ks)
DEALLOCATE (tv)
DEALLOCATE (ol)
DEALLOCATE (ei)
DEALLOCATE (t2)
DEALLOCATE (t3)
STOP
END PROGRAM QUANTUMTRANSPORT_NEGF
