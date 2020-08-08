
!====================================================================
     ! SELECT_INVENTORY !!
!====================================================================

      SUBROUTINE select_inventory(A,PBO,MeasVals_in,Nz,Np,NumZ,NumP,
     &                          chosen_prev,prev_max,H,H_length,avg_rsa,
     &                          avg_rsa_frac,avg_other,hour,Nuc,POA,SOA,
     &                                                 chosen,Asquare )

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS--------------------------------------------

      integer info, lwork, Nz, Np, NumZ, NumP, jj, ii, H, H_length
      integer PBO(Np),MeasVals(Nz),MeasVals_in(Nz),i,j,k,x
      integer Nuc, POA, SOA
      double precision s(Np), v(Np,Np), u(Nz,Nz), hour
      double precision A(Nz,Np),Adummy(Nz,Np),Ainv(Np,Nz)
      double precision rga(Nz,Np,1),rowsum(Nz),rga_dummy(Nz,Np)
      double precision rga_other(Nz,Np,1)
      double precision max_rga(Np), rga_frac(Nz,Np,1)
      integer chosen(NumP),chosen_prev(NumP),INDX(1)
      double precision next_val, last_val, Asquare(NumP,NumP)
      double precision, allocatable :: work(:)
      double precision avg_rsa(Nz,Np,H_length),meanval(Nz,Np)
      double precision avg_rsa_frac(Nz,Np,H_length),meanval_frac(Nz,Np)
      double precision avg_other(Nz,Np,H_length),meanval_other(Nz,Np)
      double precision rga_remaining, max_meas,prev_max(Np),max_values
      double precision rowsum_frac(Nz)
      integer eps
      parameter( eps=1.d-4 )
      logical mask(Nz)

      character*90 string,fullname_soacond,fname1
      character*90 runname, datname, diagname, filename
      character*90 fullname_coag, fullname_cond, fullname_param
      character*90 fullname_NucRate, fullname_measured,fullname_pbio
      common /FLNAMES/ datname, diagname, fullname_coag, fullname_cond,
     &               fullname_param,fullname_NucRate,fullname_pbio,
     &               fullname_measured, fullname_soacond, fname1      
C-----CODE-------------------------------------------------------------

cdbg      print*, '#SEL_INV: previous inventories=',chosen_prev

      print*, 'SEL_INV: soa index=', SOA
       MeasVals=MeasVals_in
        info=1
        lwork=MAX(1,3*MIN(Nz,Np)+MAX(Nz,Np),5*MIN(Nz,Np) )
        allocate( work(lwork) ) 
        

        !! calculate "relative gain array" on sensitivity matrix
        do i=1,Np
         Adummy(:,i)=A(:,i)
        enddo
        call dgesvd('A','A',Nz,Np,Adummy,Nz,s,u,Nz,v,Np,
     &                work,lwork,info)
        Ainv(:,:) = 0.d0
        do j = 1,Np
         do k = 1,Nz
          do i = 1,Np
            Ainv(j,k)=Ainv(j,k)+v(i,j)/s(i)*u(k,i)
          enddo !sum
         enddo !columns
        enddo !rows
        rga(:,:,:) = 0.d0
        rowsum(:) = 0.d0
        do k = 1,Nz
         do j = 1,Np
          rga(k,j,1) = A(k,j)*Ainv(j,k)
          rowsum(k) = rowsum(k) + abs( rga(k,j,1) )
         enddo !columns
         rga_frac(k,:,1) = rga(k,:,1)/rowsum(k)
         rowsum_frac(k)=rowsum(k)/maxval(rga(k,:,1))
         rga_other(k,:,1)=rga(k,:,1)-rowsum(k)
         enddo !rows
         meanval(:,:)=0.d0
         meanval_frac(:,:)=0.d0
         meanval_other(:,:)=0.d0
cdbg         print*, '#SEL_INV: averaging window=',H_length
         if (H.lt.H_length) then
            avg_rsa_frac(:,:,H) = rga_frac(:,:,1)
            avg_rsa(:,:,H) = rga(:,:,1)
            avg_other(:,:,H) = rga_other(:,:,1)
         else
            do i=2,H_length
              avg_rsa_frac(:,:,i-1)=avg_rsa_frac(:,:,i)
              avg_rsa(:,:,i-1)=avg_rsa(:,:,i)
              avg_other(:,:,i-1)=avg_other(:,:,i)
              meanval_frac=meanval_frac+avg_rsa_frac(:,:,i-1)/H_length
              meanval=meanval+avg_rsa(:,:,i-1)/H_length
              meanval_other=meanval_other+avg_other(:,:,i-1)/H_length
            enddo
            avg_rsa_frac(:,:,H_length)=rga_frac(:,:,1)
            avg_rsa(:,:,H_length)=rga(:,:,1)
            avg_other(:,:,H_length)=rga_other(:,:,1)
            meanval_frac=meanval_frac+avg_rsa_frac(:,:,H_length)
     &                                       /H_length
            meanval=meanval+avg_rsa(:,:,H_length)/H_length
            meanval_other=meanval_other+avg_other(:,:,H_length)/H_length
         endif
         rga_dummy = meanval
         mask(:)=.false.
      If (NumP.gt.0) Then
      IF ( (NumZ .gt. NumP).and.(H.ge.H_length) ) THEN
         max_values=0.d0
         do k=1,NumP
           mask(:)=.false.
           if (k.eq.SOA) Then
              do j=13,Nz
                mask(j)=.true.
              enddo
              print*, 'MASK:', mask
           else
             do j = 1,Nz
               if (MeasVals(j).eq.1) mask(j)=.true.
             enddo
           endif
           INDX = maxloc( rga_dummy(:,k), mask )
           chosen(k) = INDX(1)
           max_rga(k) = rga_dummy( INDX(1),k )
           rga_dummy(INDX(1),k)=0.d0
           rga_remaining=0.d0
           x=1
c           Do x=1,Nz
           Do while (x.lt.Nz)
             INDX = maxloc( rga_dummy(:,k),mask )
             max_meas=maxval( rga_dummy(:,k),mask )
c             if (max_meas.gt.max_rga(k)-5.d-2) then
             !only choose if it was measured
c             if( (MeasVals(INDX(1)).eq.1).and.(MeasVals(chosen(k))
c     &                       .eq.1) ) then
              if ( ( max_meas +1.d-4 .ge. max_rga(k) ).and.( 
     &             meanval_other( INDX(1),k ).le.
     &                    meanval_other(chosen(k),k)) )Then
                   chosen(k) = INDX(1)
                   max_rga(k) = max_meas 
                   rga_dummy(chosen(k),k) = 0.d0
                   x=x+1
              else
                  x=Nz
                  MeasVals( chosen(k) ) =0
              endif
c             elseif (MeasVals(INDX(1)).eq.0) then
c                rga_dummy( INDX(1),k)=0.d0
c                x=x+1
c             else ! MeasVals( chosen(k) ) =0
c                x=x+1
c                rga_dummy( chosen(k),k)=0.d0
c                chosen(k)=INDX(1)
c             endif
           Enddo
c           rga_remaining=rga_remaining-max_rga(k)
c           max_values=max_values+max_rga(k)
c           rga_dummy(chosen(k),:) = 0.d0
           MeasVals( chosen(k) ) = 0
         enddo

         do k=1,NumP
c            max_rga(k)=max_rga(k)+1-meanval_other( chosen(k),k)
c            if ( abs(max_rga(k)-prev_max(k)).lt.eps ) then
c               chosen(k)=chosen_prev(k)
c            endif
            max_rga(k)=max_rga(k)+1-meanval_other( chosen(k),k)
         enddo
         prev_max = max_rga
      ELSE
          chosen = chosen_prev
cdbg          print*, '#SEL_INV: use previously chosen measurements--',
cdbg     &                   chosen
      ENDIF
      IF ( mod(hour,4.).ne.0.d0 ) then
         chosen=chosen_prev
      EndIf
         do j = 1,NumP
            jj = chosen(j)
            Asquare(j,:)=A(jj,:)
         enddo
c        chosen_prev = chosen
      
         print*,'#SEL_INV: inventories',chosen
       Endif
c      print*, '#SEL_INV: MaxRS vs otherRS=',max_rga,
c     &              'Inventories=',chosen


c             endif
c             rga_dummy( INDX(1),k ) = 0.d0
c               INDX = maxloc( rga_dummy(:,k) )
c             next_val = rga_dummy( INDX(1),k )
c             last_val = meanval( chosen(k),k )
c             if ( abs( next_val-last_val) .le. 5.d-2 ) then
c              if ( (INDX(1).gt chosen(k)).and.(MeasVals(INDX(1)).eq.1) )
c     &            then
c                  chosen(k) = INDX(1)
c                  max_rga(k)=max_rga(k)-last_val+next_val
c              else
c                  x=0
c              endif
c             endif
c            rga_remaining=rga_remaining+abs(meanval(x,k))

c           if (( max_rga(k)+eps .lt. rga_remaining ).and.
c     &         (MeasVals(chosen_prev(k)).eq.1) ) THEN
c             chosen(k)=chosen_prev(k)
c           endif

c       !! Calculate condition number of LHS 3x3 matrix
c         ! Need to calculate 1-norm of matrix
c         ANORM = 0.d0
c         dummy(:) = 0.d0
c         do n = 1,NumP
c         do m = 1,NumZ
c           dummy(n) = ABS( D(m,n) ) + dummy(n)
c          enddo
c         enddo
c         ANORM = maxval( dummy )
c         print*, '##PARAM EST 1-norm of matrix:', ANORM
        !  rcond = reciprocal of condition number
c        call dgecon('1', NumZ, D, NumP, ANORM, rcond, work, iwork, info)

        open(unit=51,file=fullname_pbio,status='old',access='append')
        write(string, '("(A13,",I3,"(E15.6,","),A3)")') Nz
        write(51,string) 'RS(i,:)=[',rowsum,'];'
        write(string, '("(A13,", I3, "(E15.6,","), A4)")') Nz
        write(51,*) 'Ginv{i}=[...'
        do i=1,Np
          write(51,string) ' ',Ainv(i,:),' ...'
          if (i .ne. Np) write(51,*) '; ...'
        enddo
        write(51,*) '];'
        write(51,*) 'G{i}=[...'
        do i=1,Np
          write(51,string) ' ',A(:,i),' ...'
          if (i .ne. Np) write(51,*) '; ...'
        enddo
        write(51,*) '];'
        if (H.ge.H_length) then
        write(51,*) 'RSA{i}=[...'
        do i=1,Np
          write(51,string) ' ',meanval_frac(:,i),' ...'
          if (i .ne. Np) write(51,*) '; ...'
        enddo
        write(51,*) '];'
        endif
        close(51)
        RETURN
        END
