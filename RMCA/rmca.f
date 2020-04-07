      program rmca
      implicit none
      character*4 version
      parameter (version='3.04')
*
************************************************************************
*
*     RMCA is a general purpose RMC/MC Monte Carlo code. It will
*       i)   fit to diffraction data
*       ii)  fit to EXAFS data
*       iii) use constraints on mean coordination numbers
*       iv)  use constraints on coordination number distributions
*       v)   use a pairwise potential
*
*     For details of use see the documentation
*
*     This code is written to allow dynamic allocation of space to 
*     arrays. This is done by having one large array in the main
*     programme unit, parts of which are allocated to various 
*     purposes as data is read in and as a result it becomes known 
*     how much space is needed. Thus instead of arrays we have integer
*     variables to point to the part of the single array to be used.
*     After data is read in the remainder of the programme is contained
*     in a subroutine which allows us to refer to our data by use of
*     array names rather than pointers. This does not comply with
*     the ANSI standard but will work on most systems.
*
*-----------------------------------------------------------------------
*
*     This code originally written by Malcolm Howe 
*     (MITCHELL@UK.AC.OXFORD.VAX until Sept 92)
*
*     Corrections to coordination constraint code .......... 7/12/92
*     Corrections to multiple data sets code ......(v3.03).. 7/12/92
*     Corrections to EXAFS code ...................(v3.04).. 9/12/92
*                              Jim Wicks (MITCHELL@UK.AC.OXFORD.VAX)
*
*     Last altered 31st January 1993
*
************************************************************************
*
*     The parameter array_size determine the total amount of array
*     space available to the programme. It may need to be adjusted
*     to suit different computer systems.
*
      integer       array_size
      parameter     (array_size=6900000)
      real          array(array_size),a(array_size)
      integer       integer_array(array_size),ia(array_size)
      logical       logical_array(array_size),la(array_size)
      equivalence   (array,integer_array,logical_array,a,ia,la)
      common        /array/ array
      integer       array_pointer
*
      integer       avcoordno,avcnew,atoms,cfnew,chisq,
     -              coeff_fq,coeff_eq,coeff_gr,coeff_sq,
     -              constant,coordfrac,coordno,delta,dhistogram,dnneigh,
     -              eexpt,exafs_edge,factor,fexpt,gexpt,gnorm,gpar,gtot,
     -              histogram,ncum,ni,nneigh,nx1,nx2,nxn,
     -              offset,potentials,q,ravcoord,rcoord,rcoord_0,
     -              rcut,renorm,sexpt,sigma,sigmaavc,sigmac,sqr,
     -              sum_f,sum_ff,sum_fs,sum_q,sum_s,sum_ss,
     -              too_close,typec,typecav,typec_0,typen,typenav,
     -              typen_0,work1,work2
*
      real          vectors(3,3)
      real          axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3,
     -              d,d1,d2,d3,dr,dr_0,eunits,f,rho,rho_0,
     -              temperature,timelim,timesav,volume,weight
      integer       dim_avcoord,dim_coord,dim_eq,dim_expt,dim_fq,dim_gr,
     -              dim_q,dim_r_exafs,dim_sq,dim_work1,dim_work2,
     -              iprint,lfile,n,nacc,navcoord,ncoll,ncoord,ncoord_0,
     -              neq,nexpt,nfq,ngen,ngr,npar,nq,nr,nr_exafs,nsaved,
     -              nsq,ntried,ntypes,qweight
      logical       calchis,moveout,truncated,usepot
      character*80  title
      character*20  file
      character*20  bidon
      character*24  tempo
      logical qex
      integer       ic,i,j
      data array_pointer /0/
      rho_0=0.
      dr_0=0.
*
*-------------------------------------------------------------------------
*     Obtain name of files to be used
*
      call rmc_datafile(file)
      lfile=len(file)
      do while (file(lfile:lfile).eq.' ')
         lfile=lfile-1
      enddo
      write (*,*)
      write (*,*) '======================'
      write (*,*) 'RMCA version ',version
      write (*,*) '======================'
      write (*,*)
      write (*,*) 'Using files :',file(1:lfile)
      write (*,*)
*
*-------------------------------------------------------------------------
*     If possible read the intermediate (histogram) file
*
C      call rmca_read_his(file,lfile,array,integer_array,
C     - array_size,array_pointer,
C     - atoms,calchis,dr_0,histogram,n,nacc,ncoord_0,ngen,ni,nneigh,nr,
C     - nsaved,ntried,ntypes,rcoord_0,rho_0,title,truncated,typec_0,
C     - typen_0,vectors)
*
*-------------------------------------------------------------------------
*     If this has failed read the configuration file instead
*
200   calchis=.true.
      if (calchis) then
         typec_0=0
         typen_0=0
         rcoord_0=0
         ncoord_0=0
         call rmca_read_cfg(file,lfile,array,integer_array,
     -    array_size,array_pointer,atoms,
     -    n,nacc,ngen,ni,nsaved,ntried,ntypes,title,truncated,vectors)
      endif
      npar=ntypes*(ntypes+1)/2
*
*-------------------------------------------------------------------------
*     If cell is truncated stop as this is no longer supported
*   
      if (truncated) then
         write (*,*) 'Truncated cells are not supported'
         stop
      endif
*
*-------------------------------------------------------------------------
*     Open programme data file and read in general RMC parameters
*
      call rmca_allocate(rcut,npar,array_pointer,array_size)
      call rmca_allocate(delta,ntypes,array_pointer,array_size)
      tempo=file(1:lfile)//'.dat'
      print *,tempo
      call open_read1(3,tempo)
C      call open_read(3,file,'.dat')
      read (3,1000) title
1000  format(80A1)
      read (3,*) rho
      read (3,*) (array(rcut+ic),ic=1,npar)
      read (3,*) (array(delta+i),i=1,ntypes)
      read (3,*) dr
      read (3,*) moveout
      read (3,*) ncoll
      read (3,*) iprint
      read (3,*) timelim,timesav
*
*-------------------------------------------------------------------------
*     We must recalculate the histogram if 
*      i)  we are using the moveout option (so that we can determine 
*          which particles are too close together)
*      ii) we have changed the density or the r spacing
*
      if (moveout) then
         write (*,*) 
     -  'MOVEOUT option specified. Histogram will be recalculated'
         calchis=.true.
      endif
      if (.not.calchis) then
         if (rho.ne.rho_0) calchis=.true.
         if (dr.ne.dr_0) calchis=.true.
      endif
*
*-------------------------------------------------------------------------
*     Find volume of configuration cell
*
      volume=8.*abs(vectors(1,1)*vectors(2,2)*vectors(3,3)
     -             +vectors(2,1)*vectors(3,2)*vectors(1,3)
     -             +vectors(3,1)*vectors(1,2)*vectors(2,3)
     -             -vectors(3,1)*vectors(2,2)*vectors(1,3)
     -             -vectors(2,1)*vectors(1,2)*vectors(3,3)
     -             -vectors(1,1)*vectors(3,2)*vectors(2,3))
*
*-------------------------------------------------------------------------
*     If necessary rescale configuration cell to required density
*     and calculate the number of histogram points
*
      if (calchis) then
         f=(real(n)/(rho*volume))**(1./3.)
         do j=1,3
            do i=1,3
               vectors(i,j)=vectors(i,j)*f
            enddo
         enddo
         volume=real(n)/rho
         axb1=vectors(2,1)*vectors(3,2)-vectors(3,1)*vectors(2,2)
         axb2=vectors(3,1)*vectors(1,2)-vectors(1,1)*vectors(3,2)
         axb3=vectors(1,1)*vectors(2,2)-vectors(2,1)*vectors(1,2)
         bxc1=vectors(2,2)*vectors(3,3)-vectors(3,2)*vectors(2,3)
         bxc2=vectors(3,2)*vectors(1,3)-vectors(1,2)*vectors(3,3)
         bxc3=vectors(1,2)*vectors(2,3)-vectors(2,2)*vectors(1,3)
         cxa1=vectors(2,3)*vectors(3,1)-vectors(3,3)*vectors(2,1)
         cxa2=vectors(3,3)*vectors(1,1)-vectors(1,3)*vectors(3,1)
         cxa3=vectors(1,3)*vectors(2,1)-vectors(2,3)*vectors(1,1)
         d1=1./sqrt(axb1**2+axb2**2+axb3**2)
         d2=1./sqrt(bxc1**2+bxc2**2+bxc3**2)
         d3=1./sqrt(cxa1**2+cxa2**2+cxa3**2)
         d=volume/8.*min(d1,d2,d3)
         nr=int(d/dr)
         call rmca_allocate(histogram,nr*npar,array_pointer,
     -    array_size)
      endif
*
*-------------------------------------------------------------------------
*     Read all remaining data
*
      call rmca_read_dat(array,integer_array,logical_array,
     - array_size,array_pointer,
     - avcoordno,coeff_eq,coeff_fq,coeff_gr,coeff_sq,coordfrac,coordno,
     - dr,eexpt,
     - eunits,exafs_edge,fexpt,gexpt,navcoord,ncoord,neq,nexpt,
     - nfq,ngr,ni,npar,nq,nr,nr_exafs,ntypes,nx1,nx2,nsq,offset,
     - potentials,q,qweight,ravcoord,rcoord,renorm,
     - sexpt,sigma,sigmaavc,sigmac,temperature,
     - typec,typecav,typen,typenav,usepot,weight)
*
*-------------------------------------------------------------------------
*
*     We must recalculate the histogram if the coordination
*     constraints have been changed
*
      if (.not.calchis) then
         if (ncoord.ne.ncoord_0.and.ncoord.gt.0) then
            calchis=.true.
         else
            do i=1,ncoord
               if (integer_array(typec_0+i).ne.integer_array(typec+i)) 
     -          calchis=.true.
               if (integer_array(typen_0+i).ne.integer_array(typen+i)) 
     -          calchis=.true.
               if (array(rcoord_0+2*(i-1)+1).ne.array(rcoord+2*(i-1)+1)) 
     -          calchis=.true.
               if (array(rcoord_0+2*(i-1)+2).ne.array(rcoord+2*(i-1)+2)) 
     -          calchis=.true.
            enddo
         endif
      endif
*
*-------------------------------------------------------------------------
*     Allocate all other required array space
*
      call rmca_allocate(dhistogram,nr*npar,array_pointer,
     - array_size)
      call rmca_allocate(gnorm,nr,array_pointer,array_size)
      call rmca_allocate(gpar,nr*npar,array_pointer,array_size)
      call rmca_allocate(ncum,ntypes,array_pointer,array_size)
      call rmca_allocate(nxn,npar,array_pointer,array_size)
      too_close=0
      if (moveout) call rmca_allocate(too_close,n,array_pointer,
     - array_size)
      if (ncoord.ne.ncoord_0.and.ncoord.gt.0) then
         call rmca_allocate(nneigh,n*ncoord,array_pointer,array_size)
      endif
      call rmca_allocate(sqr,nr*nq,array_pointer,array_size)
      call rmca_allocate(gtot,nr*ngr,array_pointer,array_size)
      call rmca_allocate(chisq,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_s,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_ss,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_f,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_ff,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_fs,nexpt,array_pointer,array_size)
      call rmca_allocate(sum_q,nexpt,array_pointer,array_size)
      call rmca_allocate(constant,nexpt,array_pointer,array_size)
      call rmca_allocate(factor,nexpt,array_pointer,array_size)
      call rmca_allocate(cfnew,ncoord,array_pointer,array_size)
      call rmca_allocate(dnneigh,n*ncoord,array_pointer,
     - array_size)
      call rmca_allocate(avcnew,navcoord,array_pointer,array_size)
*
      dim_work1=max(npar*nr,npar*nq,n)
      call rmca_allocate(work1,dim_work1,array_pointer,array_size)
      dim_work2=max(nq,n)
      call rmca_allocate(work2,dim_work2,array_pointer,array_size)
*
      write (*,100) 100.*real(array_pointer)/real(array_size)
100   format(/' Fraction of array space used is ',f6.2,'%'/)
*
*-------------------------------------------------------------------------
*     Now call subroutine containing the remaining code. This allows
*     us to refer to the arrays by name rather than by pointer.
*
      dim_avcoord=max(1,navcoord)
      dim_coord=max(1,ncoord)
      dim_eq=max(1,neq)
      dim_expt=max(1,nexpt)
      dim_fq=max(1,nfq)
      dim_gr=max(1,ngr)
      dim_q=max(1,nq)
      dim_r_exafs=max(1,nr_exafs)
      dim_sq=max(1,nsq)
      call rmca_main_code(calchis,dr,eunits,file,iprint,moveout,n,
     - nacc,navcoord,ngen,ncoll,ncoord,neq,nexpt,nfq,ngr,npar,
     - nq,nr,nr_exafs,nsaved,nsq,ntried,ntypes,qweight,rho,
     - temperature,timelim,timesav,title,usepot,vectors,version,volume,
     - weight,dim_avcoord,dim_coord,dim_eq,dim_expt,dim_fq,dim_gr,
     - dim_q,dim_r_exafs,dim_sq,dim_work1,dim_work2,
     - a(atoms+1),a(avcnew+1),a(avcoordno+1),a(cfnew+1),a(chisq+1),
     - a(coeff_eq+1),a(coeff_fq+1),a(coeff_gr+1),a(coeff_sq+1),
     - a(constant+1),a(coordfrac+1),a(delta+1),a(eexpt+1),a(factor+1),
     - a(fexpt+1),a(gexpt+1),a(gnorm+1),a(gpar+1),a(gtot+1),
     - a(potentials+1),a(q+1),a(ravcoord+1),a(rcoord+1),
     - a(rcut+1),a(sexpt+1),a(sigma+1),a(sigmaavc+1),a(sigmac+1),
     - a(sqr+1),a(sum_f+1),a(sum_ff+1),a(sum_fs+1),a(sum_s+1),
     - a(sum_ss+1),a(sum_q+1),a(work1+1),a(work2+1),
     - ia(coordno+1),ia(dhistogram+1),ia(dnneigh+1),ia(exafs_edge+1),
     - ia(histogram+1),ia(ncum+1),ia(ni+1),ia(nneigh+1),
     - ia(nx1+1),ia(nx2+1),ia(nxn+1),ia(too_close+1),ia(typec+1),
     - ia(typecav+1),ia(typen+1),ia(typenav+1),
     - la(offset+1),la(renorm+1))
      close(10,status='delete')
      stop
      end




      subroutine rmca_allocate(array_name,size,array_pointer,array_size)
      implicit none
      integer       array_name,size,array_pointer,array_size
*
      array_name=array_pointer
      array_pointer=array_pointer+size
      if (array_pointer.gt.array_size) then
         write (*,*) 'Fatal error: RMCA array space exhausted'
         stop
      endif
      return
      end




      subroutine rmca_main_code(
     - calchis,dr,eunits,file,iprint,moveout,n,nacc,navcoord,ngen,
     - ncoll,ncoord,neq,nexpt,
     - nfq,ngr,npar,nq,nr,nr_exafs,nsaved,nsq,ntried,
     - ntypes,qweight,rho,temperature,timelim,timesav,
     - title,usepot,vectors,version,volume,weight,
     - dim_avcoord,dim_coord,dim_eq,dim_expt,dim_fq,dim_gr,dim_q,
     - dim_r_exafs,dim_sq,dim_work1,dim_work2,
     - atoms,avcnew,avcoordno,cfnew,chisq,coeff_eq,coeff_fq,coeff_gr,
     - coeff_sq,constant,
     - coordfrac,delta,eexpt,factor,fexpt,gexpt,gnorm,gpar,
     - gtot,potentials,q,ravcoord,rcoord,rcut,sexpt,
     - sigma,sigmaavc,sigmac,sqr,sum_f,sum_ff,sum_fs,sum_s,sum_ss,sum_q,
     - work1,work2,
     - coordno,dhistogram,dnneigh,exafs_edge,histogram,ncum,ni,nneigh,
     - nx1,nx2,nxn,too_close,
     - typec,typecav,typen,typenav,offset,renorm)
*
      implicit none
      integer       dim_avcoord,dim_coord,dim_eq,dim_expt,dim_fq,
     -              dim_gr,dim_q,dim_r_exafs,dim_sq,dim_work1,dim_work2
      integer       n,npar,nr,ntypes,lfile
      real          atoms(3,n),avcoordno(dim_avcoord),
     -              avcnew(dim_avcoord),cfnew(dim_coord),
     -              chisq(dim_expt),
     -              coeff_eq(dim_r_exafs,dim_q,ntypes,dim_eq),
     -              coeff_fq(dim_q,dim_fq,npar),coeff_gr(dim_gr,npar),
     -              coeff_sq(dim_sq,npar),constant(dim_expt),
     -              coordfrac(dim_coord),delta(ntypes),
     -              eexpt(dim_q,dim_eq),factor(dim_expt),
     -              fexpt(dim_q,dim_fq),gnorm(nr),gexpt(nr,dim_gr),
     -              gpar(nr,npar),gtot(nr,dim_gr),metric(3,3),move(6),
     -              potentials(nr,npar),q(dim_q),
     -              ravcoord(2,dim_avcoord),rcoord(2,dim_coord),
     -              rcut(npar),
     -              sexpt(dim_q,dim_sq),sigma(dim_expt),
     -              sigmaavc(dim_avcoord),sigmac(dim_coord),
     -              sqr(nr,dim_q),sum_f(dim_expt),sum_ff(dim_expt),
     -              sum_fs(dim_expt),sum_q(dim_expt),sum_s(dim_expt),
     -              sum_ss(dim_expt),vectors(3,3),
     -              work1(dim_work1),work2(dim_work2)
      integer       coordno(dim_coord),dhistogram(nr,npar),
     -              dnneigh(n,dim_coord),exafs_edge(dim_eq),
     -              histogram(nr,npar),ncum(ntypes),ni(ntypes),
     -              nneigh(n,dim_coord),
     -              nx1(dim_expt),nx2(dim_expt),nxn(npar),
     -              too_close(n),
     -              typec(dim_coord),typecav(dim_avcoord),
     -              typen(dim_coord),typenav(dim_avcoord)
      logical       offset(dim_expt),renorm(dim_expt)
      real          ran
      real          chisq0,chisq1,denergy,dr,energy,eunits,pi,prat,
     -              q2,rho,save_time,temperature,time_of_loop,
     -              timer_on,time_used,time_start,timelim,timesav,
     -              volume,weight
      integer       imove,iprint,iseed,lastprint,move_type,nacc,
     -              navcoord,ncoll,ncoord,neq,nexpt,nfq,ngen,ngen0,ngr,
     -              nq,nr_exafs,nsaved,nsq,ntoo_close,ntried,qweight
      logical       acceptable,calchis,first,looping,moveout,
     -              rerun,usepot
      character*(*) version
      character*80  title
      character*20  file
      character*24  tempo
      logical qex
      integer       i,i1,ic,icc,ieq,iexpt,ifq,igr,ir,ir1,ir2,isq,
     -              itype,iq,j,j1,jtype,k
      parameter     (pi=3.1415926)
      data          rerun /.true./
*
*-------------------------------------------------------------------------
*     Initialisation
*
      call rmc_init(iseed,time_start)
*
*-------------------------------------------------------------------------
*     List programme parameters to standard output
*
      call rmca_list(avcoordno,coeff_fq,coeff_gr,coeff_sq,
     - constant,coordno,coordfrac,delta,dim_fq,dim_gr,dim_q,dim_sq,
     - dr,eunits,exafs_edge,iprint,
     - moveout,n,navcoord,ncoll,ncoord,neq,nexpt,nfq,ngr,ni,npar,nr,
     - nr_exafs,nsq,ntypes,nx1,nx2,offset,qweight,ravcoord,rcoord,rcut,
     - renorm,rho,sigma,sigmaavc,sigmac,temperature,timelim,timesav,
     - title,typec,typecav,typen,typenav,usepot,vectors,weight)
*
*-------------------------------------------------------------------------
*     Calculate initial data sums
*
      do iexpt=1,nexpt
         sum_f(iexpt)=0.
         sum_ff(iexpt)=0.
         sum_q(iexpt)=0.
      enddo
      do igr=1,ngr
         iexpt=igr
         do ir=nx1(iexpt),nx2(iexpt)
            sum_f(iexpt)=sum_f(iexpt)+gexpt(ir,igr)
            sum_ff(iexpt)=sum_ff(iexpt)+gexpt(ir,igr)**2
            sum_q(iexpt)=sum_q(iexpt)+1.
         enddo
      enddo
      do isq=1,nsq
         iexpt=ngr+isq
         do iq=nx1(iexpt),nx2(iexpt)
            sum_f(iexpt)=sum_f(iexpt)+sexpt(iq,isq)
            sum_ff(iexpt)=sum_ff(iexpt)+sexpt(iq,isq)**2
            sum_q(iexpt)=sum_q(iexpt)+1.
         enddo
      enddo
      do ifq=1,nfq
         iexpt=ngr+nsq+ifq
         do iq=nx1(iexpt),nx2(iexpt)
            sum_f(iexpt)=sum_f(iexpt)+fexpt(iq,ifq)
            sum_ff(iexpt)=sum_ff(iexpt)+fexpt(iq,ifq)**2
            sum_q(iexpt)=sum_q(iexpt)+1.
         enddo
      enddo
      do ieq=1,neq
         iexpt=ngr+nsq+nfq+ieq
         do iq=nx1(iexpt),nx2(iexpt)
            q2=(q(iq)/q(nq))**qweight
            sum_f(iexpt)=sum_f(iexpt)+eexpt(iq,ieq)*q2
            sum_ff(iexpt)=sum_ff(iexpt)+eexpt(iq,ieq)**2*q2*q2
            sum_q(iexpt)=sum_q(iexpt)+1
         enddo
      enddo
*
*-------------------------------------------------------------------------
*     Calculate metric matrix
*
      do j=1,3
         do i=1,3
            metric(i,j)=0.
            do k=1,3
               metric(i,j)=metric(i,j)+vectors(k,i)*vectors(k,j)
            enddo
         enddo
      enddo
*
*-------------------------------------------------------------------------
*     Calculate arrays needed to convert histograms to g(r)
*
      ncum(1)=ni(1)
      do itype=2,ntypes
         ncum(itype)=ncum(itype-1)+ni(itype)
      enddo
      ic=1
      do itype=1,ntypes
         do jtype=itype,ntypes
            nxn(ic)=ni(itype)*ni(jtype)
            if (itype.ne.jtype) nxn(ic)=nxn(ic)*2
            ic=ic+1
         enddo
      enddo
      do ir=1,nr
         gnorm(ir)=(3.*ir*ir+0.25)*dr*dr*dr*2.*pi/(3.*volume)
      enddo
*
*-------------------------------------------------------------------------
*     If necessary calculate partial g(r) histograms and 
*     coordination numbers
*
      if (calchis) then
         write (*,*)
         write (*,*) 'Calculating g(r) histogram'
         call rmca_calhist(atoms,dr,histogram,metric,moveout,n,ncoord,
     -    ncum,ni,nneigh,npar,nr,ntoo_close,ntypes,rcoord,rcut,
     -    too_close,typec,typen,work1)
C         write (*,*)
C         write (*,*) 'Writing g(r) histogram to disk'
C         write (*,*) 
C     -     '*** DO NOT INTERRUPT UNTIL SAVING IS COMPLETED ***'
C         call rmca_write_his(atoms,dr,file,histogram,n,nacc,ncoord,
C     -    ngen,ni,nneigh,ntried,nr,nsaved,ntypes,rcoord,rho,title,
C     -    typec,typen,vectors)
C         write (*,*) 'Saving has been completed'
      endif
*
*-------------------------------------------------------------------------
*     Check whether the configuration satisfies the cut-offs
*
      do itype=1,ntypes
         do jtype=itype,ntypes
            ic=(itype-1)*(2*ntypes-itype)/2+jtype
            do ir=1,nint(rcut(ic)/dr)-1
               if (histogram(ir,ic).gt.0) 
     -          write (*,130) histogram(ir,ic),itype,jtype,real(ir)*dr
            enddo
         enddo
      enddo
*
*-------------------------------------------------------------------------
*     Calculate initial energy
*
      if (usepot) then
         energy=0.
         do ic=1,npar
            do ir=1,nr
               energy=energy+real(histogram(ir,ic))*potentials(ir,ic)
            enddo
         enddo
      endif
*
*-------------------------------------------------------------------------
*     Calculate the sin(Qr) table
*
      do iq=1,nq
         do ir=1,nr
            sqr(ir,iq)=4.*pi*rho*real(ir)*dr*
     -       sin(real(ir)*dr*q(iq))/q(iq)*dr
         enddo
      enddo
*
*-------------------------------------------------------------------------
*-------------------------------------------------------------------------
*     This is where the loop begins. 
*
      save_time=min(timelim,timesav)
      time_of_loop=0.
      call rmc_time(time_start,timer_on)
      timer_on=abs(timer_on)
      ngen0=ngen
      lastprint=ngen-mod(ngen,iprint)
      first=.true.
      acceptable=.true.
      looping=.true.
      ntried=ntried-1
      do while (looping)
*
*-------------------------------------------------------------------------
*        Generate a random move, (which will be zero first time
*        through the loop).
*
         if (first) then
            imove=1
            move(4)=0.
            move(5)=0.
            move(6)=0.
         else
10          continue
            if (moveout.and.ran(iseed).lt.0.1) then
               imove=int(ran(iseed)*ntoo_close)+1
               imove=too_close(imove)
            else
               imove=int(ran(iseed)*n)+1
            endif
            move_type=1
            do while (imove.gt.ncum(move_type))
               move_type=move_type+1
            enddo
            if (delta(move_type).le.0.) goto 10
            move(4)=2.*(ran(iseed)-0.5)*delta(move_type)/
     -       sqrt(metric(1,1))
            move(5)=2.*(ran(iseed)-0.5)*delta(move_type)/
     -       sqrt(metric(2,2))
            move(6)=2.*(ran(iseed)-0.5)*delta(move_type)/
     -       sqrt(metric(3,3))
         endif
         move(1)=atoms(1,imove)
         move(2)=atoms(2,imove)
         move(3)=atoms(3,imove)
         move(4)=move(1)+move(4)+3.
         move(5)=move(2)+move(5)+3.
         move(6)=move(3)+move(6)+3.
         move(4)=2.*(move(4)/2.-int(move(4)/2.))-1.
         move(5)=2.*(move(5)/2.-int(move(5)/2.))-1.
         move(6)=2.*(move(6)/2.-int(move(6)/2.))-1.
         ngen=mod(ngen+1,1000000)
         acceptable=.true.
*
*-------------------------------------------------------------------------
*        Calculate change in histogram 
*
         call rmca_changehist(acceptable,atoms,dnneigh,dhistogram,
     -    dr,first,imove,metric,move,move_type,n,ncoord,ncum,ni,npar,
     -    nr,ntypes,rcoord,rcut,typec,typen,work1,work2)
*
         if (acceptable) then
*
*----------------------------------------------------------------------
*           Calculate energy change
*
            ntried=mod(ntried+1,1000000)
            if (usepot) then
               denergy=0.
               do ic=1,npar
                  do ir=1,nr
                     denergy=denergy+
     -                real(dhistogram(ir,ic))*potentials(ir,ic)
                  enddo
               enddo
            endif
*
*----------------------------------------------------------------------
*           Calculate new partial and total g(r)'s 
* 
            if (nexpt.gt.0) then
               do ic=1,npar
                  do ir=1,nr
                     gpar(ir,ic)=real(histogram(ir,ic)+
     -                dhistogram(ir,ic))/
     -                (gnorm(ir)*real(nxn(ic)))
                  enddo
               enddo
               do igr=1,ngr
                  do ir=1,nr
                     gtot(ir,igr)=0.
                  enddo
               enddo
               do igr=1,ngr
                  do ic=1,npar
                     do ir=1,nr
                        gtot(ir,igr)=gtot(ir,igr)+coeff_gr(igr,ic)*
     -                   (gpar(ir,ic)-1.)
                     enddo
                  enddo
               enddo
            endif
*
*----------------------------------------------------------------------
*           Calculate sums for chi-squared
* 
            call rmca_sums(coeff_eq,coeff_fq,coeff_sq,dhistogram,
     -       dim_fq,dim_r_exafs,dim_q,dim_sq,exafs_edge,
     -       eexpt,fexpt,gexpt,gpar,gtot,histogram,neq,nexpt,nfq,ngr,
     -       npar,nq,nr,nr_exafs,ntypes,nx1,nx2,nsq,q,qweight,sexpt,
     -       sqr,sum_fs,sum_s,sum_ss,work1,work2) 
*
*----------------------------------------------------------------------
*           Calculate chi-squared for experimental data
*
            call rmca_chisq(nexpt,nx1,nx2,renorm,offset,
     -       factor,constant,sum_f,sum_ff,sum_fs,sum_s,sum_ss,
     -       sum_q,sigma,chisq,chisq1)
*
*----------------------------------------------------------------------
*           Add on terms for coordination constraints
*
            do icc=1,ncoord
               cfnew(icc)=0
               do i=ncum(typec(icc))-ni(typec(icc))+1,ncum(typec(icc))
                  if (nneigh(i,icc)+dnneigh(i,icc).eq.coordno(icc)) 
     -             cfnew(icc)=cfnew(icc)+1.
               enddo
               cfnew(icc)=cfnew(icc)/real(ni(typec(icc)))
               chisq1=chisq1+(cfnew(icc)-coordfrac(icc))**2/
     -          sigmac(icc)**2
            enddo
            do icc=1,navcoord
               avcnew(icc)=0.
               ir1=max(1,nint(ravcoord(1,icc)/dr))
               ir2=min(nr,nint(ravcoord(2,icc)/dr))        
               i1=min(typecav(icc),typenav(icc))
               j1=max(typecav(icc),typenav(icc))
               ic=(i1-1)*(2*ntypes-i1)/2+j1
               do ir=ir1,ir2
                  avcnew(icc)=avcnew(icc)+
     -             histogram(ir,ic)+dhistogram(ir,ic)
               enddo
               avcnew(icc)=avcnew(icc)/real(ni(typecav(icc)))
               if (typecav(icc).eq.typenav(icc)) 
     -          avcnew(icc)=avcnew(icc)*2
               chisq1=chisq1+
     -          (avcnew(icc)-avcoordno(icc))**2/sigmaavc(icc)**2
            enddo
*
*----------------------------------------------------------------------
*           Decide whether to accept move
*
            if (first) then
               acceptable=.false.
               chisq0=chisq1
            else
               prat=(chisq0-chisq1)/2.-denergy*weight
               if (prat.lt.-20) then
                  acceptable=.false.
               elseif (prat.lt.0) then
                  prat=exp(prat)
                  if (prat.lt.ran(iseed)) acceptable=.false.
               endif
            endif
            if (acceptable) nacc=mod(nacc+1,1000000)
         endif
*
*----------------------------------------------------------------------
*        Write summary
*
         if (first.or.(abs(ngen-lastprint).ge.iprint.and.acceptable)) 
     -    then
            call rmca_summary(avcnew,avcoordno,cfnew,chisq,chisq1,
     -       constant,coordfrac,(energy+denergy)/real(n),factor,nacc,
     -       navcoord,ncoord,nexpt,ngen,ntried,nx1,nx2,usepot)
            lastprint=ngen-mod(ngen,iprint)
         endif
*
         if (acceptable) then
*
*----------------------------------------------------------------------
*           Make move
*
            atoms(1,imove)=move(4)
            atoms(2,imove)=move(5)
            atoms(3,imove)=move(6)
            do ic=1,npar
               do ir=1,nr
                  histogram(ir,ic)=histogram(ir,ic)+dhistogram(ir,ic)
               enddo
            enddo
            do icc=1,ncoord
               do i=1,n
                  nneigh(i,icc)=nneigh(i,icc)+dnneigh(i,icc)
               enddo
            enddo
            chisq0=chisq1
            energy=energy+denergy

*----------------------------------------------------------------------
*           if collecting then save every nth good configuration
*
            if (ncoll.gt.0.and.nsaved.lt.ncoll.and.nacc.ge.n) then
               write (*,*)
               write (*,*) 'Saving configuration to disk'
               write (*,*)
     -          '*** DO NOT INTERRUPT UNTIL SAVING IS COMPLETED ***'
*
*
      lfile=len(file)
      do while (file(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=file(1:lfile)//'.sav'
      inquire(file=tempo,exist=qex)
      if(qex.eqv..FALSE.) go to 200
      call filedel(3,tempo)
200   call open_write1(3,tempo)
C               call open_write(3,file,'.sav')
               call write_config(title,ngen,ntried,nacc,nsaved,
     -          .false.,vectors,n,ntypes,ni,atoms,3)
               close(3)
               nsaved=nsaved+1
               write (*,*) 'Saving has been completed'
               write (*,*)
               write (*,'(a,i3,a)') 'Configuration number ',nsaved,
     -          ' has been saved'
               write (*,*)
               nacc=0
               ngen=0
               ntried=0
               if (nsaved.ge.ncoll) rerun=.false.
            endif
         endif
*
*----------------------------------------------------------------------
*        Check time
*
         call rmc_time(time_start,time_used)
         if (time_used.lt.0.or.time_used.ge.save_time) then
*
*----------------------------------------------------------------------
*           Calculate and display timing details
*
            time_of_loop=time_of_loop+(abs(time_used)-timer_on)
            write (*,*)
            write (*,*) '  Total time in loop so far :',
     -       time_of_loop*60,'s'
            write (*,*) '  Time per generated move is:',
     -       time_of_loop*60000/real(max(ngen-ngen0,3)),'ms'
*
*----------------------------------------------------------------------
*           Write out the configuration, the g(r) histogram
*           and the output file
*
            write (*,*)
            write (*,*) 'Saving configuration and results to disk'
            write (*,*) 
     -       '*** DO NOT INTERRUPT UNTIL SAVING IS COMPLETED ***'
*
*
      lfile=len(file)
      do while (file(lfile:lfile).eq.' ')
         lfile=lfile-1
      enddo
      tempo=file(1:lfile)//'.cfg'
      inquire(file=tempo,exist=qex)
      if(qex.eqv..FALSE.) go to 201
      call filedel(21,tempo)
201   call open_write1(3,tempo)
C            call open_write(3,file,'.cfg')
            call write_config(title,ngen,ntried,nacc,nsaved,
     -       .false.,vectors,n,ntypes,ni,atoms,3)
            close(3)
*
C            call rmca_write_his(atoms,dr,file,histogram,n,nacc,ncoord,
C     -       ngen,ni,nneigh,ntried,nr,nsaved,ntypes,rcoord,rho,title,
C     -       typec,typen,vectors)
*
            call rmca_output(coeff_eq,coeff_fq,coeff_gr,coeff_sq,
     -       constant,dim_fq,dim_gr,dim_q,dim_r_exafs,dim_sq,dr,
     -       eexpt,exafs_edge,factor,fexpt,file,gexpt,gnorm,
     -       gpar,gtot,histogram,neq,nexpt,nfq,ngr,npar,nq,nr,nr_exafs,
     -       nsq,ntypes,nx1,nx2,nxn,q,sexpt,sqr,title,version,work1,
     -       work2,qweight)
            write (*,*) 'Saving has been completed'
            write (*,*)
            save_time=save_time+timesav
            call rmc_time(time_start,timer_on)
         endif
*
*----------------------------------------------------------------------
*        Decide whether to continue programme
*
         if ((time_used.lt.0.).or.(time_used.ge.timelim)
     -    .or..not.rerun) looping=.false.
         first=.false.
      enddo
*
      call rmc_rerun(rerun)
      return
130   format(' WARNING: ',i3,' pairs cause g',2i1,' to be ',
     - 'non-zero at ',f4.2,'A')
      end




      subroutine rmca_read_his(file,lfile,array,integer_array,
     - array_size,array_pointer,
     - atoms,calchis,dr_0,histogram,n,nacc,ncoord_0,ngen,ni,nneigh,nr,
     - nsaved,ntried,ntypes,rcoord_0,rho_0,title,truncated,typec_0,
     - typen_0,vectors)
      implicit none
*
      integer       array_size,array_pointer
      real          array(array_size)
      integer       integer_array(array_size)
      integer       atoms,histogram,ni,nneigh,rcoord_0,typec_0,typen_0
*
      real          vectors(3,3)
      character*8   hex(9)
      real          dr_0,hex_to_real,rho_0
      integer       lfile,n,nacc,ncoord_0,ngen,npar,nr,nsaved,ntried,
     -              ntypes
      logical       calchis,error,exists,truncated
      character*80  title
      character*20  file
      character*24  tempo
      logical qex
      integer       i,j
*
C      inquire(file=file(1:lfile)//'.his',exist=exists)
C      if (.not.exists) goto 200
      go to 200
*
      write (*,*) 'Reading intermediate (histogram) file ',
     - file(1:lfile),'.his'
      tempo=file(1:lfile)//'.his'
      call open_read1(3,tempo)
      print *,tempo
C      call open_read(3,file,'.his')
*
      read (3,*,end=100,err=100) 
      read (3,1000,end=100,err=100) title
1000  format(1x,A80)
      read (3,*,end=100,err=100) ngen,ntried,nacc,nsaved,truncated
      read (3,'(1x,9a8)',end=100,err=100) (hex(i),i=1,9)
      do i=1,3
         do j=1,3
            vectors(j,i)=hex_to_real(hex(3*(i-1)+j),error)
            if (error) goto 100
         enddo
      enddo
      read (3,'(1x,a8,1x,2i10)',end=100,err=100) hex(1),n,ntypes
      rho_0=hex_to_real(hex(1),error)
      if (error) goto 100
*
      call rmca_allocate(atoms,3*n,array_pointer,array_size)
      call rmca_allocate(ni,ntypes,array_pointer,array_size)
*
      read (3,*,end=100,err=100) (integer_array(ni+i),i=1,ntypes)
      do i=1,n
         read (3,'(3(1x,a8))',end=100,err=100) (hex(j),j=1,3)
         do j=1,3
            array(atoms+3*(i-1)+j)=hex_to_real(hex(j),error)
            if (error) goto 100
         enddo
      enddo
      read (3,'(i10,1x,a8)',end=100,err=100) nr,hex(1)
      dr_0=hex_to_real(hex(1),error)
      if (error) goto 100
*
      npar=ntypes*(ntypes+1)/2
      call rmca_allocate(histogram,nr*npar,array_pointer,
     - array_size)
      do i=1,npar
         read (3,*,end=100,err=100) 
     -   (integer_array(histogram+nr*(i-1)+j),j=1,nr)
      enddo
*
      read (3,*,end=100,err=100) ncoord_0
      call rmca_allocate(typec_0,ncoord_0,array_pointer,
     - array_size)
      call rmca_allocate(typen_0,ncoord_0,array_pointer,
     - array_size)
      call rmca_allocate(rcoord_0,2*ncoord_0,array_pointer,
     - array_size)
      call rmca_allocate(nneigh,n*ncoord_0,array_pointer,
     - array_size)
      do i=1,ncoord_0
         read (3,'(2i10,1x,a8,1x,a8)',end=100,err=100) 
     -    integer_array(typec_0+i),integer_array(typen_0+i),
     -    hex(1),hex(2)
         array(rcoord_0+2*(i-1)+1)=hex_to_real(hex(1),error)
         if (error) goto 100
         array(rcoord_0+2*(i-1)+2)=hex_to_real(hex(2),error)
         if (error) goto 100
         read (3,*,end=100,err=100) 
     -    (integer_array(nneigh+(i-1)*n+j),j=1,n)
      enddo
      calchis=.false.
      close(3)
      return
100   write (*,*) 'ERROR READING FILE'
      write (*,*) 'Configuration file will be read instead'
      write (*,*)
200   calchis=.true.
      close(3)
      return
      end




      subroutine rmca_read_cfg(file,lfile,array,integer_array,
     - array_size,array_pointer,atoms,n,nacc,ngen,ni,nsaved,
     - ntried,ntypes,title,truncated,vectors)
      implicit none
*
      integer       array_size,array_pointer
      real          array(array_size)
      integer       integer_array(array_size)
      integer       atoms,ni
*
      real          vectors(3,3)
      integer       lfile,n,nacc,ngen,nsaved,ntried,ntypes
      logical       truncated,space
      character*80  title
      character*20  file
      character*24  tempo
      logical       qex
      integer       i,itype,j
*
*     Reads a configuration data file
*
      write (*,*) 'Reading configuration file ',file(1:lfile),'.cfg'
      tempo=file(1:lfile)//'.cfg'
      print *,tempo
      call open_read1(3,tempo)
C      call open_read(3,file,'.cfg')
*
*     Check format of file
*     Also have to allow for possibility of initial space
*     (Fortran carriage control) which may or may not occur at
*     beginning of each line
*
C      read (3,'(a)') title
      read (3,1000) title
      print 1000,title
1000  format(A80)
      space=.false.
      if (title(1:1).eq.' ') space=.true.
      if (space) title=title(2:len(title))
      if (title(1:10).ne.'(Version 3') then
         write (*,*) 'Configuration file has wrong format'
         stop
      endif
*
C      read (3,'(a)') title
      read (3,1000) title
      if (space) title=title(2:len(title))
      read (3,*) ngen,ntried,nacc
      read (3,*) nsaved
      read (3,*) n
      read (3,*) ntypes
      read (3,*)
      read (3,*)
      read (3,*)
      read (3,*) truncated
      read (3,*)
      read (3,*) vectors
*
      call rmca_allocate(ni,ntypes,array_pointer,array_size)
      do itype=1,ntypes
         read (3,*) integer_array(ni+itype)
         read (3,*) 
         read (3,*) 
      enddo
*
      call rmca_allocate(atoms,3*n,array_pointer,array_size)
      do i=1,n
         read (3,*) (array(atoms+3*(i-1)+j),j=1,3)
      enddo
*
      close(3)
      return
      end





      subroutine rmca_read_dat(array,integer_array,logical_array,
     - array_size,array_pointer,
     - avcoordno,coeff_eq,coeff_fq,coeff_gr,coeff_sq,coordfrac,coordno,
     - dr,eexpt,eunits,exafs_edge,fexpt,
     - gexpt,navcoord,ncoord,neq,nexpt,nfq,ngr,ni,npar,
     - nq,nr,nr_exafs,ntypes,nx1,nx2,nsq,offset,
     - potentials,q,qweight,ravcoord,rcoord,renorm,
     - sexpt,sigma,sigmaavc,sigmac,temperature,
     - typec,typecav,typen,typenav,usepot,weight)
*
      implicit none
      integer       array_size,array_pointer,lfile
      real          array(array_size)
      integer       integer_array(array_size)
      logical       logical_array(array_size)
*
      integer       avcoordno,coeff_eq,coeff_fq,coeff_gr,coeff_sq,
     -              coordfrac,coordno,eexpt,exafs_edge,fexpt,
     -              gexpt,navcoord,ncoord,neq,nexpt,nfq,ngr,ni,npar,nq,
     -              nr,nr_exafs,nx1,nx2,nsq,ntypes,offset,potentials,q,
     -              ravcoord,rcoord,renorm,sexpt,sigma,sigmaavc,sigmac,
     -              typec,typecav,typen,typenav,qweight
*
      real          boltzmannk,c,dr,eunits,r,rmax,temperature,weight
      integer       nq_co,nr_co,nrval,ntypes_co,nq1,nq2,nqsoq
      logical       usepot
      character*20  sqfile
      character*24  tempo
      logical       qex
      integer       i,ic,ieq,iexpt,ifq,igr,iq,ir,isq,itype
      parameter     (boltzmannk=1.380662e-23)
*
*     Reads all the programme parameters and any other necessary data
*     for the programme.
*
*     Read numbers of different experimental data constraints
*
      nq=0
      read (3,*) ngr,nsq,nfq,neq
      nexpt=ngr+nsq+nfq+neq
      call rmca_allocate(sigma,nexpt,array_pointer,array_size)
      call rmca_allocate(renorm,nexpt,array_pointer,array_size)
      call rmca_allocate(offset,nexpt,array_pointer,array_size)
      call rmca_allocate(nx1,nexpt,array_pointer,array_size)
      call rmca_allocate(nx2,nexpt,array_pointer,array_size)
      call rmca_allocate(gexpt,nr*ngr,array_pointer,array_size)
      call rmca_allocate(coeff_gr,ngr*npar,array_pointer,array_size)
      call rmca_allocate(coeff_sq,nsq*npar,array_pointer,array_size)
      call rmca_allocate(exafs_edge,neq,array_pointer,array_size)
      q=0
      sexpt=0
      fexpt=0
      eexpt=0
      coeff_fq=0
      coeff_eq=0
*
*     Next parameters for g(r) data constraints
*
      if (ngr.ne.0) then
         do igr=1,ngr
            iexpt=igr
            read (3,2000) sqfile
2000        format(A20)
            read (3,*) integer_array(nx1+iexpt),
     -       integer_array(nx2+iexpt)
            read (3,*) c
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C            call open_read(4,sqfile,'.dat')
            read (4,*) nrval
            integer_array(nx2+iexpt)=
     -       min(integer_array(nx2+iexpt),nrval,nr)
            read (4,*) 
            do ir=1,integer_array(nx2+iexpt)
               read (4,*) r,array(gexpt+(igr-1)*nr+ir)
               array(gexpt+(igr-1)*nr+ir)=
     -          array(gexpt+(igr-1)*nr+ir)-c
            enddo
            close(4)
            read (3,*) (array(coeff_gr+(ic-1)*ngr+igr),ic=1,npar)
            read (3,*) array(sigma+iexpt)
            read (3,*) logical_array(renorm+iexpt)
            logical_array(offset+iexpt)=.false.
         enddo
      endif
*
*     Next parameters for S(Q) data constraints
*
      if (nsq.ne.0) then
         do isq=1,nsq
            iexpt=ngr+isq
            read (3,2000) sqfile
            read (3,*) integer_array(nx1+iexpt),
     -       integer_array(nx2+iexpt)
            read (3,*) c
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C            call open_read(4,sqfile,'.dat')
            read (4,*) nq
            integer_array(nx2+iexpt)=
     -       min(integer_array(nx2+iexpt),nq)
            read (4,*) 
            if (isq.eq.1) then
               call rmca_allocate(q,nq,array_pointer,
     -          array_size)
               call rmca_allocate(sexpt,nq*nsq,array_pointer,
     -          array_size)
            endif
            do iq=1,nq
               read (4,*) array(q+iq),array(sexpt+(isq-1)*nq+iq)
               array(sexpt+(isq-1)*nq+iq)=
     -          array(sexpt+(isq-1)*nq+iq)-c
            enddo
            close(4)
            read (3,*) (array(coeff_sq+(ic-1)*nsq+isq),ic=1,npar)
            read (3,*) array(sigma+iexpt)
            read (3,*) logical_array(renorm+iexpt)
            read (3,*) logical_array(offset+iexpt)
         enddo
      endif
*
*     Next parameters for F(Q) data constraints (for the purposes of
*     this programme the difference between F(Q) and S(Q) is that 
*     F(Q) has Q dependent coefficients of the partial structure factors)
*
      if (nfq.ne.0) then
         do ifq=1,nfq
            iexpt=ngr+nsq+ifq
            read (3,2000) sqfile
            read (3,*) integer_array(nx1+iexpt),
     -       integer_array(nx2+iexpt)
            read (3,*) c
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C            call open_read(4,sqfile,'.dat')
            read (4,*) nq
            integer_array(nx2+iexpt)=
     -       min(integer_array(nx2+iexpt),nq)
            read (4,*) 
            if (ifq.eq.1) then
               if (nsq.eq.0) call rmca_allocate(q,nq,array_pointer,
     -          array_size)
               call rmca_allocate(coeff_fq,nq*nfq*npar,array_pointer,
     -          array_size)
               call rmca_allocate(fexpt,nq*nfq,array_pointer,
     -          array_size)
            endif
            do iq=1,nq
               read (4,*) array(q+iq),array(fexpt+(ifq-1)*nq+iq),
     -          (array(coeff_fq+(ic-1)*nfq*nq+(ifq-1)*nq+iq),ic=1,npar)
               array(fexpt+(ifq-1)*nq+iq)=array(fexpt+(ifq-1)*nq+iq)-c
            enddo
            close(4)
            read (3,*) array(sigma+iexpt)
            read (3,*) logical_array(renorm+iexpt)
            read (3,*) logical_array(offset+iexpt)
         enddo
      endif
*
*     Next parameters for EXAFS (here called E(Q)) data constraints
*
      nr_exafs=0
      if (neq.ne.0) then
         read (3,*) rmax
         read (3,*) qweight
         nr_exafs=min(nr,int(rmax/dr))
         do ieq=1,neq
            iexpt=ngr+nsq+nfq+ieq
            read (3,2000) sqfile
            read (3,*) nq1,nq2
            read (3,*) integer_array(nx1+iexpt),
     -       integer_array(nx2+iexpt)
            read (3,*) integer_array(exafs_edge+ieq)
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C            call open_read(4,sqfile,'.dat')
            read (4,*) nqsoq
            if (nqsoq.ne.nq2-nq1+1)
     -       stop 'EXAFS S(Q) has incorrect no of pts'
            nq=max(nq,nqsoq)
            read (4,*) 
            if (ieq.eq.1) then
               if (nsq+nfq.eq.0) call rmca_allocate(q,nq,array_pointer,
     -          array_size)
               call rmca_allocate(coeff_eq,nq*nr_exafs*ntypes*neq,
     -          array_pointer,array_size)
               call rmca_allocate(eexpt,nq*neq,array_pointer,
     -          array_size)
            endif
            do iq=nq1,nq2
               read (4,*) array(q+iq),array(eexpt+(ieq-1)*nq+iq)
            enddo
            close(4)
            read (3,2000) sqfile
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C            call open_read(4,sqfile,'.dat')
            read (4,*) ntypes_co,nq_co,nr_co
            if (ntypes_co.ne.ntypes.or.nq_co.ne.nqsoq.or.
     -       nr_co.ne.nr_exafs) then
               write (*,*) 'EXAFS coefficients incorrect'
               stop
            endif
            do itype=1,ntypes
               read (4,*)
               do iq=nq1,nq2
                  read (4,*)
                  read (4,*) (array(coeff_eq+(ieq-1)*ntypes*nq*nr_exafs
     -             +(itype-1)*nq*nr_exafs+(iq-1)*nr_exafs+ir),
     -             ir=1,nr_exafs)
                  do ir=1,nr_exafs
                     array(coeff_eq+(ieq-1)*ntypes*nq*nr_exafs
     -                +(itype-1)*nq*nr_exafs+(iq-1)*nr_exafs+ir)=
     -                array(coeff_eq+(ieq-1)*ntypes*nq*nr_exafs
     -                +(itype-1)*nq*nr_exafs+(iq-1)*nr_exafs+ir)/
     -                real(integer_array(ni+
     -                integer_array(exafs_edge+ieq)))
                  enddo
               enddo
            enddo
            close(4)
            read (3,*) array(sigma+iexpt)
            read (3,*) logical_array(renorm+iexpt)
            read (3,*) logical_array(offset+iexpt)
         enddo
      endif
*
*     Coordination constraints
*
      read (3,*) ncoord
      call rmca_allocate(typec,ncoord,array_pointer,array_size)
      call rmca_allocate(typen,ncoord,array_pointer,array_size)
      call rmca_allocate(rcoord,2*ncoord,array_pointer,array_size)
      call rmca_allocate(coordno,ncoord,array_pointer,array_size)
      call rmca_allocate(coordfrac,ncoord,array_pointer,
     - array_size)
      call rmca_allocate(sigmac,ncoord,array_pointer,array_size)
      do i=1,ncoord
         read (3,*) integer_array(typec+i),
     -    integer_array(typen+i),
     -    array(rcoord+2*(i-1)+1),array(rcoord+2*(i-1)+2),
     -    integer_array(coordno+i),
     -    array(coordfrac+i),array(sigmac+i)
      enddo
*
*     Mean coordination number constraints
*
      read (3,*) navcoord
      call rmca_allocate(typecav,navcoord,array_pointer,array_size)
      call rmca_allocate(typenav,navcoord,array_pointer,array_size)
      call rmca_allocate(ravcoord,2*navcoord,array_pointer,array_size)
      call rmca_allocate(avcoordno,navcoord,array_pointer,array_size)
      call rmca_allocate(sigmaavc,navcoord,array_pointer,array_size)
      do i=1,navcoord
         read (3,*) integer_array(typecav+i),integer_array(typenav+i),
     -    array(ravcoord+2*(i-1)+1),array(ravcoord+2*(i-1)+2),
     -    array(avcoordno+i),array(sigmaavc+i)
      enddo
*
*     Potentials
*
      read (3,*) usepot
      if (usepot) then
         read (3,*) temperature
         read (3,*) eunits
         read (3,*) weight
         read (3,2000) sqfile
      lfile=len(sqfile)
      do while (sqfile(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=sqfile(1:lfile)//'.dat'
      print *,tempo
      call open_read1(4,tempo)
C     call open_read(4,sqfile,'.dat')
         read (4,*) nrval
         read (4,*)  
         call rmca_allocate(potentials,nr*npar,array_pointer,
     -    array_size)
         do ir=1,min(nr,nrval)
            read (4,*) r,(array(potentials+(ic-1)*nr+ir),
     -       ic=1,npar)
            do ic=1,npar
               array(potentials+(ic-1)*nr+ir)=
     -          array(potentials+(ic-1)*nr+ir)
     -          *eunits/boltzmannk/temperature       
            enddo
         enddo
         close(4)
         do ir=nrval+1,nr
            do ic=1,npar
               array(potentials+(ic-1)*nr+ir)=0.
           enddo
         enddo
      else
         potentials=0
      endif
      close(3)
      return
      end







      subroutine rmca_calhist(atoms,dr,histogram,metric,moveout,n,
     - ncoord,ncum,ni,nneigh,npar,nr,ntoo_close,ntypes,rcoord,rcut,
     - too_close,typec,typen,work1)
      implicit none
      integer       n,nr
      real          atoms(3,*),metric(3,3),rcoord(2,*),rcut(*),work1(*)
      integer       histogram(nr,*),ncum(*),ni(*),
     -              nneigh(n,*),too_close(*),
     -              typec(*),typen(*)
      real          dr,x,y,z,d
      integer       ncoord,npar,ntoo_close,ntypes
      logical       moveout
      integer       i,ic,icc,ig,ir,j,type1,type2
*
      do ic=1,npar
         do ir=1,nr
            histogram(ir,ic)=0
         enddo
      enddo
      do icc=1,ncoord
         do i=1,n
            nneigh(i,icc)=0
         enddo
      enddo
      if (moveout) then
         do i=1,n
            too_close(i)=0
         enddo
      endif
      type2=1
      do i=1,n
         if (i.gt.ncum(type2)) type2=type2+1
         do type1=1,type2
            ic=(type1-1)*(2*ntypes-type1)/2+type2
            do j=ncum(type1)-ni(type1)+1,min(ncum(type1),i-1)
               x=atoms(1,i)-atoms(1,j)+3.
               y=atoms(2,i)-atoms(2,j)+3.
               z=atoms(3,i)-atoms(3,j)+3.
               x=2.*(x/2.-int(x/2.))-1.
               y=2.*(y/2.-int(y/2.))-1.
               z=2.*(z/2.-int(z/2.))-1.
               d=metric(1,1)*x*x+metric(2,2)*y*y+metric(3,3)*z*z
     -          +2.*(metric(1,2)*x*y+metric(1,3)*x*z+metric(2,3)*y*z)
               work1(j)=sqrt(d)
            enddo
            do j=ncum(type1)-ni(type1)+1,min(ncum(type1),i-1)
               ig=nint(work1(j)/dr)
               if (ig.le.nr) histogram(ig,ic)=histogram(ig,ic)+1
               if (moveout.and.work1(j).lt.rcut(ic)) then
                  too_close(i)=1
                  too_close(j)=1
               endif
            enddo
         enddo
         do icc=1,ncoord
            do type1=1,type2
               if (type1.eq.typec(icc).and.type2.eq.typen(icc)) then
                  do j=ncum(type1)-ni(type1)+1,min(ncum(type1),i-1)
                     if (work1(j).ge.rcoord(1,icc).and.
     -                work1(j).lt.rcoord(2,icc)) 
     -                nneigh(j,icc)=nneigh(j,icc)+1
                  enddo
               endif
               if (type2.eq.typec(icc).and.type1.eq.typen(icc)) then
                  do j=ncum(type1)-ni(type1)+1,min(ncum(type1),i-1)
                     if (work1(j).ge.rcoord(1,icc).and.
     -                work1(j).lt.rcoord(2,icc)) 
     -                nneigh(i,icc)=nneigh(i,icc)+1
                  enddo
               endif
            enddo
         enddo
      enddo
      if (moveout) then
         ntoo_close=0
         do i=1,n
            if (too_close(i).eq.1) then
               ntoo_close=ntoo_close+1
               too_close(ntoo_close)=i
            endif
         enddo
         if (ntoo_close.eq.0) moveout=.false.
      endif
      return
      end








      subroutine rmca_list(avcoordno,coeff_fq,coeff_gr,coeff_sq,
     - constant,coordno,coordfrac,delta,dim_fq,dim_gr,dim_q,dim_sq,
     - dr,eunits,exafs_edge,iprint,
     - moveout,n,navcoord,ncoll,ncoord,neq,nexpt,nfq,ngr,ni,npar,nr,
     - nr_exafs,nsq,ntypes,nx1,nx2,offset,qweight,ravcoord,rcoord,rcut,
     - renorm,rho,sigma,sigmaavc,sigmac,temperature,timelim,timesav,
     - title,typec,typecav,typen,typenav,usepot,vectors,weight)
      implicit none
      integer       dim_fq,dim_gr,dim_q,dim_sq
      real          avcoordno(*),coeff_fq(dim_q,dim_fq,*),
     -              coeff_gr(dim_gr,*),coeff_sq(dim_sq,*),
     -              constant(*),coordfrac(*),delta(*),
     -              ravcoord(2,*),rcoord(2,*),rcut(*),
     -              sigma(*),sigmaavc(*),sigmac(*),
     -              vectors(3,3)
      integer       coordno(*),exafs_edge(*),ni(*),nx1(*),nx2(*),
     -              typec(*),typecav(*),typen(*),typenav(*)
      logical       offset(*),renorm(*)
      real          dr,eunits,rho,temperature,timelim,timesav,weight
      integer       iprint,n,navcoord,ncoll,ncoord,neq,nexpt,nfq,
     -              ngr,npar,nr,nr_exafs,nsq,ntypes,qweight
      logical       moveout,usepot
      character*80  title
      integer       i,ic,ieq,iexpt,ifq,igr,isq,j
*
      write (*,100) title,n,ntypes
      if (ntypes.gt.1) then
         write (*,110)
         do i=0,(ntypes-1)/5
            write (*,120) (ni(j),j=i*5+1,min(ntypes,i*5+5))
         enddo
      endif
      write (*,130) vectors
*
      write (*,140) rho
      do i=1,ntypes
         do j=i,ntypes
             write (*,150) i,j,rcut((i-1)*(2*ntypes-i)/2+j)
         enddo
      enddo
      if (moveout) write (*,160)
      write (*,*)
      do i=1,ntypes
         write (*,170) delta(i),i
      enddo
      write (*,180) nr,dr,nr*dr,ncoll,iprint,timelim,timesav
*
      if (nexpt.eq.0) write (*,190)
      if (ngr.gt.0) then
         write (*,200) ngr,'g(r)''s:'
         do igr=1,ngr
            iexpt=igr
            write (*,210) igr,nx1(iexpt),nx2(iexpt),
     -       constant(iexpt),sigma(iexpt)
            if (renorm(iexpt)) write (*,220)
            write (*,240) (coeff_gr(igr,ic),ic=1,npar)
         enddo
      endif
      if (nsq.gt.0) then
         write (*,200) nsq,
     -    'S(Q)''s with constant coefficients'
         do isq=1,nsq
            iexpt=ngr+isq
            write (*,210) isq,nx1(iexpt),nx2(iexpt),
     -       constant(iexpt),sigma(iexpt)
            if (renorm(iexpt)) write (*,220)
            if (offset(iexpt)) write (*,230)
            write (*,240) (coeff_sq(isq,ic),ic=1,npar)
         enddo
      endif
      if (nfq.gt.0) then
         write (*,200) nfq,
     -    'S(Q)''s with Q dependent coefficients'
         do ifq=1,nfq
            iexpt=ngr+nsq+ifq
            write (*,210) ifq,nx1(iexpt),nx2(iexpt),
     -       constant(iexpt),sigma(iexpt)
            if (renorm(iexpt)) write (*,220)
            if (offset(iexpt)) write (*,230)
            write (*,250) (coeff_fq(1,ifq,ic),ic=1,npar)
         enddo
      endif
      if (neq.gt.0) then
         write (*,200) neq,'EXAFS spectra'
         write (*,260) nr_exafs*dr
         write (*,270) qweight
         do ieq=1,neq
            iexpt=ngr+nsq+nfq+ieq
            write (*,210) ieq,nx1(iexpt),nx2(iexpt),
     -       constant(iexpt),sigma(iexpt)
            write (*,280) exafs_edge(ieq)
            if (renorm(iexpt)) write (*,220)
            if (offset(iexpt)) write (*,230)
         enddo
      endif
*
      write (*,290) ncoord
      do i=1,ncoord
         write (*,300) typen(i),rcoord(1,i),rcoord(2,i),typec(i),
     -    coordfrac(i)*100.,coordno(i),sigmac(i)
      enddo

*
      write (*,310) navcoord
      do i=1,navcoord
         write (*,320) typenav(i),ravcoord(1,i),ravcoord(2,i),
     -    typecav(i),avcoordno(i),sigmaavc(i)
      enddo
*
      if (usepot) write (*,330) eunits,temperature,weight
     
      return
100   format(/' Title of run: ',a65//
     -        ' Configuration contains ',i5,' atoms of ',i1,' types')
110   format( ' Numbers of each type are:')
120   format(11x,5i10)
130   format(/' Configuration cell vectors are:',3(/10x,3f9.4))
140   format(/' Number density is ',f6.4,' atoms per cubic Angstrom'///
     -        ' Cut offs in partial g(r)''s are at:')
150   format(10x,i1,'-',i1,5x,f6.3,' A')
160   format(/' Particles too close to others will be preferentially ',
     - 'moved'/)
170   format(' Maximum change in any coordinate is ',f6.4,
     - 'A for particles of type ',i2)
180   format(/' Using ',i4,' r points spaced at ',f5.3,
     -          'A up to ',f6.3,'A'///
     -        ' No. of configs to save    : ',i3///
     -        ' Writing summary every ',i4,' generated moves'/
     -        ' Programme will run for ',f10.1,' events saving every ',
     -          f10.1,' events')
190   format(//' No experimental data constraints')
200   format(//' Fitting to ',i2,1x,a)
210   format(/' Data set ',i1,' uses points ',i4,' to ',i4,
     -          ' after subtraction of ',f6.3/
     -        ' Standard deviation is ',f7.4)
220   format(' Data will be renormalised')
230   format(' Data will be rezeroed')
240   format(' Coefficients of partials are:'/1000(11x,8(f7.5,1x)/))
250   format(' Coefficients of partials at minimum Q value are:'/
     -            1000(11x,8(f7.5,1x)/))
260   format(' Maximum r value for EXAFS calculation is ',f5.2)
270   format(' EXAFS fit will be weighted by Q**',i1)
280   format(' EXAFS edge is for particle type ',i2)
290   format(//' There are ',i2,' coordination constraints:')
300   format('    For atoms of type ',i2,' between ',f6.2,'A and ',
     -               f6.2,'A of atoms of type ',i2/
     -       '         constraint is ',f5.1,'% ',i2,
     -                '-fold coordination with sigma=',g10.4) 
310   format(/' There are ',i2,' coordination number constraints:')
320   format('    For atoms of type ',i2,' between ',f6.2,'A and ',
     -               f6.2,'A of atoms of type ',i2/
     -       '     constraint on average coordination number is ',
     -               f6.2,' with sigma=',g10.4) 
330   format(/' Using a potential read from file'/
     -        ' Energy units are             : ',e12.6/
     -        ' Temperature is               : ',f6.2,'K'/
     -        ' Potential weighting factor is: ',f7.4/)
      end






      subroutine rmca_changehist(acceptable,atoms,dnneigh,dhistogram,
     - dr,first,imove,metric,move,move_type,n,ncoord,ncum,ni,npar,nr,
     - ntypes,rcoord,rcut,typec,typen,work1,work2)
      implicit none
      integer       n,nr
      real          atoms(3,*),metric(3,3),move(6),rcoord(2,*),rcut(*),
     -              work1(*),work2(*)
      integer       dhistogram(nr,*),dnneigh(n,*),ncum(*),ni(*),
     -              typec(*),typen(*)
      real          dnew,dold,dr,x,y,z
      integer       imove,move_type,ncoord,npar,ntypes
      logical       acceptable,first,inthen,innow,vectorisable
      integer       i,i1,ic,icc,ig,ir,itype,j1
      parameter     (vectorisable=.false.)
*
*     Calculates the change in the histogram and check whether the move
*     is acceptable according to the cut-off restrictions
*
      do ic=1,npar
         do ir=1,nr
            dhistogram(ir,ic)=0
         enddo
      enddo
      do icc=1,ncoord
         do i=1,n
            dnneigh(i,icc)=0
         enddo
      enddo
*
      if (first) return
*
      if (vectorisable) then
         do itype=1,ntypes
            i1=min(itype,move_type)
            j1=max(itype,move_type)
            ic=(i1-1)*(2*ntypes-i1)/2+j1
            do i=ncum(itype)-ni(itype)+1,ncum(itype)
               x=atoms(1,i)-move(1)+3.
               y=atoms(2,i)-move(2)+3.
               z=atoms(3,i)-move(3)+3.
               x=2.*(x/2.-int(x/2.))-1.
               y=2.*(y/2.-int(y/2.))-1.
               z=2.*(z/2.-int(z/2.))-1.
               dold=metric(1,1)*x*x+metric(2,2)*y*y+metric(3,3)*z*z
     -          +2.*(metric(1,2)*x*y+metric(1,3)*x*z+metric(2,3)*y*z)
               work1(i)=sqrt(dold)
               x=atoms(1,i)-move(4)+3.
               y=atoms(2,i)-move(5)+3.
               z=atoms(3,i)-move(6)+3.
               x=2.*(x/2.-int(x/2.))-1.
               y=2.*(y/2.-int(y/2.))-1.
               z=2.*(z/2.-int(z/2.))-1.
               dnew=metric(1,1)*x*x+metric(2,2)*y*y+metric(3,3)*z*z
     -          +2.*(metric(1,2)*x*y+metric(1,3)*x*z+metric(2,3)*y*z)
               work2(i)=sqrt(dnew)
            enddo
            do i=ncum(itype)-ni(itype)+1,ncum(itype)
               if (work2(i).lt.rcut(ic).and.(work2(i).lt.work1(i))) 
     -          then
                   acceptable=.false.
                   return
               endif
            enddo
         enddo
      else
         do itype=1,ntypes
            i1=min(itype,move_type)
            j1=max(itype,move_type)
            ic=(i1-1)*(2*ntypes-i1)/2+j1
            do i=ncum(itype)-ni(itype)+1,ncum(itype)
               x=atoms(1,i)-move(1)+3.
               y=atoms(2,i)-move(2)+3.
               z=atoms(3,i)-move(3)+3.
               x=2.*(x/2.-int(x/2.))-1.
               y=2.*(y/2.-int(y/2.))-1.
               z=2.*(z/2.-int(z/2.))-1.
               dold=metric(1,1)*x*x+metric(2,2)*y*y+metric(3,3)*z*z
     -          +2.*(metric(1,2)*x*y+metric(1,3)*x*z+metric(2,3)*y*z)
               work1(i)=sqrt(dold)
               x=atoms(1,i)-move(4)+3.
               y=atoms(2,i)-move(5)+3.
               z=atoms(3,i)-move(6)+3.
               x=2.*(x/2.-int(x/2.))-1.
               y=2.*(y/2.-int(y/2.))-1.
               z=2.*(z/2.-int(z/2.))-1.
               dnew=metric(1,1)*x*x+metric(2,2)*y*y+metric(3,3)*z*z
     -          +2.*(metric(1,2)*x*y+metric(1,3)*x*z+metric(2,3)*y*z)
               work2(i)=sqrt(dnew)
               if (work2(i).lt.rcut(ic).and.(work2(i).lt.work1(i))) 
     -          then
                   acceptable=.false.
                   return
               endif
            enddo
         enddo
      endif
*
      itype=1
      do i=1,n
        if (imove.ne.i.and.acceptable) then
          if (i.gt.ncum(itype)) itype=itype+1
          i1=min(itype,move_type)
          j1=max(itype,move_type)
          ic=(i1-1)*(2*ntypes-i1)/2+j1
          ig=nint(work1(i)/dr)
          if (ig.le.nr) dhistogram(ig,ic)=dhistogram(ig,ic)-1
          ig=nint(work2(i)/dr)
          if (ig.le.nr) dhistogram(ig,ic)=dhistogram(ig,ic)+1
          do icc=1,ncoord
            innow=work2(i).ge.rcoord(1,icc).and.
     -          work2(i).lt.rcoord(2,icc)
            inthen=work1(i).ge.rcoord(1,icc).and.
     -          work1(i).lt.rcoord(2,icc)
            if (move_type.eq.typec(icc).and.itype.eq.typen(icc)) then
              if (innow.and..not.inthen) 
     -          dnneigh(imove,icc)=dnneigh(imove,icc)+1
              if (inthen.and..not.innow) 
     -          dnneigh(imove,icc)=dnneigh(imove,icc)-1
            endif
            if (move_type.eq.typen(icc).and.itype.eq.typec(icc)) then
              if (innow.and..not.inthen) 
     -          dnneigh(i,icc)=dnneigh(i,icc)+1
              if (inthen.and..not.innow) 
     -          dnneigh(i,icc)=dnneigh(i,icc)-1
            endif
          enddo
        endif
      enddo
      return
      end




      subroutine rmca_sums(coeff_eq,coeff_fq,coeff_sq,dhistogram,
     - dim_fq,dim_r_exafs,dim_q,dim_sq,
     - exafs_edge,eexpt,fexpt,gexpt,gpar,gtot,histogram,neq,nexpt,nfq,
     - ngr,npar,nq,nr,nr_exafs,ntypes,nx1,nx2,nsq,q,qweight,sexpt,sqr,
     - sum_fs,sum_s,sum_ss,work1,work2)
      implicit none
      integer       dim_fq,dim_r_exafs,dim_q,dim_sq,nr,ntypes
      real          coeff_eq(dim_r_exafs,dim_q,ntypes,*),
     -              coeff_fq(dim_q,dim_fq,*),coeff_sq(dim_sq,*),
     -              eexpt(dim_q,*),fexpt(dim_q,*),
     -              gexpt(nr,*),gtot(nr,*),gpar(nr,*),
     -              q(*),sexpt(dim_q,*),sqr(nr,*),
     -              sum_fs(*),sum_s(*),sum_ss(*),work1(*),work2(*)
      integer       dhistogram(nr,*),exafs_edge(*),histogram(nr,*),
     -              nx1(*),nx2(*)
      real          q2
      integer       neq,nexpt,nfq,ngr,npar,nr_exafs,nq,nsq,qweight
      integer       i,i1,ic,ieq,iexpt,ifq,igr,iq,ir,isq,itype,j1
*
      do iexpt=1,nexpt
         sum_s(iexpt)=0.
         sum_ss(iexpt)=0.
         sum_fs(iexpt)=0.
      enddo
*
*     First g(r)'s
*
      do igr=1,ngr
         iexpt=igr
         do ir=nx1(iexpt),nx2(iexpt)
            sum_s(iexpt)=sum_s(iexpt)+gtot(ir,igr)
            sum_ss(iexpt)=sum_ss(iexpt)+gtot(ir,igr)*gtot(ir,igr)
            sum_fs(iexpt)=sum_fs(iexpt)+gexpt(ir,igr)*gtot(ir,igr)
         enddo
      enddo
*
*     Next S(Q)'s and F(Q)'s. If we have F(Q)'s we must calculate the
*     partial structure factors. However if we don't and there are fewer
*     S(Q)'s than partials it is quicker to calculate total g(r)'s and
*     transform them than to calculate partials and add them.
*
      if (nsq.lt.npar.and.nfq.eq.0) then
         do isq=1,nsq
            iexpt=ngr+isq
            do i=1,nr*nsq
               work1(i)=0.
            enddo
            do ic=1,npar
               do ir=1,nr
                  work1((isq-1)*nr+ir)=work1((isq-1)*nr+ir)+
     -             coeff_sq(isq,ic)*(gpar(ir,ic)-1.)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               work2(iq)=0.
               do ir=1,nr
                  work2(iq)=work2(iq)+sqr(ir,iq)*work1((isq-1)*nr+ir)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               sum_s(iexpt)=sum_s(iexpt)+work2(iq)
               sum_ss(iexpt)=sum_ss(iexpt)+work2(iq)**2
               sum_fs(iexpt)=sum_fs(iexpt)+sexpt(iq,isq)*work2(iq)
            enddo
         enddo
      else
*
*        Calculate partial S(Q)'s
*
         do i=1,nq*npar
            work1(i)=0.
         enddo
         do iq=1,nq
            do ic=1,npar
               do ir=1,nr
                  work1((ic-1)*nq+iq)=work1((ic-1)*nq+iq)+
     -             (gpar(ir,ic)-1.)*sqr(ir,iq)
               enddo
            enddo
         enddo
*
         do isq=1,nsq
            iexpt=ngr+isq
            do iq=nx1(iexpt),nx2(iexpt)
               work2(iq)=0.
            enddo
            do ic=1,npar
               do iq=nx1(iexpt),nx2(iexpt)
                  work2(iq)=work2(iq)+
     -             coeff_sq(isq,ic)*work1((ic-1)*nq+iq)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               sum_s(iexpt)=sum_s(iexpt)+work2(iq)
               sum_ss(iexpt)=sum_ss(iexpt)+work2(iq)**2
               sum_fs(iexpt)=sum_fs(iexpt)+sexpt(iq,isq)*work2(iq)
            enddo
         enddo
*
         do ifq=1,nfq
            iexpt=ngr+nsq+ifq
            do iq=nx1(iexpt),nx2(iexpt)
               work2(iq)=0.
            enddo
            do ic=1,npar
               do iq=nx1(iexpt),nx2(iexpt)
                  work2(iq)=work2(iq)+
     -             coeff_fq(iq,ifq,ic)*work1((ic-1)*nq+iq)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               sum_s(iexpt)=sum_s(iexpt)+work2(iq)
               sum_ss(iexpt)=sum_ss(iexpt)+work2(iq)**2
               sum_fs(iexpt)=sum_fs(iexpt)+fexpt(iq,ifq)*work2(iq)
            enddo
         enddo
      endif
*
*     Finally EXAFS spectra
*
      do ieq=1,neq
         iexpt=ngr+nsq+nfq+ieq
         do iq=nx1(iexpt),nx2(iexpt)
            work2(iq)=0.
            do itype=1,ntypes
               i1=min(itype,exafs_edge(ieq))
               j1=max(itype,exafs_edge(ieq))
               ic=(i1-1)*(2*ntypes-i1)/2+j1
               do ir=1,nr_exafs
                  work2(iq)=work2(iq)+
     -             real(histogram(ir,ic)+dhistogram(ir,ic))*
     -             coeff_eq(ir,iq,itype,ieq)
               enddo
            enddo
         enddo
         do iq=nx1(iexpt),nx2(iexpt)
            q2=(q(iq)/q(nq))**qweight
            sum_s(iexpt)=sum_s(iexpt)+work2(iq)*q2
            sum_ss(iexpt)=sum_ss(iexpt)+work2(iq)**2*q2*q2
            sum_fs(iexpt)=sum_fs(iexpt)+eexpt(iq,ieq)*work2(iq)*q2*q2
         enddo
      enddo
      return
      end




      subroutine rmca_chisq(nexpt,nx1,nx2,renorm,offset,
     - factor,constant,sum_f,sum_ff,sum_fs,sum_s,sum_ss,
     - sum_q,sigma,chisq,chisq1)
      implicit none
      real          chisq(*),constant(*),factor(*),sigma(*),sum_f(*),
     -              sum_ff(*),sum_fs(*),sum_q(*),sum_s(*),sum_ss(*)
      integer       nx1(*),nx2(*)
      logical       renorm(*),offset(*)
      real          chisq1
      integer       nexpt
      integer       iexpt
*
      do iexpt=1,nexpt
         if (renorm(iexpt).and.offset(iexpt)) then
            factor(iexpt)=(sum_q(iexpt)*sum_fs(iexpt)
     -       -sum_f(iexpt)*sum_s(iexpt))/
     -       (sum_q(iexpt)*sum_ss(iexpt)-sum_s(iexpt)*sum_s(iexpt))
         elseif (renorm(iexpt)) then
            factor(iexpt)=sum_fs(iexpt)/sum_ss(iexpt)
         else
            factor(iexpt)=1.
         endif
         if (offset(iexpt)) then
            constant(iexpt)=(sum_f(iexpt)-factor(iexpt)
     -       *sum_s(iexpt))/sum_q(iexpt)
         else
            constant(iexpt)=0.
         endif
      enddo
*
      chisq1=0.
      do iexpt=1,nexpt
         chisq(iexpt)=sum_ff(iexpt)+factor(iexpt)**2*sum_ss(iexpt)
     -    +constant(iexpt)**2*sum_q(iexpt)
     -    +2.*(constant(iexpt)*factor(iexpt)*sum_s(iexpt)
     -    -constant(iexpt)*sum_f(iexpt)
     -    -factor(iexpt)*sum_fs(iexpt))
         chisq(iexpt)=chisq(iexpt)/sigma(iexpt)/sigma(iexpt)
         chisq1=chisq1+chisq(iexpt)
      enddo
      return
      end





      subroutine rmca_summary(avcnew,avcoordno,cfnew,chisq,chisq1,
     - constant,coordfrac,energy,factor,nacc,navcoord,ncoord,nexpt,
     - ngen,ntried,nx1,nx2,usepot)
      implicit none
      real          avcnew(*),avcoordno(*),cfnew(*),chisq(*),
     -              constant(*),coordfrac(*),factor(*)
      integer       nx1(*),nx2(*)
      real          chisq1,energy,points
      integer       nacc,navcoord,ncoord,ndof,nexpt,ngen,ntried
      logical       usepot
      integer       icc,iexpt
*
      ndof=ncoord+navcoord
      do iexpt=1,nexpt
         ndof=ndof+nx2(iexpt)-nx1(iexpt)+1
      enddo
      if (nexpt.eq.0.and.ncoord.eq.0) then
         write (*,'(1x,i6,a,i6,a)')
     -   nacc,' moves accepted out of ',ngen,' generated'
      else
         write (*,*)
         write (*,'(1x,i6,a,i6,a,i6,a,g10.4)') nacc,' moves acc. ',
     -    ngen,' gen. and ',ntried,
     - ' tested; Chi-squared/d.o.f.=',chisq1/real(ndof)
      endif
      do iexpt=1,nexpt
         points=nx2(iexpt)-nx1(iexpt)+1
         write (*,'(1x,a,i2,a,f6.4,a,f7.4,a,g10.4)')
     -    '     Expt ',iexpt,': Renorm.=',1./factor(iexpt),
     -    '; Constant=',constant(iexpt),'; Chi**2/nq=',
     -    chisq(iexpt)/points
      enddo
      do icc=1,ncoord
         write (*,'(1x,a,i2,a,f6.2,a,f6.2,a)')
     -    '     Coordination constraint ',icc,': Fraction=',
     -    cfnew(icc)*100.,'%; target=',coordfrac(icc)*100.,'%'
      enddo
      do icc=1,navcoord
         write (*,'(1x,a,i2,a,f5.2,a,f5.2)')
     -    '     Coordination number constraint ',icc,
     -    ' Coord. no. is ',avcnew(icc),' target is ',avcoordno(icc)
      enddo
      if (usepot) write (*,'(a,g10.4)') '      Energy/kT/n=',energy
      return
      end




      subroutine rmca_output(coeff_eq,coeff_fq,coeff_gr,coeff_sq,
     - constant,dim_fq,dim_gr,dim_q,dim_r_exafs,dim_sq,dr,eexpt,
     - exafs_edge,factor,fexpt,file,gexpt,gnorm,gpar,
     - gtot,histogram,neq,nexpt,nfq,ngr,npar,nq,nr,nr_exafs,
     - nsq,ntypes,nx1,nx2,nxn,q,sexpt,sqr,title,version,work1,work2,
     - qweight)
      implicit none
      integer       dim_fq,dim_gr,dim_q,dim_r_exafs,dim_sq,nr,ntypes
      real          coeff_eq(dim_r_exafs,dim_q,ntypes,*),
     -              coeff_fq(dim_q,dim_fq,*),
     -              coeff_gr(dim_gr,*),coeff_sq(dim_sq,*),
     -              constant(*),eexpt(dim_q,*),factor(*),
     -              fexpt(dim_q,*),gexpt(nr,*),gnorm(*),gpar(nr,*),
     -              gtot(nr,*),q(*),sexpt(dim_q,*),sqr(nr,*),
     -              work1(*),work2(*)
      integer       exafs_edge(*),histogram(nr,*),nx1(*),nx2(*),nxn(*)
      real          q2,dr,snew
      integer       neq,nexpt,nfq,ngr,npar,nq,nr_exafs,nsq,qweight
      integer       i,i1,ic,ieq,iexpt,ifq,iq,ir,isq,itype,j1,lfile
      character*(*) version
      character*80  title
      character*20  file
      character*24  tempo
      logical       qex
*
      lfile=len(file)
      do while (file(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=file(1:lfile)//'.out'
      print *,tempo
      inquire(file=tempo,exist=qex)
      if(qex.eqv..FALSE.) go to 200
      call filedel(3,tempo)
200   call open_write1(3,tempo)
C      call open_write(3,file,'.out')
      write (3,*) 'Results from RMCA version ',version
      write (3,*) 'TITLE'
      write (3,*) title
      write (3,*) nexpt,' experiments;',npar,' partials'
      write (3,*) 'Partial g(r)''s'
      write (3,*) ' r, g(r)'
      write (3,*) 'PLOTS'
      write (3,*) nr,',',npar
      do ic=1,npar
         do ir=1,nr
            gpar(ir,ic)=real(histogram(ir,ic))/
     -       (gnorm(ir)*real(nxn(ic)))
         enddo
      enddo
      do ir=1,nr
         write (3,*) real(ir)*dr,(gpar(ir,ic),ic=1,npar)
      enddo
      write (3,*) 'ENDGROUP'
      write (3,*) 'Partial S(Q)''s'
      write (3,*) ' Q, S(Q)'
      write (3,*) 'PLOTS'
      write (3,*) nq,',',npar
      do i=1,nq*npar
         work1(i)=0.
      enddo
      do ic=1,npar
         do iq=1,nq
            do ir=1,nr
               work1((ic-1)*nq+iq)=work1((ic-1)*nq+iq)+
     -          (gpar(ir,ic)-1.)*sqr(ir,iq)
            enddo
         enddo
      enddo
      do iq=1,nq
         write (3,*) q(iq),(work1((ic-1)*nq+iq),ic=1,npar)
      enddo
      if (ngr.gt.0) then
         write (3,*) 'ENDGROUP'
         write (3,*) 'Total g(r)''s'
         do iexpt=1,nexpt
            do ir=1,nr
               gtot(ir,iexpt)=0.
            enddo
         enddo
         do iexpt=1,nexpt
            do ic=1,npar
               do ir=1,nr
                  gtot(ir,iexpt)=gtot(ir,iexpt)+coeff_gr(iexpt,ic)*
     -             (gpar(ir,ic)-1.)
               enddo
            enddo
         enddo
         write (3,*) 'Comparison with input data'
         do iexpt=1,nexpt
            write (3,*) ' r, g(r) (RMC), g(r) (Expt)'
            write (3,*) 'CURVES'
            write (3,*) nx2(iexpt)-nx1(iexpt)+1,',',2
            do ir=1,nr
               write (3,*) real(ir)*dr,gtot(ir,iexpt),
     -          gexpt(ir,iexpt)/factor(iexpt)
            enddo
         enddo
      endif
      if (nsq+nfq.gt.0) then
         write (3,*) 'ENDGROUP'
         write (3,*) 'Total S(Q)''s'
         write (3,*) 'Comparison with input data'
         do isq=1,nsq
            iexpt=ngr+isq
            write (3,*) ' Q, S(Q) (RMC), S(Q) (Expt)'
            write (3,*) 'CURVES'
            write (3,*) nx2(iexpt)-nx1(iexpt)+1,',',2
            do iq=nx1(iexpt),nx2(iexpt)
               work2(iq)=0.
            enddo
            do ic=1,npar
               do iq=nx1(iexpt),nx2(iexpt)
                  work2(iq)=work2(iq)+
     -             coeff_sq(isq,ic)*work1((ic-1)*nq+iq)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               write (3,*) q(iq),work2(iq),
     -          (sexpt(iq,isq)-constant(iexpt))/factor(iexpt)
            enddo
         enddo
         do ifq=1,nfq
            iexpt=ngr+nsq+ifq
            write (3,*) ' Q, S(Q) (RMC), S(Q) (Expt)'
            write (3,*) 'CURVES'
            write (3,*) nx2(iexpt)-nx1(iexpt)+1,',',2
            do iq=nx1(iexpt),nx2(iexpt)
               work2(iq)=0.
            enddo
            do ic=1,npar
               do iq=nx1(iexpt),nx2(iexpt)
                  work2(iq)=work2(iq)+
     -             coeff_fq(iq,ifq,ic)*work1((ic-1)*nq+iq)
               enddo
            enddo
            do iq=nx1(iexpt),nx2(iexpt)
               write (3,*) q(iq),work2(iq),
     -          (fexpt(iq,ifq)-constant(iexpt))/factor(iexpt)
            enddo
         enddo
      endif
      if (neq.gt.0) then
         write (3,*) 'ENDGROUP'
         write (3,*) 'EXAFS spectra'
         write (3,*) 'Comparison with input data'
         do ieq=1,neq
            iexpt=ngr+nsq+nfq+ieq
            write (3,*) ' Q, EXAFS(Q) (RMC), EXAFS(Q) (Expt)'
            write (3,*) 'CURVES'
            write (3,*) nx2(iexpt)-nx1(iexpt)+1,',',2
            do iq=nx1(iexpt),nx2(iexpt)
               snew=0.
               do itype=1,ntypes
                  i1=min(itype,exafs_edge(ieq))
                  j1=max(itype,exafs_edge(ieq))
                  ic=(i1-1)*(2*ntypes-i1)/2+j1
                  do ir=1,nr_exafs
                     snew=snew+real(histogram(ir,ic))*
     -                coeff_eq(ir,iq,itype,ieq)
                  enddo
               enddo
               q2=q(iq)**qweight
               write (3,*) q(iq),snew*q2,
     -          (eexpt(iq,ieq)*q2-constant(iexpt))/factor(iexpt)
            enddo
         enddo
      endif
      close(3)
      return
      end




      subroutine rmca_write_his(atoms,dr,file,histogram,n,nacc,ncoord,
     - ngen,ni,nneigh,ntried,nr,nsaved,ntypes,rcoord,rho,title,
     - typec,typen,vectors)
      implicit none
      integer       n,nr,lfile
      real          atoms(3,*),rcoord(2,*),vectors(3,3)
      integer       histogram(nr,*),ni(*),nneigh(n,*),
     -              typec(*),typen(*)
      character*8   hex(9)
      real          dr,rho
      integer       nacc,ncoord,ngen,nsaved,ntried,ntypes
      character*8   real_to_hex*8
      character*80  title
      character*20  file
      character*24  tempo
      logical       qex
      integer       i,j
      external      real_to_hex
*
      lfile=len(file)
      do while (file(lfile:lfile).eq.' ')
	 lfile=lfile-1
      enddo
      tempo=file(1:lfile)//'.his'
      inquire(file=tempo,exist=qex)
      if(qex.eqv..FALSE.) go to 200
      call filedel(3,tempo)
200   call open_write1(3,tempo)
C      call open_write(3,file,'.his')
      write (3,'(1x,a)') 'RMCA (v3) intermediate (histogram) file'
      write (3,1000)title
1000  format(1x,A80)
      write (3,*) ngen,ntried,nacc,nsaved,.false.
      do i=1,3
         do j=1,3
            hex((i-1)*3+j)=real_to_hex(vectors(j,i))
         enddo
      enddo
      write (3,10) (hex(i),i=1,9)
10    format(1x,9a8)
      hex(1)=real_to_hex(rho)
      write (3,20) hex(1),n,ntypes
20    format(1x,a8,1x,2i10)
      write (3,*) (ni(i),i=1,ntypes)
      do i=1,n
         hex(1)=real_to_hex(atoms(1,i))
         hex(2)=real_to_hex(atoms(2,i))
         hex(3)=real_to_hex(atoms(3,i))
         write (3,30) (hex(j),j=1,3)
      enddo
30    format(3(1x,a8))
      hex(1)=real_to_hex(dr)
      write (3,40) nr,hex(1)
40    format(i10,1x,a8)
      do i=1,ntypes*(ntypes+1)/2
         write (3,50) (histogram(j,i),j=1,nr)
      enddo
50    format(10(1x,i7))
      write (3,*) ncoord
      do i=1,ncoord
         hex(1)=real_to_hex(rcoord(1,i))
         hex(2)=real_to_hex(rcoord(2,i))
         write (3,60) typec(i),typen(i),hex(1),hex(2)
         write (3,50) (nneigh(j,i),j=1,n)
      enddo
60    format(2i10,1x,a8,1x,a8)
      close(3)
      return
      end
*      Contains the following machine dependent routines: 
*        rmc_init        initialises clock and random number generator
*        rmc_datafile    returns name of data file to be read
*        rmc_time        calculates CPU time used
*        rmc_rerun       passes RERUN flag to operating system
*        open_read       opens a file for reading
*        open_write      opens a file for writing 
*        real_to_hex     converts real*4 to IEEE hex string
*        hex_to_real     converts IEEE hex string to real*4
*
*-----------------------------------------------------------------------
*
*     Subroutine rmc_data must return the name of the data file to
*     be read by the programme. If possible this will be supplied
*     as a parameter on the command line that invokes the programme
*
      subroutine rmc_datafile(file)
      character*20 file
      PRINT 1
1     format('  entry file (no extension) ??',$)
      READ 2,file
2     format(A20)
C      call lib$get_foreign(file)
      return
      end
*
*-----------------------------------------------------------------------
*
*     Subroutine rmc_init must set the seed for the random number 
*     generator and obtain the current cpu clock reading in minutes
*
      subroutine rmc_init(iseed,time_start)
C      use portlib
C      real tt(2)
      iseed=int(secnds(0.0))*2+1
C      time_start=etime(tt)/60.
*
*  run the random number generator 100 times
*  for avoiding affects of starting value
*
      do 10 i=1,100
      bidon=ran(iseed)
10    iseed=iseed/3
      time_start=0.
      return
      end
*
*-----------------------------------------------------------------------
*
*     Subroutine rmc_time, given the initial cpu clock reading, must 
*     return the cpu time used (in minutes). If possible should return
*     a negative value if less than 1 minute remains.
*
      subroutine rmc_time(time_start,time_used)
      real tt(2)
C      time_used=etime(tt)/60.-time_start
      time_used=time_used+1.
      return
      end
*
*-----------------------------------------------------------------------
*
*     Subroutine rmc_rerun, can be used to return the rerun flag to
*     the operating system, if desired.
*
      subroutine rmc_rerun(rerun)
      logical rerun
C      if (rerun) call exit(10)
      return
      end
*
*-----------------------------------------------------------------------
*
*     Subroutine open_read1 opens a file for reading, filling in the
*     supplied extension if none is supplied with the file name.
*     Allows use of system specific facilities of open statement.
*
*
*     file formatted
*
      subroutine open_read1(unit,file)
      integer unit
      character*(*) file
      open (unit,file=file,status='old')
      return
      end
*
*     Subroutine open_write1 opens a file for writing, filling in the
*     supplied extension if none is supplied with the file name.
*     Allows use of system specific facilities of open statement.
*
*     file formatted
*
      subroutine open_write1(unit,file)
      integer unit
      character*(*) file
      open (unit,file=file,status='new')
      return
      end
*
      subroutine filedel(unit,file)
      integer unit
      character*(*) file
      open (unit,file=file,status='old')
      close(unit,status='delete')
      return
      end
*
*
*-----------------------------------------------------------------------
*
*     Function real_to_hex converts a 32 bit real value to an
*     8 character hex string describing its internal representation.
*     If intermediate files are to be portable this must be
*     converted to IEEE form.
*
*     Assume storage is in Convex native format rather than IEEE.
*     This means we must shift exponent by 2
*
      character*8 function real_to_hex(r)
      character*8 temp
C      write (temp,'(z8)') r/4.
      write (10,*) r/4.
      real_to_hex=temp
      return
      end
*
*-----------------------------------------------------------------------
*     Function hex_to_real converts an 8 character hex string in IEEE 
*     format to a 32 bit real value.
*
      real function hex_to_real(hex,error)
      character hex*(*)
      logical error
      error=.false.
C      read (hex,'(z8)',end=10,err=10) hex_to_real
      read (10,*,end=10,err=10) hex_to_real
      hex_to_real=hex_to_real*4.
      return
10    error=.true.
      return
      end
*
*
*     Subroutines used by RMC code.
*     Last changed at version 3.00
*
      subroutine rmc_volume(vectors,truncated,volume)
      implicit none
      real     vectors(3,3),volume,triprod
      logical  truncated
*     
*     Finds the volume of the configuration cell
*
      triprod=vectors(1,1)*vectors(2,2)*vectors(3,3)
     -       +vectors(2,1)*vectors(3,2)*vectors(1,3)
     -       +vectors(3,1)*vectors(1,2)*vectors(2,3)
     -       -vectors(3,1)*vectors(2,2)*vectors(1,3)
     -       -vectors(2,1)*vectors(1,2)*vectors(3,3)
     -       -vectors(1,1)*vectors(3,2)*vectors(2,3)
      volume=8.*abs(triprod)
      if (truncated) volume=volume/2.
      return
      end




      subroutine rmc_metric(vectors,metric)
      implicit none
      real     vectors(3,3),metric(3,3)
      integer  i,j,k
*     
*     Calculates the metric matrix used for calculating distances
*
      do j=1,3
         do i=1,3
            metric(i,j)=0.
            do k=1,3
               metric(i,j)=metric(i,j)+vectors(k,i)*vectors(k,j)
            enddo
         enddo
      enddo
      return
      end



      subroutine rmc_cellsize(vectors,truncated,volume,d)
      implicit none
      real    vectors(3,3),volume,d,triprod,d1,d2,d3,d4
      real    axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3
      logical truncated
*     
*     Find shortest distance from centre of cell to a face
*
      axb1=vectors(2,1)*vectors(3,2)-vectors(3,1)*vectors(2,2)
      axb2=vectors(3,1)*vectors(1,2)-vectors(1,1)*vectors(3,2)
      axb3=vectors(1,1)*vectors(2,2)-vectors(2,1)*vectors(1,2)
      bxc1=vectors(2,2)*vectors(3,3)-vectors(3,2)*vectors(2,3)
      bxc2=vectors(3,2)*vectors(1,3)-vectors(1,2)*vectors(3,3)
      bxc3=vectors(1,2)*vectors(2,3)-vectors(2,2)*vectors(1,3)
      cxa1=vectors(2,3)*vectors(3,1)-vectors(3,3)*vectors(2,1)
      cxa2=vectors(3,3)*vectors(1,1)-vectors(1,3)*vectors(3,1)
      cxa3=vectors(1,3)*vectors(2,1)-vectors(2,3)*vectors(1,1)
      triprod=volume/8.
      if (truncated) triprod=volume/4.
      d1=triprod/sqrt(axb1**2+axb2**2+axb3**2)
      d2=triprod/sqrt(bxc1**2+bxc2**2+bxc3**2)
      d3=triprod/sqrt(cxa1**2+cxa2**2+cxa3**2)
      d=min(d1,d2,d3)
      if (truncated) then
         d1=1.5*triprod/sqrt(
     -    (axb1+bxc1+cxa1)**2+(axb2+bxc2+cxa2)**2+(axb3+bxc3+cxa3)**2)
         d2=1.5*triprod/sqrt(
     -    (axb1-bxc1+cxa1)**2+(axb2-bxc2+cxa2)**2+(axb3-bxc3+cxa3)**2)
         d3=1.5*triprod/sqrt(
     -    (axb1+bxc1-cxa1)**2+(axb2+bxc2-cxa2)**2+(axb3+bxc3-cxa3)**2)
         d4=1.5*triprod/sqrt(
     -    (axb1-bxc1-cxa1)**2+(axb2-bxc2-cxa2)**2+(axb3-bxc3-cxa3)**2)
         d=min(d,d1,d2,d3,d4)
      endif
      return
      end




      subroutine rmc_cell(d,metric,truncated,vectors,volume)
      real metric(3,3),vectors(3,3)
      logical truncated
*
      triprod=vectors(1,1)*vectors(2,2)*vectors(3,3)
     -       +vectors(2,1)*vectors(3,2)*vectors(1,3)
     -       +vectors(3,1)*vectors(1,2)*vectors(2,3)
     -       -vectors(3,1)*vectors(2,2)*vectors(1,3)
     -       -vectors(2,1)*vectors(1,2)*vectors(3,3)
     -       -vectors(1,1)*vectors(3,2)*vectors(2,3)
      volume=8.*abs(triprod)
      if (truncated) volume=volume/2.
*     
*     Find shortest distance from centre of cell to a face
*
      axb1=vectors(2,1)*vectors(3,2)-vectors(3,1)*vectors(2,2)
      axb2=vectors(3,1)*vectors(1,2)-vectors(1,1)*vectors(3,2)
      axb3=vectors(1,1)*vectors(2,2)-vectors(2,1)*vectors(1,2)
      bxc1=vectors(2,2)*vectors(3,3)-vectors(3,2)*vectors(2,3)
      bxc2=vectors(3,2)*vectors(1,3)-vectors(1,2)*vectors(3,3)
      bxc3=vectors(1,2)*vectors(2,3)-vectors(2,2)*vectors(1,3)
      cxa1=vectors(2,3)*vectors(3,1)-vectors(3,3)*vectors(2,1)
      cxa2=vectors(3,3)*vectors(1,1)-vectors(1,3)*vectors(3,1)
      cxa3=vectors(1,3)*vectors(2,1)-vectors(2,3)*vectors(1,1)
      d1=triprod/sqrt(axb1**2+axb2**2+axb3**2)
      d2=triprod/sqrt(bxc1**2+bxc2**2+bxc3**2)
      d3=triprod/sqrt(cxa1**2+cxa2**2+cxa3**2)
      d=min(d1,d2,d3)
      if (truncated) then
         d1=1.5*triprod/sqrt(
     -    (axb1+bxc1+cxa1)**2+(axb2+bxc2+cxa2)**2+(axb3+bxc3+cxa3)**2)
         d2=1.5*triprod/sqrt(
     -    (axb1-bxc1+cxa1)**2+(axb2-bxc2+cxa2)**2+(axb3-bxc3+cxa3)**2)
         d3=1.5*triprod/sqrt(
     -    (axb1+bxc1-cxa1)**2+(axb2+bxc2-cxa2)**2+(axb3+bxc3-cxa3)**2)
         d4=1.5*triprod/sqrt(
     -    (axb1-bxc1-cxa1)**2+(axb2-bxc2-cxa2)**2+(axb3-bxc3-cxa3)**2)
         d=min(d,d1,d2,d3,d4)
      endif
      do j=1,3
         do i=1,3
            metric(i,j)=0.
            do k=1,3
               metric(i,j)=metric(i,j)+vectors(k,i)*vectors(k,j)
            enddo
         enddo
      enddo
      return
      end
*     Subroutines are used for reading and writing RMC (or MD/MC) 
*     configuration files:
*
*     read_rmc_config   } general reading and writing routines
*     write_rmc_config  }
*
*     read_config       } simpler routines more appropriate
*     write_config      } for atomic systems
*
***********************************************************************

      subroutine read_rmc_config(title,ngen,ntried,nacc,nsaved,
     - max_molecules,max_moltypes,max_sites,nvar,
     - truncated,vectors,
     - nmol,nmoltypes,ni,nsites,sites,
     - neuler,centres,ichan)

      real          centres(nvar,max_molecules)
      real          sites(3,max_sites,max_moltypes)
      real          vectors(3,3)
      integer       ni(max_moltypes),nsites(max_moltypes)
      character*80  title
      logical truncated,space
*
*     Reads a configuration data file
*
*     Check format of file
*     Also have to allow for possibility of initial space
*     (Fortran carriage control) which may or may not occur at
*     beginning of each line
*
      read (ichan,1000) title
1000  format(A80)
      print 1000,title
      space=.false.
      if (title(1:1).eq.' ') space=.true.
      if (space) title=title(2:len(title))
      if (title(1:10).ne.'(Version 3') then
         write (*,*) 'Configuration file has wrong format'
         stop
      endif
*
      read (ichan,1000) title
      if (space) title=title(2:len(title))
      read (ichan,*) ngen,ntried,nacc
      read (ichan,*) nsaved
      read (ichan,*) nmol
      read (ichan,*) nmoltypes
      read (ichan,*) nsites_max
      read (ichan,*) neuler
*
      if (nmol.gt.max_molecules) then
         write (*,*) 'Configuration contains ',nmol,' particles'
         write (*,*) 'Arrays are only large enough for ',max_molecules
         stop
      endif
      if (nmoltypes.gt.max_moltypes) then
         write (*,*) 'Configuration contains ',nmoltypes,
     -     ' types of particle'
         write (*,*) 'Arrays are only large enough for ',max_moltypes
         stop
      endif
      if (nsites_max.gt.max_sites) then
         write (*,*) 'Configuration contains molecules with',nsites_max,
     -    'atomic sites'
         write (*,*) 'Arrays are only large enough for ',max_sites
         stop
      endif
      if (neuler+3.gt.nvar) then
         write (*,*) 'Configuration specifies ',neuler,' Euler angles'
         write (*,*) 'Arrays are only large enough for ',nvar-3
         stop
      endif
* 
      read (ichan,*)
      read (ichan,*) truncated
      read (ichan,*)
      read (ichan,*) vectors
*
      do itype=1,nmoltypes
         read (ichan,*) ni(itype)
         read (ichan,*) nsites(itype)
         read (ichan,*) ((sites(i,isite,itype),i=1,3),
     -    isite=1,nsites(itype))
      enddo
*
      do i=1,nmol
         read (ichan,*) (centres(j,i),j=1,3+neuler)
      enddo
      return
      end

***********************************************************************

      subroutine write_rmc_config(title,ngen,ntried,nacc,nsaved,
     - max_sites,nvar,
     - truncated,vectors,
     - nmol,nmoltypes,ni,nsites,sites,
     - neuler,centres,ichan)
      real          centres(nvar,*)
      real          sites(3,max_sites,*)
      real          vectors(3,3)
      integer       ni(*),nsites(*)
      character*80  title
      logical truncated
*
*     Writes a configuration data file
*
      write (ichan,*) '(Version 3 format configuration file)'
      write (ichan,1000) title
1000  format(1X,A80)
*
      write (ichan,*) 
      write (ichan,'(1x,3i10,a)') ngen,ntried,nacc,
     - ' moves generated, tried, accepted'
      write (ichan,'(1x,i10,20x,a)') nsaved,' configurations saved'
      write (ichan,*) 
      write (ichan,'(1x,i10,a)') nmol,' molecules (of all types)'
      write (ichan,'(1x,i10,a)') nmoltypes,' types of molecule'
      nsites_max=1
      do i=1,nmoltypes
         nsites_max=max(nsites_max,nsites(i))
      enddo
      write (ichan,'(1x,i10,a)') nsites_max,
     - ' is the largest number of atoms in a molecule'
      write (ichan,'(1x,i10,a)') neuler,' Euler angles are provided'
*
      write (ichan,*)
      if (truncated) then
         write (ichan,'(1x,l10,a)') truncated,
     -    ' (Box is truncated octahedral)'
      else
         write (ichan,'(1x,l10,a)') truncated,
     -    ' (Box is not truncated octahedral)'
      endif
      write (ichan,'(11x,a)') ' Defining vectors are:'
      write (ichan,'(3(11x,3(1x,f10.6)/))') vectors
*
      do itype=1,nmoltypes
         write (ichan,'(1x,i10,a,i2)') ni(itype),
     -     ' molecules of type ',itype
         write (ichan,'(1x,i10,a)') nsites(itype),' atomic sites'
         write (ichan,'(100(11x,3(1x,f10.6)/))')
     -    ((sites(i,isite,itype),i=1,3),isite=1,nsites(itype))
      enddo
*
      do i=1,nmol
         write (ichan,'(3(1x,f10.7),3(1x,f10.8))') 
     -    (centres(j,i),j=1,3+neuler)
      enddo
      return
      end

************************************************************************

      subroutine read_config(title,ngen,ntried,nacc,nsaved,
     - max_atoms,max_types,
     - truncated,vectors,
     - n,ntypes,ni,
     - atoms,ichan)
      real          atoms(3,max_atoms)
      real          vectors(3,3)
      integer       ni(*)
      real          sites(3,1,100)
      integer       nsites(100)
      character*80  title
      logical truncated
*
      call read_rmc_config(title,ngen,ntried,nacc,nsaved,
     - max_atoms,max_types,1,3,truncated,vectors,
     - n,ntypes,ni,nsites,sites,neuler,atoms,ichan)
*
      return
      end

***********************************************************************

      subroutine write_config(title,ngen,ntried,nacc,nsaved,
     - truncated,vectors,n,ntypes,ni,atoms,ichan)
      real          atoms(3,*)
      real          vectors(3,3)
      integer       ni(*)
      real          sites(3,1,100)
      integer       nsites(100)
      character*80  title
      logical truncated
*
      do i=1,ntypes
         nsites(i)=1
         sites(1,1,i)=0.
         sites(2,1,i)=0.
         sites(3,1,i)=0.
      enddo
      call write_rmc_config(title,ngen,ntried,nacc,nsaved,
     - 1,3,truncated,vectors,n,ntypes,ni,nsites,sites,
     - 0,atoms,ichan)
*
      return
      end
