module alloc_mod
!!
!!  allocates multiple variables used in est_noise. These are essentially global variables
!   that are accessible from the main program and its subroutines. This allows the contents
!   of these variable to be passed with invoking their names in the arguments of the subroutine
!~
! This module also allocates memory dependent upon the type of noise model being used.
!
  use iso_fortran_env
  implicit none
  private
  !  stuff for least squares inversion
  real(kind=real64), allocatable, public :: covar(:,:), covinv(:,:),A(:,:),t(:),x(:),ex(:), t_year(:),xsave(:)
  real(kind=real32), allocatable, public :: d(:),res(:)
  integer, public :: nmod, ic     !!  number of model parameters, number of data  

!  parameters used to define the least-squares model
  real(kind=real32), public :: bexp(10),expmax
  character(len=7), public :: exp_choice(10),exp_type(10)
  character(len=3), public :: net
  real(kind=real64), public :: rat_chg1(82),rat_chg2(82),off(82),per(82), &
      texp(10), rat_chng_norm(82),dec_time1,dec_time2,aux_norm(82),rate_norm, &
      t_start,t_stop,dt_sam,Anor(82)
  integer, public :: n_file_press,n_rat_chg,nbt,noff,nper,nrate,n_exp,nexp_fix, iseed, &
       nmiss,irowmiss(10000),ipl_flag_1,ipl_flag_2,ibp_flag,kmax,jmax,maxnmod,mdmax, irowOffset,nmod_orig
  real(kind=real64), allocatable, public :: filtpl1(:),filtpl2(:),filtbp(:),covarpl1(:,:),covarpl2(:,:),covarbp(:,:), &
      Finv(:)
      
!  variables used when ModType=f --  additive noise with "fast" inversion
!!  this set is used in the mle subroutine found in mle_mod
    real(kind=real64), public, allocatable :: f(:),E(:,:),AAA(:,:),AAAA(:,:),B(:,:)
    integer, public :: mmax,max_time
!    integer, allocatable :: irow(:) 
 
 ! this set found in calcres subroutine found in est_noise_mod
   real (kind=real64), public, allocatable :: a1(:,:),a2(:,:),evd(:),ai(:,:),ainv(:,:), &
      Aw(:,:),Awt(:,:),a3(:,:)
    real(kind=real64), public, allocatable :: dw(:),din(:),dr(:)
    real(kind=real64), public, allocatable :: dxx(:),corrxx(:),Axx(:,:)

    public :: CreateModTypeF, DestructModTypeF, CreateModTypeQ,DestructModTypeQ
    public :: CreateModTypeC, DestructModTypeC, CreateGlobal, DestructGlobal      
  contains
  
  subroutine CreateModTypeQ
! invoked for ModType=Q
    use iso_fortran_env
    implicit none

    mdmax=ic
    mdmax=size(t_year)

    allocate(a1(maxnmod,mdmax))
    allocate(a2(maxnmod,maxnmod))
    allocate(evd(maxnmod))
    allocate(ai(maxnmod,maxnmod))
    allocate(ainv(maxnmod,mdmax))

  end subroutine CreateModTypeQ
  subroutine DestructModTypeQ
! invoked for ModType=Q
    use iso_fortran_env
    implicit none
    if (allocated(a1)) deallocate(a1)
    if (allocated(a2)) deallocate(a2)
    if (allocated(evd))  deallocate(evd)
    if (allocated(ai)) deallocate(ai)
    if (allocated(ainv)) deallocate(ainv)
  end subroutine DestructModTypeQ


  
  subroutine CreateModTypeC
!  allocates arrays for ModType=c
    use iso_fortran_env
    implicit none
 !! this set found in mle subroutine      
       mdmax=size(t_year)
       max_time=size(t_year)

       
       if (.not. allocated(f)) allocate(f(mdmax+irowOffset))
       allocate(dxx(mdmax))
       if (.not. allocated(dr)) allocate(dr(mdmax))
       allocate(Axx(mdmax,maxnmod))
       allocate(corrxx(mdmax))
       if (.not. allocated(a1)) allocate(a1(maxnmod,mdmax))
       if (.not. allocated(a2)) allocate(a2(maxnmod,maxnmod))
       if (.not. allocated(evd)) allocate(evd(maxnmod))
       if (.not. allocated(ai)) allocate(ai(maxnmod,maxnmod))
       if (.not. allocated(ainv)) allocate(ainv(maxnmod,mdmax))
  end subroutine CreateModTypeC
  
  subroutine DestructModTypeC
!  deallocates arrays for ModType=c
    use iso_fortran_env
    implicit none
 !! this set found in mle subroutine      

    if (allocated(f)) deallocate(f)
    if (allocated(dxx)) deallocate(dxx)
    if (allocated(dr)) deallocate(dr)
    if (allocated(corrxx)) deallocate(corrxx)
    if (allocated(Axx)) deallocate(Axx)
    if (allocated(a1)) deallocate(a1)
    if (allocated(a2)) deallocate(a2)
    if (allocated(evd)) deallocate(evd)
    if (allocated(ai)) deallocate(ai)
    if (allocated(ainv)) deallocate(ainv)

  end subroutine DestructModTypeC


  
  subroutine CreateModTypeF
!  allocates arrays for ModType=f
    use iso_fortran_env
    implicit none
 !! this set found in mle subroutine      
       mdmax=size(t_year)
       max_time=size(t_year)
       allocate(f(mdmax))
       allocate(finv(mdmax))
!       allocate(irow(kmax))
       allocate(E(mdmax,nmiss))
       allocate(AAA(nmiss,nmiss))
       allocate(AAAA(mdmax,nmiss))
       allocate(B(mdmax,mdmax))
       
! this set found in calcres subroutine
    allocate(a1(maxnmod,mdmax))
    allocate(a2(maxnmod,maxnmod))
    allocate(evd(maxnmod))
    allocate(ai(maxnmod,maxnmod))
    allocate(ainv(maxnmod,mdmax))
     allocate(dw(mdmax))
      allocate(din(mdmax))
     allocate(dr(mdmax))
      allocate(Aw(mdmax,nmod))
      allocate(Awt(maxnmod,mdmax))
      allocate(a3(maxnmod,mdmax))


  end subroutine CreateModTypeF
  subroutine DestructModTypeF
    use iso_fortran_env
    implicit none

      if (allocated(E)) deallocate(E)
      if (allocated(AAA)) deallocate(AAA)
      if (allocated(AAAA)) deallocate(AAAA)
      if (allocated(B)) deallocate(B)
      if (allocated(f)) deallocate(f)
       if (allocated(finv)) deallocate(finv)
!  found in calcres subroutine    
      if (allocated(dw)) deallocate(dw)
      if (allocated(din)) deallocate(din)
      if (allocated(dr)) deallocate(dr)
      if (allocated(Aw)) deallocate(Aw)
      if (allocated(Awt)) deallocate(Awt)
      if (allocated(a3))  deallocate(a3)   
      if (allocated(a1)) deallocate(a1)
      if (allocated(a2)) deallocate(a2)
      if (allocated(evd)) deallocate(evd)
      if (allocated(ai)) deallocate(ai)
      if (allocated(ainv)) deallocate(ainv)
  end subroutine DestructModTypeF


  
  subroutine CreateGlobal
!  allocate and make public the following variables
!  These variables are the least squares problem, d=Ax
!    subject to the covariance and its inverse
    use iso_fortran_env
    implicit none
    maxnmod=nmod+n_exp
    ic=ic+5
    allocate(covinv(ic,ic))    !! data covariance matrix
    allocate(covar(ic,ic))     !! inverser covariance
    allocate(A(ic,maxnmod))       !! design or trajectory matrix
    allocate(t(ic))            !!  time since t_start in days
    allocate(t_year(ic))            !!  time since t_start in yearss
    allocate(x(maxnmod))          !! the estimated model
    allocate(xsave(maxnmod))          !! the estimated model
    if (maxnmod .lt. 7) then
      allocate(ex(7))
    else
      allocate(ex(maxnmod))          !!  the error bars for the model
    end if
    allocate(d(ic))            !!  the data
    allocate(res(ic))          !! residuals, difference between obs and predicted
   
  end subroutine CreateGlobal
  subroutine DestructGlobal
!  deallocates the global variable
    use iso_fortran_env
    implicit none  
    if (allocated(covinv)) deallocate(covinv)
    if (allocated(covar)) deallocate(covar)
    if (allocated(A)) deallocate(A)
    if (allocated(t)) deallocate(t)
    if (allocated(t_year)) deallocate(t_year)
    if (allocated(x)) deallocate(x)
    if (allocated(ex)) deallocate(ex)
    if (allocated(d)) deallocate(d)
    if (allocated(res)) deallocate(res)
    return
  end subroutine DestructGlobal  
  
  
end module alloc_mod

  
