module tblite_scf_mixer_diis
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_lapack, only : getrf, getri, gesv
   use tblite_blas, only: gemm
   private

   public :: new_diis

   type, public :: diis_input
      !> Number of steps to keep in memory
      integer :: memory
   end type diis_input

   type, public, extends(mixer_type) :: diis_mixer
      integer :: ndim
      integer :: memory
      integer :: iter
      real(wp) :: e_max
      real(wp), allocatable :: err_v(:)
      real(wp), allocatable :: F(:)
      real(wp), allocatable :: D(:)
      real(wp), allocatable :: S(:,:)
   contains
      !> Apply mixing
      procedure :: next
      !> Set Fock-matrix from 1D array
      procedure :: set_1d_F
      !> Set density-Matrix from 1D array
      procedure :: set_1d
      !> Get Fock-Matrix as 1D array
      procedure :: get_1d_F
      !> Get error metric from mixing
      procedure :: get_error
      !> Not used, only to satisfy the deferred statement
      procedure :: diff_1d
      !> Not used, only to satisfy the deferred statement
      procedure :: get_1d      
   end type diis_mixer

   contains

   !> Initialize DIIS mixer
   subroutine new_diis(self, S, input)
      class(diis_mixer), intent(out) :: self
      integer :: ndim
      real(wp), intent(in) :: S(:,:)
      type(diis_input), intent(in) :: input
      ndim = size(S, 1)
      self%ndim = ndim
      self%iter = 0
      self%memory = input%memory
      self%e_max = 0
      allocate(self%S(ndim,ndim))
      self%S = S
      allocate(self%err_v(ndim*(ndim+1)/2))
      allocate(self%D(ndim*ndim))
      allocate(self%F(ndim*(ndim+1)/2))
   end subroutine new_diis

   !> Routines for accessing D and F

   !> Set new density-matrix from 1D array
   subroutine set_1d(self, qvec)
      !> Instance of the mixer
      class(diis_mixer), intent(inout) :: self
      !> Density vector
      real(wp), intent(in) :: qvec(:)
      self%D = qvec
   end subroutine set_1d

   !> Set new Fock-matrix from 1D array
   subroutine set_1d_F(self, f_1d)
      !> Instance of the mixer
      class(diis_mixer), intent(inout) :: self
      !> Fock-matrix
      real(wp), intent(in) :: f_1d(:)
      integer :: indx
      indx = 0
      do i=1,self%ndim
         do j=i,self%ndim
            indx = indx + 1
            self%F(indx) = f_1d(j+self%ndim*(i-1))
         end do
      end do
   end subroutine set_1d_F
   
   !> Apply mixing 
   subroutine next(self, error)
      class(diis_mixer), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: info

      self%iter = self%iter + 1
      call get_emax(self)
      call diis(self%ndim, self%iter, self%F, self%e_max, self%err_v, self%memory, info)
      if (info /= 0) call fatal_error(error, "DIIS failed to obtain next iteration")
   end subroutine next

   !> Get Fock-matrix as 1D array
   subroutine get_1d_F(self, coeff)
      class(diis_mixer), intent(inout) :: self
      real(wp), intent(out) :: coeff(:)
      integer :: indx
      indx = 0
      do i=1,self%ndim
         do j=i,self%ndim
            indx = indx + 1
            coeff(i+self%ndim*(j-1)) = self%F(indx)
            coeff(j+self%ndim*(i-1)) = self%F(indx)
         end do
      end do
   end subroutine get_1d_F

   !> Not used, only to satisfy the deferred statement
   subroutine diff_1d(self, qvec)
      class(diis_mixer), intent(inout) :: self
      real(wp), intent(in) :: qvec(:)
   end subroutine diff_1d
   
   !> Not used, only to satisfy the deferred statement
   subroutine get_1d(self, qvec)
      class(diis_mixer), intent(inout) :: self
      real(wp), intent(out) :: qvec(:)
   end subroutine get_1d

   !> Get error for output
   pure function get_error(self) result(error)
      class(diis_mixer), intent(in) :: self
      real(wp) :: error
      error = self%e_max
   end function get_error

   !> Calculate error for DIIS
   subroutine get_emax(self)
      class(diis_mixer), intent(inout) :: self
      integer :: indx
      real(wp) :: F_2d(self%ndim, self%ndim), D_2d(self%ndim, self%ndim)
      real(wp) :: tmp(self%ndim,self%ndim), tmp1(self%ndim,self%ndim)
      real(wp) :: tmp2(self%ndim,self%ndim)
      real(wp) :: tmp3(self%ndim,self%ndim), tmp4(self%ndim,self%ndim)
      indx = 0
      do i=1,self%ndim
         do j=i,self%ndim
            indx = indx + 1
            F_2d(i,j) = self%F(indx)
            F_2d(j,i) = self%F(indx)
         end do
      end do
      D_2d = reshape(self%D, shape(D_2d))
      tmp = 0.0_wp
      call gemm(D_2d, F_2d, tmp1, beta=0.0_wp)
      call gemm(D_2d, self%S, tmp2, beta=0.0_wp)
      call gemm(self%S, tmp1, tmp3, beta=0.0_wp)
      call gemm(F_2d, tmp2, tmp4, beta=0.0_wp)
      tmp = tmp3 - tmp4
      self%e_max = abs(tmp(1,1))
      indx = 1
      do i=1,self%ndim
         do j=i,self%ndim
            self%err_v(indx) = tmp(j,i)
            if (abs(self%err_v(indx)) > self%e_max) self%e_max = abs(self%err_v(indx))
            indx=indx+1
         end do
      end do
   end subroutine get_emax

   !> General linear algebra
   subroutine hadamard(vect1,vect2,n)
      integer, intent(in) :: n
      real(wp), intent(in) :: vect1(n)
      real(wp) :: prod_vect(n,1)
      real(wp), intent(inout) :: vect2(n,1)
      do i=1,n
         prod_vect(i,1) = vect1(i) * vect2(i,1)
      end do
      vect2 = prod_vect
   end subroutine hadamard

   subroutine frobenius(amat,bmat,prod)
      integer :: indx, cnt
      real(wp),intent(in) :: amat(:), bmat(:)
      real(wp),intent(out) :: prod
      prod = 0.0_wp
      cnt = 0
      indx = 1
      do i=1,size(amat)
         prod = prod + amat(i)*bmat(i)
         if (i/=indx) then
            prod = prod + (amat(i)*bmat(i))
         else
            indx = indx + size(amat)-cnt
            cnt = cnt+1
         end if
      end do
   end subroutine frobenius
   
   subroutine norm_i(amat, norm)
      real(wp),intent(in) :: amat(:,:)
      real(wp) :: norm, tmp
      norm = 0
      do i=1,size(amat,1)
         tmp = 0
         do j=1,size(amat,1)
            tmp = tmp + abs(amat(j,i))
         end do
         if (tmp>norm) norm = tmp
      end do
   end subroutine norm_i
   
   subroutine build_scale_matrix(amat, scale_matr)
      real(wp), intent(in) :: amat(:,:)
      real(wp) :: tmp
      real(wp), intent(out) :: scale_matr(:)
      do i=1,size(amat,1)
         tmp = 0
         do j=1,size(amat,1)
            tmp = tmp + abs(amat(i,j))
         end do
         scale_matr(i) = 1/tmp
      end do
   end subroutine build_scale_matrix
   
   !> Specific to DIIS 
   subroutine build_F(B, n, it, omit, saved_F, sol, F, info)
      integer, intent(in) :: n, it, omit
      integer, intent(out) :: info
      integer :: perm(it-omit)
      real(wp), intent(inout) :: B(it-omit,it-omit), sol(it-omit,1)
      real(wp), intent(in) :: saved_F(:)
      real(wp), intent(out) :: F(:)

      call call_lapack_gesv(B, sol, info, perm)
      if (info==0) then
         F = 0.0_wp
         do i=1,it-omit-1
            F = F + sol(i,1)*saved_F((i+omit-1)*n*(n+1)/2+1:(i+omit)*n*(n+1)/2)
         end do
      end if
   end subroutine build_F

   subroutine call_lapack_gesv(a, b, info, ipiv)
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(inout) :: b(:,:)
      integer, intent(out) :: info
      integer, intent(out) :: ipiv(:)
      
      call gesv(a, b, info, ipiv)
   end subroutine call_lapack_gesv

   subroutine build_B(err_v, sv_err, n, it, omit, old_B, B)
      integer, intent(in) :: n, it, omit
      real(wp), intent(in) :: sv_err(:), err_v(:), old_B(:, :)
      real(wp), intent(out) :: B(:, :)
      real(wp) :: Tr
      B=1.0_wp
      B(it-omit+1, it-omit+1) = 0.0_wp
      if (it /= 1) then
         B(1:it-omit-1, 1:it-omit-1) = old_B
         do i=1,it-omit-1
            call frobenius(sv_err((i+omit-1)*n*(n+1)/2+1:(i+omit)*n*(n+1)/2), err_v, Tr)
            B(it-omit, i) = Tr
            B(i, it-omit) = Tr
         end do
      end if
      call frobenius(err_v, err_v, Tr)
      B(it-omit, it-omit) = Tr
   end subroutine build_B

   subroutine save_vectors(F, err_v, n, it, old_sv_F, old_sv_err, sv_F, sv_err)
      integer, intent(in) :: n, it
      real(wp), intent(in) :: F(:), err_v(:)
      real(wp), intent(in) :: old_sv_F(:), old_sv_err(:)
      real(wp), intent(out) :: sv_F(it*n*(n+1)/2), sv_err(it*n*(n+1)/2)
      if (it/=1) then
         sv_err(1:(it-1)*n*(n+1)/2) = old_sv_err
         sv_F(1:(it-1)*n*(n+1)/2) = old_sv_F
      end if        
      sv_err((it-1)*n*(n+1)/2+1:it*n*(n+1)/2) = err_v
      sv_F((it-1)*n*(n+1)/2+1:it*n*(n+1)/2) = F
   end subroutine save_vectors

   subroutine rescale_B(k,it,omit, B, B_scaled, sol)
      integer, intent(in) :: it, omit
      integer :: info, perm(it-omit)
      real(wp), intent(in) :: B(:,:)
      real(wp) :: scale_matr(it-omit)
      real(wp) :: tmp_matr(it-omit,it-omit), norm1, norm2
      real(wp), intent(out) :: sol(it-omit,1), B_scaled(it-omit,it-omit), k

      call build_scale_matrix(B,scale_matr)
      do i=1,it-omit
         do j=1,it-omit
            tmp_matr(j,i) = scale_matr(j)*B(j,i)
         end do
      end do
      sol = 0d0
      sol(it-omit,1) = 1d0
      call hadamard(scale_matr, sol, it-omit)
      B_scaled = tmp_matr
      call call_lapack_getrx(tmp_matr, perm, info)
      call norm_i(B_scaled, norm1)      
      call norm_i(tmp_matr, norm2)
      k = norm1*norm2
   end subroutine rescale_B

   !> Calls the Lapack routines with the correct intent
   subroutine call_lapack_getrx(a, ipiv, info)
      real(wp), intent(inout) :: a(:,:)
      integer, intent(out) :: ipiv(:)
      integer, intent(out) :: info

      call getrf(a, ipiv, info)
      call getri(a, ipiv, info)
   end subroutine call_lapack_getrx

   subroutine diis(n, it, F, e_max, error_v, memory, info)
      !parameters
      integer :: max_k = 100000, min_dim = 5
      real(wp) :: diis_thresh = 1.0d0
      logical :: start_diis

      integer, intent(in) :: n, it, memory
      integer, intent(out) :: info
      real(wp), intent(inout) :: F(:)
      real(wp), intent(in) :: error_v(:), e_max
      real(wp) :: k, e_3
      real(wp), allocatable :: old_sv_err(:), old_sv_F(:), B_scaled(:,:), old_B(:,:), sol(:,:)
      real(wp), allocatable, save :: B(:,:), saved_err(:), saved_F(:)
      real(wp), save :: save_emax(3)
      integer, save :: mov, var
      
      !Getting e_3 and checking if in range to start
      if (e_max<diis_thresh) then
         start_diis = .TRUE.
      else
         start_diis = .FALSE.
         info = 0
      end if
      if (it == 1) then
         mov = 0
         var = 1
      end if
      if (var == 1) then
         save_emax(1) = e_max
         var = 2
         if (it>2) then
            e_3 = (save_emax(2)-save_emax(3))+(save_emax(3)-save_emax(1))
         else
            e_3 = 1
         end if
      else if (var == 2) then
         save_emax(2) = e_max
         var = 3
         if (it>2) then
            e_3 = (save_emax(3)-save_emax(1))+(save_emax(1)-save_emax(2))
         else
            e_3 = 1
         end if
      else
         save_emax(3) = e_max
         var = 1
         e_3 = (save_emax(1)-save_emax(2))+(save_emax(2)-save_emax(3))
      end if

      !Preparing for B construction
      allocate(old_sv_err((it-1)*n*(n+1)/2))
      allocate(old_sv_F((it-1)*n*(n+1)/2))
      allocate(old_B(it-1-mov,it-1-mov))
      if (it/=1) then
         old_sv_err = saved_err
         old_sv_F = saved_F
         old_B = B(1:it-mov-1,1:it-mov-1)
         deallocate(saved_err)
         deallocate(saved_F)
         deallocate(B)
      end if
      allocate(B(it-mov+1, it-mov+1))
      allocate(saved_err(it*n*(n+1)/2))
      allocate(saved_F(it*n*(n+1)/2))
      allocate(sol(it-mov+1,1))
      allocate(B_scaled(it-mov+1, it-mov+1))
      call save_vectors(F, error_v, n, it, old_sv_F, old_sv_err, saved_F, saved_err)
      deallocate(old_sv_err)
      deallocate(old_sv_F)

      !Constructing B (and reducing k via rescaling + omitting older itations)
      call build_B(error_v, saved_err, n, it, mov, old_B, B)
      call rescale_B(k,it+1,mov, B, B_scaled, sol)
      if (e_3 <= 0 .OR. size(B, dim=1) >= memory) then
         do
            if (e_3 <= 0) then
               if (it-mov <=min_dim) exit
               if (k < max_k) exit
            else
               if (size(B, dim=1) <= memory) exit
            end if
            if (allocated(old_B) .eqv. .TRUE.) deallocate(old_B)
            allocate(old_B(it-mov,it-mov))
            old_B = B(2:it+1-mov,2:it+1-mov)
            deallocate(B)
            deallocate(B_scaled)
            deallocate(sol)
            allocate(B(it-mov,it-mov))
            allocate(B_scaled(it-mov,it-mov))
            allocate(sol(it-mov,1))
            B = old_B
            deallocate(old_B)
            call rescale_B(k,it,mov, B,B_scaled, sol)
            mov = mov+1
         end do
      else
         deallocate(old_B)
      end if
      
      !Constructing F if below threshold
      if (start_diis .eqv. .TRUE.) call build_F(B_scaled,n,it+1,mov,saved_F,sol, F, info)
      deallocate(sol)
      deallocate(B_scaled)
      !Use for next round: saved_D, saved_err, B, mov, var, save_emax
   end subroutine diis
end module tblite_scf_mixer_diis
