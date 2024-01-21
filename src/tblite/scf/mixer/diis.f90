module tblite_scf_mixer_diis
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_scf_mixer_type, only : mixer_type
   use tblite_lapack, only : getrf, getri, gesv
   use tblite_blas, only: gemm
   private

   public :: new_diis

   integer, parameter :: max_k = 100000, min_dim = 5
   real(wp), parameter :: diis_thresh = 0.5_wp

   type, public :: diis_input
      !> Number of steps to keep in memory
      integer :: memory
      integer :: mode
   end type diis_input

   type, public, extends(mixer_type) :: diis_mixer
      integer :: ndim
      integer :: memory
      integer :: mode
      integer :: iter
      logical :: first_cycle
      logical :: first_occurrence
      logical :: start_diis
      real(wp) :: e_max
      real(wp), allocatable :: err_v(:)
      real(wp), allocatable :: F(:)
      real(wp), allocatable :: D(:)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: S(:,:)
      !> For main DIIS routine
      integer :: mov
      integer :: last_mov
      integer :: var
      real(wp) :: save_emax(3)
      real(wp), allocatable :: saved_err(:)
      real(wp), allocatable :: saved_mat(:)
      real(wp), allocatable :: B(:,:)
   contains
      !> Apply mixing
      procedure :: next
      !> Set Fock-matrix from 1D array
      procedure :: set_1d_F
      !> Set density matrix from 1D array
      procedure :: set_1d_D
      !> Set density from 1D array (only for diis_d)
      procedure :: set_1d
      !> Get Fock-Matrix as 1D array (only for diis_f)
      procedure :: get_1d_F
      !> Get density from 1D array (only for diis_d)
      procedure :: get_1d
      !> Get error metric from mixing
      procedure :: get_error
      !> Get error for first run
      procedure :: diff_1d
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
      self%e_max = 0
      self%first_cycle = .true.
      self%first_occurrence = .false.
      self%start_diis = .false.
      self%memory = input%memory
      self%mode = input%mode
      allocate(self%S(ndim,ndim))
      self%S = S
      allocate(self%err_v(ndim*(ndim+1)/2))
      allocate(self%D(ndim*ndim))
      allocate(self%F(ndim*(ndim+1)/2))
      !> For main DIIS routine
      self%mov = 0
      self%last_mov = 0
      self%var = 1
   end subroutine new_diis

   !> Routines for accessing D and F

   !> Set new density-matrix from 1D array
   subroutine set_1d_D(self, d_1d)
      !> Instance of the mixer
      class(diis_mixer), intent(inout) :: self
      !> Density vector
      real(wp), intent(in) :: d_1d(:)
      if ((self%first_cycle .eqv. .false.) &
         & .and. (self%first_occurrence .eqv. .false.)) self%D = d_1d
   end subroutine set_1d_D

   !> Set new Fock-matrix from 1D array
   subroutine set_1d_F(self, f_1d)
      !> Instance of the mixer
      class(diis_mixer), intent(inout) :: self
      !> Fock-matrix
      real(wp), intent(in) :: f_1d(:)
      integer :: indx
      indx = 0
      if ((self%first_cycle .eqv. .false.) &
         & .and. (self%first_occurrence .eqv. .false.)) then
         do i=1,self%ndim
            do j=i,self%ndim
               indx = indx + 1
               self%F(indx) = f_1d(j+self%ndim*(i-1))
            end do
         end do
         call get_emax(self)
      end if
   end subroutine set_1d_F

   !> Set new density from 1D array (only for diis_d)
   subroutine set_1d(self, qvec)
      !> Instance of the mixer
      class(diis_mixer), intent(inout) :: self
      real(wp), allocatable :: lastq(:)
      !> Density vector
      real(wp), intent(in) :: qvec(:)
      if ((self%first_cycle .eqv. .true.) &
         & .or. (self%first_occurrence .eqv. .true.)) then
         if (.not. allocated(self%q)) then
            allocate(self%q(size(qvec)))
            self%q = qvec
         else
            allocate(lastq(size(self%q)))
            lastq = self%q
            deallocate(self%q)
            allocate(self%q(size(lastq)+size(qvec)))
            self%q(1:size(lastq)) = lastq
            deallocate(lastq)
            self%q(size(lastq)+1:size(lastq)+size(qvec)) = qvec
         end if
      end if
   end subroutine set_1d

   !> Apply mixing 
   subroutine next(self, error)
      class(diis_mixer), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: info
      if (self%start_diis .eqv. .true.) then
         self%iter = self%iter + 1
         select case(self%mode)
         case(1)
            call diis(self, self%F, info)
         case(2)
            call diis(self, self%q, info)
            self%first_occurrence = .false.
         end select
         if (info /= 0) call fatal_error(error, "DIIS failed to obtain next iteration")
      end if
   end subroutine next

   !> Get Fock-matrix as 1D array (for diis_f)
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
   
   !> Get density as 1D array (for diis_d)
   subroutine get_1d(self, qvec)
      class(diis_mixer), intent(inout) :: self
      real(wp), intent(out) :: qvec(:)
      real(wp), allocatable :: truncq(:)
      if (size(qvec) < size(self%q)) then
         allocate(truncq(size(self%q)-size(qvec)))
         truncq = self%q(size(qvec)+1:size(self%q))
         qvec = self%q(1:size(qvec))
         deallocate(self%q)
         allocate(self%q(size(truncq)))
         self%q = truncq
         deallocate(truncq)
      else
         qvec = self%q
         deallocate(self%q)
      end if
   end subroutine get_1d

   !> Get error for first run
   subroutine diff_1d(self, qvec)
      class(diis_mixer), intent(inout) :: self
      real(wp), intent(in) :: qvec(:)
      real(wp), allocatable :: lastq(:)
      real(wp), allocatable, save :: old_q(:)
      if (self%first_cycle .eqv. .true.) then
         if (.not. allocated(old_q)) then
            allocate(old_q(size(self%q)))
            old_q = self%q
            deallocate(self%q)
            allocate(self%q(size(qvec)))
            self%q = qvec
         else
            allocate(lastq(size(self%q)))
            lastq = self%q
            deallocate(self%q)
            allocate(self%q(size(lastq)+size(qvec)))
            self%q(1:size(lastq)) = lastq
            deallocate(lastq)
            self%q(size(lastq)+1:size(lastq)+size(qvec)) = qvec
            if (size(self%q) == size(old_q)) then
               self%e_max = 0.0_wp
               do i = 1, size(self%q)
                  self%e_max = self%e_max + (self%q(i) - old_q(i))**2 / size(self%q)
               end do
               self%e_max = sqrt(self%e_max)
               deallocate(self%q)
               deallocate(old_q)
               self%first_cycle = .false.
            end if   
         end if
      else if (self%mode == 2) then
         self%first_occurrence = .true.
      end if
   end subroutine diff_1d

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
      real(wp) :: tmp3(self%ndim,self%ndim)
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
      call gemm(self%S, tmp1, D_2d, beta=0.0_wp)
      call gemm(F_2d, tmp2, tmp1, beta=0.0_wp)
      tmp = D_2d - tmp1
      self%e_max = abs(tmp(1,1))
      indx = 1
      do i=1,self%ndim
         do j=i,self%ndim
            self%err_v(indx) = tmp(j,i)
            if (abs(self%err_v(indx)) > self%e_max) self%e_max = abs(self%err_v(indx))
            indx=indx+1
         end do
      end do
      if (self%e_max < diis_thresh) self%start_diis = .true.
   end subroutine get_emax

   !> General linear algebra
   subroutine hadamard(vect1,vect2)
      real(wp), intent(in) :: vect1(:)
      real(wp), allocatable :: prod_vect(:,:)
      real(wp), intent(inout) :: vect2(:,:)
      allocate(prod_vect(size(vect1),1))
      do i=1,size(vect1)
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
   subroutine build_mat(B, it, omit, saved_mat, sol, mat, info)
      integer, intent(in) :: it, omit
      integer, intent(out) :: info
      integer :: perm(it-omit)
      real(wp), intent(inout) :: B(:,:), sol(:,:)
      real(wp), intent(in) :: saved_mat(:)
      real(wp), intent(out) :: mat(:)

      call call_lapack_gesv(B, sol, info, perm)
      if (info==0) then
         mat = 0.0_wp
         do i=1,it-omit-1
            mat = mat + sol(i,1)*saved_mat((i-1)*size(mat)+1:i*size(mat))
         end do
      end if
   end subroutine build_mat

   subroutine build_B(err_v, sv_err, it, omit, old_B, B)
      integer, intent(in) :: it, omit
      real(wp), intent(in) :: sv_err(:), err_v(:), old_B(:, :)
      real(wp), intent(out) :: B(:, :)
      real(wp) :: Tr
      B=1.0_wp
      B(it-omit+1, it-omit+1) = 0.0_wp
      if (it /= 1) then
         B(1:it-omit-1, 1:it-omit-1) = old_B
         do i=1,it-omit-1
            call frobenius(sv_err((i-1)*size(err_v)+1:i*size(err_v)), err_v, Tr)
            B(it-omit, i) = Tr
            B(i, it-omit) = Tr
         end do
      end if
      call frobenius(err_v, err_v, Tr)
      B(it-omit, it-omit) = Tr
   end subroutine build_B

   subroutine save_vectors(mat, err_v, it, last_sv_mat, last_sv_err, sv_mat, sv_err, omit)
      integer, intent(in) :: it, omit
      real(wp), intent(in) :: mat(:), err_v(:), last_sv_mat(:), last_sv_err(:)
      real(wp), intent(out) :: sv_mat(:), sv_err(:)
      if (it/=1) then
         sv_err(1:(it-omit-1)*size(err_v)) = last_sv_err
         sv_mat(1:(it-omit-1)*size(mat)) = last_sv_mat
      end if        
      sv_err((it-omit-1)*size(err_v)+1:(it-omit)*size(err_v)) = err_v
      sv_mat((it-omit-1)*size(mat)+1:(it-omit)*size(mat)) = mat
   end subroutine save_vectors

   subroutine rescale_B(k,it,omit, B, B_scaled, sol)
      integer, intent(in) :: it, omit
      integer :: info, perm(it-omit)
      real(wp), intent(in) :: B(:,:)
      real(wp) :: scale_matr(it-omit)
      real(wp) :: tmp_matr(it-omit,it-omit), norm1, norm2
      real(wp), intent(out) :: sol(:,:), B_scaled(:,:), k

      call build_scale_matrix(B,scale_matr)
      do i=1,it-omit
         do j=1,it-omit
            tmp_matr(j,i) = scale_matr(j)*B(j,i)
         end do
      end do
      sol = 0.0_wp
      sol(it-omit,1) = 1.0_wp
      call hadamard(scale_matr, sol)
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

   subroutine call_lapack_gesv(a, b, info, ipiv)
      real(wp), intent(inout) :: a(:,:)
      real(wp), intent(inout) :: b(:,:)
      integer, intent(out) :: info
      integer, intent(out) :: ipiv(:)
      
      call gesv(a, b, info, ipiv)
   end subroutine call_lapack_gesv

   subroutine diis(self, mat, info)
      class(diis_mixer), intent(inout) :: self
      integer, intent(out) :: info
      real(wp), intent(inout) :: mat(:)
      real(wp) :: k, e_3
      real(wp), allocatable :: last_sv_err(:), last_sv_mat(:), B_scaled(:, :), old_B(:,:), sol(:,:)
      
      info = 0
      !> Calculating e3
      select case(self%var)
      case(1)
         self%save_emax(1) = self%e_max
         self%var = 2
         if (self%iter>2) then
            e_3 = self%save_emax(2)-self%save_emax(1)
         else
            e_3 = 1
         end if
      case(2)
         self%save_emax(2) = self%e_max
         self%var = 3
         if (self%iter>2) then
            e_3 = self%save_emax(3)-self%save_emax(2)
         else
            e_3 = 1
         end if
      case(3)
         self%save_emax(3) = self%e_max
         self%var = 1
         e_3 = self%save_emax(1)-self%save_emax(3)
      end select

      !Preparing for B construction
      allocate(last_sv_mat((self%iter-1-self%mov)*size(mat)))
      allocate(last_sv_err((self%iter-1-self%mov)*size(self%err_v)))
      allocate(old_B(self%iter-1-self%mov,self%iter-1-self%mov))
      if (self%iter/=1) then
         last_sv_err = self%saved_err(1+(self%mov-self%last_mov)*size(self%err_v) &
            & :(self%iter-1-self%last_mov)*size(self%err_v))
         last_sv_mat = self%saved_mat(1+(self%mov-self%last_mov)*size(mat) &
            & :(self%iter-1-self%last_mov)*size(mat))
         old_B = self%B(1:self%iter-self%mov-1,1:self%iter-self%mov-1)
         deallocate(self%saved_err)
         deallocate(self%saved_mat)
         deallocate(self%B)
      end if
      allocate(self%B(self%iter-self%mov+1, self%iter-self%mov+1))
      allocate(self%saved_err((self%iter-self%mov)*size(self%err_v)))
      allocate(self%saved_mat((self%iter-self%mov)*size(mat)))
      allocate(sol(self%iter-self%mov+1,1))
      allocate(B_scaled(self%iter-self%mov+1, self%iter-self%mov+1))
      call save_vectors(mat, self%err_v, self%iter, last_sv_mat, last_sv_err, &
         & self%saved_mat, self%saved_err, self%mov)
      deallocate(last_sv_err)
      deallocate(last_sv_mat)

      !Constructing B (and reducing k via rescaling + omitting older itations)
      call build_B(self%err_v, self%saved_err, self%iter, self%mov, old_B, self%B)
      call rescale_B(k,self%iter+1,self%mov, self%B, B_scaled, sol)
      self%last_mov = self%mov
      if (e_3 <= 0 .or. size(self%B, 1) >= self%memory) then
         do
            if (e_3 <= 0) then
               if (self%iter-self%mov <= min_dim .or. k < max_k) exit
            else
               if (size(self%B, 1) <= self%memory) exit
            end if
            if (allocated(old_B) .eqv. .true.) deallocate(old_B)
            allocate(old_B(self%iter-self%mov,self%iter-self%mov))
            old_B = self%B(2:self%iter+1-self%mov,2:self%iter+1-self%mov)
            deallocate(self%B)
            deallocate(B_scaled)
            deallocate(sol)
            allocate(self%B(self%iter-self%mov,self%iter-self%mov))
            allocate(B_scaled(self%iter-self%mov,self%iter-self%mov))
            allocate(sol(self%iter-self%mov,1))
            self%B = old_B
            deallocate(old_B)
            call rescale_B(k,self%iter,self%mov, self%B,B_scaled, sol)
            self%mov = self%mov+1
         end do
      else
         deallocate(old_B)
      end if
      
      !Mixing mat
      call build_mat(B_scaled,self%iter+1,self%mov,self%saved_mat,sol, mat, info)
      deallocate(sol)
      deallocate(B_scaled)

   end subroutine diis
end module tblite_scf_mixer_diis
