module tblite_scf_mixer_diis
    use mctc_env, only : wp, error_type, fatal_error
    use tblite_lapack, only : getrf, getri
    use tblite_scf_mixer_type, only : mixer_type
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
        !> Apply mixing to the density
        procedure :: next
        !> Set Fock-matrix from 1D array
        procedure :: set_1d_F
        !> Set Density-Matrix from 1D array
        procedure :: set_1d
        !> Get Density-matrix as 1D array
        procedure :: get_1d
        !> Get error metric from mixing
        procedure :: get_error
    end type diis_mixer

    contains

    !> Initialize DIIS mixer
    subroutine new_diis(self, ndim, S, input)
        type(diis_mixer), intent(out) :: self
        integer, intent(in) :: ndim
        real(wp), intent(in) :: S(:,:)
        type(diis_input), intent(in) :: input

        self%ndim = ndim
        self%iter = 0
        self%memory = input%memory
        self%e_max = 0
        allocate(self%S(ndim,ndim))
        self%S = S
        allocate(self%err_v((ndim**2+ndim)/2))
        allocate(self%D(ndim*ndim))
        allocate(self%F(ndim*ndim))
    end subroutine new_diis

    !> Routines for accessing D and F

    !> Set new density from 1D array
    subroutine set_1d(self, d_1d)
        !> Instance of the mixer
        class(diis_mixer), intent(inout) :: self
        !> Density vector
        real(wp), intent(in) :: d_1d(:)
        self%D = d_1d
    end subroutine set_1d

    subroutine set_1d_F(self, f_1d)
        class(diis_mixer), intent(inout) :: self
        real(wp) :: f_1d(:)

        self%F = f_1d
    end subroutine set_1d_F
    
    !> Apply mixing to the density
    subroutine next(self, error)
        class(diis_mixer), intent(inout) :: self
        type(error_type), allocatable, intent(out) :: error
        integer :: info

        self%iter = self%iter + 1
        call diis(self%ndim, self%iter, self%D, self%e_max, self%err_v, self%memory, info)
        if (info /= 0) call fatal_error(error, "DIIS failed to obtain next iteration")
    end subroutine next
    
    !> Get density as 1D array
    subroutine get_1d(self, d_1d)
        class(diis_mixer), intent(inout) :: self
        real(wp), intent(out) :: d_1d(:)

        d_1d(:) = self%D
    end subroutine get_1d

    pure function get_error(self) result(error)
        class(diis_mixer), intent(inout) :: self
        integer :: indx
        real(wp) :: F_2d(self%ndim, self%ndim), D_2d(self%ndim, self%ndim)
        real(wp) :: tmp(self%ndim,self%ndim), tmp1(self%ndim,self%ndim)
        real(wp) :: tmp2(self%ndim,self%ndim), error
        real(wp) :: tmp3(self%ndim,self%ndim), tmp4(self%ndim,self%ndim)
        F_2d = reshape(self%F, shape(F_2d))
        D_2d = reshape(self%D, shape(D_2d))
        tmp = 0.0_wp
        tmp1 = 0.0_wp
        tmp2 = 0.0_wp
        tmp3 = 0.0_wp
        tmp4 = 0.0_wp
        !call dgemm('N', 'N', self%ndim, self%ndim, self%ndim, 1, D_2d, self%ndim, &
        !& F_2d, self%ndim, 0d0, tmp1, self%ndim)
        !call dgemm('N', 'N', self%ndim, self%ndim, self%ndim, 1, D_2d, self%ndim, &
        !& self%S, self%ndim, 0d0, tmp2, self%ndim)
        !call dgemm('N', 'N', self%ndim, self%ndim, self%ndim, 1, self%S, self%ndim, &
        !& tmp1, self%ndim, 0d0, tmp3, self%ndim)
        !call dgemm('N', 'N', self%ndim, self%ndim, self%ndim, 1, F_2d, self%ndim, &
        !& tmp2, self%ndim, 0d0, tmp4, self%ndim)
        !tmp = tmp3 - tmp4
        tmp=MATMUL(self%S,MATMUL(D_2d,F_2d))-MATMUL(F_2d,MATMUL(D_2d,self%S))
        self%e_max = abs(tmp(1,1))
        indx = 1
        do i=1,self%ndim
            do j=i,self%ndim
                self%err_v(indx) = tmp(j,i)
                if (abs(self%err_v(indx)) > self%e_max) self%e_max = abs(self%err_v(indx))
                indx=indx+1
            end do
        end do
        error = self%e_max
    end function get_error

    !> General linear algebra
    subroutine hadamard(vect1,vect2,n)
        integer, intent(in) :: n
        real(wp), intent(in) :: vect1(n)
        real(wp) :: prod_vect(n)
        real(wp), intent(inout) :: vect2(n)
        do i=1,n
            prod_vect(i) = vect1(i) * vect2(i)
        end do
        vect2 = prod_vect
    end subroutine

    subroutine frobenius(MA,MB,prod,n)
        integer, intent(in) :: n
        integer :: indx, cnt
        real(wp),intent(in), dimension((n**2+n)/2) :: MA,MB
        real(wp),intent(out) :: prod
        prod = 0
        cnt = 0
        indx = 1
        do i=1,n*(n+1)/2
            prod = prod + MA(i)*MB(i)
            if (i/=indx) then
                prod = prod + (MA(i)*MB(i))
            else
                indx = indx + n-cnt
                cnt = cnt+1
            end if
        end do
    end subroutine
    
    subroutine norm_i(MA, norm, n)
        integer :: n
        real(wp),intent(in), dimension(n, n) :: MA
        real(wp) :: norm, tmp
        norm = 0
        do i=1,n
            tmp = 0
            do j=1,n
                tmp = tmp + abs(MA(j,i))
            end do
            if (tmp>norm) norm = tmp
        end do
    end subroutine
    
    subroutine build_scale_matrix(MA, scale_matr, n)
        integer, intent(in) :: n
        real(wp), intent(in) :: MA(n,n)
        real(wp) :: tmp
        real(wp), intent(out) :: scale_matr(n)
        do i=1,n
            tmp = 0
            do j=1,n
                tmp = tmp + abs(MA(i,j))
            end do
            scale_matr(i) = 1/tmp
        end do
    end subroutine
    
    !> Specific to DIIS 
    subroutine build_D(B, n, it, omit, saved_D, sol, D, info)
        integer, intent(in) :: n, it, omit
        integer, intent(out) :: info
        real(wp), intent(in) :: B(it-omit,it-omit), saved_D(it*n**2)
        real(wp), intent(out) :: D(n*n)
        real(wp), dimension(it-omit) :: perm, sol

        call lapack_gesv(it-omit,1,B,it-omit, perm, sol,it-omit, info)
        if (info==0) then
            D = 0
            do i=1,it-omit-1
                D = D + sol(i)*saved_D((i+omit-1)*n**2+1:(i+omit)*n**2)
            end do
        end if
    end subroutine

    subroutine build_B(err_v, sv_err, n, it, omit, old_B, B)
        integer, intent(in) :: n, it, omit
        real(wp), intent(in) :: sv_err(it*(n**2+n)/2), err_v((n**2+n)/2), old_B(it-omit-1, it-omit-1)
        real(wp), intent(out) :: B(it-omit+1, it-omit+1)
        real(wp) :: Tr
        B=1d0
        B(it-omit+1, it-omit+1) = 0d0
        if (it /= 1) then
            B(1:it-omit-1, 1:it-omit-1) = old_B
            do i=1,it-omit-1
                call frobenius(sv_err((i+omit-1)*(n**2+n)/2+1:(i+omit)*(n**2+n)/2), err_v, Tr, n)
                B(it-omit, i) = Tr
                B(i, it-omit) = Tr
            end do
        end if
        call frobenius(err_v, err_v, Tr, n)
        B(it-omit, it-omit) = Tr
    end subroutine

    subroutine save_vectors(D, err_v, n, it, old_sv_D, old_sv_err, sv_D, sv_err)
        integer, intent(in) :: n, it
        real(wp), intent(in) :: D(n*n), err_v((n**2+n)/2)
        real(wp), intent(in) :: old_sv_D((it-1)*n**2), old_sv_err((it-1)*(n**2+n)/2)
        real(wp), intent(out) :: sv_D(it*n**2), sv_err(it*(n**2+n)/2)
        if (it/=1) then
            sv_err(1:(it-1)*(n**2+n)/2) = old_sv_err
            sv_D(1:(it-1)*n**2) = old_sv_D
        end if          
        sv_err((it-1)*(n**2+n)/2+1:it*(n**2+n)/2) = err_v
        sv_D((it-1)*n**2+1:it*n**2) = D
    end subroutine

    subroutine rescale_B(k,it,omit, B, B_scaled, sol)
        integer, intent(in) :: it, omit
        integer :: info
        real(wp), intent(in) :: B(it-omit,it-omit)
        real(wp) :: perm(it-omit), scale_matr(it-omit)
        real(wp) :: tmp_matr(it-omit,it-omit), work(3*(it-omit)-1), norm1, norm2
        real(wp), intent(out) :: sol(it-omit), B_scaled(it-omit,it-omit), k

        call build_scale_matrix(B,scale_matr, it-omit)
        do i=1,it-omit
            do j=1,it-omit
                tmp_matr(j,i) = scale_matr(j)*B(j,i)
            end do
        end do
        sol = 0d0
        sol(it-omit) = 1d0
        call hadamard(scale_matr, sol, it-omit)
        B_scaled = tmp_matr
        call getrf(tmp_matr, perm, info)
        call getri(tmp_matr, perm, info)
        call norm_i(B_scaled, norm1,it-omit)        
        call norm_i(tmp_matr, norm2, it-omit)
        k = norm1*norm2
    end subroutine

    subroutine diis(n, it, D, e_max, error_v, memory, info)
        !parameters
        integer :: max_k = 100000, min_dim = 5
        real(wp) :: diis_thresh = 1.0d0
        logical :: start_diis

        integer, intent(in) :: n, it, memory
        integer, intent(out) :: info
        real(wp), intent(inout) :: D(n,n), e_max
        real(wp) :: k, error_v((n**2+n)/2), e_3
        real(wp), allocatable :: old_sv_err(:), old_sv_D(:), B_scaled(:,:), old_B(:,:), sol(:)
        real(wp), allocatable, save :: B(:,:), saved_err(:), saved_D(:)
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

        print *,start_diis
        print *,e_3
        print *,e_max

        !Preparing for B construction
        allocate(old_sv_err((it-1)*(n**2+n)/2))
        allocate(old_sv_D((it-1)*n**2))
        allocate(old_B(it-1-mov,it-1-mov))
        if (it/=1) then
            old_sv_err = saved_err
            old_sv_D = saved_D
            old_B = B(1:it-mov-1,1:it-mov-1)
            deallocate(saved_err)
            deallocate(saved_D)
            deallocate(B)
        end if
        allocate(B(it-mov+1, it-mov+1))
        allocate(saved_err(it*(n**2+n)/2))
        allocate(saved_D(it*n**2))
        allocate(sol(it-mov+1))
        allocate(B_scaled(it-mov+1, it-mov+1))
        call save_vectors(D, error_v, n, it, old_sv_D, old_sv_err, saved_D, saved_err)
        deallocate(old_sv_err)
        deallocate(old_sv_D)

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
                allocate(sol(it-mov))
                B = old_B
                deallocate(old_B)
                call rescale_B(k,it,mov, B,B_scaled, sol)
                mov = mov+1
            end do
        else
            deallocate(old_B)
        end if
        
        !Constructing F if below threshold
        if (start_diis .eqv. .TRUE.) call build_D(B_scaled,n,it+1,mov,saved_D,sol, D, info)
        deallocate(sol)
        deallocate(B_scaled)
        !Use for next round: saved_D, saved_err, B, mov, var, save_emax
    end subroutine
end module tblite_scf_mixer_diis