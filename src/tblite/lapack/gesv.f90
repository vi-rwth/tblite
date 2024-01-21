module tblite_lapack_gesv
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: wrap_gesv

   interface wrap_gesv
      module procedure :: wrap_sgesv
      module procedure :: wrap_dgesv
   end interface wrap_gesv

   interface lapack_gesv
      pure subroutine sgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(inout) :: a(lda,n)
         real(sp), intent(inout) :: b(lda,nrhs)
         integer, intent(out) :: ipiv(n)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine sgesv
      pure subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(inout) :: a(lda,n)
         real(dp), intent(inout) :: b(lda,nrhs)
         integer, intent(out) :: ipiv(n)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dgesv
   end interface lapack_gesv

contains

pure subroutine wrap_sgesv(amat, bmat, info, ipiv)
   real(sp), intent(inout) :: amat(:,:)
   real(sp), intent(inout) :: bmat(:,:)
   integer, intent(out) :: info
   integer, intent(out) :: ipiv(:)
   integer :: ldb, lda, n, nrhs

   n = size(amat, 1)
   lda = max(1, n)
   ldb = max(1, n)
   nrhs = size(bmat, 2)
   call lapack_gesv(n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_sgesv

pure subroutine wrap_dgesv(amat, bmat, info, ipiv)
   real(dp), intent(inout) :: amat(:,:)
   real(dp), intent(inout) :: bmat(:,:)
   integer, intent(out) :: info
   integer, intent(out) :: ipiv(:)
   integer :: ldb, lda, n, nrhs

   n = size(amat, 1)
   lda = max(1, n)
   ldb = max(1, n)
   nrhs = size(bmat, 2)
   call lapack_gesv(n, nrhs, amat, lda, ipiv, bmat, ldb, info)
end subroutine wrap_dgesv

end module tblite_lapack_gesv
   
