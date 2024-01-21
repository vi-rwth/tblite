module test_double_dictionary
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_double_dictionary, only : double_dictionary_type
   implicit none
   private

   public :: collect_double_dictionary

contains

subroutine collect_double_dictionary(testsuite)

   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("add entries", test_add_entries), &
      new_unittest("get entries by label", test_get_entries_by_label), &
      new_unittest("get entries by index", test_get_entries_by_index), &
      new_unittest("get label from index", test_get_label_from_index), &
      new_unittest("test = assigment operator", test_assigment_operator), &
      new_unittest("test + operator", test_addition_operator), &
      new_unittest("test remove functionality by index", test_removal_index), &
      new_unittest("test remove functionality by label", test_removal_label), &
      new_unittest("access entries with invalid label or index", test_invalid_label_index), &
      new_unittest("access valid entries with wrong size of array as return", test_invalid_array_size), &
      new_unittest("compare index and label lookup", test_equivalence_index_label_lookup), &
      new_unittest("initialize labels", test_initialize_labels), &
      new_unittest("update entries", test_update_entries_label) &
      ]
end subroutine collect_double_dictionary

subroutine test_removal_index(error)
   type(error_type), allocatable, intent(out) :: error
   character(len=:), allocatable :: label1, label2
   real(wp), allocatable :: array1(:), array2(:, :)

   type(double_dictionary_type) :: dict, ini_dict

   call fill_test_dict(dict)
   ini_dict = dict
   call dict%remove_entry(3)
   call check(error, dict%get_n_entries(), 2)

   call dict%remove_entry(3)
   call check(error, dict%get_n_entries(), 2)
   call dict%get_label(1, label1)
   call dict%get_label(2,label2)

   call check(error, (label1 == "test1"))
   call check(error, (label1 == "test2"))
   call dict%get_entry("test1", array1)
   call dict%get_entry("test2", array2)

   call check(error, sum(ini_dict%record(1)%array1 - array1), 0.0_wp)
   call check(error, sum(ini_dict%record(2)%array2 - array2), 0.0_wp)

   call dict%remove_entry(1)
   call check(error, dict%get_n_entries(), 1)

   call dict%get_entry("test2", array2)
   call check(error, sum(ini_dict%record(2)%array2 - array2), 0.0_wp) 

end subroutine

subroutine test_removal_label(error)
   type(error_type), allocatable, intent(out) :: error
   character(len=:), allocatable :: label1, label2
   real(wp), allocatable :: array1(:), array2(:, :)

   type(double_dictionary_type) :: dict, ini_dict

   call fill_test_dict(dict)
   ini_dict = dict
   call dict%remove_entry("test3")
   call check(error, dict%get_n_entries(), 2)

   call dict%remove_entry("test3")
   call check(error, dict%get_n_entries(), 2)
   call dict%get_label(1, label1)
   call dict%get_label(2,label2)

   call check(error, (label1 == "test1"))
   call check(error, (label1 == "test2"))
   call dict%get_entry("test1", array1)
   call dict%get_entry("test2", array2)

   call check(error, sum(ini_dict%record(1)%array1 - array1), 0.0_wp)
   call check(error, sum(ini_dict%record(2)%array2 - array2), 0.0_wp)

   call dict%remove_entry("test1")
   call check(error, dict%get_n_entries(), 1)

   call dict%get_entry("test2", array2)
   call check(error, sum(ini_dict%record(2)%array2 - array2), 0.0_wp) 

end subroutine

subroutine test_initialize_labels(error)

   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict

   call dict%initialize_entry("test")

   call dict%initialize_entry("test1", 4)

   call dict%initialize_entry("test2", 4, 6)

   call dict%initialize_entry("test3", 4, 6, 9)

   call check(error, size(dict%record), 4)
   !Checking records directly since lookup will be tested in another test
   call check(error, (.not.(allocated(dict%record(1)%array1))))
   call check(error, (.not.(allocated(dict%record(1)%array2))))
   call check(error, (.not.(allocated(dict%record(1)%array3))))

   call check(error, actual = size(dict%record(2)%array1), expected = 4)
   call check(error, actual = sum(dict%record(2)%array1), expected = 0.0_wp)

   call check(error, actual = size(dict%record(3)%array2, dim = 1), expected = 4)
   call check(error, actual = size(dict%record(3)%array2, dim = 2), expected = 6)
   call check(error, actual = sum(dict%record(3)%array2), expected = 0.0_wp)

   call check(error, actual = size(dict%record(4)%array3, dim = 1), expected = 4)
   call check(error, actual = size(dict%record(4)%array3, dim = 2), expected = 6)
   call check(error, actual = size(dict%record(4)%array3, dim = 2), expected = 9)
   call check(error, actual = sum(dict%record(4)%array3), expected = 0.0_wp)

end subroutine

subroutine test_add_entries(error)

   type(error_type), allocatable, intent(out) :: error
   type(double_dictionary_type) :: dict

   call fill_test_dict(dict)

   call check(error, size(dict%record), 3)
   call check(error, actual = size(dict%record(1)%array1), expected = 4)
   call check(error, actual = sum(dict%record(1)%array1), expected = 4.0_wp)

   call check(error, actual = size(dict%record(2)%array2, dim = 1), expected = 4)
   call check(error, actual = size(dict%record(2)%array2, dim = 2), expected = 6)
   call check(error, actual = sum(dict%record(2)%array2), expected = (4*6*2.0_wp))

   call check(error, actual = size(dict%record(3)%array3, dim = 1), expected = 4)
   call check(error, actual = size(dict%record(3)%array3, dim = 2), expected = 6)
   call check(error, actual = size(dict%record(3)%array3, dim = 2), expected = 9)
   call check(error, actual = sum(dict%record(3)%array3), expected = 0.0_wp)

end subroutine

subroutine test_get_entries_by_label(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict

   real(wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   real(wp), allocatable :: array1_in(:), array2_in(:, :), array3_in(:, :, :)

   allocate(array1(4), source = 1.0_wp)
   allocate(array2(4, 6), source = 2.0_wp)
   allocate(array3(4, 6, 9), source = 0.0_wp)

   array1_in = array1
   array2_in = array2
   array3_in = array3

   call dict%add_entry("test1", array1)

   call dict%add_entry("test2", array2)

   call dict%add_entry("test3", array3)

   call dict%get_entry("test1", array1)
   call dict%get_entry("test2", array2)
   call dict%get_entry("test3", array3)

   call check(error, actual = sum(array1 - array1_in), expected = 0.0_wp)
   call check(error, actual = sum(array2 - array2_in), expected = 0.0_wp)
   call check(error, actual = sum(array3 - array3_in), expected = 0.0_wp)

end subroutine

subroutine test_get_entries_by_index(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict

   real(wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   real(wp), allocatable :: array1_in(:), array2_in(:, :), array3_in(:, :, :)

   allocate(array1(4), source = 1.0_wp)
   allocate(array2(4, 6), source = 2.0_wp)
   allocate(array3(4, 6, 9), source = 0.0_wp)

   array1_in = array1
   array2_in = array2
   array3_in = array3

   call dict%add_entry("test1", array1)

   call dict%add_entry("test2", array2)

   call dict%add_entry("test3", array3)

   call dict%get_entry(1, array1)
   call dict%get_entry(2, array2)
   call dict%get_entry(3, array3)

   call check(error, actual = sum(array1 - array1_in), expected = 0.0_wp)
   call check(error, actual = sum(array2 - array2_in), expected = 0.0_wp)
   call check(error, actual = sum(array3 - array3_in), expected = 0.0_wp)

end subroutine

subroutine test_equivalence_index_label_lookup(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict

   real(wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)

   real(wp), allocatable :: array1_in(:), array2_in(:, :), array3_in(:, :, :)

   call fill_test_dict(dict)

   call dict%get_entry("test1", array1)
   call dict%get_entry("test2", array2)
   call dict%get_entry("test3", array3)

   call dict%get_entry(1, array1_in)
   call dict%get_entry(2, array2_in)
   call dict%get_entry(3, array3_in)

   call check(error, actual = sum(array1 - array1_in), expected = 0.0_wp)
   call check(error, actual = sum(array2 - array2_in), expected = 0.0_wp)
   call check(error, actual = sum(array3 - array3_in), expected = 0.0_wp)

end subroutine

subroutine test_invalid_label_index(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict
   real(wp), allocatable :: array(:), array2(:, :), array3(:, :, :)

   call fill_test_dict(dict)

   call dict%get_entry("label4", array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(4, array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(0, array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(-1, array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry("label4", array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry(4, array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry(0, array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry(-1, array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry("label4", array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry(4, array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry(0, array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry(-1, array3)
   call check(error, (.not.allocated(array3)))

end subroutine

subroutine test_invalid_array_size(error)
   type(error_type), allocatable, intent(out) :: error

   type(double_dictionary_type) :: dict
   real(wp), allocatable :: array(:), array2(:, :), array3(:, :, :)

   call fill_test_dict(dict)

   call dict%get_entry("test1", array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry(1, array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry(1, array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry("test1", array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry("test2", array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(2, array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(2, array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry("test2", array3)
   call check(error, (.not.allocated(array3)))

   call dict%get_entry("test3", array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(3, array)
   call check(error, (.not.allocated(array)))

   call dict%get_entry(3, array2)
   call check(error, (.not.allocated(array2)))

   call dict%get_entry("test3", array2)
   call check(error, (.not.allocated(array2)))

end subroutine

subroutine test_get_label_from_index(error)
   type(error_type), allocatable, intent(out) :: error
   type(double_dictionary_type) :: dict
   character(len=:), allocatable :: name1

   call fill_test_dict(dict)
   call dict%get_label(1, name1)
   call check(error, (name1 == "test1"))
   call dict%get_label(2, name1)
   call check(error, (name1 == "test2"))
   call dict%get_label(3, name1)
   call check(error, (name1 == "test3"))
end subroutine

subroutine fill_test_dict(dict)
   type(double_dictionary_type), intent(inout) :: dict
   real(wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   allocate(array1(4), source = 1.0_wp)
   allocate(array2(4, 6), source = 2.0_wp)
   allocate(array3(4, 6, 9), source = 0.0_wp)

   call dict%add_entry("test1", array1)

   call dict%add_entry("test2", array2)

   call dict%add_entry("test3", array3)
end subroutine

subroutine test_assigment_operator(error)
   type(error_type), allocatable, intent(out) :: error
   type(double_dictionary_type) :: dict1, dict2
   call fill_test_dict(dict1)
   write(*,*) dict1%get_n_entries()
   write(*,*) size(dict1%record, dim=1)
   dict2 = dict1

   call check(error, dict1%n , dict2%n)

   call check(error, sum(dict1%record(1)%array1 - dict2%record(1)%array1), 0.0_wp)

   call check(error, sum(dict1%record(2)%array2 - dict2%record(2)%array2), 0.0_wp)

   call check(error, sum(dict1%record(3)%array3 - dict2%record(3)%array3), 0.0_wp)

end subroutine

subroutine test_addition_operator(error)
   type(error_type), allocatable, intent(out) :: error
   type(double_dictionary_type) :: dict1, dict2, dict3, dict4 
   character(len=:), allocatable :: label1, label2
   real(wp), allocatable :: array(:), array1(:), array2(:)
   call fill_test_dict(dict1)
   array = [0.0_wp, 0.0_wp]
   call dict2%add_entry("dict2", array)
   call dict2%add_entry("dict3", array)

   dict3 = dict1 + dict2
   !Check that total number of entries matches up
   call check(error, actual = dict3%get_n_entries(), expected = (dict1%get_n_entries()+ dict2%get_n_entries()))
   !Check Ordering LHS should vom before RHS
   call dict3%get_label(2, label1)
   call dict1%get_label(2, label2)
   call check(error, (label1 == label2))

   call dict3%get_label(5, label1)
   call dict2%get_label(2, label2)
   call check(error, (label1 == label2))

   call dict3%get_entry(1, array1)
   call dict1%get_entry(1, array2)
   call check(error, (sum(array1) == sum(array2)))

   call dict3%get_entry("dict3", array1)
   call dict2%get_entry(2, array2)
   call check(error, (sum(array1) == sum(array2)))
   
   dict3 = dict4 + dict1

   call check(error, (dict3%get_n_entries() == dict1%get_n_entries()))

   call dict4%concatenate(dict1)

   call check(error, (dict4%get_n_entries() == dict1%get_n_entries()))
   call dict4%get_entry(1, array1)
   call dict1%get_entry(1, array2)
   call check(error, (sum(array1) == sum(array2)))
end subroutine



subroutine test_update_entries_label(error)
   type(error_type), allocatable, intent(out) :: error
   type(double_dictionary_type) :: dict1, dict2, dict3 
   character(len=:), allocatable :: label1, label2
   real(wp), allocatable :: array(:), array1(:), array2(:, :), array3(:, :, :)
   call fill_test_dict(dict1)
   allocate(array(1), source= 0.0_wp)
   call dict1%update("test1", array)
   call dict1%get_entry("test1", array1)

   call check(error, (sum(array) == sum(array1)))

   call dict1%update("test2", array)
   call dict1%get_entry("test2", array1)

   call check(error, (sum(array) == sum(array1)))

   call dict1%update("test3", array)
   call dict1%get_entry("test3", array1)

   call check(error, (sum(array) == sum(array1)))

   allocate(array2(2, 2), source= 0.0_wp)

   call fill_test_dict(dict2)

   call dict2%update("test1", array2)
   call dict2%update("test2", array2)
   call dict2%update("test3", array2)

   deallocate(array2)

   call dict2%get_entry("test1", array2)
   call check(error, (size(array2, dim=1) == 2))
   call dict2%get_entry("test2", array2)
   call check(error, (size(array2, dim=1) == 2))
   call dict2%get_entry("test3", array2)
   call check(error, (size(array2, dim=1) == 2))

   allocate(array3(3, 2, 1), source= 1.0_wp)
   call fill_test_dict(dict3)

   call dict3%update("test1", array3)
   call dict3%update("test2", array3)
   call dict3%update("test3", array3)

   deallocate(array2)

   call dict3%get_entry("test1", array3)
   call check(error, (size(array3, dim=1) == 3))
   call dict3%get_entry("test2", array3)
   call check(error, (size(array3, dim=1) == 3))
   call dict3%get_entry("test3", array2)
   call check(error, (size(array3, dim=1) == 3))

end subroutine

end module test_double_dictionary
