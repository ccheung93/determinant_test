program detstr
    Use, Intrinsic :: iso_fortran_env, Only : dp => real64, int64
    implicit none

    integer :: n_orbitals
    integer, parameter :: bits_per_int = 32   ! bits per integer
    integer :: n_ints
    integer, dimension(:), allocatable :: bit_rep !array to store the packed bitstring
    integer :: i
    integer, dimension(:), allocatable :: bdet1, bdet2
    Integer, allocatable, dimension(:)   :: idet1, idet2, nfs

    integer :: Ne, Nd, Nsu, nf, nf2, j
    integer, dimension(:,:), allocatable :: Iarr, Iarr2, Barr ! array of int_rep determinants, array of bit_rep determinants

    integer :: start_time, end_time, clock_rate
    real :: total_time, total_time2
    logical :: testing = .false.

    Call system_clock(count_rate=clock_rate)

    ! Read list of determinants in integer representation
    call system_clock(start_time)
    print*, 'Enter number of electrons: '
    read(*,*) Ne
    Open (16,file='CONF.DET',status='UNKNOWN',form='UNFORMATTED')
    Read(16) Nd, n_orbitals
    print*, 'Ns=', n_orbitals
    allocate(Iarr(Ne,Nd), Iarr2(Ne,Nd))
    Iarr=0
    do i=1,Nd
       read(16) Iarr(1:Ne,i)
       if (maxval(Iarr(1:Ne,i)) > n_orbitals) n_orbitals = maxval(Iarr(1:Ne,i))
    end do
    call system_clock(end_time)

    print*, 'number of electrons: ', Ne
    print*, 'number of determinants: ', Nd
    print*, 'number of orbitals: ', n_orbitals
    print*, 'reading Iarr took', end_time - start_time, 'ms'

    ! number of integers needed to store all orbitals in bit representation
    n_ints = (n_orbitals + bits_per_int - 1) / bits_per_int
    allocate(Barr(n_ints,Nd))
    Barr=0

    ! compare integer representation to bit representation
    print*, 'number of ints used in integer representation:', Ne, ',', Ne*4, 'bytes'
    print*, 'number of ints used in bit representation:', n_ints, ',', n_ints*4, 'bytes'
    print*, 'reduction of ', real(Ne)/n_ints,'x'

    ! allocate arrays
    allocate(idet1(Ne), idet2(Ne))
    allocate(bdet1(n_ints), bdet2(n_ints))
    
    ! convert integer representation to bit representation
    call system_clock(start_time)
    do i=1,Nd
        idet1 = Iarr(1:Ne,i)
        call convert_int_rep_to_bit_rep(idet1, bdet1)
        Barr(1:n_ints, i) = bdet1
    end do
    call system_clock(end_time)
    print*, 'converting integer representation to bit representation took', end_time - start_time, 'ms'

    ! convert bit representation to integer representation
    call system_clock(start_time)
    do i=1,Nd
        bdet1 = Barr(1:n_ints, i)
        call convert_bit_rep_to_int_rep(bdet1, idet1)
        Iarr2(1:Ne,i) = idet1
    end do
    call system_clock(end_time)
    print*, 'converting bit representation to integer representation took', end_time - start_time, 'ms'

    ! check if converted integer representation is the same as the original
    if (all(Iarr == Iarr2)) then
        print*, 'success! Iarr == Iarr2'
    else
        print*, 'error! Iarr /= Iarr2'
    end if

    ! Comparison of determinants in integer representation
    call system_clock(start_time)
    idet1 = Iarr(1:Ne,1)
    allocate(nfs(Nd))
    print*, '===== comparison stage ======'

    do i=1,Nd
        idet2 = Iarr(1:Ne,i)
        call CompD(idet1, idet2, nf)
        nfs(i) = nf
    end do
    call system_clock(end_time)
    total_time = end_time - start_time
    print*, 'comparing integer representation took: ', total_time, 'ms'

    ! Comparison of determinants in bit representation
    call system_clock(start_time)
    bdet1 = Barr(1:n_ints,1)
    do j=1,Nd
        
        bdet2 = Barr(1:n_ints,j)
        nf = 0
        do i = 1, n_ints
            nf = nf + popcnt(ieor(bdet1(i), bdet2(i)))
        end do
        ! convert hamming distance to number of swaps/differences
        nf2 = nf/2 + mod(nf, 2)
        if (testing == .true. .and. nfs(j) /= nf2) then
            print*, '==============================='
            print*, j, nfs(j), nf, nf2, nfs(j) == nf
            idet2 = Iarr(1:Ne,j)
            call print_bit_rep(idet1, bdet1)
            call print_bit_rep(idet2, bdet2)
            print*, '==============================='
        end if
    end do
    call system_clock(end_time)
    total_time2 = end_time - start_time
    print*, 'comparing bit representation took: ', total_time2, 'ms'
    print*, 'speedup:', total_time/total_time2,'x'

contains

    subroutine print_bit_rep(int_rep, bit_rep)
        implicit none

        integer, dimension(:), allocatable :: int_rep, bit_rep
        integer :: i

        ! print resulting packed bit representation for each integer
        print*, 'int rep of determinant: '
        print*, int_rep
        do i=1, n_ints
            print*, 'packed bit rep in integer ', i, ': ', bit_rep(i)
            print*, 'binary rep of integer ', i, ': ', print_bits(bit_rep(i), bits_per_int)
        end do

    end subroutine print_bit_rep

    function print_bits(num, n_bits) result(bitstr)
        ! print bit string from bit representation of determinant
        implicit none
        integer, intent(in) :: num, n_bits
        character(len=n_bits) :: bitstr
        integer :: i

        do i=1,n_bits
            if (btest(num, i-1)) then
                bitstr(i:i) = '1'
            else
                bitstr(i:i) = '0'
            end if
        end do

    end function print_bits

    subroutine CompD(id1,id2,nf)
        ! this subroutine compares determinants and counts number of differences in orbitals
        ! return the number of differences nf
        implicit None
        integer  :: ni, nj, nf, i, j, l1, l2
        integer, allocatable, dimension(:)   :: id1, id2

        ni=0
        nj=0
        nf=3
        i=1
        j=1
  
        do while (i <= Ne .and. j <= Ne)
            l1=id1(i) ! id1(i) is the i-th element of determinant 1 (number of orbital occupied)
            l2=id2(j) ! id1(j) is the j-th element of determinant 2 (number of orbital occupied)
            if (l1 == l2) then
                i=i+1
                j=j+1
            else if (l1 > l2) then
                nj=nj+1
                if (nj > 2) return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                j=j+1
            else
                ni=ni+1
                if (ni > 2) return ! If difference > 2 Then matrix element between id1 and id2 will be 0
                i=i+1
            end if
        end do
        nf = max(ni,nj)
        return
    end subroutine CompD

    subroutine convert_int_rep_to_bit_rep(int_rep, bit_rep)
        implicit none
        integer, dimension(:), allocatable, intent(in) :: int_rep ! integer representation
        integer, dimension(:), allocatable :: bit_rep ! bit representation
        integer :: i, index, bit_position ! loop index, integer index, bit position

        ! initialize bit representation array to 0 (all orbitals unoccupied)
        bit_rep = 0

        ! loop over occupied orbitals and set the corresponding bits
        do i=1, size(int_rep)
            ! find which integer the orbital falls in
            index = (int_rep(i) - 1) / bits_per_int + 1
    
            ! find which bit position within the integer the orbital falls in
            bit_position = mod(int_rep(i) - 1, bits_per_int)
    
            ! set the bit corresponding to the orbital
            bit_rep(index) = ibset(bit_rep(index), bit_position)
        end do

    end subroutine convert_int_rep_to_bit_rep

    subroutine convert_bit_rep_to_int_rep(bit_rep, int_rep)
        implicit none
        integer, dimension(:), allocatable, intent(in) :: bit_rep ! bit representation
        integer, dimension(:), allocatable :: int_rep ! integer representation
        integer :: i, j, cnt, mask, index, bit_position ! loop index, integer index, bit position
        character(len=bits_per_int) :: bit_int

        ! initialize int representation array to 0
        int_rep = 0

        ! loop over occupied orbitals and set the corresponding integer positions
        cnt = 1
        do i=1, size(bit_rep)
            ! we can skip '0' elements of the bit_rep since they correspond to no occupancy 
            if (bit_rep(i) == 0) cycle
            bit_int = print_bits(bit_rep(i), bits_per_int)
            do j = 1, bits_per_int
                if (bit_int(j:j) == '1') then
                    int_rep(cnt) = j + (i-1)*bits_per_int
                    cnt = cnt + 1
                end if
            end do
        end do

    end subroutine convert_bit_rep_to_int_rep

end program detstr