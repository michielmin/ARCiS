module mod_references
  implicit none

  !-------------------------------------------------------------------
  ! Reference type constants
  !-------------------------------------------------------------------
  integer, parameter :: REF_OPACITY    = 1
  integer, parameter :: REF_REFRACTIVE = 2
  integer, parameter :: REF_SURFACE    = 3
  integer, parameter :: REF_METHOD     = 4
  integer, parameter :: REF_CHEMISTRY  = 5

  integer, parameter :: NREF_TYPES = 5

  character(len=24), parameter :: REF_TYPE_LABEL(NREF_TYPES) = [ &
    'Molecular opacities     ', &
    'Refractive indices      ', &
    'Surface reflectance     ', &
    'Methods                 ', &
    'Chemistry               '  &
  ]

  character(len=12), parameter :: REF_TYPE_KEY(NREF_TYPES) = [ &
    'opacity     ', &
    'refrind     ', &
    'surface     ', &
    'method      ', &
    'chemistry   '  &
  ]

  !-------------------------------------------------------------------
  ! Single reference record
  !-------------------------------------------------------------------
  type :: t_reference
    character(len=48)  :: key        ! internal lookup key, e.g. "H2O_Polyansky2018"
    integer            :: ref_type   ! one of REF_* constants above
    character(len=192)  :: adslink
    character(len=64)  :: bibtex_key ! used in \citealt{}
  end type t_reference

  !-------------------------------------------------------------------
  ! Dynamic list of references
  !-------------------------------------------------------------------
  type :: t_reflist
    type(t_reference), allocatable :: refs(:)
    integer :: n = 0
  contains
    procedure :: init    => reflist_init
    procedure :: add     => reflist_add
    procedure :: has_key => reflist_has_key
    procedure :: clear   => reflist_clear
  end type t_reflist

  type(t_reflist) :: master_refs
  type(t_reflist) :: active_refs

contains

  !clean helper routine to call from ARCiS
  subroutine references_init(database_file)
    character(len=*), intent(in) :: database_file
    call load_reference_database(database_file, master_refs)
    call active_refs%init()
  end subroutine references_init



  !===================================================================
  ! t_reflist type-bound procedures
  !===================================================================

  subroutine reflist_init(self, capacity)
    class(t_reflist), intent(inout) :: self
    integer, intent(in), optional   :: capacity
    integer :: cap
    cap = 64
    if (present(capacity)) cap = capacity
    if (allocated(self%refs)) deallocate(self%refs)
    allocate(self%refs(cap))
    self%n = 0
  end subroutine reflist_init

  !-------------------------------------------------------------------

  subroutine reflist_add(self, ref)
    class(t_reflist), intent(inout) :: self
    type(t_reference), intent(in)   :: ref
    type(t_reference), allocatable  :: tmp(:)
    integer :: newcap

    ! Silently skip duplicates
    if (self%has_key(ref%key)) return

    ! Grow backing array if full
    if (self%n == size(self%refs)) then
      newcap = 2 * size(self%refs)
      allocate(tmp(newcap))
      tmp(1:self%n) = self%refs(1:self%n)
      call move_alloc(tmp, self%refs)
    end if

    self%n = self%n + 1
    self%refs(self%n) = ref
  end subroutine reflist_add

  !-------------------------------------------------------------------

  logical function reflist_has_key(self, key)
    class(t_reflist), intent(in) :: self
    character(len=*), intent(in) :: key
    integer :: i
    reflist_has_key = .false.
    do i = 1, self%n
      if (trim(self%refs(i)%key) == trim(key)) then
        reflist_has_key = .true.
        return
      end if
    end do
  end function reflist_has_key

  !-------------------------------------------------------------------

  subroutine reflist_clear(self)
    class(t_reflist), intent(inout) :: self
    self%n = 0
    ! Keep allocation; just reset counter
  end subroutine reflist_clear

  !===================================================================
  ! Load master database from pipe-delimited text file
  !
  ! File format (one reference per line, # lines are comments):
  !
  !   key | type | bibtex_key | adslink
  !
  ! Example:
  !   H2O_Polyansky2018 | opacity | Polyansky2018 | Polyansky et al. | 2018 | MNRAS | 10.1093/mnras/sty1877
  !===================================================================

  subroutine load_reference_database(filename, master)
    character(len=*),  intent(in)    :: filename
    type(t_reflist),   intent(inout) :: master

    integer            :: iunit, ios, itype, j
    character(len=512) :: line
    character(len=48)  :: f_key
    character(len=12)  :: f_typestr
    character(len=64)  :: f_bibtex
    character(len=192) :: f_adslink
    type(t_reference)  :: ref

    call master%init(128)

    open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(2a)') 'ERROR mod_references: cannot open ', trim(filename)
      return
    end if

    do
      read(iunit, '(a)', iostat=ios) line
      if (ios /= 0) exit

      ! Skip blank lines and comment lines
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#')    cycle

      ! Parse 8 pipe-delimited fields
      call parse_pipe_field(line, 1, f_key)
      call parse_pipe_field(line, 2, f_typestr)
      call parse_pipe_field(line, 3, f_bibtex)
      call parse_pipe_field(line, 4, f_adslink)

      ! Map type string to integer constant
      itype = 0
      do j = 1, NREF_TYPES
        if (trim(adjustl(f_typestr)) == trim(REF_TYPE_KEY(j))) then
          itype = j
          exit
        end if
      end do
      if (itype == 0) then
        write(*,'(3a)') 'WARNING mod_references: unknown type "', &
          trim(f_typestr), '" — skipping entry'
        cycle
      end if

      ref%key        = trim(strip_ctrl(adjustl(f_key)))
      ref%ref_type   = itype
      ref%bibtex_key = trim(strip_ctrl(adjustl(f_bibtex)))
      ref%adslink    = trim(adjustl(f_adslink))

      call master%add(ref)
    end do

    close(iunit)
    write(*,'(a,i4,2a)') 'mod_references: loaded ', master%n, &
      ' references from ', trim(filename)
  end subroutine load_reference_database

  !===================================================================
  ! Register a reference into the active list by key.
  ! Looks it up in master; warns if key is not found.
  !===================================================================

  subroutine register_ref(key)		!active, master, key)
!    type(t_reflist),   intent(inout) :: active
!    type(t_reflist),   intent(in)    :: master
    character(len=*),  intent(in)    :: key
    integer :: i

    do i = 1, master_refs%n
      if (trim(master_refs%refs(i)%key) == trim(key)) then
        call active_refs%add(master_refs%refs(i))
        return
      end if
    end do

    write(*,'(2a)') 'WARNING mod_references: key not found in database: ', trim(key)
  end subroutine register_ref

  !===================================================================
  ! Write LaTeX snippet: one booktabs table per reference type,
  ! only for types that have at least one active reference.
  ! Intended for \input{} in your paper.
  !===================================================================

  subroutine write_latex_reftable(filename,datadir)		!active, filename)
!    type(t_reflist),  intent(in) :: active
    character(len=*), intent(in) :: filename,datadir

    integer :: iunit, ios, itype, i
    logical :: found

    open(newunit=iunit, file=trim(filename), status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(2a)') 'ERROR mod_references: cannot write ', trim(filename)
      return
    end if

    write(iunit,'(a)') '\documentclass[12pt]{report}'
    write(iunit,'(a)') '\usepackage{hyperref}'
    write(iunit,'(a)') '\include{' // trim(datadir) // '/general.tex}'
    write(iunit,'(a)') '\begin{document}'
    write(iunit,'(a)') '\title{ARCiS auto generated references}'
    write(iunit,'(a)') '\author{}'
    write(iunit,'(a)') '\date{\today}'
    write(iunit,'(a)') '\maketitle'

    write(iunit,'(a)') '% ARCiS auto-generated reference tables'
    write(iunit,'(a)') '% \usepackage{booktabs} required in preamble'
    write(iunit,'(a)') ''

    do itype = 1, NREF_TYPES

      ! Check whether any active refs belong to this type
      found = .false.
      do i = 1, active_refs%n
        if (active_refs%refs(i)%ref_type == itype) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) cycle

      write(iunit,'(a)')    '\begin{longtable}{l|l|l}'
      write(iunit,'(3a)')   '  \caption{', trim(REF_TYPE_LABEL(itype)), '}\\'
      write(iunit,'(a)')    '    \hline'
      write(iunit,'(a)')    '  \centering'
      write(iunit,'(a)')    '    Name & citation & link \\'
      write(iunit,'(a)')    '    \hline'

      do i = 1, active_refs%n
        if (active_refs%refs(i)%ref_type /= itype) cycle
        ! Build the row as a single string to avoid line-length wrapping
        write(iunit,'(a)') &
          '    \texttt{' // trim(active_refs%refs(i)%key) // '} & ' // &
          '    {\cite{' // trim(active_refs%refs(i)%bibtex_key) // '}} & ' // &
          '\href{' // trim(strip_ctrl(active_refs%refs(i)%adslink)) // '}{ADS link} \\'
      end do

      write(iunit,'(a)') '    \hline'
      write(iunit,'(a)') '\end{longtable}'
      write(iunit,'(a)') ''

    end do

    write(iunit,'(a)') '\bibliographystyle{' // trim(datadir) // '/aa.bst}'
    write(iunit,'(a)') '\bibliography{' // trim(datadir) // '/biblist.bib}'

    write(iunit,'(a)') '\end{document}'
    close(iunit)
    write(*,'(2a)') 'mod_references: wrote reference tables to ', trim(filename)
  end subroutine write_latex_reftable

  !===================================================================
  ! Internal helper: extract the n-th pipe-delimited field from a line
  !===================================================================

  subroutine parse_pipe_field(line, n, field)
    character(len=*), intent(in)  :: line
    integer,          intent(in)  :: n
    character(len=*), intent(out) :: field

    integer :: pos, count, start, finish, linelen

    field    = ''
    linelen  = len_trim(line)
    count    = 0
    start    = 1

    do pos = 1, linelen
      if (line(pos:pos) == '|') then
        count = count + 1
        if (count == n) then
          ! Field n runs from start to pos-1
          finish = pos - 1
          if (finish >= start) then
            field = adjustl(line(start:finish))
          end if
          return
        end if
        start = pos + 1
      end if
    end do

    ! Last field (no trailing pipe needed)
    if (count == n - 1 .and. start <= linelen) then
      field = adjustl(line(start:linelen))
    end if
  end subroutine parse_pipe_field


  !===================================================================
  ! Replace bare & with \& so author strings are safe in LaTeX tables
  !===================================================================

  function latex_escape_amp(str) result(out)
    character(len=*), intent(in)  :: str
    character(len=256)            :: out
    integer :: i, j
    j = 0
    out = ''
    do i = 1, len_trim(str)
      if (str(i:i) == '&') then
        j = j + 1; out(j:j) = '\'
        j = j + 1; out(j:j) = '&'
      else
        j = j + 1; out(j:j) = str(i:i)
      end if
    end do
  end function latex_escape_amp

  !===================================================================
  ! Strip trailing control characters (CR, LF, tab) from a string
  !===================================================================

  function strip_ctrl(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: out
    integer :: i, last
    out  = str
    last = len_trim(str)
    do while (last > 0)
      i = iachar(out(last:last))
      if (i == 10 .or. i == 13 .or. i == 9) then
        last = last - 1
      else
        exit
      end if
    end do
    if (last < len(out)) out(last+1:) = ' '
  end function strip_ctrl


end module mod_references

