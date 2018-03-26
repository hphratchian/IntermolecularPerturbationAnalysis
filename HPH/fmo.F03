      module fmo_mod
      use mqc_general
      use mqc_algebra2
      use mqc_gaussian
      use iso_fortran_env

      contains

      function countNBasis(iAtom,atom2bfMap) Result(nBasis)
!
!     This function returns an (intrinsic) integer giving the number of basis
!     functions, <nBasis>, belonging to atomic center number <iAtom> according
!     to the atom --> basis function map <atom2bfMap).
!
!
      implicit none
      integer(kind=int64),intent(in)::iAtom
      type(MQC_Variable)::atom2bfMap
      integer(kind=int64)::nBasis
      integer(kind=int64)::bf1,bf2
!
      bf1 = INT(MQC_Variable_get_MQC(atom2bfMap,[iAtom]))
      bf2 = INT(MQC_Variable_get_MQC(atom2bfMap,[iAtom+1]))
      nBasis = bf2-bf1
!
      return
      end function countNBasis


      function build_atom2bfMap(nAtoms,bf2AtomMap)  &
        Result(atom2bfMap)
!
!     This function produces an atom --> basis function map. The output array is
!     dimensioned (nAtoms+1) long. Basis functions on atom center number i runs
!     from atom2bfMap(i) to atom2bfMap(i+1)-1.
!
!
      implicit none
      integer(kind=int64),intent(in)::nAtoms
      class(MQC_Variable)::bf2AtomMap
      type(MQC_Variable)::temp
      integer(kind=int64),dimension(:),allocatable::atom2bfMap
      integer(kind=int64),dimension(1)::tempArray
      integer(kind=int64)::i,currentAtom
!
!
!     Allocate memory for atom2bfMap and then fill it.
!
      allocate(atom2bfMap(nAtoms+1))
      atom2bfMap = 0
      atom2bfMap(nAtoms+1) = SIZE(bf2AtomMap)+1
      call bf2AtomMap%print(iout=6,header='bf2AtomMap')
      tempArray(1) = 1
      do i = 1,SIZE(bf2AtomMap)
        tempArray(1) = i
        temp = MQC_Variable_get_MQC(bf2AtomMap,tempArray)
        if(atom2bfMap(INT(temp)).eq.0) atom2bfMap(INT(temp)) = i
      endDo
!
!     In case there are atomic centers without any basis functions -- as there
!     will be in the fragment calculations -- we have to step through
!     <atom2bfMap> and replace zeros.
!
      do i = nAtoms,1,-1
        if(atom2bfMap(i).eq.0) atom2bfMap(i) = atom2bfMap(i+1)
      endDo
!
      return
      end function build_atom2bfMap


      function build_full2fragmentMap(iFragment,atom2fragment,  &
        atom2bfFull,bf2AtomFull,atom2bfFragment,bf2AtomFragment)  &
        Result(bfFull2FragmentMap)
!
!     This function produces a fragment --> full basis set map for a particular
!     fragment (number iFragment). To do thise work, the function needs the
!     fragment number being considered <iFragment>, the fragment list over atoms
!     <atom2fragment>, the full molecule's basis function --> atomic center
!     map <bf2AtomFull>, and the current fragment job's basis function -->
!     atomic center map <bf2AtomFragment>.
!
!
      implicit none
      integer(kind=int64),intent(in)::iFragment
      type(mqc_variable)::atom2fragment,atom2bfFull,bf2AtomFull
      type(mqc_variable)::atom2bfFragment,bf2AtomFragment
      type(mqc_variable)::bfFull2FragmentMap
!
      integer(kind=int64)::i,j,myAtom,myFragment
      integer(kind=int64),dimension(:,:),allocatable::bfFull2FragmentMapTemp
!
!     Do the map work.
!
      allocate(bfFull2FragmentMapTemp(SIZE(bf2AtomFull,1),2))
!hph      bfFull2FragmentMapTemp(:,:) = int(0.0)
      do i = 1,SIZE(bf2AtomFull,1)
        do j = 1,2
          bfFull2FragmentMapTemp(i,j) = 0
        endDo
      endDo
!      write(*,*)' Hrant: ',bfFull2FragmentMapTemp(1,1)
!      call mqc_print(6,bfFull2FragmentMapTemp,header='Hrant: Flag 1')
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - Building frag->full map...'
      write(*,*)'         Fragment = ',iFragment
      write(*,*)
      do myAtom = 1,SIZE(atom2fragment)
        myFragment = INT(MQC_Variable_get_MQC(atom2fragment,[myAtom]))
        write(*,*)'         myAtom = ',myAtom
        write(*,*)'         myFrag = ',myFragment
        if(myFragment.eq.iFragment) then
          j = INT(MQC_Variable_get_MQC(atom2bfFull,[myAtom]))
          write(*,*)'           MATCH! j = ',j,  &
            ' nBasis = ',countNBasis(myAtom,atom2bfFull)
          do i = INT(MQC_Variable_get_MQC(atom2bfFragment,[myAtom])),  &
            INT(MQC_Variable_get_MQC(atom2bfFragment,[myAtom+1]))-1
            bfFull2FragmentMapTemp(j,1) = i
            bfFull2FragmentMapTemp(j,2) = iFragment
            j = j+1
          endDo
        endIf
      endDo
      bfFull2FragmentMap = bfFull2FragmentMapTemp
!
      return
      end function build_full2fragmentMap


      function build_blockFragmentMatrixAO(iFragment,bfFull2Fragment,  &
        fragmentMatrix) Result(blockMatrix)
!
!     This function pushes a matrix in a fragment AO basis (and sized NFragBasis
!     x NFragBasis) out to a block in a full molecule AO basis.
!
!
      implicit none
      integer(kind=int64),intent(in)::iFragment
      type(mqc_variable)::bfFull2Fragment,fragmentMatrix
      type(mqc_variable)::blockMatrix
      real(kind=real64),dimension(:,:),allocatable::blockMatrixTemp
!
      integer(kind=int64)::i,j,nBasisFull,iFrag,jFrag,iFragBF,jFragBF
      real(kind=real64)::temp
!
!     Do the work...
!
      nBasisFull = SIZE(bfFull2Fragment,1)
      allocate(blockMatrixTemp(nBasisFull,nBasisFull))
!hph      blockMatrixTemp = 0.0
      do i = 1,nBasisFull
        do j = 1,nBasisFull
          blockMatrixTemp(i,j) = dfloat(0)
        endDo
      endDo

      write(*,*)
      write(*,*)
      write(*,*)' Hrant - block matrix work:'
      do i = 1,nBasisFull
        iFrag = INT(MQC_Variable_get_MQC(bfFull2Fragment,[i,2]))
        if(iFrag.eq.iFragment) then
          iFragBF = INT(MQC_Variable_get_MQC(bfFull2Fragment,[i,1]))
          do j = 1,nBasisFull
            jFrag = INT(MQC_Variable_get_MQC(bfFull2Fragment,[j,2]))
            if(jFrag.eq.iFragment) then
              jFragBF = INT(MQC_Variable_get_MQC(bfFull2Fragment,[j,1]))
              temp = float(MQC_Variable_get_MQC(fragmentMatrix,[iFragBF,jFragBF]))
              write(*,'(3x,''i='',I3,'' j='',I3,'' iFrag='',I3,'' jFrag='',I3,'' val='',f10.6)') i,j,iFrag,jFrag,temp
              blockMatrixTemp(i,j) = float(MQC_Variable_get_MQC(fragmentMatrix,[iFragBF,jFragBF]))
            endIf
          endDo
        endIf
      endDo
      write(*,*)
      write(*,*)
      blockMatrix = blockMatrixTemp
!
      return
      end function build_blockFragmentMatrixAO





      end module fmo_mod




      program fmo
!
!     This program carries out a simple fragment-MO calculation. The program
!     requires 2 command line arguments: (1) the name of a Gaussian matrix file
!     containing results from the full-molecule calculation; and (2) the name of
!     a Gaussian FChk file containing results from the fragment calculation.
!
!
      use mqc_general
      use mqc_algebra2
      use mqc_gaussian
      use iso_fortran_env
      use fmo_mod
!
      implicit none
      integer(kind=int64),parameter::IOut=6
      integer(kind=int64)::i,j,nAtoms,nFragments,nBasis,nBasisUse,nBasisFragment
      integer(kind=int64)::intTmp
      character(len=512)::filenameFullMolecule,filenameFragment
      logical::DEBUG=.false.,ok
      type(mqc_gaussian_unformatted_matrix_file)::fileFullMolecule,fileFragment
      type(MQC_Variable)::atom2fragment,atom2bfFull,atom2bfFragment,  &
        bf2AtomFull,bf2AtomFragment,bfFull2FragmentMap
      type(MQC_Variable)::fullMoleculeSao,fullMoleculeFaoAlpha
      type(MQC_Variable)::neAlphaFragment,neBetaFragment,  &
        fragmentMatrixTemp,fragmentSaoZero,fragmentFaoAlphaZero
!
!     Format statements
!
 1000 Format(1x,'Fragment Molecular Orbital Program',/,1x,'Number of Fragments=',I2)
 1010 Format(1x,'Full-Molecule Matrix Filename: ',A)
 1020 Format(1x,'Fragment ',I1,' FChk Filename: ',A)
 2000 Format(1x,'NBasis=',I4,' NBasisUse=',I4)
 3900 Format(1x,'Done getting data from the full molecule matrix file.')
!
!
!     Begin by getting the user-provided data file names from the command line.
!
      nFragments = command_argument_count() - 1
      write(IOut,1000) nFragments
      call get_command_argument(1,filenameFullMolecule)
      write(IOut,1010) TRIM(filenameFullMolecule)
      call neAlphaFragment%clear(0,dimensions=[nFragments+1])
      call neBetaFragment%clear(0,dimensions=[nFragments+1])
!
!     Load the full-molecule info.
!
      call fileFullMolecule%load(filenameFullMolecule)
      call fileFullMolecule%getArray('INTEGER FRAGMENT',mqcVarOut=atom2fragment)
      call atom2fragment%print(header='Fragment List')
      nAtoms = fileFullMolecule%getVal('nAtoms')
      nBasis = fileFullMolecule%getVal('nBasis')
      nBasisUse = fileFullMolecule%getVal('nBasisUse')
      call bf2AtomFull%clear(0,dimensions=[nBasis])
      call bfFull2FragmentMap%clear(0,dimensions=[nBasis,2])
      do i = 1,nBasis
        call bf2AtomFull%put(fileFullMolecule%getBasisInfo(i,'basis2Atom'),[i])
      endDo
      call bf2AtomFull%print(header='bf2AtomFull')
      write(IOut,2000) nBasis,nBasisUse
      atom2bfFull = build_atom2bfMap(nAtoms,bf2AtomFull)
      call mqc_print(atom2bfFull,header='atom --> basis function map (FULL)')
!
!     Now, get the overlap, fock, and MO coefficient matrices for the full
!     molecule..
!
      call fileFullMolecule%getArray('OVERLAP',mqcVarOut=fullMoleculeSao)
      call fileFullMolecule%getArray('ALPHA FOCK MATRIX',mqcVarOut=fullMoleculeFaoAlpha)
      call fileFullMolecule%closeFile()
!
!     Initialize/clear a few MQC variables that we will use below for the
!     zero-order matrices.
!
      call fragmentSaoZero%clear(0.0,[nBasis,nBasis])
      call fragmentFaoAlphaZero%clear(0.0,[nBasis,nBasis])
!
!     Load the fragment data.
!
      do i = 1,nFragments
        call get_command_argument(i+1,filenameFragment)
        call fileFragment%load(filenameFragment)
        write(IOut,1020) i,TRIM(filenameFragment)
        intTmp = fileFragment%getVal('nAtoms')
        write(IOut,*)' Fragment nAtoms    = ',intTmp
        nBasisFragment = fileFragment%getVal('nBasis')
        write(IOut,*)' Fragment nBasis    = ',nBasisFragment
        call bf2AtomFragment%init('int',dimensions=[nBasisFragment])
        intTmp = fileFragment%getVal('nBasisUse')
        write(IOut,*)' Fragment nBasisUse = ',intTmp
        do j = 1,nBasisFragment
          call bf2AtomFragment%put(fileFragment%getBasisInfo(i,'basis2Atom'),[j])
        endDo
        call bf2AtomFragment%print(header='bf2AtomFragment')
        atom2bfFragment = build_atom2bfMap(nAtoms,bf2AtomFragment)
        call mqc_print(atom2bfFragment,header='atom --> basis function map (Fragment)')
        bfFull2FragmentMap =  &
          bfFull2FragmentMap + build_full2fragmentMap(i,atom2fragment,  &
          atom2bfFull,bf2AtomFull,atom2bfFragment,bf2AtomFragment)
        call bfFull2FragmentMap%print(header='bfFull2FragmentMap')
!
!       Get the overlap matrix from the fragment file.
        call fileFragment%getArray('OVERLAP',mqcVarOut=fragmentMatrixTemp)
        if(DEBUG) call fragmentMatrixTemp%print(header='fragment OVERLAP')
        if(DEBUG) call mqc_print(  &
          mqcVariable=build_blockFragmentMatrixAO(i,bfFull2FragmentMap,  &
          fragmentMatrixTemp),header='Hrant - Temporary OVERLAP pushed out to FULL.')
        fragmentSaoZero =  &
          fragmentSaoZero + build_blockFragmentMatrixAO(i,bfFull2FragmentMap,  &
          fragmentMatrixTemp)
        call fragmentMatrixTemp%print(header='S for fragment')
        if(DEBUG) call fragmentSaoZero%print( &
          header='S for fragment after pushing to full')
!
!       Get the Fock matrix from the fragment file.
        call fileFragment%getArray('ALPHA FOCK MATRIX',mqcVarOut=fragmentMatrixTemp)
        if(DEBUG) call fragmentMatrixTemp%print(header='fragment FOCK(alpha)')
        if(DEBUG) call mqc_print(  &
          mqcVariable=build_blockFragmentMatrixAO(i,bfFull2FragmentMap,  &
          fragmentMatrixTemp),header='Hrant - Temporary Fock(alpha) matrix pushed out to FULL.')
        fragmentFaoAlphaZero =  &
          fragmentFaoAlphaZero + build_blockFragmentMatrixAO(i,bfFull2FragmentMap,  &
          fragmentMatrixTemp)
        call fragmentMatrixTemp%print(header='F(alpha) for fragment')
        if(.TRUE.) call fragmentFaoAlphaZero%print( &
          header='F(alpha) for fragment after pushing to full')
!
!       Done with the full molecule file. Close it.
        call fileFragment%closeFile()
        write(IOut,3900)
      endDo
!
!     Print fragment matrix results.
!
      call bfFull2FragmentMap%print(header='FINAL: bfFull2FragmentMap')
      call fragmentSaoZero%print(header='FINAL: BLOCK Overlap Matrix')
      call fragmentFaoAlphaZero%print(header='FINAL: BLOCK Fock(alpha) Matrix')
!
!     Print full-molecule overlap and Fock minus fragment-blocked overlap and
!     Fock matrices.
!
      call mqc_print(fullMoleculeSao-fragmentSaoZero,header='S Difference')
      call mqc_print(fullMoleculeFaoAlpha-fragmentFaoAlphaZero,header='F Difference')


!
!     


      write(*,*)
      write(*,*)
      write(*,*)' DONE'
      end program fmo
