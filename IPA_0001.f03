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
      do i = 1,nBasisFull
        iFrag = INT(MQC_Variable_get_MQC(bfFull2Fragment,[i,2]))
        if(iFrag.eq.iFragment) then
          iFragBF = INT(MQC_Variable_get_MQC(bfFull2Fragment,[i,1]))
          do j = 1,nBasisFull
            jFrag = INT(MQC_Variable_get_MQC(bfFull2Fragment,[j,2]))
            if(jFrag.eq.iFragment) then
              jFragBF = INT(MQC_Variable_get_MQC(bfFull2Fragment,[j,1]))
              temp = float(MQC_Variable_get_MQC(fragmentMatrix, &
                          [iFragBF,jFragBF]))
              write(*,'(3x,''i='',I3,'' j='',I3,'' iFrag='',I3, &
                        '' jFrag='',I3,'' val='',f10.6)') i, &
                        j,iFrag,jFrag,temp
              blockMatrixTemp(i,j) = &
                float(MQC_Variable_get_MQC(fragmentMatrix, &
                [iFragBF,jFragBF]))
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

!
!     *********************************************************************
!     *********************************************************************
!

        Program IPA_0001
        Use MQC_Gaussian
        Use MQC_Algebra2
        Use MQC_Files
        Use MQC_General
        Use iso_fortran_env

!       Written By: Samantha Bidwell 
!       Code Version: 0001
!       Date: 2/25/20

!       Purpose
!               The program is developed to calculate the energy
!               contribution from the mixing of the MO's of three fragments.
!               This is done using the 1st, 2nd and 3rd order energy corrections
!               of perturbation theory.

!       General Information
!               - This code ising gaussian matrix files as its input
!               - The MQC package developed by the Hratchian group is
!                 used to analyze the matrix file.
!               - To compile the code run the associated makefile.

!       General Instructions for Using the Code
!               1. Run an optimization and frequency of the molecule you
!               would like to study with the keyword output=matrix to
!               create the matrix file.
!               2. Run individual jobs for each of fragments with the
!               remaining molecules as ghost atoms also using the
!               keyword output=matrix.
!               3. To complile the code...
!               4. To run the code (after it is compliled) use ./IPA_100
!               molecule.mat frag1.mat frag2.mat frag3.mat

!       Code Outline
!               -- Overall variable declaration
!               -- Read in the matrix files
!               -- Save the fock, overlap, and coefficent matricies and
!               the energy vectors
!               -- Run the basic_information subroutine to calculate the
!               HF_E, the IPM_E, the 2e_E, the 1e_E, and print the
!               number of atoms, basis functions and electrons for each
!               of the files.
!               -- Check the the sum of the fragment jobs is equivilant
!               to the overall molecule.
!               -- Form the composite matricies using the proper basis
!               function mapping technique.
!               -- Transform all of the matricies into the FMO basis.
!               -- Calculate the 1st, 2nd and 3rd order mixing for the
!               fragments
!               -- Print each of the energy corrections and emphasize
!               the orbitals with the largest contributions.

!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------

!       Variable Declarations
        Implicit None
       
        Character :: M, F1, F2, F3 
        Integer :: iout=6,p,q,r,s,i
        Integer :: NAtom,NBasis,NEle
        Integer :: NAtomTotal,NBasisTotal,NEleTotal
        Integer :: M_NAtom,M_NBasis,M_NEle
        Integer :: F1_NAtom,F1_NBasis,F1_NEle
        Integer :: F2_NAtom,F2_NBasis,F2_NEle
        Integer :: F3_NAtom,F3_NBasis,F3_NEle
        Real :: N
        Type(MQC_Variable) :: temp,temp1,temp2,temp3
        Type(MQC_Variable) :: F_Comp,F_FMO_Comp,F_FMO_Molecule
        Type(MQC_Variable) :: C_Comp,E_Comp
        Type(MQC_Variable) :: First_Order_E,Second_Order_E,Third_Order_E
        Type(MQC_Variable) :: DeltaE,VonE,VVonE,V,VV,VVV,TermOne,TermTwo
        Type(MQC_Variable) :: DoubleDeltaE,EoneEtwo,term1minusterm2
        Type(MQC_Variable) :: M_F,M_S,M_D,M_C,M_H,M_E
        Type(MQC_Variable) :: F1_F,F1_S,F1_D,F1_C,F1_H,F1_E
        Type(MQC_Variable) :: F2_F,F2_S,F2_D,F2_C,F2_H,F2_E
        Type(MQC_Variable) :: F3_F,F3_S,F3_D,F3_C,F3_H,F3_E
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: M_Mat
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F1_Mat
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F2_Mat
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F3_Mat


!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------

!       Read in each of the command line arguments.
!       Save the Fock, Overlap, Energy, Density and Coefficent matricies
!       from each of the matrix files using the Basic_In formation
!       subroutine.

        Write(*,*) 'STEP: Calling each of the command line arguments & 
        and reading in all of the information'
        Write(*,*) 'Reading in Molecule information.'
        Call Get_Command_Argument(1,M)
        Call Basic_Information(M,M_Mat,M_NAtom,M_NBasis,M_NEle, &
             M_F,M,M_D,M_C,M_H,M_E)
        Write(*,*) 'Reading in Fragment_1 information.'
        Call Get_Command_Argument(1,F1)
        Call Basic_Information(F1,F1_Mat,F1_NAtom,F1_NBasis,F1_NEle, &
             F1_F,F1_S,F1_D,F1_C,F1_H,F1_E)
        Write(*,*) 'Reading in Fragment_2 information.'
        Call Get_Command_Argument(1,F2)
        Call Basic_Information(F2,F2_Mat,F2_NAtom,F2_NBasis,F2_NEle, &
             F2_F,F2_S,F2_D,F2_C,F2_H,F2_E)
        Write(*,*) 'Reading in Fragment_3 information.'
        Call Get_Command_Argument(1,F3)
        Call Basic_Information(F3,F3_Mat,F3_NAtom,F3_NBasis,F3_NEle, &
             F3_F,F3_S,F3_D,F3_C,F3_H,F3_E)

!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------

!       Check that the number of basis functions, atoms and electrons
!       are equal to their expected values.

        Write(*,*) 'STEP: Checking that the number of atoms in each &
        fragment input are equal:'
        
        NAtomTotal=F1_NAtom+F2_NAtom+F3_NAtom
        
        If (NAtomTotal.eq.3*M_NAtom) then
         Write(iout,*) 'Atoms are equal. Continuing...'
        Else
         Write(iout,*) 'Atoms not equal. Program terminating.'
         STOP
        EndIf

        NAtom=M_NAtom

        Write(*,*) 'STEP: Checking that the number of basisfunctions &
        in each fragment input are equal:'
        
        NBasisTotal=F1_NBasis+F2_NBasis+F3_NBasis
        
        If (NBasisTotal.eq.3*M_NBasis) then
         Write(iout,*) 'Basis functions are equal. Continuing...'
        Else
         Write(iout,*) 'Basis functions not equal. Program terminating.'
         STOP
        EndIf

        NBasis=M_NBasis

        Write(*,*) 'STEP: Checking that the number of electrons &
        in each fragment input are equal:'
        
        NEleTotal=F1_NEle+F2_NEle+F3_NEle
        
        If (NEleTotal.eq.3*M_NEle) then
         Write(iout,*) 'Electrons are equal. Continuing...'
        Else
         Write(iout,*) 'Electrons not equal. Program terminating.'
         STOP
        EndIf

        NEle=M_NEle
        
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       Form the composite matricies.
        Write(*,*) 'STEP: Forming the composite matricies.'
!       The composite matrices that are formed are: Fock_Composite,
!       Coefficent_Composite, Overlap_Composite and Energy_Composite.
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       Transform each of the matricies into the FMO basis.
        Write(*,*) 'STEP: Transforming each of the matricies into the &
        fragment molecular orbital (FMO) basis.'
!       These matricies are created to simplify the following energy
!       corrections. 
        F_FMO_M = 
        F_FMO_C =
        S_FMO_M =
        S_FMO_C = 
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       First Order Energy Correction
        Write(*,*) 'STEP: Calculating the first order energy &
        correction.'
        Call First_Order_E%init(NBasisTotal,0)
        Do p=1,NBasis
          Call First_Order_E%put(F_FMO_Molecule%at(p,p)-F_FMO_Comp%at(p,p),p)
        EndDo
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        write(iout,*) 'The First Order Energy Correction'
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call First_Order_E%print(iout,'First order:')
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       Second Order Energy Correction
        Write(*,*) 'STEP: Calculating the second order energy &
        correction.'
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       Third Order Energy Correction       
        Write(*,*) 'STEP: Calculating the third order energy &
        correction.'
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
!       Print the energy corrections and emphasize the mixing that leads
!       to the largest contributions.
!       --------------------------------------------------------------------------- 
!       ---------------------------------------------------------------------------
        End Program IPA_0001


!       ****************************************************************************
!       ****************************************************************************
        Subroutine Basic_Information(Commandline,Mat,NAtoms,NBasis, &
        Nele,F,S,D,C,H,E)
        
        Use MQC_Gaussian
        Use MQC_Algebra2
        Use MQC_Files
        Use MQC_General
        Use iso_fortran_env

!       This subroutine is written to define:
!               1. Restricted or Unrestricted?
!               2. Print natoms, nbasis and nelectrons
!               3. The hartree fock energy
!               4. The IPM energy
!               5. The 1 and 2 electron energy

        Implicit None
        Character(len=256) :: Commandline
        Integer :: iout=6, NAtoms, NBasis, Nele, i
        Type(MQC_Variable) :: HFE, IPME, OneEle, TwoEle
        Type(MQC_Variable) :: F, S, D, C, H, E
        Type(MQC_Variable) :: FD, HD, tempIPM 
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: Mat 

        Call Mat%load(Commandline)
        Write(*,*) '  '
        Write(*,*) 'SUB-STEP: Loading information from the *.mat &
                    file', Commandline 
        
!       Testing and printing if the job is restricted or unrestricted
        Write(iout,*) 'Is restricted?', Mat%isrestricted()
        Write(iout,*) 'Is unrestricted?', Mat%isunrestricted()

!       Printing the number of ele, atoms and basis from the mat file
!       and saving them.
        Nele=Mat%getval('NElectrons')
        Write(iout,*) 'NElectrons=', Nele
        NAtoms=Mat%getval('NAtoms')
        Write(iout,*) 'NAtoms=', NAtoms
        Nele=Mat%getval('NBasis')
        Write(iout,*) 'NBasis=', NBasis
       

!       Load and save all need matricies/vectors (F,S,D,C,H,E)
        Call Mat%getarray('alpha fock matrix', &
        mqcVarOut=F)
        Call Mat%getarray('overlap',mqcVarOut=S)
        Call Mat%getarray('alpha density matrix', &
        mqcVarOut=D)
        Call Mat%getarray('alpha mo coefficents', &
        mqcVarOut=C)
        Call Mat%getarray('core hamiltonian alpha', &
        mqcVarOut=H)
        Call Mat%getarray('alpha orbital energies', &
        mqcVarOut=E)
        
        Write(*,*) '  '
        Write(*,*) 'SUB-STEP: Calculating basic information using &
        information for the *.mat file', Commandline

!       Calculate the Hartree Fock Energy (HF=1/2(<FD>+1/2<HD>))
        Call MQC_Variable_initialize(HFE,'real')
        Call MQC_Variable_initialize(HD,'real')
        Call MQC_Variable_initialize(FD,'real')
        FD=contraction(F,D)
        HD=contraction(H,D)
        HFE=(1/2)*(float(FD)+((1/2)*float(HD)))
        Call HFE%print(iout,'The HF energy is: ')

!       Calculate the IPM Energy
        Call MQC_Variable_initialize(IPME,'real')
        Call MQC_Variable_initialize(tempIPM,'real')
        
        Do i=1,Nele/2
!hph           tempIPM=E%at(i)
           IPME=float(IPME)+(2*float(tempIPM))
        EndDo
        Call IPME%print(iout,'The IPM energy is: ')

!       Calculate the one-electron energy
        Call MQC_Variable_initialize(OneEle,'real')
        OneEle=contraction(F,D)-contraction(H,D)
        Call OneEle%print(iout,'The One Electron Energy is: ')

!       Calculate the two-electron energy
        Call MQC_Variable_initialize(TwoEle,'real')
        TwoEle=contraction(D,H)
        Call TwoEle%print(iout,'The Two Electron Energy is: ')
        
        write(*,*) '  '
        write(*,*) 'SUB-STEP: Finishing basic information subroutine.'  
        End Subroutine Basic_Information
 
