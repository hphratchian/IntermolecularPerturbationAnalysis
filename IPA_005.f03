        Program IPA_004
        Use MQC_Gaussian
        Use MQC_Algebra2
        Use MQC_Files
        Use MQC_General
        use iso_fortran_env

!       This program is designed to calculate the orbital mixing between
!       three fragments of a molecule and determine the types of bonding 
!       occuring between the fragments.This will be done using first,
!       second and third order perturbation theory. 

!       In order to run the program you must first run four gaussian
!       jobs. The jobs will be:
!       First --> The entire molecule
!       Second --> Fragment One, with the remaining atoms ghosted
!       Third --> Fragment Two, with the remaining atoms ghosted
!       Fourth --> Fragment Three, with the remaining atoms ghosted

!       Each of the gaussian jobs need to form the associated matrix
!       file using the output=matrix keyword. The program will then 
!       read in each of the matrix files and run.

!       How to run:
!       ./CodeName molecule.mat frag1.mat frag2.mat frag3.mat 

!       ______________________________________________________________________
!       ______________________________________________________________________
!       General Varible Declarations
        Implicit None
       
        Integer :: iout=6, p, q, r, i
        Real(kind=real64) :: N
        Type(MQC_Variable) :: temp, temp1 
        Integer(kind=int64) :: NAtoms,NBasis,NElectrons
        Integer(kind=int64) :: NAtomsTotal,NBasisTotal,NElectronsTotal
        Type(MQC_Variable) :: Fock_AO_Composite,Fock_FMO_Composite
        Type(MQC_Variable) :: Coefficent_Composite,Fock_FMO_Molecule
        Type(MQC_Variable) :: Energy_Composite 
        Type(MQC_Variable) :: First_Order_Energy
        Type(MQC_Variable) :: Second_Order_Energy
        Type(MQC_Variable) :: DeltaE, VonE, VVonE
        Type(MQC_Variable) :: V
        Type(MQC_Variable) :: Third_Order_Energy
        Type(MQC_Variable) :: DoubleDeltaE, VV, VVV, TermOne, TermTwo
        Type(MQC_Variable) :: TermOneminusTermTwo, EoneEtwo

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Open each of the command arguments and use the Basic_Information
!       subroutine to calculate/collect the basis information need for
!       the remainder of the program and provide basic information the
!       user may want to use.

!       Varible declaration for subroutines
        Character(len=256) :: CL_1, CL_2, CL_3, CL_4
        
        Type(MQC_Variable) :: M_Fock, M_Overlap, M_Density, M_Coeff
        Type(MQC_Variable) :: M_Hamiltonian, M_Energy
        Integer(kind=int64) :: M_NAtoms, M_NBasis, M_NElectrons
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: M_Matrix_File 
        
        Type(MQC_Variable) :: F1_Fock, F1_Overlap, F1_Density, F1_Coeff
        Type(MQC_Variable) :: F1_Hamiltonian, F1_Energy
        Integer(kind=int64) :: F1_NAtoms, F1_NBasis, F1_NElectrons
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F1_Matrix_File 
        
        Type(MQC_Variable) :: F2_Fock, F2_Overlap, F2_Density, F2_Coeff
        Type(MQC_Variable) :: F2_Hamiltonian, F2_Energy
        Integer(kind=int64) :: F2_NAtoms, F2_NBasis, F2_NElectrons
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F2_Matrix_File 
        
        Type(MQC_Variable) :: F3_Fock, F3_Overlap, F3_Density, F3_Coeff
        Type(MQC_Variable) :: F3_Hamiltonian, F3_Energy
        Integer(kind=int64) :: F3_NAtoms, F3_NBasis, F3_NElectrons
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F3_Matrix_File 
        
        Call Get_Command_Argument(1, CL_1)
        Write(iout,*) 'The basic information for the molecule is:'
        Call Basic_Information(CL_1,M_Matrix_File,M_NAtoms, &
             M_NBasis,M_NElectrons,M_Fock,M_Overlap,M_Density, &
             M_Coeff,M_Hamiltonian,M_Energy)
        
        Call Get_Command_Argument(1, CL_2)
        Write(iout,*) 'The basic information for fragment_1 is:'
        Call Basic_Information(CL_2,F1_Matrix_File,F1_NAtoms, &
             F1_NBasis,F1_NElectrons,F1_Fock,F1_Overlap,F1_Density, &
             F1_Coeff,F1_Hamiltonian,F1_Energy)
        
        Call Get_Command_Argument(1, CL_3)
        Write(iout,*) 'The basic information for fragment_2 is:'
        Call Basic_Information(CL_3,F2_Matrix_File,F2_NAtoms, &
             F2_NBasis,F2_NElectrons,F2_Fock,F2_Overlap,F2_Density, &
             F2_Coeff,F2_Hamiltonian,F2_Energy)

        Call Get_Command_Argument(1, CL_4)
        Write(iout,*) 'The basic information for fragment_3 is:'
        Call Basic_Information(CL_4,F3_Matrix_File,F3_NAtoms, &
             F3_NBasis,F3_NElectrons,F3_Fock,F3_Overlap,F3_Density, &
             F3_Coeff,F3_Hamiltonian,F3_Energy)

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Double check that the fragments and molecule infomation matches.
!       This section checks that the number of basis functions, atoms
!       and electrons are the same in each of the files. If they are not
!       the same an error message will be printed and the program will
!       end. Once this section normally terminates the number of basis
!       functions, atoms and electrons will be printed for the system. 

!       Variable Declarations        
!^        Integer(kind=int64) :: NAtoms,NBasis,NElectrons
!^        Integer(kind=int64) :: NAtomsTotal,NBasisTotal,NElectronsTotal
        
        Write(iout,*) 'Checking that the basic information matches...'
        
!       Make sure that the number of atoms is equivelent. 
        NAtomsTotal=F1_NAtoms+F2_NAtoms+F3_NAtoms
        If (NAtomsTotal.eq.3*M_NAtoms) then
         Write(iout,*) 'Atoms are equal'
        Else
         Write(iout,*) 'Atoms not equal. Program terminating.'
         STOP
        EndIf
        NAtoms=M_NAtoms

!       Make sure that the number of electrons is equivelent. 
        NElectronsTotal=F1_NElectrons+F2_NElectrons+F3_NElectrons
        If (NElectronsTotal.eq.3*M_NElectrons) then
         Write(iout,*) 'Electrons are equal'
        Else
         Write(iout,*) 'Electrons not equal. Program terminating.'
         STOP
        EndIf
        NElectrons=M_NElectrons

!       Make sure that the number of basis functions is equivelent. 
        NBasisTotal=F1_NBasis+F2_NBasis+F3_NBasis
        If (NBasisTotal.eq.3*M_NBasis) then
         Write(iout,*) 'Basis functions are equal'
        Else
         Write(iout,*) 'Basis functions not equal. Program terminating.'
         STOP
        EndIf
        NBasis=M_NBasis

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Composite Matricies
!       ______________________________________________________________________
!       Make the composite coefficent matrix and 
!       fock matrix. To do this we will add the three fragment matricies
!       together. This assumes that the matrix will be ~0 in all of the
!       columns/rows associated with the ghosted atoms. These matricies
!       will then be transformed into the Fragment Molecular Orbital
!       (FMO) basis.

!       Variable Declarations and initializing
!^        Type(MQC_Variable) :: Fock_AO_Composite,Fock_FMO_Composite
!^        Type(MQC_Variable) :: Coefficent_Composite,Fock_FMO_Molecule
!^        Type(MQC_Variable) :: Energy_Composite 

        Call MQC_Variable_initialize(Coefficent_Composite,'real', &
             rank=1,dimensions=[NBasis])
        Call MQC_Variable_initialize(Fock_AO_Composite,'real', & 
             rank=2,dimensions=[NBasis,NBasis])
        Call MQC_Variable_initialize(Fock_FMO_Composite,'real', &       
             rank=2,dimensions=[NBasis,NBasis])
        Call MQC_Variable_initialize(Fock_FMO_Molecule,'real', & 
             rank=2,dimensions=[NBasis,NBasis])

        Energy_Composite=F1_Energy+F2_Energy+F3_Energy        
        Write(*,*) 'F1_Coeff Rank:', RANK(F1_Coeff)
        Write(*,*) 'F2_Coeff Rank:', RANK(F2_Coeff)
        Write(*,*) 'F3_Coeff Rank:', RANK(F3_Coeff)
        Write(*,*) 'F1_Coeff Size:', SIZE(F1_Coeff)
        Write(*,*) 'F2_Coeff Size:', SIZE(F2_Coeff)
        Write(*,*) 'F3_Coeff Size:', SIZE(F3_Coeff)
        Coefficent_Composite=(F1_Coeff+F2_Coeff)+F3_Coeff 
!*        Fock_AO_Composite=F1_Fock+F2_Fock+F3_Fock
        
!*        Fock_FMO_Composite=transpose(Coefficent_Composite)* &
!*                           Fock_AO_Composite*Coefficent_Composite
!*        Fock_FMO_Molecule=transpose(Coefficent_Composite)* &
!*                          M_Fock*Coefficent_Composite

!       ______________________________________________________________________
!       ______________________________________________________________________
!       First Order Energy Correction
!       ______________________________________________________________________

!       Variable declarations and initializing
!^        Type(MQC_Variable) :: First_Order_Energy
     
        Call MQC_Variable_initialize(First_Order_Energy,'real', &
             rank=1,dimensions=[NBasis])
      
        Do p=1,NBasis       
!*         Call First_Order_Energy%put(Fock_FMO_Molecule%at(p,p) &
!*              -Fock_FMO_Composite(p,p),p)
        EndDo

        Call First_Order_Energy%print(iout,'The first order energy &
             correction is: ')

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Second Order Energy Correction
!       ______________________________________________________________________
!       Here two intermediates are formed to solve for the final energy
!       correction. This is due to the fact that one value must be
!       excluded in the final summation.
!       ______________________________________________________________________

!       Variable Declarations
!^        Type(MQC_Variable) :: Second_Order_Energy
!^        Type(MQC_Variable) :: DeltaE, VonE, VVonE
!^        Type(MQC_Variable) :: V
        
        Call MQC_Variable_initialize(V,'real',rank=2, &
             dimensions=[NBasis,NBasis])
        Call MQC_Variable_initialize(Second_Order_Energy,'real',rank=1)
        Call MQC_Variable_initialize(DeltaE,'real',rank=2, &
             dimensions=[NBasis,NBasis])
        Call MQC_Variable_initialize(VonE,'real',rank=2)
        Call MQC_Variable_initialize(VVonE,'real',rank=2)

!       Denominator
        Do p=1,NBasis
         Do q=1,NBasis
!*          Call DeltaE%put(Energy_Composite(p)- &
!*              Energy_Composite(q),p)
         EndDo
        EndDo

!       The diagonal of DeltaE will be zero, so this will be set to one
!       now for divison then changed later.
        Do p=1,NBasis
!*         DeltaE%at(p,p)=1
        EndDo

!       Now calculate V on E, or the amplitude, then multiply it by the
!       second V in the expression.
        Do p=1,NBasis
         Do q=1,NBasis
!*          Call V%put(Fock_FMO_Molecule(p,q)-Fock_FMO_Composite(p,q),p,q)
         EndDo
        EndDo

        Do p=1,NBasis
         Do q=1,NBasis
!*          Call VonE%put(V(p,q)/DeltaE(p,q),p,q)
         EndDo
        EndDo

        Do p=1,NBasis
         Do q=1,NBasis
!*          Call VVonE%put(VonE(p,q)*V(p,q),p,q)
         EndDo
        EndDo

!       Set the diagonal back to zero.
        
        Do p=1,NBasis
!*         VVonE%at(p,p)=0
        EndDo

!       Calculate the second order correction by summing each of the
!       columns in the VVonE matrix and putting them into a vector that
!       is NBasis long.
        N=0.0

        Do p=1,NBasis
         Do q=1,NBasis
!*          N=N+VVonE%at(p,q)
         EndDo
!*         Call Second_Order_Energy%put(N,p)
        EndDo

        Call Second_Order_Energy%print(iout,'The second order energy &
             correction is: ')
        
!       ______________________________________________________________________
!       ______________________________________________________________________
!       Third Order Energy Correction 
!       ______________________________________________________________________

!       Variable declarations        
!^        Type(MQC_Variable) :: Third_Order_Energy
!^        Type(MQC_Variable) :: DoubleDeltaE, VV, VVV, TermOne, TermTwo
!^        Type(MQC_Variable) :: TermOneminusTermTwo, EoneEtwo
        
        Call MQC_Variable_initialize(Third_Order_Energy,'real',rank=1)
        Call MQC_Variable_initialize(EoneEtwo,'real',rank=1)
        Call MQC_Variable_initialize(DoubleDeltaE,'real',rank=2, &
             dimensions=[NBasis,NBasis])
        Call MQC_Variable_initialize(VV,'real',rank=2)
        Call MQC_Variable_initialize(VVV,'real',rank=2)
        Call MQC_Variable_initialize(TermOne,'real',rank=2)
        Call MQC_Variable_initialize(TermTwo,'real',rank=2)
        Call MQC_Variable_initialize(TermOneminusTermTwo,'real',rank=2)

!       Denominator #1
        Do p=1,NBasisTotal
          Do q=1,NBasisTotal
            Do r=1,NBasisTotal
!*              Call DoubleDeltaE%put(((Energy_Composite(p)- &
!*              Energy_Composite(q))*(Energy_Composite(p)- &
!*              Energy_Composite(r))),p,q)
            EndDo
          EndDo
        EndDo

        Do p=1,NBasisTotal
!*          DoubleDeltaE%at(p,p)=1
        EndDo

!       Denominator #2 will be the same as DeltaE in the second order
!       energy correction

!       Numerator #1: Three V multiplied by eachother. This will be done
!       by first multiplyinh two together to make VV, then multiplying
!       that by another V.

        VV=(Fock_FMO_Molecule-Fock_FMO_Composite)*(Fock_FMO_Molecule- &
           Fock_FMO_Composite)                
        VVV=VV*(Fock_FMO_Molecule-Fock_FMO_Composite)

!       Numerator #2 is just the first order energy correction
!       multiplied by the second order energy correction.
        
        Do p=1,NBasisTotal
!*          Call EoneEtwo%put(First_Order_Energy_Correction%at(p)*  &
!*               Second_Order_Energy_Correction%at(p),p)
        EndDo
        
!       Term #1: This will be numerator one divided by denominator
!       number one

        TermOne=VVV/DoubleDeltaE        
        Do p=1,NBasisTotal
!*          TermOne%at(p,p)=0
        EndDo
        
!       Term #2: This will be numerator two divided by denominator
!       number two
        TermTwo=EoneEtwo/DeltaE
        Do p=1,NBasisTotal
!*          TermTwo%at(p,p)=0
        EndDo
        
!       Subtract the two terms to get the final third order energy
!       corrrection.

        TermOneminusTermTwo=TermOne-TermTwo
        
        N=0
        Do p=1,NBasisTotal
          Do q=1,NBasisTotal
!*            N=N+TermOneminusTermTwo%at(p,q)
          EndDo
!*          Call Third_Order_Energy_Correction%put(N,p)
        EndDo 

        Call Third_Order_Energy%print(iout,'The third order energy &
             correction is: ')
                
!       ______________________________________________________________________
!       ______________________________________________________________________

        End Program IPA_004

!       ********************************************************************
!       ********************************************************************
        Subroutine Basic_Information(CommandLine,MatrixFile,NAtoms,  &
        NBasis,NElectrons,Fock,Overlap,Density,Coeff,Hamiltonian,Energy)
        Use MQC_Gaussian
        Use MQC_Algebra2
        Use MQC_Files
        Use MQC_General
        use iso_fortran_env
!
!       This subroutine is designed to load all necessary information
!       and calculate/print basic information that comes from
!       the unformated matrix file. The info includes:
!       1. If the calculation is restricted or unrestricted? 
!       2. Print the number of atoms, electrons, and basis functions.
!       3. Calulate the number of electrons in the system with <PS>
!       4. Calculate the Hartree Fock Energy.
!       5. Calculate the IPM energy. 
!       6. Calculate the one and two electron energy 

!       Variable Declarations
        Implicit None
        Character(len=256) :: CommandLine
        Integer(kind=int64) :: iout=6, i
        Integer(kind=int64) :: NAtoms, NBasis, NElectrons
        Type(MQC_Variable) :: HFEnergy, IPMEnergy, OneEnergy, TwoEnergy
        Type(MQC_Variable) :: FD, HD, Half, temp
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: MatrixFile
        Type(MQC_Variable) :: Fock, Overlap, Density, Coeff, Hamiltonian
        Type(MQC_Variable) Energy
        Half=0.5
        Call MatrixFile%load(CommandLine)

!       Test if the job was restricted or unrestricted.
        Write(iout,*) 'Is restricted? ', MatrixFile%isrestricted()
        Write(iout,*) 'Is unrestricted? ', MatrixFile%isunrestricted()

!       Print the number of electrons, atoms and basis functions from
!       the matrix file.
        NElectrons=MatrixFile%getval('NElectrons')
        Write(iout,*) 'NElectrons=', NElectrons
        NAtoms=MatrixFile%getval('NAtoms')
        Write(iout,*) 'NAtoms=', NAtoms
        NBasis=MatrixFile%getval('NBasis')
        Write(iout,*) 'NBasis=', NBasis

!       Load and save all of the important matricies
        Call MatrixFile%getarray('overlap',mqcVarOut=Overlap)
!!        Call Overlap%print(iout,'The Overlap is: ') 
        Call MatrixFile%getarray('alpha fock matrix',mqcVarOut=Fock)
        Call MatrixFile%getarray('alpha density matrix', &
        mqcVarOut=Density)
        Density=Density+Density
        Call MatrixFile%getarray('alpha mo coefficients', &
        mqcVarOut=Coeff)
        Call MatrixFile%getarray('core hamiltonian alpha', &
        mqcVarOut=Hamiltonian)
        Call MatrixFile%getarray('alpha orbital energies', &
        mqcVarOut=Energy)

!       Calculate and print the number of electrons in the system.
!*        Call MQC_Print(contraction(Density,Overlap),iout,'<PS>= ')

!       Calculate and print the Hartree Fock Energy
!*        Call MQC_Variable_initialize(HFEnergy,'real')
!*        Call MQC_Variable_initialize(FD,'real')
!*        Call MQC_Variable_initialize(HD,'real')
        
!*        FD=contraction(Fock,Density)
!*        HD=contraction(Hamiltonian,Density)
!*        HFEnergy=Half*(FD+(Half*HD))
!*        Call HFEnergy%print(iout,'The Hartree Fock Energy is:')        

!       Calculate and print the energy from the IPM

!*        Call MQC_Variable_initialize(IPMEnergy,'real') 

!*        Do i=1,NElectrons/2
!*          temp=Energy%at(i)
!*          IPMEnergy=IPMEnergy+(temp+temp)
!*        EndDo

!*        Call IPMEnergy%print(iout,'The IPM Energy is: ')

!       Calculate and print the one electron energy
!*        Call MQC_Variable_initialize(OneEnergy,'real')
         
!*        OneEnergy=contraction(Fock,Density)-  &
!*        contraction(Hamiltonian,Density)

!*        Call OneEnergy%print(iout,'The One Electron Energy is: ') 

!       Calculate and print the two electron energy
!*        Call MQC_Variable_initialize(TwoEnergy,'real')
        
!*        TwoEnergy=contraction(Density,Hamiltonian)
!*        Call TwoEnergy%print(iout,'The Two Electron Energy is: ')

        Call MatrixFile%closefile()

        End Subroutine Basic_Information
!       ********************************************************************
!       ********************************************************************

