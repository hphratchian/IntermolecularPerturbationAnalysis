        Program PT_Bonding_Analysis
        Use MQC_Gaussian
        Use MQC_Algebra
        Use MQC_Files
        Use MQC_General
 
!       This program is designed to calculate the orbital mixing between
!       three fragments of a molecule and to determine the types of bonding 
!       occuring between the fragments.This will be done using up to
!       third order perturbation theory. 

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
!       Varible Declarations
        Implicit None
        Integer :: iout=6, iunit_matrixfile

!       Input Characters
        Character(len=256) :: CommandLine1,CommandLine2
        Character(len=256) :: CommandLine3,CommandLine4
                
!       Varibles needed for basic information.
        Integer :: M_NAtoms, M_NBasis, M_NBasisUse
        Integer :: M_NElectrons,F1_NElectrons
        Integer :: F2_NElectrons,F3_NElectrons
        Integer :: F1_NAtoms,  F1_NBasis, F1_NBasisUse
        Integer :: F2_NAtoms,  F2_NBasis, F2_NBasisUse
        Integer :: F3_NAtoms,  F3_NBasis, F3_NBasisUse
        Integer :: NAtomsTotal, NBasisTotal, NElectronsTotal

!       Input matrix files
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: M_Matrix_File 
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F1_Matrix_File 
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F2_Matrix_File 
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: F3_Matrix_File 

!       Matricies pulled from the .mat files      
        Type(MQC_Matrix) :: M_Fock, M_Overlap, M_Density
        Type(MQC_Matrix) :: M_Coeff, M_Hamiltonian
        Type(MQC_Vector) :: M_Energy        
        Type(MQC_Matrix) :: F1_Fock, F1_Overlap, F1_Density
        Type(MQC_Matrix) :: F1_Coeff, F1_Hamiltonian
        Type(MQC_Vector) :: F1_Energy        
        Type(MQC_Matrix) :: F2_Fock, F2_Overlap, F2_Density
        Type(MQC_Matrix) :: F2_Coeff, F2_Hamiltonian
        Type(MQC_Vector) :: F2_Energy        
        Type(MQC_Matrix) :: F3_Fock, F3_Overlap, F3_Density
        Type(MQC_Matrix) :: F3_Coeff, F3_Hamiltonian
        Type(MQC_Vector) :: F3_Energy        

!       Random 
        Integer :: p,q,i,r
        Type(MQC_Matrix) :: temp, temp1 
        Real :: N

!       Terms needed to form composite matricies and put everything in
!       the FMO basis
        Type(MQC_Matrix) :: Fock_AO_Composite
        Type(MQC_Matrix) :: Fock_FMO_Composite
        Type(MQC_Matrix) :: Coefficent_Composite
        Type(MQC_Matrix) :: Fock_FMO_Molecule 
        
!       Vectors containing the values for each of the energy corrections
        Type(MQC_Vector) :: First_Order_Energy_Correction
        Type(MQC_Vector) :: Second_Order_Energy_Correction
        Type(MQC_Vector) :: Third_Order_Energy_Correction

!       Terms for second order energy correction
        Type(MQC_Matrix) :: DeltaE, VonE, VVonE
 
!       Terms for the third order energy correction
        Type(MQC_Matrix) :: DoubleDeltaE
        Type(MQC_Matrix) :: VV
        Type(MQC_Matrix) :: VVV
        Type(MQC_Matrix) :: TermOne
        Type(MQC_Matrix) :: TermTwo
        Type(MQC_Matrix) :: TermOneminusTermTwo
        Type(MQC_Vector) :: EoneEtwo


!       ______________________________________________________________________
!       Get all of the basic information for each of the matrix files
!       that are given to the program.
 
        Call Get_Command_Argument(1,CommandLine1)
        Write(iout,*) '**********************************************&
        ************************************************************'
        Write(iout,*) 'THE BASIC INFORMATION FOR THE MOLECULE IS: '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Basic_Information(Commandline1, M_Matrix_File,M_NAtoms,  &
        M_NBasis,M_NBasisUse,M_NElectrons,M_Fock,M_Overlap,  &
        M_Density,M_Coeff,M_Hamiltonian,M_Energy)
!*        Call M_Fock%print(iout, 'Molecular Fock Matrix')       
!*        Call M_Coeff%print(iout, 'Molecular Coefficent Matrix')       
 
        Call Get_Command_Argument(2,CommandLine2)
        Write(iout,*) ' '
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Write(iout,*) 'THE BASIC INFORMATION FOR FRAGMENT 1 IS: '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Basic_Information(Commandline2, F1_Matrix_File,F1_NAtoms, &
        F1_NBasis,F1_NBasisUse,F1_NElectrons,F1_Fock,F1_Overlap,  &
        F1_Density,F1_Coeff,F1_Hamiltonian,F1_Energy)
!*        Call F1_Fock%print(iout, 'F1 Fock Matrix')       
!*        Call F1_Coeff%print(iout, 'F1 Coefficent Matrix')       
            
        Call Get_Command_Argument(3,CommandLine3)
        Write(iout,*) ' '
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Write(iout,*) 'THE BASIC INFORMATION FOR FRAGMENT 2 IS: '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Basic_Information(Commandline3, F2_Matrix_File,F2_NAtoms, &
        F2_NBasis,F2_NBasisUse,F2_NElectrons,F2_Fock,F2_Overlap,  &
        F2_Density,F2_Coeff,F2_Hamiltonian,F2_Energy)
!*        Call F2_Fock%print(iout, 'F2 Fock Matrix')       
!*        Call F2_Coeff%print(iout, 'F2 Coefficent Matrix')       
 
        Call Get_Command_Argument(4,CommandLine4)
        Write(iout,*) ' '
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Write(iout,*) 'THE BASIC INFORMATION FOR FRAGMENT 3 IS: '
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Basic_Information(Commandline4, F3_Matrix_File,F3_NAtoms, &
        F3_NBasis,F3_NBasisUse,F3_NElectrons,F3_Fock,F3_Overlap,  &
        F3_Density,F3_Coeff,F3_Hamiltonian,F3_Energy)
!*        Call F3_Fock%print(iout, 'F3 Fock Matrix')       
!*        Call F3_Coeff%print(iout, 'F3 Coefficent Matrix')       
 
!       ______________________________________________________________________
!       Double check that the fragments and molecule infomation matches

       
        Write(iout,*) ' '
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        write(iout,*) 'Checking that the number of basis functions, &
        atoms and eletrons in each of the fragments and molecule  &
        are the same.' 
        Write(iout,*) '**********************************************&
        ************************************************************'
        Write(iout,*) ' '
        NAtomsTotal=F1_Natoms+F2_Natoms+F3_Natoms
        If (NAtomsTotal==3*M_NAtoms) then
          write(iout,*) 'The number of atoms in the the fragments is &
        equal to the number of atoms in the molecule.'
        Else
          write(iout,*) 'ERROR: The total number of atoms in the  &
        fragments is not the same as the total number in the   &
        molecule file. Please check the input files to correct this  &
        error.'
        EndIf
        
        Write(iout,*) ' '
        NBasisTotal=(F1_NBasis+F2_NBasis+F3_NBasis)/3
        If (NBasisTotal==M_NBasis) then
          write(iout,*) 'The number of basis functions in  &
        the fragments is equal to the number of atoms in the molecule.'
        Else
          write(iout,*) 'ERROR: The total number of basis in the  &
        fragments is not the same as the total number in the   &
        molecule file. Please check the input files to correct this  &
        error.'
        EndIf

        Write(iout,*) ' '
        NElectronsTotal=F1_NElectrons+F2_NElectrons+F3_NElectrons
        If (NElectronsTotal==M_NElectrons) then
          write(iout,*) 'The number of electrons in the  &      
         fragments is equal to the number of atoms in the molecule.'
        Else
          write(iout,*) 'ERROR: The total number of electrons in the  &
        fragments is not the same as the total number in the   &
        molecule file. Please check the input files to correct this &
        error.'
        EndIf
        
        Write(iout,*) ' '
        Write(iout,*) 'The total number of atoms is:',M_NAtoms 
        Write(iout,*) 'The total number of basis functions is:',  &
        M_NBasis 
        Write(iout,*) 'The total number of electrons is:',  &
        NElectronsTotal 

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Composite Matricies: Make the composite coefficent matrix and 
!       fock matrix. 

!       Add the three fragment matricies together and print them.
        Coefficent_Composite=F1_Coeff+F2_Coeff+F3_Coeff
!*        Call Coefficent_Composite%print(iout, 'Composite Matrix')

        Fock_AO_Composite=F1_Fock+F2_Fock+F3_Fock
!*        Call Fock_AO_Composite%print(iout, 'Fock Composite AO')

!       Transform everything into the FMO basis

        Fock_FMO_Composite=transpose(Coefficent_Composite).dot. &
        Fock_AO_Composite.dot.Coefficent_Composite
!*        Call Fock_FMO_Composite%print(iout, 'Fock Composite FMO')

        Fock_FMO_Molecule=transpose(Coefficent_Composite).dot. &
        M_Fock.dot.Coefficent_Composite
!*        Call Fock_FMO_Molecule%print(iout, 'Fock Molecule FMO')
         
!       ______________________________________________________________________
!       ______________________________________________________________________
!       First Order Energy Correction
        Call First_Order_Energy_Correction%init(NBasisTotal,0)
        Do p=1,M_NBasis
          Call First_Order_Energy_Correction%put(Fock_FMO_Molecule%at(p,p)- &
          Fock_FMO_Composite%at(p,p),p)
        EndDo       
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        write(iout,*) 'The First Order Energy Correction'
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call First_Order_Energy_Correction%print(iout,'First order:')

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Second Order Energy Correction
!       Here we will form two intermediates to find the final energy
!       correction. This is because the sum in the final expression
!       excludes one value and we want to avoid this in the code.

!       Find the denominator of the energy correction        
        
        Call DeltaE%init(NBasisTotal,NBasisTotal,0)
        Do p=1,NBasisTotal 
          Do q=1,NBasisTotal
            Call DeltaE%put(Fock_FMO_Composite%at(p,p)- &
            Fock_FMO_Composite%at(q,q),p,q)
          EndDo
        EndDo

!       The diagonal of DeltaE will be zero.. since we cannot divide by
!       zero we will set the zeros to one and take care of the incorrect
!       values in a later step.

        Do p=1,NBasisTotal        
          DeltaE%at(p,p)=1
        EndDo

!       Now calcualte V on E. This will also be called the amplitude.
!       Then multiply it by the second V in the expression. 

        VonE=((Fock_FMO_Molecule-Fock_FMO_Composite)).ewd.DeltaE
        VVonE=VonE.ewp.(Fock_FMO_Molecule-Fock_FMO_Composite)

!      Set the diagonals to zero.
 
       Do p=1,NBasisTotal
          VVonE%at(p,p)=0
       EndDo

!       Calculate the second order correction by summing each of the
!       columns in the VVonE matrix and putting them into a vector that
!       is NBasis long 
        
        Call Second_Order_Energy_Correction%init(NBasisTotal,0)
!*        N=0
!*        Do p=1,NBasisTotal
!*          Do q=1,NBasisTotal
!*            N=N+VVonE%at(p,q)
!*          EndDo
!*          Call Second_Order_Energy_Correction%put(N,p)
!*        EndDo 
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        write(iout,*) 'The Second Order Energy Correction'
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Second_Order_Energy_Correction%print(iout,'Second order:')

!       ______________________________________________________________________
!       ______________________________________________________________________
!       Third Order Energy Correction

!       Define Denominator #1
        Call DoubleDeltaE%init(NBasisTotal,NBasisTotal,0)
        Do p=1,NBasisTotal
          Do q=1,NBasisTotal
            Do r=1,NBasisTotal
              Call DoubleDeltaE%put(((Fock_FMO_Composite%at(p,p)- &
              Fock_FMO_Composite%at(q,q))* &
              (Fock_FMO_Composite%at(p,p)- &
              Fock_FMO_Composite%at(r,r))),p,q)
            EndDo
          EndDo
        EndDo

        Do p=1,NBasisTotal
          DoubleDeltaE%at(p,p)=1
        EndDo

!       Denominator #2 will be the same as DeltaE in the second order
!       energy correction

!       Numerator #1: Three V multiplied by eachother. This will be done
!       by first multiplyinh two together to make VV, then multiplying
!       that by another V.
        VV=(Fock_FMO_Molecule-Fock_FMO_Composite)* & 
        (Fock_FMO_Molecule-Fock_FMO_Composite)                
        VVV=VV*(Fock_FMO_Molecule-Fock_FMO_Composite)

!       Numerator #2 is just the first order energy correction
!       multiplied by the second order energy correction.
        Call EoneEtwo%init(NBasisTotal,0)
        
        Do p=1,NBasisTotal
          Call EoneEtwo%put(First_Order_Energy_Correction%at(p)*  & 
          Second_Order_Energy_Correction%at(p),p)
        EndDo
!        Call EoneEtwo%print(iout,'EoneEtwo')
        
!       Term #1: This will be numerator one divided by denominator
!       number one
        TermOne=VVV.ewd.DoubleDeltaE        
        Do p=1,NBasisTotal
          TermOne%at(p,p)=0
        EndDo
!        Call TermOne%print(iout,'E3Term1')
        
!       Term #2: This will be numerator two divided by denominator
!       number two
!*        TermTwo=EoneEtwo/DeltaE
!*        Do p=1,NBasisTotal
!*          TermTwo%at(p,p)=0
!*        EndDo
!*        Call TermTwo%print(iout,'E3Term2')
        
!       Subtract the two terms to get the final third order energy
!       corrrection.

!*        TermOneminusTermTwo=TermOne-TermTwo
        
        Call Third_Order_Energy_Correction%init(NBasisTotal,0)
!*        N=0
!*        Do p=1,NBasisTotal
!*          Do q=1,NBasisTotal
!*            N=N+TermOneminusTermTwo%at(p,q)
!*          EndDo
!*          Call Third_Order_Energy_Correction%put(N,p)
!*        EndDo 
        Write(iout,*) ' '
        Write(iout,*) '**********************************************&
        ************************************************************'
        write(iout,*) 'The Third Order Energy Correction'
        Write(iout,*) '**********************************************&
        ************************************************************'
        Call Third_Order_Energy_Correction%print(iout,'Third order:')

!       ______________________________________________________________________
!       ______________________________________________________________________

        End Program PT_Bonding_Analysis

!       ********************************************************************
!       ********************************************************************
        Subroutine Basic_Information(CommandLine,MatrixFile,NAtoms,  &
        NBasis,NBasisUse,NElectrons,Fock,Overlap,Density,  &
        Coeff,Hamiltonian,Energy)
        Use MQC_Gaussian
        Use MQC_Algebra
        Use MQC_Files
        Use MQC_General
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
 
!       Varible Declaration
        Implicit None
        Type(MQC_Matrix) :: Sam
        Integer :: i,iout=6,NAtoms,NBasis,NBasisUse,NElectrons        
        Type(MQC_Scalar) :: HFEnergy,IPMEnergy,OneEnergy,TwoEnergy,FD,HD
        Type(MQC_Scalar) :: Half,temp
        Type(MQC_Gaussian_Unformatted_Matrix_File) :: MatrixFile
        Type(MQC_Matrix) :: Fock,Overlap,Density,Coeff,Hamiltonian
        Type(MQC_Vector) :: Energy
        Character(len=256) :: CommandLine
        Half=0.5
        
!       Load the unformatted matrix file
        Call MatrixFile%load(CommandLine)
!
!       Test if the job was restricted or unrestricted.
        Write(iout,*) 'Is restricted? ', MatrixFile%isrestricted()
        Write(iout,*) 'Is unrestricted? ', MatrixFile%isunrestricted()
        
!       Find and print the number of atoms, basis functions and used
!       basis functions.
 
        NElectrons=MatrixFile%getval('NElectrons')
        NAtoms=MatrixFile%getval('NAtoms')
        NBasis=MatrixFile%getval('NBasis')
        NBasisUse=MatrixFile%getval('NBasisUse')
        
        Write(iout,*) 'NAtom= ', NAtoms
        Write(iout,*) 'NBasis= ', NBasis
        Write(iout,*) 'NElectrons=', NElectrons
        Write(iout,*) 'NBasisUse= ', NBasisUse

!       Load and save all of the matricies
        Call MatrixFile%getarray('overlap',Overlap)
!        Call Overlap%print(iout,'Overlap')
        Call MatrixFile%getarray('alpha fock matrix',Fock)
!        Call Fock%print(iout,'Fock')
        Call MatrixFile%getarray('alpha density matrix',Density)
        Density=Density+Density
        Call MatrixFile%getarray('alpha MO coefficients',Coeff)
        Call MatrixFile%getarray('core hamiltonian alpha',Hamiltonian)
        Call MatrixFile%getarray('alpha orbital energies',  &
        vectorout=Energy)
!
!       Calculate the number of electrons in the system       
        Call MQC_Print(contraction(Density,Overlap),iout,'<PS>= ') 
!
!       Calculate the Hartree Fock Energy
        FD=contraction(Fock,Density)
        HD=contraction(Hamiltonian,Density)
        HFEnergy=Half*(FD+(Half*HD))
        Call HFEnergy%print(iout,'The Hartree Fock Energy is: ')
!
!       Calculate the independent particle model energy
        IPMEnergy=0.0
        Do i=1,NElectrons/2
          temp=Energy%at(i)
          IPMEnergy=IPMEnergy+(temp+temp)
        EndDo
        Call IPMEnergy%print(iout,'The Independent Particle Model  &
        Energy is: ')
!
!       Calculate the one electron energy
        OneEnergy=contraction(Fock,Density)-  &
        contraction(Hamiltonian,Density)
        Call OneEnergy%print(iout,'The One Electron Energy is: ') 

!       Calculate the two electron energy
        TwoEnergy=contraction(Density,Hamiltonian)
        Call TwoEnergy%print(iout,'The Two Electron Energy is: ')
!
!       Close the matrix file
        Call MatrixFile%closefile()
!
        End Subroutine Basic_Information
!       ********************************************************************
!       ********************************************************************
       
