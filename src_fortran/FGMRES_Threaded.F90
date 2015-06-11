!----------------------------------------------------------------------------------------------------
! Method_FGMRES.F90
! Driver of the FGMRES Routine to solve a linear system Ax = b
! To implement GMRES I followed the algorithm proposed by Y. Saad and M Schultz in GMRES: A generalized minimal residual 
! algorithm for solving nonsymmetric linear systems, SIAM J. Sci. Stat. Comput., Vol 7, No 3, 
! July 1986
! Then I edited it to add a varying right preconditioner (Flexible GMRES) following the paper 
! 'A set of flexible GMRES routines for real and complex arithmetics on high Performance computers', 
! V. Fraysse, L. Giraud, S. Gratton, CERFACS Technical Report TR/PA/06/09
!----------------------------------------------------------------------------------------------------
!  Copyright(c) 2005 Argonne National Laboratory
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! Externally threaded version of the subroutine
!----------------------------------------------------------------------------------------------------
! Barriers if needed must be erected around all user subroutines
! Incoming vectors are shared and users must put barriers in place to prevent issues
! This routine supports combined MPI and OpenMP so long as MyThreadID=0 matches the MPI rank
! 
! Because this subroutine is called by each thread, all local variables will be unique to each thread
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
RECURSIVE SUBROUTINE FGMRES_Threaded(User_Krylov,Solution,RightHandSide,                    &

                     GuessIsNonZero,ReasonForConvergence,IterationCount,  &
                     ResidualNorm,VectorNorm, &
                     Apply_A,Apply_PC)
#ifdef WITHOMP
  USE OMP_LIB
#endif
USE CommonBlock
USE Method_Krylov
IMPLICIT NONE
#ifdef USEMPI
#include "mpif.h"
#endif
!#define Local_Debug_FGMRES_Driver
! Passed In
TYPE (Method_Krylov_Type)  User_Krylov      ! The allocated structure
PROTEUS_Real Solution(User_Krylov%Local_Owned)      ! In: an initial guess/ Out: the updated approximation of the solution
PROTEUS_Real RightHandSide(User_Krylov%Local_Owned) ! Right hand side of the equation (= b vector)
! Thread shared variables
PROTEUS_Log  GuessIsNonZero             ! If true the vector stored in Solution is used as the initial guess, otherwise a null vector is used
PROTEUS_Int  ReasonForConvergence       ! Divergence (<0), MaxIterationCount (=0), Convergence (>0)
PROTEUS_Int  IterationCount             ! Count the total number of inner iteration (number of time we apply A and orthogonalize)
PROTEUS_Real ResidualNorm,VectorNorm      ! Norm of the residual that is returned (a shared variable between the threads)
PROTEUS_Real HessenNorm_Shared(User_Krylov%BackVectors)
! Subroutines that are called
EXTERNAL  Apply_A                       ! The subroutine which applies the A matrix
EXTERNAL  Apply_PC                      ! The subroutine which applies the preconditioner

! Local Thread specific variables
PROTEUS_Int  Maximum_Outer              ! Maximum of outer iteration to stop even if we do not get convergence 
PROTEUS_Real LocalConst,Relative_Stop,Divergence_Stop ! These could be thread shared
PROTEUS_Real Cosinus,Sinus,aconst,bconst ! These are only needed on the root thread
PROTEUS_Int  I,J,K,Outer,Inner            ! Thes must be thread specific to avoid unnecessary barriers

!$OMP  PARALLEL &
!$OMP& default(none) &
!$OMP& shared(NumVertices, NumAngles, User_Krylov, Solution, RightHandSide, &
!$OMP&   NumThreads, GuessIsNonZero, &
!$OMP&   ReasonForConvergence, IterationCount, ParallelComm, ParallelRank, &
!$OMP&   ResidualNorm, VectorNorm, &
!$OMP&   HessenNorm_Shared) &
!$OMP& private(Maximum_Outer,LocalConst,Relative_Stop,Divergence_Stop, &
!$OMP&   Cosinus,Sinus,aconst,bconst,I,J,K,Outer,Inner)

! Grab the value from the structure
Maximum_Outer = User_Krylov%Maximum_Iterations/User_Krylov%BackVectors
IF (Maximum_Outer*User_Krylov%BackVectors .NE. User_Krylov%Maximum_Iterations) Maximum_Outer = Maximum_Outer + 1

#ifdef Local_Debug_FGMRES_Driver
   WRITE(Output_Unit,*)'NumThreads = ',NumThreads
!$OMP SINGLE
      WRITE(Output_Unit,'("[GMRES]...Initial guess ",I9)') User_Krylov%Local_Owned
      DO I = 1,User_Krylov%Local_Owned
         WRITE(Output_Unit,'("[GMRES]...solution(",I4,") = ",1PE13.6," RHS=",1PE13.6)') I, Solution(I),RightHandSide(I)
      END DO
!$OMP END SINGLE
#endif


IF (GuessIsNonZero) THEN ! Apply A to the initial guess
!$OMP DO
   DO I = 1, NumVertices*NumAngles
      User_Krylov%Basis(I,1) = 0.0d0
   END DO
!$OMP END DO
   ! Start: Compute the initial residual r0 = Ax0 - b and the first basis vector v1 = r0/||r0||
   CALL Apply_A(Solution,User_Krylov%Basis) ! (1,1)
!$OMP DO
   DO I = 1, NumVertices*NumAngles
      User_Krylov%Basis(I,1) = RightHandSide(I) - User_Krylov%Basis(I,1)
   END DO
!$OMP END DO
ELSE ! Zero solution if GuessIsNonZero .EQ. False
!========================================================================
!$OMP DO
   DO I = 1, NumVertices*NumAngles
      Solution(I) = 0.0d0
      User_Krylov%Basis(I,1) = RightHandSide(I)
   END DO
!$OMP END DO
!========================================================================
! An alternative is to use the WORKSHARE construct, but in this case, the
! performance improvements are trivial.
! !$OMP WORKSHARE
!       Solution = 0.0d0
!       User_Krylov%Basis(:,1) = RightHandSide(:)
! !$OMP END WORKSHARE
!========================================================================
END IF

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
      WRITE(Output_Unit,'("[GMRES]...Initial residual")')
      DO I = 1,User_Krylov%Local_Owned
         WRITE(Output_Unit,'("[GMRES]...Basis(",I4,") = ",1PE13.6)') I, User_Krylov%Basis(I,1)
      END DO
!$OMP END SINGLE
#endif

!$OMP SINGLE
ReasonForConvergence = 0
IterationCount = 0
ResidualNorm = 0.0d0
!$OMP END SINGLE

!$OMP DO REDUCTION(+:ResidualNorm)
DO I = 1, NumVertices*NumAngles
   ResidualNorm = ResidualNorm + User_Krylov%Basis(I,1) * User_Krylov%Basis(I,1)
END DO
!$OMP END DO

!$OMP SINGLE
ResidualNorm = DSQRT(ResidualNorm)
!$OMP END SINGLE

!Setup the relative tolerance limit
Relative_Stop = ResidualNorm * User_Krylov%Relative_Tolerance 
! Setup the divergence tolerance limit
Divergence_Stop = ResidualNorm *(1.0d0 + User_Krylov%Divergence_Tolerance)

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
      WRITE(Output_Unit,'("[GMRES]...Max Outer        = ",I14)') Maximum_Outer
      WRITE(Output_Unit,'("[GMRES]...Max Inner        = ",I14)') User_Krylov%BackVectors
      WRITE(Output_Unit,'("[GMRES]...Total Max Iter   = ",I14)') User_Krylov%Maximum_Iterations
      WRITE(Output_Unit,'("[GMRES]...Absolute target  = ",1PE13.6)') User_Krylov%Absolute_Tolerance
      WRITE(Output_Unit,'("[GMRES]...Relative target  = ",1PE13.6)') Relative_Stop
      WRITE(Output_Unit,'("[GMRES]...Divergence Stop  = ",1PE13.6)') Divergence_Stop
!$OMP END SINGLE
#endif

IF (ResidualNorm .EQ. 0.0d0) GOTO 1000 ! Yes, there is nothing to do because we have the exact answer

DO Outer = 1, Maximum_Outer
! It might be possible to workshare this in parallel...
!$OMP SINGLE
      DO J = 1,User_Krylov%BackVectors+1
         DO I = 1,User_Krylov%BackVectors
            User_Krylov%Hessenberg(I,J) = 0.0d0
         END DO
      END DO
      DO J = 1,User_Krylov%BackVectors
         User_Krylov%Modified_RHS(J) = 0.0d0
      END DO
      User_Krylov%Modified_RHS(1) = ResidualNorm ! No barrier is needed for this as I assume we need one just before and after Apply_PC and Apply_A
!$OMP END SINGLE

   ! For the initialization the Residual is in Basis(*,1) and its norm in ResidualNorm    
   LocalConst = 1.0d0 / ResidualNorm
!$OMP DO
   DO I = 1, NumVertices*NumAngles
      User_Krylov%Basis(I,1)    =  User_Krylov%Basis(I,1) * LocalConst
      User_Krylov%Basis(I,2)    = 0.0d0
      User_Krylov%PC_Basis(I,1) = 0.0d0
   END DO
!$OMP END DO

   ! Iterate: For j =1,2,...m do:
   !   h(i,j) = (Avj,vi), i = 1,...,j
   !   v*(j+1) = Av(j) - sum(h(i,j)*v(i),i=1,j)
   !   h(j+1,j) = ||v*(j+1)||
   !   v(j+1) = v*(j+1)/ h(j+1,j)
   DO Inner = 1, User_Krylov%BackVectors
      ! Apply the preconditionner to the basis vector and store it in a intermediate vector
      CALL Apply_PC(User_Krylov%Basis(1,Inner),User_Krylov%PC_Basis(1,Inner))
#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
         WRITE(Output_Unit,'("[GMRES]...At Outer ",I4," and Inner ",I4," some vectors:")') Outer,Inner
         DO I = 1, User_Krylov%Local_Owned
            WRITE(Output_Unit,'("[GMRES]...residual/norm=V(",I4,") = ",1PE13.6," PC*V(",I4,") = ",1PE13.6)') &
                 I, User_Krylov%Basis(I,Inner),I,User_Krylov%PC_Basis(I,Inner)
         END DO
!$OMP END SINGLE
#endif
      ! Apply A to the intermediate vector and store it in the next basis point before orthonormalization
      CALL Apply_A(User_Krylov%PC_Basis(1,Inner),User_Krylov%Basis(1,Inner+1))
#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
         WRITE(Output_Unit,'("[GMRES]...At Outer ",I4," and Inner ",I4," intermediate vectors:")') Outer,Inner
         DO I = 1, User_Krylov%Local_Owned
            WRITE(Output_Unit,'("[GMRES]...A*Z=V(",I4,") = ",1PE13.6)') I, User_Krylov%Basis(I,Inner+1)
         END DO
!$OMP END SINGLE
#endif

!$OMP SINGLE
HessenNorm_Shared = 0.0d0
!$OMP END SINGLE

!$OMP DO REDUCTION(+:HessenNorm_Shared)
   DO J=1, NumVertices*NumAngles
      DO I=1, Inner
         HessenNorm_Shared(I) = HessenNorm_Shared(I) +User_Krylov%Basis(J,Inner+1)*User_Krylov%Basis(J,I)
      END DO
   END DO
!$OMP END DO

!$OMP SINGLE
   DO I = 1, Inner
#ifdef USEMPI
         User_Krylov%Hessenberg(I,User_Krylov%BackVectors+1) = HessenNorm_Shared(I)
#else
         User_Krylov%Hessenberg(I,Inner) = HessenNorm_Shared(I)
#endif
   END DO
#ifdef USEMPI
         ! Reduce the dot product on the whole communicator space
         CALL MPI_ALLREDUCE(User_Krylov%Hessenberg(1,User_Krylov%BackVectors+1),User_Krylov%Hessenberg(1,Inner),  &
                            User_Krylov%Local_Owned,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
         CALL Basic_CheckError(Output_Unit,J,'MPI_ALLREDUCE Hessenberg_Column in (Basic_FillHessenberg)')
#endif
!$OMP END SINGLE

!=====================================================================================================

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
         WRITE(Output_Unit,'("[GMRES]...Hessenberg matrix K=",I3)') 
         DO I = 1,Inner
            WRITE(Output_Unit,'("[GMRES]...Hessenberg(",I3,") = ",100(1PE13.6,1X))') I, (User_Krylov%Hessenberg(I,J),J=1,Inner)
         END DO
!$OMP END SINGLE
#endif
      ! Compute an new orthogonal vector
      DO J = 1,Inner
!$OMP DO
         DO I = 1, NumVertices*NumAngles
            User_Krylov%Basis(I,Inner+1) = User_Krylov%Basis(I,Inner+1) - User_Krylov%Hessenberg(J,Inner)*User_Krylov%Basis(I,J)
         END DO
!$OMP END DO
      END DO

!========================================================================
!$OMP SINGLE
   VectorNorm = 0.0d0
!$OMP END SINGLE

!$OMP DO REDUCTION(+:VectorNorm)
      DO I = 1, NumVertices*NumAngles
         VectorNorm = VectorNorm + User_Krylov%Basis(I,Inner+1)*User_Krylov%Basis(I,Inner+1)
      END DO
!$OMP END DO

!$OMP SINGLE
#ifdef USEMPI
         LocalConst = VectorNorm
         CALL MPI_ALLREDUCE(LocalConst,VectorNorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
         CALL Basic_CheckError(Output_Unit,J,"MPI_ALLREDUCE for ResidualNorm in (Method_FGMRES)")
#endif
         VectorNorm = DSQRT(VectorNorm)
         User_Krylov%Hessenberg(Inner+1,Inner) = VectorNorm
         IF (VectorNorm .EQ. 0.0d0) THEN
            ResidualNorm = VectorNorm ! This was moved from below and the check by all below will cause an exit
         ELSE
            VectorNorm = 1.0d0 / VectorNorm
            ! An efficient way to compute the norm of the residual is to use a QR decomposition of the Hessenberg matrix
            ! We also solve the minimization problem
            ! Apply previous rotation to the last column of Hessenberg
            DO I = 1,Inner-1,1
               aconst  = User_Krylov%Hessenberg(I,Inner)
               bconst  = User_Krylov%Hessenberg(I+1,Inner) ! We cannot remove this because of the redefinition in the next two lines
               User_Krylov%Hessenberg(I,Inner)   = User_Krylov%Givens(I,1)*aconst - User_Krylov%Givens(I,2)*bconst
               User_Krylov%Hessenberg(I+1,Inner) = User_Krylov%Givens(I,2)*aconst + User_Krylov%Givens(I,1)*bconst
            END DO
            ! Compute the givens coefficients for this iteration
            ResidualNorm = User_Krylov%Hessenberg(Inner,Inner)
            LocalConst = User_Krylov%Hessenberg(Inner+1,Inner)
            ResidualNorm = DSQRT(ResidualNorm*ResidualNorm + LocalConst*LocalConst)
            LocalConst = 1.0d0 / ResidualNorm
            User_Krylov%Givens(Inner,1) =  User_Krylov%Hessenberg(Inner,Inner)   * LocalConst
            User_Krylov%Givens(Inner,2) = -User_Krylov%Hessenberg(Inner+1,Inner) * LocalConst
   
            ! Apply this rotation to the Hessenberg matrix and to the modified right hand side
            aconst  = User_Krylov%Hessenberg(Inner,Inner)
            bconst  = User_Krylov%Hessenberg(Inner+1,Inner) ! We cannot remove this because of the redefinition in the next two lines
            User_Krylov%Hessenberg(Inner,Inner)   = User_Krylov%Givens(Inner,1)*aconst - User_Krylov%Givens(Inner,2)*bconst
            User_Krylov%Hessenberg(Inner+1,Inner) = User_Krylov%Givens(Inner,2)*aconst + User_Krylov%Givens(Inner,1)*bconst
            aconst  = User_Krylov%Modified_RHS(Inner)
            bconst  = User_Krylov%Modified_RHS(Inner+1) ! We cannot remove this because of the redefinition in the next two lines
            User_Krylov%Modified_RHS(Inner)   = User_Krylov%Givens(Inner,1)*aconst - User_Krylov%Givens(Inner,2)*bconst
            User_Krylov%Modified_RHS(Inner+1) = User_Krylov%Givens(Inner,2)*aconst + User_Krylov%Givens(Inner,1)*bconst
            ! Now we have the residual norm at no extra cost
            ResidualNorm = DABS(User_Krylov%Modified_RHS(Inner+1))
            ! update the number of iteration before exit
            IterationCount         = IterationCount + 1
            User_Krylov%Iterations = User_Krylov%Iterations + 1
         END IF ! Vector norm = 0
!$OMP END SINGLE

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
         WRITE(Output_Unit,'("[GMRES]...VectorNorm=",1PE13.6)') VectorNorm
         DO I = 1, User_Krylov%Local_Owned
            WRITE(Output_Unit,'("[GMRES]...New XX(",I4,") = ",1PE13.6)') I, User_Krylov%Basis(I,Inner+1) * VectorNorm
         END DO
         DO I = 1,Inner+1
            WRITE(Output_Unit,'("[GMRES]...Updated Hessenberg(",I3,") = ",100(1PE13.6,1X))') I, (User_Krylov%Hessenberg(I,J),J=1,Inner)
         END DO
         DO I = 1,Inner+1
            WRITE(Output_Unit,'("[GMRES]...Rotated Hessenberg(",I3,") = ",100(1PE13.6,1X))') I, (User_Krylov%Hessenberg(I,J),J=1,Inner)
         END DO
         DO I = 1,Inner
            WRITE(Output_Unit,'("[GMRES]...Givens(",I3,") = ",100(1PE13.6,1X))') I,User_Krylov%Givens(I,1),User_Krylov%Givens(I,2)
         END DO
   
         DO I = 1,Inner+1
            WRITE(Output_Unit,'("[GMRES]......RHS(",I4,") = ",1PE13.6)') I,User_Krylov%Modified_RHS(I)
         END DO
         WRITE(Output_Unit,'("[GMRES]...At Outer ",I4," and Inner ",I4," residual = ",1PE13.6)') Outer,Inner,ResidualNorm
!$OMP END SINGLE
#endif
      ! if Vector Norm == 0, we reach the happy break down.
      IF (VectorNorm .EQ. 0.0d0) EXIT ! Inner iteration loop

      IF (Inner .EQ. User_Krylov%BackVectors) THEN ! We have nothing left to initialize so...
!$OMP DO
         DO I = 1, NumVertices*NumAngles
            User_Krylov%Basis(I,Inner+1) = User_Krylov%Basis(I,Inner+1) * VectorNorm
         END DO
!$OMP END DO
      ELSE
!$OMP DO
         DO I = 1, NumVertices*NumAngles
            User_Krylov%Basis(I,Inner+1)    = User_Krylov%Basis(I,Inner+1) * VectorNorm
            User_Krylov%Basis(I,Inner+2)    = 0.0d0
            User_Krylov%PC_Basis(I,Inner+1) = 0.0d0
         END DO
!$OMP END DO
      END IF

      ! Test if the Residual meets the stopping criterion
      IF ((ResidualNorm .LT. User_Krylov%Absolute_Tolerance) .OR. (ResidualNorm .LT. Relative_Stop)) EXIT

      ! Test if we diverged
      IF (ResidualNorm .GT. Divergence_Stop) EXIT

      ! Test if the iteration limit was hit
      IF (IterationCount .GE. User_Krylov%Maximum_Iterations) EXIT
   END DO ! Inner
   IF (Inner .GT. User_Krylov%BackVectors) Inner = User_Krylov%BackVectors
   ! Form the approximate solution x_m = x_0 + V_m*y_m
   ! where y_m minimizes ||Beta.e_1 - H_m.y||, y in R^m
   ! All we need is the matrix H_m and the vectors (v1,...vm)
   ! We only need to compute h(i,m) = (A*vm, vi)
#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
      WRITE(Output_Unit,'("[GMRES]...Exit inner loop with inner = ",I4)') Inner
!$OMP END SINGLE
#endif
   IF (ResidualNorm .GT. Divergence_Stop) THEN
      ReasonForConvergence = -4 ! KSP_DIVERGED_DTOL
      EXIT
   END IF

   ! Compute Ym by solving Ym = R-1 * G
   ! R is an upper triangular matrix, the modified Hessenberg matrix
! OMP BARRIER  ! This barrier is not needed as both instances of exiting out of the inner loop will have consistent data
!$OMP SINGLE
      DO I = Inner, 1, -1
         DO J = I+1,Inner
            User_Krylov%Modified_RHS(I) = User_Krylov%Modified_RHS(I) - User_Krylov%Hessenberg(I,J)*User_Krylov%Modified_RHS(J)
         END DO
         User_Krylov%Modified_RHS(I) = User_Krylov%Modified_RHS(I) / User_Krylov%Hessenberg(I,I) ! This cannot be moved
      END DO
#ifdef Local_Debug_FGMRES_Driver
      DO I = 1,User_Krylov%Local_Owned
         WRITE(Output_Unit,'("[GMRES]...Exiting PC Basis(",I3,") = ",100(1PE13.6,1X))') I, (User_Krylov%PC_Basis(I,J),J=1,Inner)
      END DO
      DO I = 1,User_Krylov%Local_Owned
         WRITE(Output_Unit,'("[GMRES]...Exiting Basis(",I3,") = ",100(1PE13.6,1X))') I, (User_Krylov%Basis(I,J),J=1,Inner+1)
      END DO
      DO I = 1,Inner+1
         WRITE(Output_Unit,'("[GMRES]...Exiting RHS(",I4,") = ",1PE13.6)') I,User_Krylov%Modified_RHS(I)
      END DO
#endif
!$OMP END SINGLE

   ! Update the solution Xm = X0 + Zm*Y
   DO J = 1,Inner
   ! WORKSHARE also works here, but I found it was slower
!$OMP DO
      DO I = 1, NumVertices*NumAngles
         Solution(I) = Solution(I) + User_Krylov%PC_Basis(I,J)*User_Krylov%Modified_RHS(J)
      END DO
!$OMP END DO
   END DO

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
      WRITE(Output_Unit,'("[GMRES]...Solution after outer iteration ",I4)') Outer
      DO I = 1,User_Krylov%Local_Owned
         WRITE(Output_Unit,'("[GMRES]...X(",I4,") = " 1PE13.6)') I,Solution(I)
      END DO
!$OMP END SINGLE
#endif

   ! Initialize the basis so the user does not have to do it
!$OMP DO
   DO I = 1, NumVertices*NumAngles
      User_Krylov%Basis(I,1) = 0.0d0
   END DO
!$OMP END DO

   ! Compute r_m = f- A x_m and store it in Scratch_Basis 
   CALL Apply_A(Solution,User_Krylov%Basis)
   
!============================================================================================
!$OMP SINGLE
   ResidualNorm = 0.0d0
!$OMP END SINGLE

!$OMP DO REDUCTION(+:ResidualNorm)
   DO I = 1, NumVertices*NumAngles
      User_Krylov%Basis(I,1) = RightHandSide(I) - User_Krylov%Basis(I,1)
      ResidualNorm = ResidualNorm + User_Krylov%Basis(I,1)*User_Krylov%Basis(I,1)
   END DO
!$OMP END DO

!$OMP SINGLE
#ifdef USEMPI
      LocalConst = ResidualNorm
      CALL MPI_ALLREDUCE(LocalConst,VectorNorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
      CALL Basic_CheckError(Output_Unit,J,"MPI_ALLREDUCE for ResidualNorm in (Method_FGMRES)")
#endif
      ResidualNorm = DSQRT(ResidualNorm)
!$OMP END SINGLE

#ifdef Local_Debug_FGMRES_Driver
!$OMP SINGLE
   WRITE(Output_Unit,'("[GMRES]...At FinalCheck residual = ",1PE13.6)') ResidualNorm
!$OMP END SINGLE
#endif
       ! Test if the ResidualNorm meet the stopping criterion
    IF (ResidualNorm .LE. 0.0D0) THEN ! Happy breakdown convergence
       ReasonForConvergence = 5 ! KSP_CONVERGED_BREAKDOWN
       EXIT
    ELSE IF (ResidualNorm .LT. User_Krylov%Absolute_Tolerance) THEN
       ReasonForConvergence = 3 ! KSP_CONVERGED_ATOL 
       EXIT
    ELSE IF (ResidualNorm .LT. Relative_Stop) THEN
       ReasonForConvergence = 2 ! KSP_CONVERGED_RTOL
       EXIT
    ELSE IF (IterationCount .GE. User_Krylov%Maximum_Iterations) THEN
       ReasonForConvergence = 1 ! KSP_HIT_ITERATION_LIMIT
       EXIT
    ELSE
      IF (ParallelRank .EQ. 0) THEN
!$OMP SINGLE
         WRITE(Output_Unit,'("[SN-KERNEL]...Outer ",I6," after ",I6," inners has residual = ",1PE13.6)') Outer,Inner,ResidualNorm
!$OMP END SINGLE
      END IF
    END IF
END DO ! Outer

1000 CONTINUE

!$OMP END PARALLEL

END SUBROUTINE FGMRES_Threaded
