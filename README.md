# SNbone
This application targets the primary computational solve burden of a SN,
continuous finite element based transport equation solver

## Branch Features

In the master branch, worksharing is explicitly coded in SPMD. and work is distributed using thread IDs.  In this branch, the worksharing is handled through the OpenMP runtime, using OMP REDUCTION and OMP SINGLE pragmas; the thread ID is not explicitly needed for worksharing.  

This has been implemented for the FGMRES_Threaded subroutine and the AVE3 method.  
