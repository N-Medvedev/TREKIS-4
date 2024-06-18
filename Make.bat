:: This file was written by N.Medvedev 
:: in 2018-2022
:: -----------------------------------------------------
:: Brief notes on cmd programming :  https://ss64.com/nt/if.html
:: -----------------------------------------------------

@echo off
setlocal EnableDelayedExpansion

:: Go into the directory with source files:
cd Source_files

:: read argument from user
   SET arg1=%1
   
   SET "Starline=************************************************************************************"
   echo %Starline%
   echo Started compilation: %date% %time%
   
:: in case of an empty argument, assume: no debug
   IF [%1] == [] (
      SET arg1=NODEBUG
   )
   
   :: shorthand expressions for debug options (convert them into DEBUGOMP option):
   IF /I %arg1%==DBG (
      SET arg1=DEBUGOMP
   )
   IF /I %arg1%==DB (
      SET arg1=DEBUGOMP
   )
:: shorthand expressions for no-optimization and no debug options (fast compiling):
   IF /I %arg1%==SLOW (
      SET arg1=FAST
   )


:: create a list of all files to be compiled
   SET "List_of_files=Universal_constants.f90 Geometries.f90 Little_subroutines.f90 Objects.f90 Variables.f90 Gnuplotting.f90 Periodic_table.f90 Dealing_with_files.f90 Dealing_with_EADL.f90 Dealing_with_cdf_files.f90 Read_numerical_parameters.f90 Relativity.f90 SHI_charge_state.f90 CS_general_tools.f90 CS_integration_limits.f90 CDF_Ritchi.f90 CDF_Mermin.f90 CDF_delta.f90 Dealing_with_XYZ_files.f90 MD_general_tools.f90 Dealing_with_LAMMPS.f90 Read_MD_parameters.f90 Read_input_data.f90 Dealing_with_DOS.f90 Initial_conditions.f90 CDF_get_from_data.f90 CS_photons_pair_creation.f90 CS_photons_Compton.f90 CS_photons_Rayleigh.f90 CS_electrons_inelastic.f90 CS_electrons_elastic.f90 CS_electrons_Bremsstrahlung.f90 CS_ion_inelastic.f90 CS_positrons_inelastic.f90 CS_positrons_elastic.f90 CS_positrons_Bremsstrahlung.f90 CS_positrons_annihilation.f90 CS_holes_elastic.f90 CS_holes_inelastic.f90 Output.f90 MC_general_tools.f90 MC_data_analysis.f90 MC_photon.f90 MC_electron.f90 MC_positron.f90 MC_hole.f90 MC_SHI.f90 MC.f90 MD_data_analysis.f90 MD_Pot_Simple.f90 MD_Pot_ZBL.f90 MD_Pot_Buck.f90 MD_Pot_SW.f90 MD_Pot_Coulomb.f90 MD.f90 TREKIS_main_file.f90"
   
   
:: Compile all those listed files with options provided by the user   
   IF /I "%arg1%"=="debug" (
       SET Program_name=TREKIS_DEBUG.exe
       
       echo %Starline%
       echo Compiling with DEBUG option, no OpenMP or optimizations are included
       echo Started at: %date% %time%
       echo %Starline%
       
       ifx.exe -c  /F9999999999 /QxHost /QaxAVX /Qmkl /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  %List_of_files% /standard-semantics
       
       echo %Starline%
       echo Assembling compiled modules into an executable: %Program_name%
       echo Started at: %date% %time%
       echo %Starline%
       
       ifx.exe /F9999999999 /QxHost /QaxAVX /Qmkl /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  *.obj   /exe:!Program_name! /standard-semantics

   ) ELSE (
      IF  /I "%arg1%"=="debugomp" (
         
         SET Program_name=TREKIS_DEBUG_OMP.exe
         
         echo %Starline%
         echo Compiling with DEBUGOMP option, OpenMP but no optimizations are included
         echo Started at: %date% %time%
         echo %Starline%
         
         ifx.exe -c  /F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  %List_of_files% /standard-semantics

         echo %Starline%
         echo Assembling compiled modules into an executable: %Program_name%
         echo Started at: %date% %time%
         echo %Starline%
           
         ifx.exe /F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  *.obj   /exe:!Program_name! /standard-semantics

      ) ELSE (
         IF /I %arg1%==FAST (
            SET Program_name=TREKIS_OMP.exe

            echo %Starline%
            echo Compiling with FAST option, OpenMP, no optimizations, no debug
            echo Started at: %date% %time%
            echo %Starline%

            :: List compiler options
            ifx.exe -c /F9999999999 /fpp /Qopenmp /Qmkl=parallel /real-size:64 /Od /fpe:0 /fp:precise /Qvec %List_of_files% /standard-semantics

            echo %Starline%
            echo Assembling compiled modules into an executable: %Program_name%
            echo Started at: %date% %time%
            echo %Starline%

            ifx.exe /F9999999999 /fpp /Qopenmp /Qmkl=parallel /real-size:64 /Od /fpe:0 /fp:precise /Qvec *.obj   /exe:!Program_name! /standard-semantics

            del *.pdb
         ) ELSE (
            SET Program_name=TREKIS.exe
      
            echo %Starline%
            echo Compiling for release, OpenMP and optimizations are included
            echo Started at: %date% %time%
            echo %Starline%
         
            ifx.exe -c /F9999999999 /fpp /Qopenmp /Qmkl=parallel  /Ot /O3 /Qvec /Qipo /real-size:64 /QxHost /QaxAVX  %List_of_files% /standard-semantics
          
            echo %Starline%
            echo Assembling compiled modules into an executable: %Program_name%
            echo Started at: %date% %time%
            echo %Starline%
          
            ifx.exe /F9999999999 /fpp /Qopenmp /Qmkl=parallel  /Ot /O3 /Qvec /Qipo /real-size:64 /MP:4 /QxHost /QaxAVX   *.obj   /exe:!Program_name! /standard-semantics
            del *.pdb
         )
      )
   )
   
echo %Starline%
echo The program %Program_name% was created at %date% %time%
echo %Starline%
   
:: Remove files that are no longer needed
del *.obj *.mod

:: Go back into the parent directory from the source files:
cd ..\

:: Copy exe file from the source directory into the parent directory:
xcopy source_files\%Program_name% %Program_name%* /Y /Q

:: Delete the exe file from the soure directory:
del source_files\%Program_name%


:: %--------------------------------------
:: NOTES:
:: 1) Options  /inline or /Qparallel  give *Internal compiler error*
:: 2) Option /MP[:n] can only be used for compiling objects, not f90 files, because for those the order matters (links to each other)
:: 
:: %--------------------------------------
::  OTHER OPTIONS:

:: ifx.exe -c /F9999999999 /fpp /Qopenmp /D OMP_inside /Qmkl=parallel /O3 /Qvec /Qipo /real-size:64

:: ifx.exe -c /F9999999999 /fpp /Qparallel /Qopt-prefetch=2 /Qopenmp /D OMP_inside /Qmkl=parallel /O3 /Qvec /Qipo /libs:static /real-size:64 /threads

:: ifx.exe -c /F9999999999 /fpp /Qmkl=parallel /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec
