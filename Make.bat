:: 0000000000000000000000000000000000000000000000000000000000000
:: This file is part of TREKIS-4
:: available at: https://github.com/N-Medvedev/TREKIS-4
:: This file was written by N.Medvedev
:: in 2018-2025
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


:: Read the ordered source list shared with the GNU Make build.
   SET "List_of_files="
   FOR /F "usebackq eol=# tokens=*" %%f IN ("sources.list") DO SET "List_of_files=!List_of_files! %%f"
   
   
:: Compile all those listed files with options provided by the user   
   IF /I "%arg1%"=="debug" (
       SET Program_name=TREKIS_DEBUG.exe
       
       echo %Starline%
       echo Compiling with DEBUG option, no OpenMP or optimizations are included
       echo Started at: %date% %time%
       echo %Starline%
       
       ifx.exe -c  /F9999999999 /QxHost /QaxAVX  /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  %List_of_files% /standard-semantics
       
       echo %Starline%
       echo Assembling compiled modules into an executable: %Program_name%
       echo Started at: %date% %time%
       echo %Starline%
       
       ifx.exe /F9999999999 /QxHost /QaxAVX  /fpp /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  *.obj   /exe:!Program_name! /standard-semantics

   ) ELSE (
      IF  /I "%arg1%"=="debugomp" (
         
         SET Program_name=TREKIS_DEBUG_OMP.exe
         
         echo %Starline%
         echo Compiling with DEBUGOMP option, OpenMP but no optimizations are included
         echo Started at: %date% %time%
         echo %Starline%
         
         ifx.exe -c  /F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp  /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  %List_of_files% /standard-semantics

         echo %Starline%
         echo Assembling compiled modules into an executable: %Program_name%
         echo Started at: %date% %time%
         echo %Starline%
           
         ifx.exe /F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp  /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec  *.obj   /exe:!Program_name! /standard-semantics

      ) ELSE (
         IF /I %arg1%==FAST (
            SET Program_name=TREKIS_OMP.exe

            echo %Starline%
            echo Compiling with FAST option, OpenMP, no optimizations, no debug
            echo Started at: %date% %time%
            echo %Starline%

            :: List compiler options
            ifx.exe -c /F9999999999 /fpp /Qopenmp  /real-size:64 /Od /fpe:0 /fp:precise /Qvec %List_of_files% /standard-semantics

            echo %Starline%
            echo Assembling compiled modules into an executable: %Program_name%
            echo Started at: %date% %time%
            echo %Starline%

            ifx.exe /F9999999999 /fpp /Qopenmp  /real-size:64 /Od /fpe:0 /fp:precise /Qvec *.obj   /exe:!Program_name! /standard-semantics

            del *.pdb
         ) ELSE (
            SET Program_name=TREKIS.exe
      
            echo %Starline%
            echo Compiling for release, OpenMP and optimizations are included
            echo Started at: %date% %time%
            echo %Starline%
         
            ifx.exe -c /F9999999999 /fpp /Qopenmp   /Ot /O3 /Qvec /Qipo /real-size:64 /QxHost /QaxAVX  %List_of_files% /standard-semantics
          
            echo %Starline%
            echo Assembling compiled modules into an executable: %Program_name%
            echo Started at: %date% %time%
            echo %Starline%
          
            ifx.exe /F9999999999 /fpp /Qopenmp   /Ot /O3 /Qvec /Qipo /real-size:64 /MP:4 /QxHost /QaxAVX   *.obj   /exe:!Program_name! /standard-semantics
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

:: ifx.exe -c /F9999999999 /fpp /Qopenmp /D OMP_inside  /O3 /Qvec /Qipo /real-size:64

:: ifx.exe -c /F9999999999 /fpp /Qparallel /Qopt-prefetch=2 /Qopenmp /D OMP_inside  /O3 /Qvec /Qipo /libs:static /real-size:64 /threads

:: ifx.exe -c /F9999999999 /fpp  /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec
