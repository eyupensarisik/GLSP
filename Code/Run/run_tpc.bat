@ECHO OFF

set BINDIR=x64\Debug
set DATADIR=..\..\Data\%1
set PARAMFILE=..\Run\Parameters\%2.txt
set RESULTDIR=Results\%1_%2_%3_%4
set SUMMARYFILE=LSSP_%1_%2_%3_%4.csv
set EXECUTABLE=GLSP.exe

set NEW_P=%3
set NEW_T=%4

IF EXIST %RESULTDIR% rmdir /S /Q %RESULTDIR%
md %RESULTDIR%

ECHO Name,^
P,^
T,^
TwoPhase_Callback_Cut,^
TwoPhase_CPU,^
TwoPhase_Obj,^
TwoPhase_Callback_CPU,^
TwoPhase_Callback_Call,^
 > %RESULTDIR%\%SUMMARYFILE%

FOR %%f IN (%DATADIR%\*.dat) DO (
	ECHO Started solving %%f with %PARAMFILE%
	REM ECHO %BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_P% %NEW_T% 
	%BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_P% %NEW_T% > %RESULTDIR%\%%~nxf.log & IF ERRORLEVEL 1 GOTO:EOF
)