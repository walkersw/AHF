@echo off
rem %comspec% /k ""C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat"" x86
call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86
set include=%include%;C:\FILES\FELICITY\Static_Codes\AHF\src_code\
set include=%include%;C:\FILES\C++\eigen-3.3.7\
cd "C:\FILES\FELICITY\Static_Codes\AHF\Unit_Test_src"
start "" %windir%\system32\cmd.exe
exit
