rem script to compile and link a generic unit test
rem %1 comes from the command line

rem @echo off
cl.exe /nologo /MTd /W3 /EHsc /Od /D_WIN32 /D_DEBUG /D_CONSOLE /D_MBCS /RTC1 /c %1.cpp
link.exe /nologo /subsystem:console /debug /out:"%1.exe" %1.obj
del %1.obj
del %1.ilk
del %1.pdb
