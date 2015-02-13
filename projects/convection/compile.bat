@echo off

takeown /f a.exe
rm a.exe
gfortran config.f90 ../../modules/file_io.f90 main.f90
a.exe