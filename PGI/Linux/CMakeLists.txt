### nTRACER CMake Script

set (CMAKE_Fortran_COMPILER mpif90)

project (
	nTRACER
)

enable_language(Fortran)
cmake_minimum_required(VERSION 3.12)

set (MKL_DIR "/proj/intel/compilers_and_libraries/linux/mkl")
set (MKL_LIB "-L${MKL_DIR}/lib/intel64 -R${MKL_DIR}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_vml_avx")

set (CMAKE_Fortran_FLAGS "-r8 -fast -Mfma -Mpreprocess -mp -acc -Mcuda=fastmath,ptxinfo -ta=tesla,cuda10.2,cc75 -Mcudalib=cublas,cusparse")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DVAR=2")

add_definitions(-D__PGI)

OPTION(MKL "MKL Specification" OFF)
IF(MKL)
    ADD_DEFINITIONS(-D__INTEL_MKL)
	message("==== Make with MKL")
else()
    message("==== Make without MKL")
ENDIF(MKL)

if(MKL)
    include_directories(
	    ../../src
		${MKL_DIR}/include
	)
else()
	include_directories (
		../../src
	)
endif()

file (GLOB nTRACER_FILES ../../src/*.*)

add_executable (
	nTRACER
	${nTRACER_FILES}
)

if(MKL)
	target_link_libraries (
		nTRACER
		${MKL_LIB}
	)
else()
	target_link_libraries (
		nTRACER
		lapack
		blas
	)
endif()