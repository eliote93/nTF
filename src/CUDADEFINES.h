#define		GPU_MAX_GROUP			47

#define		GPU_MAX_ORDER			9
#define		GPU_MAX_AZI				32
#define		GPU_MAX_POLAR 			4
#define		GPU_MAX_POLAR_1D		10

#define		GPU_PRECISION			    4
#define   GPU_XS_PRECISION      8
#define   GPU_SM_PRECISION      8
#define   GPU_PNUM_PRECISION    8
#define   GPU_RES_PRECISION     4
#define		GPU_FLUX_PRECISION		4
#define		GPU_SOURCE_PRECISION	4
#define		GPU_CMFD_PRECISION		4
!#define		GPU_CMFD_PRECISION		8
#define     GPU_NODAL_PRECISION     4
#define		P0_BLOCK_SIZE			16
#define		PN_BLOCK_SIZE			8
#define   SG_BLOCK_SIZE     4
#define   XS_NUM_BLOCK     4
#define   XS_BLOCK_SIZE     8

#define   XS_BLOCK_GDIM     16
#define   XS_BLOCK_RDIM     4

#define   FSP_BLOCK_GDIM    4
#define   FSP_BLOCK_RDIM    16

#define   EFF_BLOCK_RDIM    8
#define   EFF_BLOCK_IDIM    4
#define   EFF_BLOCK_GDIM    4

#if (GPU_CMFD_PRECISION == 4)

#define		CSR_CMFD_PRECISION			CSR_MIXED
#define		BSR_CMFD_PRECISION			BSR_MIXED
#define		cublasPcopy_v2				cublasScopy_v2
#define		cublasPaxpy_v2				cublasSaxpy_v2
#define		cublasPdot_v2				cublasSdot_v2
#define		cublasPnrm2_v2				cublasSnrm2_v2
#define		cublasPscal_v2				cublasSscal_v2
#define		cublasPgemv_v2				cublasSgemv_v2
#define		cublasPgemm_v2				cublasSgemm_v2
#define		cublasPsymv_v2				cublasSsymv_v2
#define		cublasPcsrsv_solve			cublasScsrsv_solve
#define		cublasPgetrfBatched			cublasSgetrfBatched
#define		cublasPgetriBatched			cublasSgetriBatched
#define		cusparsePcsrmv				cusparseScsrmv
#define		cusparsePcsrsv_analysis		cusparseScsrsv_analysis
#define		cusparsePcsrsv_solve		cusparseScsrsv_solve
#define		cusparsePcsrilu0			cusparseScsrilu0
#define		cusparsePbsrsv2_bufferSize	cusparseSbsrsv2_bufferSize
#define		cusparsePbsrsv2_analysis	cusparseSbsrsv2_analysis
#define		cusparsePbsrsv2_solve		cusparseSbsrsv2_solve

#elif (GPU_CMFD_PRECISION == 8)

#define		CSR_CMFD_PRECISION			CSR_DOUBLE
#define		BSR_CMFD_PRECISION			BSR_DOUBLE
#define		cublasPcopy_v2				cublasDcopy_v2
#define		cublasPaxpy_v2				cublasDaxpy_v2
#define		cublasPdot_v2				cublasDdot_v2
#define		cublasPnrm2_v2				cublasDnrm2_v2
#define		cublasPscal_v2				cublasDscal_v2
#define		cublasPgemv_v2				cublasDgemv_v2
#define		cublasPgemm_v2				cublasDgemm_v2
#define		cublasPsymv_v2				cublasDsymv_v2
#define		cublasPcsrsv_solve			cublasDcsrsv_solve
#define		cublasPgetrfBatched			cublasDgetrfBatched
#define		cublasPgetriBatched			cublasDgetriBatched
#define		cusparsePcsrmv				cusparseDcsrmv
#define		cusparsePcsrsv_analysis		cusparseDcsrsv_analysis
#define		cusparsePcsrsv_solve		cusparseDcsrsv_solve
#define		cusparsePcsrilu0			cusparseDcsrilu0
#define		cusparsePbsrsv2_bufferSize	cusparseDbsrsv2_bufferSize
#define		cusparsePbsrsv2_analysis	cusparseDbsrsv2_analysis
#define		cusparsePbsrsv2_solve		cusparseDbsrsv2_solve

#endif

#if (GPU_NODAL_PRECISION == 4)

#define		CSR_NODAL_PRECISION			CSR_FLOAT
#define		cusparseCsrmv				cusparseScsrmv
#define		cublasScal				  cublasSscal_v2
#define		cublasAxpy				  cublasSaxpy_v2
#define		cublasCopy				  cublasScopy_v2

#elif (GPU_NODAL_PRECISION == 8)

#define		CSR_NODAL_PRECISION			CSR_DOUBLE
#define		cusparseCsrmv				cusparseDcsrmv
#define		cublasScal				  cublasDscal_v2
#define		cublasAxpy				  cublasDaxpy_v2
#define		cublasCopy				  cublasDcopy_v2

#endif
