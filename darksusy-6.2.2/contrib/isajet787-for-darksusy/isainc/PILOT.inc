#define DOUBLE_X
#define STDIO_X
#define MOVEFTN_X
#define RANLUX_X
#define NOCERN_X
#define ISATOOLS_X
#define ISAFLAVR_X
#undef  PRTEESIG_X

#ifdef  LINUX_X
#define IMPNONE_X
#define IDATE_X
#define ETIME_X
#endif

#ifdef  MACOS_X
#define IMPNONE_X
#define IDATE_X
#define ETIME_X
#endif

#ifdef  CDC_X
#define SINGLE_X
#undef  DOUBLE_X
#define LEVEL2_X
#define CDCPACK_X
#undef  STDIO_X
#undef  RANFCALL_X
#undef  RANLUX_X
#endif

#ifdef  SUN_X
#define IMPNONE_X
#define IDATE_X
#define ETIME_X
#endif 

#ifdef  VAX_X
#define IMPNONE_X
#endif

#ifdef  RANFTN_X
#undef  RANLUX_X
#define RANFCALL_X
#endif

#ifdef  CERN_X
#undef  NOCERN_X
#endif

#ifdef  NOISATOOLS_X
#undef ISATOOLS_X
#endif
