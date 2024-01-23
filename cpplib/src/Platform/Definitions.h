#pragma once

/*  Standart definitions
	__linux__       Defined on Linux
	__sun           Defined on Solaris
	__FreeBSD__     Defined on FreeBSD
	__NetBSD__      Defined on NetBSD
	__OpenBSD__     Defined on OpenBSD
	__APPLE__       Defined on Mac OS X
	__hpux          Defined on HP-UX
	__osf__         Defined on Tru64 UNIX (formerly DEC OSF1)
	__sgi           Defined on Irix
	_AIX            Defined on AIX
	_WIN32          Defined on Windows
*/

// API implementation for Windows and Linux (not UNIX)
#if defined(_WIN32)
#define EXPORT __declspec(dllexport)
#define IMPORT __declspec(dllimport)
#elif defined(__linux__)
#define EXPORT __attribute__((visibility("default")))
#define IMPORT
#else
#define EXPORT
#define IMPORT
#endif // OS specific

#if defined(cpplib_EXPORTS)
#define API extern "C" EXPORT
#else
#define API extern "C"
#endif

#if !defined(NDEBUG) || defined(CPPLIBLIB)
#undef API 
#define API extern "C"
#endif

API int* SearchMain(const char* search, const char** data, const int data_s, const int np, const bool exact);
API bool CompareGraph(const char* search1, const char* search2, const bool exact);
API const char* FindMoleculesInCell(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s);
API const char* FindMoleculesWithoutCell(const int* types, const float* xyz, const int types_s);

API const char* FindDistanceWC(const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value);
API const char* FindDistanceIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s, const int type1, const int type2, const float min_value, const float max_value);
API const char* FindAngleWC(const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123);
API const char* FindAngleIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const float min12, const float max12, const float min23, const float max23, const float min123, const float max123);
API const char* FindTorsionWC(const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
	const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234);
API const char* FindTorsionIC(const float* unit_cell, const char** symm, const int symm_s, const int* types, const float* xyz, const int types_s,
	const int type1, const int type2, const int type3, const int type4, const float min12, const float max12, const float min23, const float max23, const float min34, const float max34,
	const float min123, const float max123, const float min234, const float max234, const float min1234, const float max1234);



