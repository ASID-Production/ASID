#pragma once

#ifdef DEBUG_MESS_ON
#include <iostream>
inline void deb_write(const char* mes) {
	std::cerr << mes << std::endl;
}
template <class T>
inline void deb_write(const char* mes, const T& val) {
	std::cerr << mes  << val << std::endl;
}
#else
inline void deb_write(const char*) { /*empty debug function*/ }
template <class T>
inline void deb_write(const char*, const T&) { /*empty debug function*/ }

#endif // DEBUG_MESS_ON
#ifndef WIN32
#define _ASSERT(a) (0);
#endif