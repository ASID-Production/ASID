#pragma once

#ifdef DEBUG_MESS_ON
#include <iostream>
inline void deb_write(const char* mes) {
	std::cout << mes << std::endl;
}
template <class T>
inline void deb_write(const char* mes, const T& val) {
	std::cout << mes  << val << std::endl;
}
#else
inline void deb_write(const char*) { /*empty debug function*/ }
template <class T>
inline void deb_write(const char*, const T&) { /*empty debug function*/ }

#endif // DEBUG_MESS_ON