#ifndef TO_SIZE_T
#define TO_SIZE_T
#include <sstream>

// function which converts string to size_t
template<typename stringtype>
std::size_t to_size_t(stringtype number)
{
	std::stringstream sstream(number);
	std::size_t result;
	sstream >> result;
	return result;
}

#endif