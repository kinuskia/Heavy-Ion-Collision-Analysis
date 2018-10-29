#ifndef LAST_CHARACTER
#define LAST_CHARACTER

// get the k-th last character of a string
char last_character(std::string str, int k = 0)
{
	const int len = str.length();
	assert(len >=k);
	return str[len-1-k];
}

#endif