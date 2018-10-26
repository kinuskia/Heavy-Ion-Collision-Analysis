#ifndef POLAR_HPP
#define POLAR_HPP


// class to save values that depend on polar coordinates
template<typename REAL>
class Polar
{
public:
	typedef REAL number_type;
	typedef std::size_t size_type;
	Polar(number_type r, number_type phi, number_type value)
	: phi_(phi)
	, r_(r)
	, value_(value)
	{}

	number_type r() const
	{
		return r_;
	}

	number_type phi() const
	{
		return phi_;
	}

	number_type value() const
	{
		return value_;
	}

private:
	number_type r_;
	number_type phi_;
	number_type value_;
};

#endif