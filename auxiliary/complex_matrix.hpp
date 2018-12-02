#ifndef COMPLEX_MATRIX_HPP
#define COMPLEX_MATRIX_HPP

// Complex-matrix headers
#include <vector>


template <class REAL>
class complex_matrix
{
	typedef std::size_t size_type;
	typedef REAL number_type;
private:
	std::vector<number_type> real_data_; // real part of matrix data
	std::vector<number_type> imag_data_; // imaginary part of matrix data
	size_type m_rows_; // number of matrix rows
	size_type m_cols_; // number of matrix columns

public:
	// constructor for empty matrix
	complex_matrix()
	: real_data_(0, 0)
	, imag_data_(0, 0)
	, m_rows_(0)
	, m_cols_(0)
	{}

	// constructor 
	complex_matrix(const size_type rows, const size_type cols)
	: real_data_(rows*cols, 0)
	, imag_data_(rows*cols, 0)
	, m_rows_(rows)
	, m_cols_(cols)
	{}

	// assignment operator
	complex_matrix & operator= (const complex_matrix & A)
	{
		real_data_ = A.real_data_;
		imag_data_ = A.imag_data_;
		m_rows_ = A.m_rows_;
		m_cols_ = A.m_cols_;
		return *this;
	}

	// Assignment operator for scalar
	complex_matrix & operator= (number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				this->set_entry(i, j, s, 0);
			}
		}
		return *this;
	}

	// getter variables
	size_type rowsize() const
	{
		return m_rows_;
	}

	size_type colsize() const
	{
		return m_cols_;
	}

	// return proper index for entry (i, j)
	size_type index (const size_type i, const size_type j) const
	{
		return i*m_cols_ + j;
	}

	bool index_in_range(const size_type i, const size_type j) const
	{
		bool in_range;
		if (i < m_rows_ && j < m_cols_)
		{
			in_range = true;
		}
		else 
		{
			in_range = false;
		}
		return in_range;
	}

	// Set entry
	void set_entry (const size_type row, const size_type col, const number_type real, const number_type imag)
	{
		assert( index_in_range(row, col) );
		real_data_[index(row, col)] = real;
		imag_data_[index(row, col)] = imag;
	}

	void increase_entry (const size_type row, const size_type col, const number_type change_real, const number_type change_imag)
	{
		assert( index_in_range(row, col) );
		real_data_[index(row, col)] += change_real;
		imag_data_[index(row, col)] += change_imag;		
	}

	// Get entry
	number_type get_real(const size_type row, const size_type col) const
	{
		assert( index_in_range(row, col) );
		return real_data_[index(row, col)];
	}

	number_type get_imag(const size_type row, const size_type col) const
	{
		assert( index_in_range(row, col) );
		return imag_data_[index(row, col)];
	}


	// Addition assignment
	complex_matrix & operator+= (const complex_matrix & B)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				this->increase_entry(i, j, B.get_real(i, j), B.get_imag(i, j));
			}
		}
		return *this;
	}

	// sum of two matrices
	complex_matrix operator+ (const complex_matrix & x) const
	{
		assert(this->colsize() == x.colsize());
		assert(this->rowsize() == x.rowsize());
		complex_matrix y(rowsize(), colsize());
		y = *this;
		y += x;

		return y;
	}

	// Multiplication assignment
	complex_matrix & operator*= (const number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				this->set_entry(i, j, s*(this->get_real(i, j)), s*(this->get_imag(i, j)));
			}
		}
		return *this;
	}

	// Division assignment
	complex_matrix & operator/= (const number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				this->set_entry(i, j, (this->get_real(i, j))/s, (this->get_imag(i, j))/s);
			}
		}
		return *this;
	}


	// Print entries to console
	void print () const
	{
		std::cout << "\n";
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				std::cout << "(" << get_real(i, j) << ", " << get_imag(i, j) << ")" << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
};


#endif