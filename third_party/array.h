// 2024-2-26-09点18分

#ifndef __ARRAY_H_
#define __ARRAY_H_

// #define NDEBUG  //Uncomment this macro definition when boundary check becomes unnecessary.

#ifndef NDEBUG
#include <cassert>
#endif
#include <iostream>
#include <cmath>
#include <fstream>

/************************************************************************/
/*                            Array 1D                                  */
/************************************************************************/
template <typename T>
class Array1D
{
public:
	Array1D();
	Array1D(unsigned int row);
	Array1D(unsigned int row, const T &val);

	Array1D(const Array1D<T> &A);
	inline Array1D<T> &operator=(const Array1D<T> &A);

	~Array1D();

public:
	inline void Resize(unsigned int row);
	inline void Resize(unsigned int row, const T &val);

	inline unsigned int dim() const;
	inline void setZero();
	inline void abs();

	inline void output(std::string filename = "output.txt");

	inline double sum();
	inline double max();
	inline double min();

	inline double norm(std::string type = "l2");

	inline T &operator()(unsigned int row);
	inline const T &operator()(unsigned int row) const;

	inline T &operator[](unsigned int row);
	inline const T &operator[](unsigned int row) const;

private:
	T *m_A;
	unsigned int m_Row;
};
///////////////////////////////////////////////////////////
//////////////
template <typename T>
Array1D<T>::Array1D(void)
	: m_A(nullptr),
	  m_Row(0)
{
}
template <typename T>
Array1D<T>::Array1D(unsigned int row)
	: m_Row(row)
{
	m_A = new T[m_Row];
}
template <typename T>
Array1D<T>::Array1D(unsigned int row, const T &val)
	: m_Row(row)
{
	m_A = new T[m_Row];

	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = val;
	}
}
template <typename T>
Array1D<T>::Array1D(const Array1D &A)
	: m_Row(A.dim())
{
	m_A = new T[m_Row];
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = A[i];
	}
}
template <typename T>
Array1D<T>::~Array1D(void)
{
	if (m_A != nullptr)
	{
		delete[] m_A;
		m_A = nullptr;
	}
}
template <typename T>
void Array1D<T>::Resize(unsigned int row)
{
	if (row == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A;
	}

	m_Row = row;
	m_A = new T[m_Row];
}
template <typename T>
void Array1D<T>::Resize(unsigned int row, const T &val)
{
	if (row == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A;
	}

	m_Row = row;
	m_A = new T[m_Row];
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = val;
	}
}
template <typename T>
unsigned int Array1D<T>::dim() const
{
	return m_Row;
}
template <typename T>
void Array1D<T>::setZero()
{
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = 0;
	}
}
template <typename T>
void Array1D<T>::abs()
{
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = fabs(m_A[i]);
	}
}
template <typename T>
void Array1D<T>::output(std::string filename)
{
	std::ofstream outputFile(filename);
	if (outputFile.is_open())
	{
		for (unsigned int i = 0; i != m_Row; ++i)
		{
			outputFile << m_A[i] << " ";
		}
		outputFile.close();
	}
	else
	{
	}
}
template <typename T>
double Array1D<T>::sum()
{
	double s = 0;
	for (unsigned int i = 0; i != m_Row; ++i)
		s += m_A[i];
	return s;
}
template <typename T>
double Array1D<T>::max()
{
	double maxVal = m_A[0];
	for (unsigned int i = 1; i != m_Row; ++i)
	{
		if (maxVal < m_A[i])
			maxVal = m_A[i];
	}
	return maxVal;
}
template <typename T>
double Array1D<T>::min()
{
	double minVal = m_A[0];
	for (unsigned int i = 1; i != m_Row; ++i)
	{
		if (minVal > m_A[i])
			minVal = m_A[i];
	}
	return minVal;
}
template <typename T>
double Array1D<T>::norm(std::string type)
{
	Array1D<T> temp(m_Row);
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		temp[i] = m_A[i];
	}

	if (type == "l1")
	{
		temp.abs();
		T val = 0;
		for (unsigned int i = 0; i != m_Row; ++i)
		{
			val += temp[i];
		}
		return val;
	}
	else if (type == "l2")
	{
		T s = 0;
		for (unsigned int i = 0; i != m_Row; ++i)
		{
			s += temp[i] * temp[i];
		}
		return sqrt(s);
	}
	else if (type == "inf")
	{
		temp.abs();
		double maxVal = 0;
		for (unsigned int i = 0; i != m_Row; ++i)
		{
			// std::cout << m_A[i] << std::endl;
			if (temp[i] > maxVal)
				maxVal = temp[i];
		}
		return maxVal;
	}
	else
	{
		// If the type is not recognized, throw an exception
		throw std::invalid_argument("Invalid norm type: " + type);
	}
}
template <typename T>
Array1D<T> &Array1D<T>::operator=(const Array1D &A)
{
	if (this == &A) // 自赋值检查
	{
		return *this;
	}

	if (m_A != nullptr)
	{
		delete[] m_A;
	}

	m_Row = A.dim();
	m_A = new T[m_Row];
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		m_A[i] = A[i];
	}

	return *this;
}
template <typename T>
T &Array1D<T>::operator()(unsigned int row)
{
#ifndef NDEBUG
	assert(row < m_Row);
#endif

	return m_A[row];
}
template <typename T>
const T &Array1D<T>::operator()(unsigned int row) const
{
#ifndef NDEBUG
	assert(row < m_Row);
#endif

	return m_A[row];
}
template <typename T>
T &Array1D<T>::operator[](unsigned int row)
{
#ifndef NDEBUG
	assert(row < m_Row);
#endif

	return m_A[row];
}
template <typename T>
const T &Array1D<T>::operator[](unsigned int row) const
{
#ifndef NDEBUG
	assert(row < m_Row);
#endif
	return m_A[row];
}
template <typename T>
std::ostream &operator<<(std::ostream &os, const Array1D<T> &A)
{
	for (unsigned int i = 0; i < A.dim(); ++i)
	{
		os << A[i] << " ";
	}
	return os;
}
/************************************************************************/
/*                            Array  2D                                 */
/************************************************************************/
template <typename T>
class Array2D
{
public:
	Array2D();
	Array2D(unsigned int row, unsigned int col);
	Array2D(unsigned int row, unsigned int col, const T &val);

	Array2D(const Array2D<T> &A);
	inline Array2D<T> &operator=(const Array2D<T> &A);

	~Array2D();

public:
	void Resize(unsigned int row, unsigned int col);
	void Resize(unsigned int row, unsigned int col, const T &val);

	inline unsigned int row() const;
	inline unsigned int col() const;
	inline void setZero();
	inline void transpose(Array2D<T> &At);

	inline T &operator()(unsigned int row, unsigned int col);
	inline const T &operator()(unsigned int row, unsigned int col) const;

	inline T *operator[](unsigned int row);
	inline const T *operator[](unsigned int row) const;

private:
	T **m_A;
	unsigned int m_Row;
	unsigned int m_Col;
};
///////////////////////////////////////////////////////////
//////////////
template <typename T>
Array2D<T>::Array2D(void)
	: m_A(nullptr),
	  m_Row(0),
	  m_Col(0)
{
}
template <typename T>
Array2D<T>::Array2D(unsigned int row, unsigned int col)
	: m_Row(row),
	  m_Col(col)
{
	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col]; // 使用连续的内存存储二维数组

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}
}
template <typename T>
Array2D<T>::Array2D(unsigned int row, unsigned int col, const T &val)
	: m_Row(row),
	  m_Col(col)
{
	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col]; // 使用连续的内存存储二维数组

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}

	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			m_A[i][j] = val;
		}
	}
}
template <typename T>
Array2D<T>::Array2D(const Array2D<T> &A)
	: m_Row(A.row()),
	  m_Col(A.col())
{
	// 分配新内存
	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col];

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}

	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			m_A[i][j] = A[i][j];
		}
	}
}
template <typename T>
Array2D<T> &Array2D<T>::operator=(const Array2D<T> &A)
{
	if (this == &A) // 自赋值检查
	{
		return *this;
	}

	m_Row = A.row();
	m_Col = A.col();

	// 释放原有内存
	if (m_A != nullptr)
	{
		delete[] m_A[0];
		delete[] m_A;
	}

	// 分配新内存
	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col];

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}

	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			m_A[i][j] = A[i][j];
		}
	}

	return *this;
}
template <typename T>
Array2D<T>::~Array2D(void)
{
	if (m_A != nullptr)
	{
		delete[] m_A[0];
		delete[] m_A;

		m_A = nullptr;
	}
}
template <typename T>
void Array2D<T>::Resize(unsigned int row, unsigned int col)
{
	if (row == 0 || col == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A[0];
		delete[] m_A;
	}

	m_Row = row;
	m_Col = col;

	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col];

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}
}
template <typename T>
void Array2D<T>::Resize(unsigned int row, unsigned int col, const T &val)
{
	if (row == 0 || col == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A[0];
		delete[] m_A;
	}

	m_Row = row;
	m_Col = col;

	m_A = new T *[m_Row];
	m_A[0] = new T[m_Row * m_Col];

	for (unsigned int i = 1; i != m_Row; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Col;
	}

	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			m_A[i][j] = val;
		}
	}
}
template <typename T>
unsigned int Array2D<T>::row() const
{
	return m_Row;
}
template <typename T>
unsigned int Array2D<T>::col() const
{
	return m_Col;
}
template <typename T>
void Array2D<T>::setZero()
{
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			m_A[i][j] = 0;
		}
	}
}
template <typename T>
void Array2D<T>::transpose(Array2D<T> &At)
{
	for (unsigned int i = 0; i != m_Row; ++i)
	{
		for (unsigned int j = 0; j != m_Col; ++j)
		{
			At[j][i] = m_A[i][j];
		}
	}
}
template <typename T>
T &Array2D<T>::operator()(unsigned int row, unsigned int col)
{
#ifndef NDEBUG
	assert(row < m_Row && col < m_Col);
#endif

	return m_A[row][col];
}
template <typename T>
const T &Array2D<T>::operator()(unsigned int row, unsigned int col) const
{
#ifndef NDEBUG
	assert(row < m_Row && col < m_Col);
#endif

	return m_A[row][col];
}
template <typename T>
T *Array2D<T>::operator[](unsigned int row)
{
#ifndef NDEBUG
	assert(row < m_Row); // 注意：此时,只能检查第一维是否存在数组越界!
#endif

	return m_A[row];
}
template <typename T>
const T *Array2D<T>::operator[](unsigned int row) const
{
#ifndef NDEBUG
	assert(row < m_Row);
#endif
	return m_A[row];
}

template <typename T>
Array1D<T> operator*(const Array2D<T> &A, const Array1D<T> &B)
{
	unsigned int rowA = A.row();
	unsigned int colA = A.col();
	unsigned int sizeB = B.dim();

#ifndef NDEBUG
	assert(colA == sizeB); // Matrix multiplication requires the number of columns in A to be equal to the size of B
#endif

	Array1D<T> C(rowA);
	C.setZero();
	for (unsigned int i = 0; i < rowA; ++i)
	{
		for (unsigned int j = 0; j < colA; ++j)
		{
			C[i] += A[i][j] * B[j];
		}
	}

	return C;
}

template <typename T>
Array2D<T> operator*(const Array2D<T> &A, const Array2D<T> &B)
{
	unsigned int rowA = A.row();
	unsigned int colA = A.col();
	unsigned int rowB = B.row();
	unsigned int colB = B.col();

	assert(colA == rowB); // 矩阵乘法要求A的列数等于B的行数

	Array2D<T> C(rowA, colB);
	C.setZero();
	for (unsigned int i = 0; i < rowA; ++i)
	{
		for (unsigned int j = 0; j < colB; ++j)
		{
			for (unsigned int k = 0; k < colA; ++k)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return C;
}

/************************************************************************/
/*                              Array  3D                               */
/************************************************************************/
template <typename T>
class Array3D
{
public:
	Array3D();
	Array3D(unsigned int nx, unsigned int ny, unsigned int nz);
	Array3D(unsigned int nx, unsigned int ny, unsigned int nz, const T &val);

	Array3D(const Array3D<T> &A);
	inline Array3D<T> &operator=(const Array3D<T> &A);

	~Array3D();

public:
	void Resize(unsigned int nx, unsigned int ny, unsigned int nz);
	void Resize(unsigned int nx, unsigned int ny, unsigned int nz, const T &val);

	inline unsigned int dim1() const;
	inline unsigned int dim2() const;
	inline unsigned int dim3() const;
	inline void setZero();

	inline T &operator()(unsigned int nx, unsigned int ny, unsigned int nz);
	inline const T &operator()(unsigned int nx, unsigned int ny, unsigned int nz) const;

	inline T **operator[](unsigned int nx);
	inline const T *const *operator[](unsigned int nx) const;

private:
	T ***m_A;
	unsigned int m_Nx;
	unsigned int m_Ny;
	unsigned int m_Nz;
};
///////////////////////////////////////////////////////////
//////////////
template <typename T>
Array3D<T>::Array3D(void)
	: m_A(nullptr),
	  m_Nx(0),
	  m_Ny(0),
	  m_Nz(0)
{
}
template <typename T>
Array3D<T>::Array3D(unsigned int nx, unsigned int ny, unsigned int nz)
	: m_Nx(nx),
	  m_Ny(ny),
	  m_Nz(nz)
{
	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}
}
template <typename T>
Array3D<T>::Array3D(unsigned int nx, unsigned int ny, unsigned int nz, const T &val)
	: m_Nx(nx),
	  m_Ny(ny),
	  m_Nz(nz)
{
	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}

	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 0; j != m_Ny; ++j)
		{
			for (unsigned int k = 0; k != m_Nz; ++k)
			{
				m_A[i][j][k] = val;
			}
		}
	}
}
template <typename T>
Array3D<T>::Array3D(const Array3D<T> &A)
	: m_Nx(A.dim1()),
	  m_Ny(A.dim2()),
	  m_Nz(A.dim3())
{
	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}

	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 0; j != m_Ny; ++j)
		{
			for (unsigned int k = 0; k != m_Nz; ++k)
			{
				m_A[i][j][k] = A[i][j][k];
			}
		}
	}
}
template <typename T>
Array3D<T> &Array3D<T>::operator=(const Array3D<T> &A)
{
	if (this == &A) // 自赋值检查
	{
		return *this;
	}

	m_Nx = A.dim1();
	m_Ny = A.dim2();
	m_Nz = A.dim3();

	if (m_A != nullptr)
	{
		delete[] m_A[0][0];
		delete[] m_A[0];
		delete[] m_A;
	}

	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}

	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 0; j != m_Ny; ++j)
		{
			for (unsigned int k = 0; k != m_Nz; ++k)
			{
				m_A[i][j][k] = A[i][j][k];
			}
		}
	}

	return *this;
}
template <typename T>
Array3D<T>::~Array3D(void)
{
	if (m_A != nullptr)
	{
		delete[] m_A[0][0];
		delete[] m_A[0];
		delete[] m_A;

		m_A = nullptr;
	}
}
template <typename T>
void Array3D<T>::Resize(unsigned int nx, unsigned int ny, unsigned int nz)
{
	if (nx == 0 || ny == 0 || nz == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A[0][0];
		delete[] m_A[0];
		delete[] m_A;
	}

	m_Nx = nx;
	m_Ny = ny;
	m_Nz = nz;

	// 分配内存
	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}
}
template <typename T>
void Array3D<T>::Resize(unsigned int nx, unsigned int ny, unsigned int nz, const T &val)
{
	if (nx == 0 || ny == 0 || nz == 0)
	{
		return;
	}
	if (m_A != nullptr)
	{
		delete[] m_A[0][0];
		delete[] m_A[0];
		delete[] m_A;
	}

	m_Nx = nx;
	m_Ny = ny;
	m_Nz = nz;

	// 分配内存
	m_A = new T **[m_Nx];

	// 分配内存
	m_A[0] = new T *[m_Nx * m_Ny];
	// 为其余的m_A[i]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i] = m_A[i - 1] + m_Ny;
	}

	// 分配内存
	m_A[0][0] = new T[m_Nx * m_Ny * m_Nz];
	// 为其余的m_A[i][j]赋值
	for (unsigned int i = 1; i != m_Nx; ++i)
	{
		m_A[i][0] = m_A[i - 1][0] + m_Ny * m_Nz;
	}
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 1; j != m_Ny; ++j)
		{
			m_A[i][j] = m_A[i][j - 1] + m_Nz;
		}
	}

	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 0; j != m_Ny; ++j)
		{
			for (unsigned int k = 0; k != m_Nz; ++k)
			{
				m_A[i][j][k] = val;
			}
		}
	}
}
template <typename T>
unsigned int Array3D<T>::dim1() const
{
	return m_Nx;
}
template <typename T>
unsigned int Array3D<T>::dim2() const
{
	return m_Ny;
}
template <typename T>
unsigned int Array3D<T>::dim3() const
{
	return m_Nz;
}
template <typename T>
void Array3D<T>::setZero()
{
	for (unsigned int i = 0; i != m_Nx; ++i)
	{
		for (unsigned int j = 0; j != m_Ny; ++j)
		{
			for (unsigned int k = 0; k != m_Nz; ++k)
			{
				m_A[i][j][k] = 0;
			}
		}
	}
}
template <typename T>
T &Array3D<T>::operator()(unsigned int nx, unsigned int ny, unsigned int nz)
{
#ifndef NDEBUG
	assert(nx < m_Nx && ny < m_Ny && nz < m_Nz);
#endif

	return m_A[nx][ny][nz];
}
template <typename T>
const T &Array3D<T>::operator()(unsigned int nx, unsigned int ny, unsigned int nz) const
{
#ifndef NDEBUG
	assert(nx < m_Nx && ny < m_Ny && nz < m_Nz);
#endif

	return m_A[nx][ny][nz];
}
template <typename T>
T **Array3D<T>::operator[](unsigned int nx)
{
#ifndef NDEBUG
	assert(nx < m_Nx);
#endif

	return m_A[nx];
}
template <typename T>
const T *const *Array3D<T>::operator[](unsigned int nx) const
{
#ifndef NDEBUG
	assert(nx < m_Nx);
#endif
	return m_A[nx];
}
#endif
