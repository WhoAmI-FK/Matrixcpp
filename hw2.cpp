#include <iostream>
#include <iomanip>
struct matrixException : std::exception {
	const char* what() const noexcept { return "Matrix rows or columns are defined in wrong way."; }
};
struct customException : std::exception {
	const char* what() const noexcept { return "You cannot assign same matrix..."; }
};
struct additionException : std::exception {
	const char* what() const noexcept { return "Matrices should have the same size.."; }
};
struct deterException : std::exception {
	const char* what() const noexcept { return "Only square matrices have determinant..."; }
};

template<class T>
class Matrix {
private:
	size_t _n;
	size_t _m;
	T** _matrix;
	void defineMatrix() {
		if (_n <= 0 || _m <= 0) throw matrixException();
		_matrix = new T * [_n];
		for (int i = 0; i < _n; ++i)
			_matrix[i] = new T[_m];
	}
public:
	Matrix(size_t n, size_t m) : _n(n),_m(m) {
		defineMatrix();
	}
	Matrix(size_t n = 1) : _n(n),_m(n) {
		defineMatrix();
	}
	~Matrix() {
		for (size_t i = 0; i < _n; ++i)
			delete[] _matrix[i];
		delete[] _matrix;
		_matrix = nullptr;
		_n = _m = 0;
	}
	Matrix(const Matrix& mtr) {
		this->_m = mtr._m;
		this->_n = mtr._n;
		this->defineMatrix();
		for (size_t i = 0; i < this->_n; ++i) {
			for (size_t j = 0; j < this->_m; ++j)
				this->_matrix[i][j] = mtr._matrix[i][j];
		}
	}
	Matrix& operator=(const Matrix& mtr) {
		if (this!=&mtr) {
			if (this->_matrix != nullptr) {
					for (size_t i = 0; i < _n; ++i)
						delete[] this->_matrix[i];
				delete[] this->_matrix;
			}
			this->_m = mtr._m;
			this->_n = mtr._n;
			this->defineMatrix();
			for (size_t i = 0; i < this->_n; ++i) {
				for (size_t j = 0; j < this->_m; ++j)
					this->_matrix[i][j] = mtr._matrix[i][j];
			}
		}
		
			return *this;
	}
	 T** getMatrix() const {
		return this->_matrix;
	}
	 int getNumOfRows() const {
		return this->_n;
	}
	 int getNumOfColumns() const {
		return this->_m;
	}
	T operator()(int r, int c) const {
		return this->_matrix[r][c]; // to get cell
	}
	T& operator()(int r, int c){
		return this->_matrix[r][c];
	}

	Matrix& operator+(const Matrix& mtr) {
		if (this->_m == mtr._m && this->_n == mtr._n) {
			for (size_t i = 0; i < this->_n; ++i)
				for (size_t j = 0; j < this->_m; ++j)
					this->_matrix[i][j] += mtr._matrix[i][j];
			return *this;
		}
		else throw additionException();
	}
	Matrix& operator*(T mult) {
		for (size_t i = 0; i < this->_n; ++i)
			for (size_t j = 0; j < this->_m; ++j)
				this->_matrix[i][j] *= mult;
		return *this;
	}
	Matrix& operator*(const Matrix& mtr) {
		std::string message = "Number of columns in the first matrix must be equal to the number of rows in the second matrix";
		if (this->_m == mtr._n) {
			Matrix new_matrix(this->_n, mtr._m);
			for (size_t i = 0; i < new_matrix._n; ++i) {
				for (size_t j = 0; j < new_matrix._m; ++j) {
					int value = 0;
					for (size_t z = 0; z < mtr._n; z++) {
						value += (this->_matrix[i][z] * mtr._matrix[z][j]);
					}
					new_matrix._matrix[i][j] = value;
				}
			}
			*this = new_matrix;
		}
		else std::invalid_argument(message.c_str());
		return *this;
	}
	void displayMatrix() {
		std::cout << "\n\n";
		for (size_t i = 0; i < this->_n; ++i) {
			std::cout << "|";
			for (size_t j = 0; j < this->_m; ++j) {
				std::cout << std::setw(3) << this->_matrix[i][j];
			}
			std::cout << " |\n";
		}
	}
	void fillMatrix() {
		for (size_t i = 0; i < this->_n; ++i)
			for (size_t j = 0; j < this->_m; ++j) {
				std::cout << "Row " << i + 1 << " Cell " << j + 1 << ": ";
				std::cin >> this->_matrix[i][j];
			}
	}
	T getDeterminant() {
		if ((this->_n / this->_m) == 1) {
			T _det = 0;
			T mult;
			if (this->_m == 2 && this->_n == 2) {
				return (this->_matrix[0][0] * this->_matrix[1][1] - this->_matrix[0][1] * this->_matrix[1][0]);
			}
			else {
				for (size_t i = 0; i < this->_m; ++i) {
					mult = i%2 == 0 ? this->_matrix[0][i] : (-1)*this->_matrix[0][i];
					Matrix tmp(this->_n - 1, this->_m - 1);
					for (size_t j = 1,n=0; j < this->_n; ++j) {
						for (size_t k = 0,m=0; k < this->_m; ++k) {
							if (k != i) {
								tmp._matrix[n][m++] = this->_matrix[j][k];
							}
						}
						n++;
					}
					_det += mult * (tmp.getDeterminant());
				}
				return _det;
			}
		}
		else {
			throw deterException();
		}
	}
	Matrix getInverseMatrix() {
		T _deter = this->getDeterminant();
		if (this->getDeterminant() != 0) {
			Matrix<T> inverse(this->_n, this->_m);
			T mult;
			for (size_t i = 0; i < this->_n; ++i) {
				for (size_t j = 0; j < this->_m; ++j) {
					Matrix<T> tmp(this->_n - 1, this->_m - 1);
					int cc, rr;
					rr = cc = 0;
					mult = (i + j + 2) % 2 == 0 ? 1 : -1;
					for (size_t r = 0; r < this->_n; ++r) {
						for (size_t c = 0; c < this->_m; ++c) {
							if ((r != i) && (c != j)) {
								tmp._matrix[rr][cc++] = this->_matrix[r][c];

								if (cc == (this->_m - 1)) {
									cc = 0;
									rr++;
								}
							}
						}
					}
					std::cout << "\n\n";
					tmp.displayMatrix();
					std::cout << "\n\n";
					system("pause");
					inverse._matrix[i][j] = (mult * tmp.getDeterminant()) / this->getDeterminant();
				}
			}

			inverse.transpose();

			return inverse;
		}
		else {
			throw std::invalid_argument("Inverse of this matrix does not exist...");
		}
	}
	void transpose() {
		Matrix<T>* tmp;
		tmp= new Matrix<T>(this->_m, this->_n);
		for (size_t i = 0, m = 0; i < this->_n; ++i, ++m) {
			for (size_t j = 0, n = 0; j < this->_m; ++j, ++n) {
				tmp->_matrix[n][m] = this->_matrix[i][j];
			}
		}
		*this = *tmp;
	}
};

int main()
{
	Matrix<double> mtr(3,3);
	try {
		// here are the operations...
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}
	return 0;
}