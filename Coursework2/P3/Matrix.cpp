#include <cmath>
#include <cassert>
#include <iomanip>
#include "Matrix.hpp"
#include "Vector.hpp"


// Overwritten copy constructor
// Allocate memory for new matrix, and copy
// entries into this matrix
Matrix::Matrix(const Matrix& otherMatrix)
{
   mNumRows = otherMatrix.mNumRows;
   mNumCols = otherMatrix.mNumCols;
   mData = new double* [mNumRows];
   for (int i=0; i<mNumRows; i++)
   {
      mData[i] = new double [mNumCols];
   }
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = otherMatrix.mData[i][j];
      }
   }
}

// Constructor for vector of a given length
// Allocates memory, and initialises entries
// to zero
Matrix::Matrix(int numRows, int numCols)
{
   assert(numRows > 0);
   assert(numCols > 0);
   mNumRows = numRows;
   mNumCols = numCols;
   mData = new double* [mNumRows];
   for (int i=0; i<mNumRows; i++)
   {
      mData[i] = new double [mNumCols];
   }
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = 0.0;
      }
   }
}

// Overwritten destructor to correctly free memory
Matrix::~Matrix()
{
   for (int i=0; i<mNumRows; i++)
   {
      delete[] mData[i];
   }
   delete[] mData;
}

// Method to get number of rows of matrix
int Matrix::GetNumberOfRows() const
{
   return mNumRows;
}

// Method to get number of columns of matrix
int Matrix::GetNumberOfColumns() const
{
   return mNumCols;
}

// Overloading the round brackets
// Note that this uses `one-based' indexing,
// and a check on the validity of the index
double& Matrix::operator()(int i, int j)
{
   assert(i > 0);
   assert(i < mNumRows+1);
   assert(j > 0);
   assert(j < mNumCols+1);
   return mData[i-1][j-1];
}

// Overloading the assignment operator
Matrix& Matrix::operator=(const Matrix& otherMatrix)
{
   assert(mNumRows = otherMatrix.mNumRows);
   assert(mNumCols = otherMatrix.mNumCols);

   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mData[i][j] = otherMatrix.mData[i][j];
      }
   }
   return *this;
}

// Overloading the unary + operator
Matrix Matrix::operator+() const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j];
      }
   }
   return mat;
}

// Overloading the unary - operator
Matrix Matrix::operator-() const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = -mData[i][j];
      }
   }
   return mat;
}

// Overloading the binary + operator
Matrix Matrix::operator+(const Matrix& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j] + m1.mData[i][j];
      }
   }
   return mat;
}

// Overloading the binary - operator
Matrix Matrix::operator-(const Matrix& m1) const
{
   assert(mNumRows == m1.mNumRows);
   assert(mNumCols == m1.mNumCols);
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = mData[i][j] - m1.mData[i][j];
      }
   }
   return mat;
}

// Overloading scalar multiplication
Matrix Matrix::operator*(double a) const
{
   Matrix mat(mNumRows, mNumCols);
   for (int i=0; i<mNumRows; i++)
   {
      for (int j=0; j<mNumCols; j++)
      {
         mat(i+1,j+1) = a*mData[i][j];
      }
   }
   return mat;
}

void Matrix::overwriteRow(int row_index,
                          const Vector& content){
    assert(row_index > 0);
    assert(row_index < mNumRows+1);
    assert(content.GetSize() == mNumCols);
    for (int i=1; i<=mNumCols; i++){
        mData[row_index][i] = content.Read(i-1);
    }
}

void Matrix::overwriteColumn(int col_index,
                             const Vector& content){
    assert(col_index > 0);
    assert(col_index < mNumCols+1);
    assert(content.GetSize() == mNumRows);
    for (int i=1; i<=mNumRows; i++){
        mData[i][col_index] = content.Read(i-1);
    }
}

// Function Kronecker Product
/*      Added by Giordi Azonos      */
Matrix Matrix::KroneckerProduct( const Matrix& m1) const{
    int m1_n_rows = m1.GetNumberOfRows();
    int m1_n_cols = m1.GetNumberOfColumns();
    Matrix result(mNumRows*m1_n_rows, mNumCols*m1_n_cols);

    /*"cur_row" and "cur_col are the indices of the
    result matrix.*/
    int cur_row = 0;
    int cur_col = 0;

    // Iterate through rows of this matrix
    for(int r_a=0; r_a<mNumRows; r_a++){
        // Iterate through rows of the given matrix
        for (int r_b=0; r_b<m1_n_rows; r_b++){

            /* We will first fill the entire columns of the
             * first row, of the result matrix, and then change row,
             * as if we were writing.
            */
            // Iterate through columns of this matrix.
            for (int c_a=0; c_a<mNumCols; c_a ++){
                // Iterate through columns of given matrix.
                for (int c_b=0; c_b<m1_n_cols; c_b++){
                    result(cur_row+1, cur_col+1) = mData[r_a][c_a] * m1.mData[r_b][c_b];
                    cur_col += 1;
                }
            }
            // Restart column of result matrix,
            cur_col = 0;
            // Change row of the result matrix
            cur_row += 1;
            // Equivalent to a change of line when writing.
        }
    }
    return result;
}

// Overloading matrix multiplied by another matrix
/*      Added by Giordi Azonos        */
Matrix Matrix::operator*(const Matrix& m1) const{
    assert(mNumCols == m1.GetNumberOfRows());
    Matrix result(mNumRows, mNumCols);
    for (int r=0; r < mNumRows; r++){
        for (int c=0; c<mNumCols; c++){
            double new_entry = 0;
            for (int entry = 0; entry < mNumRows; entry++){
                new_entry += mData[r][entry] * m1.mData[entry][c];
            }
            result(r+1, c+1) = new_entry;
        }
    }
    return result;
}

// Overloading matrix multiplied by a vector
Vector operator*(const Matrix& m, const Vector& v)
{
   int original_vector_size = v.GetSize();
   assert(m.GetNumberOfColumns() == original_vector_size);
   int new_vector_length = m.GetNumberOfRows();
   Vector new_vector(new_vector_length);

   for (int i=0; i<new_vector_length; i++)
   {
      for (int j=0; j<original_vector_size; j++)
      {
         new_vector[i] += m.mData[i][j]*v.Read(j);
      }
   }

   return new_vector;
}

// Overloading vector multiplied by a matrix
Vector operator*(const Vector& v, const Matrix& m)
{
   int original_vector_size = v.GetSize();
   assert(m.GetNumberOfRows() == original_vector_size);
   int new_vector_length = m.GetNumberOfColumns();
   Vector new_vector(new_vector_length);

   for (int i=0; i<new_vector_length; i++)
   {
      for (int j=0; j<original_vector_size; j++)
      {
         new_vector[i] += v.Read(j)*m.mData[j][i];
      }
   }

   return new_vector;
}

// Calculate determinant of square matrix recursively
double Matrix::CalculateDeterminant() const
{
   assert(mNumRows == mNumCols);
   double determinant = 0.0;

   if (mNumRows == 1)
   {
      determinant = mData[0][0];
   }
   else
   {
      // More than one entry of matrix
      for (int i_outer=0; i_outer<mNumRows; i_outer++)
      {
         Matrix sub_matrix(mNumRows-1, 
                             mNumRows-1);
         for (int i=0; i<mNumRows-1; i++)
         {
            for (int j=0; j<i_outer; j++)
            {
               sub_matrix(i+1,j+1) = mData[i+1][j];
            }
            for (int j=i_outer; j<mNumRows-1; j++)
            {
               sub_matrix(i+1,j+1) = mData[i+1][j+1];
            }
         }
         double sub_matrix_determinant = 
                  sub_matrix.CalculateDeterminant();

         determinant += pow(-1.0, i_outer)*
                  mData[0][i_outer]*sub_matrix_determinant;
      }
   }
   return determinant;
}

// Added by Daniele Avitabile
std::ostream& operator<<(std::ostream& output,
                        const Matrix& m)
{

  // Print formatted output
  for (int i=0; i< m.GetNumberOfRows(); i++)
  {
    for (int j=0; j< m.GetNumberOfColumns(); j++)
    {
      output << std::setw(14)
             << std::setprecision(5)
	     << std::scientific 
	     << m.mData[i][j];
    }
    output << std::endl;
  }
  output << std::endl;

  return output;

}
