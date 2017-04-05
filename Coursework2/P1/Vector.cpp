#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "Vector.hpp"

// Overridden copy constructor
// Allocates memory for new vector, and copies
// entries of other vector into it
Vector::Vector(const Vector& otherVector)
{
   mSize = otherVector.GetSize();
   mData = new double [mSize];
   for (int i=0; i<mSize; i++)
   {
      mData[i] = otherVector.mData[i];
   }
}

// Constructor for vector of a given size
// Allocates memory, and initialises entries
// to zero
Vector::Vector(int size)
{
   if(size < 0)
       throw std::runtime_error("Can't initialize Vector of negative size.");
   mSize = size;
   mData = new double [mSize];
   for (int i=0; i<mSize; i++)
   {
      mData[i] = 0.0;
   }
}

// Overridden destructor to correctly free memory
Vector::~Vector()
{
   delete[] mData;
}

// Method to get the size of a vector
int Vector::GetSize() const
{
   return mSize;
}

// Overloading square brackets
// Note that this uses `zero-based' indexing,
// and a check on the validity of the index
double& Vector::operator[](int i)
{
   checkIndex(i, 0, mSize-1, "operator []");
   return mData[i];
}

// Read-only variant of []
// Note that this uses `zero-based' indexing,
// and a check on the validity of the index
double Vector::Read(int i) const
{
   checkIndex(i, 0, mSize-1, "Read()");
   return mData[i];
}

// Overloading round brackets
// Note that this uses `one-based' indexing,
// and a check on the validity of the index
double& Vector::operator()(int i)
{
   checkIndex(i, 0, mSize-1, "operator()");
   return mData[i-1];
}

// Overloading the assignment operator
Vector& Vector::operator=(const Vector& otherVector)
{
   if(mSize != otherVector.mSize)
       throw std::runtime_error("Error on Vector::operator =: Different sizes.");
   for (int i=0; i<mSize; i++)
   {
      mData[i] = otherVector.mData[i];
   }
   return *this;
}

// Overloading the unary + operator
Vector Vector::operator+() const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i];
   }
   return v;
}

// Overloading the unary - operator
Vector Vector::operator-() const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = -mData[i];
   }
   return v;
}

// Overloading the binary + operator
Vector Vector::operator+(const Vector& v1) const
{
    if(mSize != v1.mSize)
        throw std::runtime_error("Error on Vector::operator +: Different sizes.");
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i] + v1.mData[i];
   }
   return v;
}

// Overloading the binary - operator
Vector Vector::operator-(const Vector& v1) const
{
    if(mSize != v1.mSize)
        throw std::runtime_error("Error on Vector::operator -: Different sizes.");
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = mData[i] - v1.mData[i];
   }
   return v;
}

// Overloading scalar multiplication
Vector Vector::operator*(double a) const
{
   Vector v(mSize);
   for (int i=0; i<mSize; i++)
   {
      v[i] = a*mData[i];
   }
   return v;
}

// Method to calculate norm (with default value p=2)
// corresponding to the Euclidean norm
double Vector::Norm(int p) const
{
   double norm_val, sum = 0.0;
   for (int i=0; i<mSize; i++)
   {
      sum += pow(fabs(mData[i]), p);
   }
   norm_val = pow(sum, 1.0/((double)(p)));
   return norm_val;
}

// Method to calculate the infinity norm
double Vector::InfinityNorm() const
{
   double norm_val = 0;
   for (int i=0; i<mSize; i++)
   {
     double abs_val = fabs(mData[i]);
     if ( norm_val < abs_val )
     {
       norm_val = abs_val;
     }
   }
   return norm_val;
}

// Method to calculate scalar product with another vector
double Vector::dot(const Vector& v) const
{
   double scalar_product = 0.0;
   if(mSize != v.mSize)
       throw std::runtime_error("Error on Vector::dot: Different sizes.");
   for (int i=0; i<mSize; i++)
   {
      scalar_product += mData[i]*v.Read(i);
   }
   return scalar_product;
}

// MATLAB style friend to get the size of a vector
int length(const Vector& v)
{
   return v.mSize;
}
//Code from Chapter10.tex line 60 save as Vector.cpp

// Added by Daniele Avitabile
std::ostream& operator<<(std::ostream& output,
                        const Vector& v)
{
   for (int i=0; i<v.mSize; i++)
   {
      output << std::setw(14)
             << std::setprecision(5)
	     << std::scientific 
	     << v.Read(i)
	     << std::endl;
   }
   output << std::endl;

   return output;
}

// Added by Giordi Azonos
void Vector::checkIndex(int index, int min, int max, std::string prefix) const {
    if (index < min || index > max) {
        std::ostringstream out;
        out << "Vector::" << prefix << ": index " << index
            << " is outside of valid range ";

        out << "[";
        if (min < max) {
            out << min << ".." << max;
        } else if (min == max) {
            out << min;
        } // else min > max, no range, empty vector
        out << "]";

        throw std::runtime_error(out.str());
    }
}

// Added by Giordi Azonos
std::string Vector::toString()const{

    std::string str = "[";
    for (int i=0; i<mSize-1; i++){
        str += std::to_string(mData[i]) + ", ";
    }
    str += std::to_string(mData[mSize-1]) + "]";

    return str;
}
