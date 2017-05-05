#include "Helpers.hpp"

Helpers::Helpers(){}

double Helpers::infinityNorm(const Vector& v){
    int n = v.GetSize();
    double norm_val = 0;
    for (int i=0; i<n; i++){
        double abs_val = fabs(v.Read(i));
        if ( norm_val < abs_val ){
            norm_val = abs_val;
        }
    }
    return norm_val;
}

void Helpers::printV(const Vector& v){
    int size = v.GetSize();
    std::cout << "[ ";
    for (int i=0; i<size-1; i++){
        if ( isAlmostZero(v.Read(i)) ) std::cout << 0 << ", ";
        else std::cout << v.Read(i) << ", ";
    }
    std::cout << v.Read(size-1) << " ]" << std::endl;
}

bool Helpers::isAlmostZero(double number, double error_tol){
    if (number < 0){ // I make cases for negative number because abs(number) rounds the number.
        return (number >= -error_tol);
    }else{
        return (number <= error_tol);
    }
}

void Helpers::saveVectorToFile(const Vector& v, std::string fileName){

    int n = v.GetSize();

    // Format output
    std::ofstream output_file;
    output_file.setf(std::ios::scientific,std::ios::floatfield);
    output_file.precision(6);

    // Open file (and perform a check)
    output_file.open(fileName);
    if( !output_file.is_open() ){
        throw std::runtime_error("Error opening file.");
    }

    // Write data
    for (int i=0; i<n; i++){
        output_file << std::setw(15) << v.Read(i) << std::endl;
    }

    // Close file
    output_file.close();

}

Vector Helpers::linspace(double first, double last, int n){
    double h = (last - first)/(n-1);
    Vector mesh(n);
    mesh[0] = first;
    for(int i=1; i<n; i++){
        mesh[i] = mesh[i-1] + h;
    }
    mesh[n-1] = last;
    return mesh;
}

void Helpers::solveTridiagonal(const Vector& diag, const Vector& offd, Vector& x){

    int n = diag.GetSize();

    Vector a(n);
    Vector b(n-1);

    /* LU decomposition */
    a[0] = diag.Read(0);
    for (int k=0; k<n-1; k++) {
        b[k] = offd.Read(k)/a[k];
        a[k+1] = diag.Read(k+1) - b[k]*offd.Read(k);
    }

    /* Inverting L */
    for (int k=1; k<n; k++) {
        x[k] = x[k] - b[k-1]*x[k-1];
    }

    /* Inverting U */
    x[n-1] /= a[n-1];
    for (int k=n-2; k>=0; k--) {
        x[k] = x[k]/a[k] - b[k]*x[k+1];
    }

}

Vector Helpers::tridiagTimesVec(const Vector& ld,
                                const Vector& d,
                                const Vector& ud,
                                const Vector& v){
    int n = d.GetSize();
    Vector result(n);
    result[0] = d.Read(0)*v.Read(0) + ud.Read(0)*v.Read(1);
    for (int i=0; i<ld.GetSize()-1; i++){
        result[i+1] = ld.Read(i)*v.Read(i) + d.Read(i+1)*v.Read(i+1) + ud.Read(i)*v.Read(i+2);
    }
    result[n-1] = ld.Read(n-2)*v.Read(n-2) + d.Read(n-1)*v.Read(n-1);
    return result;
}

void Helpers::openOutputFile(std::string file_name, std::ofstream& output_file){

    output_file.setf(std::ios::scientific,std::ios::floatfield);
    output_file.precision(7);

    output_file.open(file_name);

    if(!output_file.is_open()){
        throw std::runtime_error("Error when opening the file to write the solution.");
    }
}

void Helpers::saveSolution(const double t, const double l_BC,
                           const Vector& us, const double r_BC,
                           std::ofstream& output_file){
    // Save time on first column.
    output_file << std::setw(15) << t;
    // Save the first component, which is the value
    // of the left Boundary Condition.
    output_file << std::setw(15) << l_BC;
    // Save the interior solution.
    for (int i=0; i<us.GetSize(); i++){
        output_file << std::setw(15) << us.Read(i);
    }
    // Save the last component, which is the value
    // of the right Boundary Condition.
    output_file << std::setw(15) << r_BC;
    output_file << "\n"; // End of Line.
}

void Helpers::saveData(const double h, double norm, std::ofstream& output_file){
    output_file << std::setw(15) << h;
    output_file << std::setw(15) << norm << std::endl;
}

void Helpers::saveData(const double t, const Vector& us, std::ofstream& output_file){
    output_file << std::setw(15) << t;
    for (int i=0; i<us.GetSize(); i++){
        output_file << std::setw(15) << us.Read(i);
    }
    output_file << "\n"; // End of Line.
}

void Helpers::closeOutputFile(std::ofstream& output_file){
    output_file.close();
}

void Helpers::solveTridiagonal(const Vector& ld, const Vector& d, const Vector& ud, Vector& rhs){
    int n =d.GetSize();
    Vector delta(d);
    Vector f(rhs);
    // Elimination Stage.
    for (int i=1; i<n; i++){
        delta[i] = delta[i] - ud.Read(i-1)*ld.Read(i-1)/delta[i-1];
        f[i] = f[i] - f[i-1]*ld.Read(i-1)/delta[i-1];
    }

    // Backsolve stage
    rhs[n-1] = f[n-1]/delta[n-1];
    for(int i=n-2; i>=0; i--){
        rhs[i] = (f[i]  - ud.Read(i)*rhs[i+1])/delta[i];
    }

}
