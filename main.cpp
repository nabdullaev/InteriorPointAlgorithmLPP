#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
const double eps = 1e-5;

class Matrix {
public:
    int n, m;
    vector<vector<double>> matrix;
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        matrix.resize(n);
        for (int i = 0; i < n; i++) {
            matrix[i].resize(m);
        }
    }
    Matrix operator+(Matrix B) {
        Matrix *C = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                C->matrix[i][j] = matrix[i][j] + B.matrix[i][j];
            }
        }
        return *C;
    }
    Matrix operator-(Matrix B) {
        Matrix *C = new Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                C->matrix[i][j] = matrix[i][j] - B.matrix[i][j];
            }
        }
        return *C;
    }
    Matrix operator*(Matrix B) {
        Matrix *C = new Matrix(n, B.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < B.m; j++) {
                for (int k = 0; k < m; k++) {
                    C->matrix[i][j] += matrix[i][k] * B.matrix[k][j];
                }
            }
        }
        return *C;
    }
    Matrix operator=(Matrix B) {
        n = B.n;
        m = B.m;
        matrix.resize(n);
        for (int i = 0; i < n; i++) {
            matrix[i].resize(m);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                matrix[i][j] = B.matrix[i][j];
            }
        }
        return *this;
    }
    // overloading operator <<
    friend ostream& operator<<(ostream& os, const Matrix& B) {
        for (int i = 0; i < B.n; i++) {
            for (int j = 0; j < B.m; j++) {
                if (abs(B.matrix[i][j]) < eps)
                    os << fixed << setprecision(4) << 0.00 << " ";
                else
                    os << fixed << setprecision(4) << B.matrix[i][j] << " ";
            }
            os << endl;
        }
        return os;
    }
    // overloading operator >>
    friend istream& operator>>(istream& is, Matrix& B) {
        for (int i = 0; i < B.n; i++) {
            for (int j = 0; j < B.m; j++) {
                is >> B.matrix[i][j];
            }
        }
        return is;
    }
    Matrix T() {
        Matrix *C = new Matrix(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                C->matrix[j][i] = matrix[i][j];
            }
        }
        return *C;
    }

    double norm(int a) {
        double sum = 0;
        for (int i = 0; i < a; i++) {
            double temp = 0;
            for (int j = 0; j < m; j++) {
                temp += matrix[i][j] * matrix[i][j];
            }
            sum = max(sum, temp);
        }
        return sqrt(sum);
    }
};

Matrix projection(Matrix A_hat, int m) {
    Matrix I(m, m);
    for (int i = 0; i < m; i++)
        I.matrix[i][i] = 1;
    Matrix A_hat_T = A_hat.T();
    Matrix A_hat_A_hat_T = A_hat * A_hat_T;

    Matrix inv = A_hat_A_hat_T;
    Matrix inverted(inv.n, 2 * inv.m);
    for (int i = 0; i < inv.n; i++) {
        for (int j = 0; j < inv.m; j++) {
            inverted.matrix[i][j] = inv.matrix[i][j];
        }
    }
    for (int i = 0; i < inv.n; i++) {
        for (int j = inv.m; j < 2 * inv.m; j++) {
            if (i == j - inv.m)
                inverted.matrix[i][j] = 1;
            else
                inverted.matrix[i][j] = 0;
        }
    }
    for (int i = 0; i < inv.n; i++) {
        // Swap the row with the maximum element row
        double max = inverted.matrix[i][i];
        int max_row = i;
        for (int j = i + 1; j < inv.n; j++) {
            if (abs(inverted.matrix[j][i]) > abs(max)) {
                max = inverted.matrix[j][i];
                max_row = j;
            }
        }
        if (max_row != i) {
            swap(inverted.matrix[i], inverted.matrix[max_row]);
        }
        // Make all the elements below the pivot equal to zero
        for (int j = i + 1; j < inv.n; j++) {
            if (abs(inverted.matrix[j][i]) < eps)
                continue;
            double c = inverted.matrix[j][i] / inverted.matrix[i][i];
            for (int k = i; k < 2 * inv.m; k++) {
                inverted.matrix[j][k] -= inverted.matrix[i][k] * c;
            }
        }
    }

    // Make all the elements above the pivot equal to zero
    for (int i = inv.n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (abs(inverted.matrix[j][i]) < eps)
                continue;
            double c = inverted.matrix[j][i] / inverted.matrix[i][i];
            for (int k = 2 * inv.m - 1; k >= i; k--) {
                inverted.matrix[j][k] -= inverted.matrix[i][k] * c;
            }
        }
    }
    // Divide the row by the pivot element
    for (int i = 0; i < inv.n; i++) {
        double c = inverted.matrix[i][i];
        for (int j = i; j < 2 * inv.m; j++) {
            inverted.matrix[i][j] /= c;
        }
    }
    Matrix invA(inv.n, inv.m);
    for (int i = 0; i < inv.n; i++) {
        for (int j = inv.m; j < 2 * inv.m; j++) {
            invA.matrix[i][j - inv.m] = inverted.matrix[i][j];
        }
    }
    Matrix P = I - A_hat_T * invA * A_hat;
    return P;
}

//check if the solution is out of original bounds
bool out_of_boundaries(Matrix& x, Matrix& A, Matrix& b){
    int n = A.n;
    int m = A.m-n;
    for(int i = 0; i < n; ++i){
        double left_hand_side = 0;
        for(int j = 0; j < m; ++j){
            left_hand_side += x.matrix[j][0]*A.matrix[i][j];
        }
        if(left_hand_side > b.matrix[i][0]){
            return true;
        }
    }
    return false;
}


int main() {
    int max_iterations = 1000;
    int iteration_count = 0;
    // This is the program for interior point method to solve linear programming problem
    // Create a matrix A and c and assign values to them
    cout << "Enter the number of rows and columns of matrix A (n x m): ";
    int n, m;
    cin >> n >> m;
    cout << "Enter the matrix A (matrix of coefficients of constraint function) (n x m): ";
    Matrix A(n, m);
    cin >> A;
    cout << "Enter the matrix c (vector of coefficients of objective function): ";
    Matrix c(m, 1);
    cin >> c;
    // Create a diagonal matrix D which represents the initial feasible solution
    cout << "Enter the vector of right-hand side numbers b: ";
    Matrix b(n, 1);
    cin >> b;
    vector<bool> need_superplus_variable(n, false);
    for(int i = 0; i < n; ++i){
        if(b.matrix[i][0] < 0){
            b.matrix[i][0] *= -1;
            for(int j=0; j < m; ++j){
                A.matrix[i][j] *= -1;
            }
            need_superplus_variable[i] = true;
        }
    }
    // Read accuracy
    double accuracy;
    cout << "Enter the approximation accuracy: ";
    cin >> accuracy;
    Matrix D(m+n, m+n);
    cout << "Enter initial solution with " << n+m << " variables\n";
    for(int i = 0; i < m+n; ++i){
        cout << "x" << i+1 << ": ";
        cin >> D.matrix[i][i];

    }
    Matrix x_init(n+m, 1);
    for(int i = 0; i < m+n; ++i){
        x_init.matrix[i][0] = D.matrix[i][i];
    }
    if(out_of_boundaries(x_init, A, b)){
        cout << "Method is not applicable\n";
        return 0;
    }


    Matrix A_new(n, m+n);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            A_new.matrix[i][j] = A.matrix[i][j];
        }
    }
    for(int j = m; j < m+n; ++j){
        A_new.matrix[j-m][j] = 1;
    }

    A = A_new;

    Matrix c_new(m+n, 1);
    for(int i = 0; i < m; ++i){
        c_new.matrix[i][0] = c.matrix[i][0];
    }
    c = c_new;

    Matrix x_prev(m+n, 1);
    for (int i = 0; i < m+n; i++)
        x_prev.matrix[i][0] = 0;

    while (iteration_count <= max_iterations) {
        Matrix D_inverse(m+n, m+n);
        for(int i = 0; i < n+m; ++i){
            D_inverse.matrix[i][i] = 1/D.matrix[i][i];
        }

        // Calculate A_hat = A * D
        Matrix A_hat = A * D;
        // Calculate c_hat = D * c
        Matrix c_hat = D * c;
        // Calculate P = I - A_hat^T * (A_hat * A_hat^T)^(-1) * A_hat
        Matrix P = projection(A_hat, m+n);

        // Calculate c_p = P * c_hat
        Matrix c_p = P * c_hat;
        // Find number v such that v = min(c_p[i]), for all c_p[i] < 0
        double v = 0;
        for (int i = 0; i < c_p.n; i++) {
            if (c_p.matrix[i][0] < v)
                v = c_p.matrix[i][0];
        }
        // Find x_hat = O + (alpha / v) * c_p
        double alpha = 0.5;
        for (int i = 0; i < c_p.n; i++)
            c_p.matrix[i][0] = 1 + (alpha / abs(v)) * c_p.matrix[i][0];

        Matrix x = D * c_p;
        for (int i = 0; i < D.m; i++)
            D.matrix[i][i] = x.matrix[i][0];

        Matrix diff = x - x_prev;
        if (diff.norm(n) < accuracy || out_of_boundaries(x, A, b)) {
            x_prev = x;
            break;
        }

        x_prev = x;
        iteration_count++;
    }
    if(iteration_count > max_iterations){
        cout << "The problem does not have solution!\n";
    } else {
        cout << "The optimal solution is the vector x with the corresponding elements: " << endl;
        double result = 0;
        for(int i = 0; i < m; ++i){
            if (abs(x_prev.matrix[i][0]) < eps)
                cout << fixed << setprecision(5)  << 0.00<<'\n';
            else
                cout << fixed << setprecision(5)  << x_prev.matrix[i][0]<<'\n';
            result+=c.matrix[i][0]*x_prev.matrix[i][0];
        }
        cout << "Objective function value is " << result;
    }
    return 0;
}
