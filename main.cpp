#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
const double eps = 1e-9;

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

    double norm() {
        double sum = 0;
        for (int i = 0; i < n; i++) {
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



int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
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
    Matrix D(m, m);
    for (int i = 0; i < m; i++)
        cin >> D.matrix[i][i];

    // Read accuracy
    double accuracy;
    cout << "Enter the approximation accuracy: ";
    cin >> accuracy;

    Matrix x_prev(m, 1);
    for (int i = 0; i < m; i++)
        x_prev.matrix[i][0] = 0;
    while (true) {
        // Calculate A_hat = A * D
        Matrix A_hat = A * D;
        // Calculate c_hat = D * c
        Matrix c_hat = D * c;
        // Calculate P = I - A_hat^T * (A_hat * A_hat^T)^(-1) * A_hat
        Matrix P = projection(A_hat, m);
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
        if (diff.norm() < accuracy) {
            x_prev = x;
            break;
        }

        x_prev = x;
    }
    cout << "The optimal solution is the vector x with the corresponding elements: " << endl;
    cout << x_prev << endl;
    return 0;
}
