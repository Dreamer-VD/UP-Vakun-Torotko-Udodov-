#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

void printMatrix(const Matrix& mat) // Функция для вывода матрицы на консоль.
{
    if (mat.empty()) {
        cout << "Матрица пуста для вывода.\n";
        return;
    }
    for (const auto& row : mat)
    {
        for (double val : row)
        {
            cout << fixed << setprecision(2) << setw(8) << val << " ";
        }
        cout << endl;
    }
}

Matrix getMinor(const Matrix& mat, int skip_row, int skip_col) // Функция для получения минора матрицы.
{
    int n = mat.size();
    if (n <= 1)
    {
        return Matrix();
    }
    Matrix minor(n - 1, vector<double>(n - 1));
    int r = 0;
    for (int i = 0; i < n; ++i)
    {
        if (i == skip_row) continue;
        int c = 0;
        for (int j = 0; j < n; ++j)
        {
            if (j == skip_col) continue;
            minor[r][c] = mat[i][j];
            c++;
        }
        r++;
    }
    return minor;
}

double calculateDeterminant(const Matrix& mat) // Рекурсивная функция для вычисления определителя матрицы.
{
    int n = mat.size();

    if (n == 1)
    {
        return mat[0][0];
    }

    if (n == 2)
    {
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    double determinant = 0;
    for (int k = 0; k < n; ++k)
    {
        Matrix minor = getMinor(mat, 0, k);
        determinant += pow(-1, k) * mat[0][k] * calculateDeterminant(minor);
    }
    return determinant;
}

Matrix createMatrixForKramer(const Matrix& A, const Vector& b, int columnIndexToReplace) // Функция для создания вспомогательной матрицы для метода Крамера.
{
    int n = A.size();
    if (n == 0 || A[0].size() != n || b.size() != n)
    {
        cerr << "Ошибка: Некорректные размеры матрицы A или вектора b для createMatrixForKramer.\n";
        return Matrix();
    }

    Matrix Ak = A;
    for (int i = 0; i < n; ++i)
    {
        Ak[i][columnIndexToReplace] = b[i];
    }
    return Ak;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Программа для решения СЛАУ методом Крамера.\n";

    int n;
    while (true)
    {
        cout << "Введите порядок системы (количество неизвестных n не должно превышать 3): ";
        if (!(cin >> n))
        {
            cout << "Ошибка ввода. Пожалуйста, введите целое число.\n";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        else if (n <= 0 || n >= 4)
        {
            cout << "Порядок системы должен быть 1, 2 или 3.\n";
        }
        else
        {
            break;
        }
    }

    Matrix A(n, Vector(n));
    Vector b(n);
    Vector x(n);

    cout << "Введите элементы матрицы коэффициентов A (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            while (true)
            {
                cout << "A[" << i + 1 << "][" << j + 1 << "]: ";
                if (!(cin >> A[i][j]))
                {
                    cout << "Ошибка ввода. Пожалуйста, введите число.\n";
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                else
                {
                    break;
                }
            }
        }
    }

    cout << "Введите элементы столбца свободных членов b (" << n << "x1):\n";
    for (int i = 0; i < n; ++i)
    {
        while (true)
        {
            cout << "b[" << i + 1 << "]: ";
            if (!(cin >> b[i]))
            {
                cout << "Ошибка ввода. Пожалуйста, введите число.\n";
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
            else
            {
                break;
            }
        }
    }

    cout << "\nВведенная матрица A:\n";
    printMatrix(A);
    cout << "\nВведенный столбец свободных членов b:\n";
    for (int i = 0; i < n; ++i)
    {
        cout << "b[" << i + 1 << "] = " << fixed << setprecision(2) << b[i] << endl;
    }

    double detA = calculateDeterminant(A);
    cout << "\nГлавный определитель системы det(A) = " << fixed << setprecision(4) << detA << endl;

    const double epsilon = 1e-9;
    if (abs(detA) < epsilon)
    {
        bool all_det_i_zero = true;
        for (int j = 0; j < n; ++j)
        {
            Matrix Aj = createMatrixForKramer(A, b, j);
            if (abs(calculateDeterminant(Aj)) > epsilon)
            {
                all_det_i_zero = false;
                break;
            }
        }
        if (all_det_i_zero)
        {
            cout << "Система имеет бесконечно много решений (det(A)=0 и все det(Aj)=0)." << endl;
        }
        else
        {
            cout << "Система не имеет решений (det(A)=0 и хотя бы один det(Aj)!=0)." << endl;
        }
    }
    else
    {
        cout << "\nРешения системы (x_i = det(A_i) / det(A)):\n";
        for (int j = 0; j < n; ++j)
        {
            Matrix Aj = createMatrixForKramer(A, b, j);

            double detAj = calculateDeterminant(Aj);
            x[j] = detAj / detA;
            cout << "x[" << j + 1 << "] = " << fixed << setprecision(4) << detAj << " / " << detA << " = " << x[j] << endl;
        }

        cout << "\nВектор решений x:\n";
        for (int i = 0; i < n; ++i)
        {
            cout << "x[" << i + 1 << "] = " << fixed << setprecision(4) << x[i] << endl;
        }
    }

    return 0;
}
