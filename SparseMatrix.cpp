#include "SparseMatrix.h"
#include <vector>

SparseMatrix::SparseMatrix() : n(0), ig(nullptr), jg(nullptr), ggl(nullptr), ggu(nullptr),
                                di(nullptr), b(nullptr),
                                Atb(nullptr), r(nullptr), z(nullptr), p(nullptr),
                                Ap(nullptr), tmp(nullptr), invDiag(nullptr)
{ }
 
bool SparseMatrix::AllocateWorkVectors()
{
    if (n <= 0) return false;

    if (Atb) delete[] Atb;
    if (r)   delete[] r;
    if (z)   delete[] z;
    if (p)   delete[] p;
    if (Ap)  delete[] Ap;
    if (tmp) delete[] tmp;
    if (tmp2) delete[] tmp2;
    if (invDiag) delete[] invDiag;

    if (ggl_row_start) delete[] ggl_row_start;
    if (ggu_row_start) delete[] ggu_row_start;

    if (ilu_di) delete[] ilu_di;
    if (ilu_ggl) delete[] ilu_ggl;
    if (ilu_ggu) delete[] ilu_ggu;

    ilu_di = new double[n];
    ilu_ggl = new double[sizeGGL];
    ilu_ggu = new double[sizeGGU];

    Atb = new double[n];
    r = new double[n];
    z = new double[n];
    p = new double[n];
    Ap = new double[n];
    tmp = new double[n];
    tmp2 = new double[n];
    invDiag = new double[n];
    ggl_row_start = new int[n];
    ggu_row_start = new int[n];

    for (int i = 0; i < n; ++i) 

    {
        Atb[i] = r[i] = z[i] = p[i] = Ap[i] = tmp[i] = tmp2[i] = tmp2[i]= 0.0;
        ggl_row_start[i] = ggu_row_start[i] = 0;
    }

    return true;
}

SparseMatrix::~SparseMatrix()
{
    if (ig)   delete[] ig;
    if (jg)   delete[] jg;
    if (ggl)  delete[] ggl;
    if (ggu)  delete[] ggu;
    if (di)   delete[] di;
    if (b)    delete[] b;

    if (Atb) delete[] Atb;
    if (r)   delete[] r;
    if (z)   delete[] z;
    if (p)   delete[] p;
    if (Ap)  delete[] Ap;
    if (tmp) delete[] tmp;
    if (invDiag) delete[] invDiag;

    if (ilu_di && ilu_di != reinterpret_cast<double*>(1)) delete[] ilu_di;
    if (ilu_ggl) delete[] ilu_ggl;
    if (ilu_ggu) delete[] ilu_ggu;

    ig = jg = nullptr;
    ggl = ggu = nullptr;
    di = b = nullptr;
    Atb = r = z = p = Ap = tmp = invDiag = nullptr;
    ilu_di = ilu_ggl = ilu_ggu = nullptr;
}
bool SparseMatrix::ReadIntArray(const string& file, int*& masData, int size)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new int[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            cerr << "Ошибка чтения: " << file << endl;
            delete[] masData;
            masData = nullptr;
            return false;
        }
    }
    return true;
}
bool SparseMatrix::ReadDoubleArray(const string& file, double*& masData, int size)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new double[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            cerr << "Ошибка чтения: " << file << endl;
            delete[] masData;
            masData = nullptr;
            return false;
        }
    }
    return true;
}
bool SparseMatrix::ReadDoubleArray(const string& file, double*& masData, int size, int& GetRealSize)
{
    if (size < 0)
    {
        cerr << "Размерность массива не действиетльна: " << endl;
        return false;
    }
    ifstream fmas(file);
    if (!fmas)
    {
        cerr << "Невозможно открыть: " << file << endl;
        return false;
    }
    masData = new double[size];
    for (int i = 0; i < size; i++)
    {
        if (!(fmas >> masData[i]))
        {
            if (fmas.eof())
            {
                break;
            }
            else
            {
                cerr << "Ошибка чтения: " << file << " на элементе " << i << endl;
                delete[] masData;
                masData = nullptr;
                return false;
            }
        }
        GetRealSize++;
    }
    return true;
}
bool SparseMatrix::ReadFromFiles(const string& kuslauFile)
{
    ifstream fkus(kuslauFile.c_str());
    if (!fkus)
    {
        cerr << "Невозможно открыть " << kuslauFile << "\n";
        return false;
    }
    fkus >> n >> maxIter >> eps;
    fkus.close();

    if (!ReadIntArray("ig.txt", ig, n + 1)) return false;
    int jgSize = ig[n];
    if (jgSize <= 0)
    {
        cerr << "Неверный jgSize (ig[n]) = " << jgSize << endl;
        return false;
    }
    if (!ReadIntArray("jg.txt", jg, jgSize)) return false;
    if (!ReadDoubleArray("ggl.txt", ggl, jgSize, sizeGGL)) return false;
    if (!ReadDoubleArray("ggu.txt", ggu, jgSize, sizeGGU)) return false;
    if (!ReadDoubleArray("di.txt", di, n)) return false;
    if (!ReadDoubleArray("pr.txt", b, n)) return false;
    AllocateWorkVectors();
    return true;
}
void SparseMatrix::Multiply(const double* x, double* y) const
{// y = A*x
    for (int i = 0; i < n; i++)
    {
        y[i] = di[i] * x[i];
    }

    int ggl_k = 0;
    int ggu_k = 0;

    for (int i = 0; i < n; i++)
    {
        int i0 = ig[i], i1 = ig[i + 1];
        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            {
                if (ggl_k < sizeGGL)
                {
                    y[i] += ggl[ggl_k++] * x[j];
                }
            }
            else if (j > i)
            {
                if (ggu_k < sizeGGU)
                {
                    y[i] += ggu[ggu_k++] * x[j];
                }
            }
        }
    }
}

void SparseMatrix::MultiplyTranspose(const double* v, double* y) const
{
    for (int i = 0; i < n; i++)
        y[i] = di[i] * v[i];

    int kgl = 0;
    int kgu = 0;

    for (int i = 0; i < n; i++)
    {
        int start = ig[i];
        int end = ig[i + 1];

        for (int p = start; p < end; p++)
        {
            int j = jg[p];

            if (j < i)
            {
                y[j] += ggl[kgl] * v[i];
                kgl++;
            }
            else if (j > i)
            {
                y[j] += ggu[kgu] * v[i];
                kgu++;
            }
        }
    }
}

void SparseMatrix::PrintDense() const
{
    cout << "Плотная матрица " << n << "x" << n << ":\n";
    int lowerUpperSize = (ig[n] - n) / 2;

    double** dense = new double* [n];

    for (int i = 0; i < n; i++)
    {
        dense[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            dense[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        dense[i][i] = di[i];
    }

    int ggl_index = 0;
    int ggu_index = 0;

    for (int i = 0; i < n; i++)
    {
        int start = ig[i];
        int end = ig[i + 1];

        for (int p = start; p < end; p++)
        {
            int j = jg[p];

            if (i > j)
            { 
                if (ggl_index < sizeGGL)
                {
                    dense[i][j] = ggl[ggl_index++];
                }
            }
            else if (i < j)
            {
                if (ggu_index < sizeGGU)
                {
                    dense[i][j] = ggu[ggu_index++];
                }
            }
        }
    }

    // Печать матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(8) << setprecision(3) << dense[i][j] << " ";
        }
        cout << endl;
    }
}
double SparseMatrix::Dot(const double* a, const double* b) const
{
    double result = 0.0;
    for (int i = 0; i < n; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

double SparseMatrix::Norm2(const double* a) const
{
    return sqrt(Dot(a, a));
}


bool SparseMatrix::SolveMsgLUFactorization(double* x)
{
    BuildILU();

    for (int i = 0; i < n; ++i) p[i] = 0.0;

    Multiply(x, tmp);
    for (int i = 0; i < n; i++) r[i] = b[i] - tmp[i];

    Solve_Ly_b(r, tmp);          // tmp = L^-1 r
    SolveLTranspose(tmp, Ap);    // Ap = L^-T tmp
    MultiplyTranspose(Ap, tmp2); // tmp2 = A^T Ap
    SolveUTranspose(tmp2, r);    // r = U^-T tmp2 -> r_tilde

    for (int i = 0; i < n; ++i) z[i] = r[i];
    double rho_prev = Dot(r, r);
    double initial_norm = sqrt(rho_prev);
    cout << "Начальная норма невязки: " << initial_norm << endl;

    if (initial_norm < eps)
    {
        cout << "Решение найдено на начальном шаге." << endl;
        return true;
	}

    for (int k = 1; k <= maxIter; k++)
    {
        Solve_Ux_y(z, tmp);
        Multiply(tmp, Ap);
        Solve_Ly_b(Ap, tmp2);
        SolveLTranspose(tmp2, tmp);
        MultiplyTranspose(tmp, Ap);
        SolveUTranspose(Ap, Atb);

        double denom = Dot(Atb, z);
        if (fabs(denom) < eps) break;

		double alpha = rho_prev / denom;

        for (int i = 0; i < n; ++i)
        {
            p[i] += alpha * z[i];
            r[i] -= alpha * Atb[i];
        }

        double rho_new = Dot(r, r);
        double current_norm = sqrt(rho_new);

        if (current_norm < eps) break;

        double beta = rho_new / rho_prev;
        rho_prev = rho_new;

        for (int i = 0; i < n; ++i) z[i] = r[i] + beta * z[i];
    }
    Solve_Ux_y(p, x);

    return true;
}

void SparseMatrix::BuildILU()
{
    for (int i = 0; i < n; ++i) ilu_di[i] = di[i];
    for (int i = 0; i < sizeGGL; ++i) ilu_ggl[i] = ggl[i];
    for (int i = 0; i < sizeGGU; ++i) ilu_ggu[i] = ggu[i];

    int k_l = 0;
    int k_u = 0;
    for (int i = 0; i < n; ++i)
    {
        ggl_row_start[i] = k_l;
        ggu_row_start[i] = k_u;

        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int p = i0; p < i1; ++p)
        {
            int j = jg[p];
            if (j < i) k_l++;
            else if (j > i) k_u++;
        }
    }
}

void SparseMatrix::Solve_Ly_b(const double* b, double* y) const
{
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;

        int i0 = ig[i];
        int i1 = ig[i + 1];
        int k = ggl_row_start[i];

        for (int p = i0; p < i1; ++p)
        {
            int j = jg[p];
            if (j < i)
            {
                sum += ilu_ggl[k] * y[j];
                k++;
            }
        }
        y[i] = (b[i] - sum) / ilu_di[i];
    }
}
void SparseMatrix::Solve_Ux_y(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0.0;

        int i0 = ig[i];
        int i1 = ig[i + 1];
        int k = ggu_row_start[i];

        int p = i0;
        while (p < i1 && jg[p] < i) p++;

        for (; p < i1; ++p) 
        {
            int j = jg[p];
            if (j > i)
            {
                sum += ilu_ggu[k] * x[j];
                k++;
            }
        }

        x[i] -= sum;
    }
}
void SparseMatrix::SolveLTranspose(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = n - 1; i >= 0; i--)
    {
        int i0 = ig[i];
        int i1 = ig[i + 1];
        int k = ggl_row_start[i];
        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            {
                x[j] -= ilu_ggl[k] * x[i];
                k++;
            }
        }
		x[i] /= ilu_di[i];
    }
}
void SparseMatrix::SolveUTranspose(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = 0; i < n; i++)
    {
        int i0 = ig[i];
        int i1 = ig[i + 1];
        int k = ggu_row_start[i];
        int p = i0;
        while (p < i1 && jg[p] < i) p++;
        for (; p < i1; p++)
        {
            int j = jg[p];
            if (j > i)
            {
                x[j] -= ilu_ggu[k] * x[i];
                k++;
            }
        }
    }
}