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

    delete[] Atb; delete[] r; delete[] z; delete[] p; delete[] Ap; delete[] tmp; delete[] invDiag; delete[] diagAtA;

    Atb = new double[n];
    r = new double[n];
    z = new double[n];
    p = new double[n];
    Ap = new double[n];
    tmp = new double[n];
    invDiag = new double[n];
    diagAtA = new double[n];
    inxL = new double[n];
    inxU = new double[n];
    for (int i = 0; i < n; ++i) 
    {
        Atb[i] = r[i] = z[i] = p[i] = Ap[i] = tmp[i] = invDiag[i] = inxU[i] = inxL[i] = 0.0;
    }

    return true;
}

SparseMatrix::~SparseMatrix()
{
    delete[] ig;  ig = nullptr;
    delete[] jg;  jg = nullptr;
    delete[] ggl; ggl = nullptr;
    delete[] ggu; ggu = nullptr;
    delete[] di;  di = nullptr;
    delete[] b;   b = nullptr;
	delete[] Atb; Atb = nullptr;
    delete[] r;   r = nullptr;
    delete[] z;   z = nullptr;
    delete[] p;   p = nullptr;
    delete[] Ap;  Ap = nullptr;
    delete[] tmp; tmp = nullptr;
    delete[] invDiag; invDiag = nullptr;
	delete[] diagAtA; diagAtA = nullptr;
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
{
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
    {
        y[i] = di[i] * v[i];
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
            { // A[i][j]
                if (ggl_k < sizeGGL)
                {
                    y[j] += ggl[ggl_k++] * v[i];
                }
            }
            else if (j > i)
            { // A[i][j]
                if (ggu_k < sizeGGU)
                {
                    y[j] += ggu[ggu_k++] * v[i];
                }
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
void SparseMatrix::PrintResults(int iter, double normR, const double* x, const double* xTrue, int n) const
{
    cout << ", Iter = " << iter;
    cout << scientific << setprecision(15) << ", normR = " << normR << endl;

    cout << "i\t x_i" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i] << endl;
    }

    cout << "i\t x* _i - x_i" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << xTrue[i] - x[i] << endl;
    }
    cout << endl << endl;
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
bool SparseMatrix::InitialSolution(double* x) const
{
    Multiply(x, tmp);
    for (int i = 0; i < n; ++i) tmp[i] = b[i] - tmp[i];
    MultiplyTranspose(tmp, r);

    double norm0 = Norm2(r);
    if (norm0 < eps)
    {
        cout << "MSG: начальное приближение является решением, 0 итераций.\n";
        return true;
    }
    return false;
}

bool SparseMatrix::SolveMSG(double* x,double*xTrue, const string& prec) const
{
    bool firstSolve = InitialSolution(x);
    if (firstSolve)
    {
        return true;
    }

    MultiplyTranspose(b, Atb);

    for (int i = 0; i < n; ++i)
    {
        //x[i] = 0.0;
        r[i] = Atb[i];
    }
    
    if (prec == "diag")
    {
        ComputeDiagAtA(diagAtA);

        for (int i = 0; i < n; ++i)
        {
            double v = diagAtA[i];
            if (fabs(v) < 1e-30) throw runtime_error("\nОШИБКА!!! Диагональный элемент равен нулю. Метод диальногальной предобусловленности не дал результат\n");
            invDiag[i] = 1.0 / v;
            z[i] = r[i] * invDiag[i];
        }
    }
    else if (prec == "ilu")
    {
        SolveUTranspose(r, p);   // p = U^-T r
        SolveLTranspose(p, Ap);  // Ap = L^-T p
        SolveL(Ap, p);           // p = L^-1 Ap
        SolveU(p, z);            // z = U^-1 p
    }
    else
    {
        for (int i = 0; i < n; ++i) z[i] = r[i];
    }


    for (int i = 0; i < n; ++i) p[i] = z[i];

    Multiply(p, tmp);           // tmp = A * p
    MultiplyTranspose(tmp, Ap); // Ap = A^T * (A * p)

    double rz_old = Dot(r, z);
    double normAtb = Norm2(Atb);
    if (normAtb < 1e-16) normAtb = 1.0;
    
    for (int iter = 0; iter < maxIter; ++iter)
    {
        Multiply(p, tmp);
        MultiplyTranspose(tmp, Ap);

        double pAp = Dot(p, Ap);
        if (fabs(pAp) < 1e-30)
        {
            cout << "Решение сошлось (pAp ~ 0) на итерации " << iter;
            return true;
        }

        double alpha = rz_old / pAp;

        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double relRes = Norm2(r) / normAtb;
        if (relRes < eps)
        {
            cout << "MSG сошёлся на итерации ";
            PrintResults(iter, relRes, x, xTrue, n);
            return true;
        }

        if (prec == "ilu")
        {
            SolveUTranspose(r, tmp);
            SolveLTranspose(tmp, Ap);
            
            SolveUTranspose(r, tmp);
            SolveLTranspose(tmp, Ap);
            SolveL(Ap, tmp);
            SolveU(tmp, z);
        }
        else if (prec == "diag")
        {
            for (int i = 0; i < n; ++i) z[i] = r[i] * invDiag[i];
        }
        else
        {
            for (int i = 0; i < n; ++i) z[i] = r[i];
        }

        double rz_new = Dot(r, z);
        double beta = rz_new / rz_old;
        rz_old = rz_new;

        for (int i = 0; i < n; ++i)
        {
            p[i] = z[i] + beta * p[i];
        }
        Multiply(p, tmp);
        MultiplyTranspose(tmp, Ap);
    }

    cerr << "Решение не сошлось за " << maxIter << " итераций.\n";
    return false;
}
bool SparseMatrix::SolveLOS(double* x, double* xTrue, const string& prec) const
{
    Multiply(x, r);
	for (int i = 0; i < n; ++i) r[i] = b[i] - r[i];

    if (prec == "diag")
    {
        ApplyDiagPreconditioner(r, z);
    }
    else if (prec == "ilu") 
    {
        ApplyILUPreconditioner(r, z);
    }
    else
    {
        for (int i = 0; i < n; ++i) 
        {
            z[i] = r[i];
        }
    }

	Multiply(z, p);

    double normB = Norm2(b);
    double normR = Norm2(r);
    double relRes = normR / normB;

    cout << "LOS начальная невязка: " << relRes << endl;

    if (relRes < eps)
    {
        cout << "LOS сошёлся на 0 итерации" << endl;
        return true;
    }

    for (int iter = 0; iter < maxIter; ++iter)
    {
        double pDotR = Dot(p, r);
        double pDotP = Dot(p, p);

        if (fabs(pDotP) < 1e-30)
        {
            cout << "LOS: (p,p) ~ 0 на итерации " << iter << endl;
            return true;
        }

		double alpha = pDotR / pDotP;

        for (int i = 0; i < n; i++)
        {
            x[i] += alpha * z[i];
        }
        for (int i = 0; i < n; i++)
        {
            r[i] -= alpha * p[i];
		}

        normR = Norm2(r);
        relRes = normR / normB;

        if (relRes < eps)
        {
            cout << "LOS сошёлся на итерации ";
            PrintResults(iter, relRes, x, xTrue, n);
            return true;
        }

        if (prec == "diag") 
        {
            ApplyDiagPreconditioner(r, tmp);
        }
        else if (prec == "ilu")
        {
            ApplyILUPreconditioner(r, tmp);
        }
        else 
        {
            for (int i = 0; i < n; ++i) 
            {
                tmp[i] = r[i];
            }
        }

        Multiply(r, Ap);
        double pDotAp = Dot(p, Ap);
		double beta = -pDotAp / pDotP;

        for (int i = 0; i < n; i++)
        {
            z[i] = r[i] + beta * z[i];
            p[i] = Ap[i] + beta * p[i];
		}
        // Я пишу этот код где-то ночью
        // Ветер шёпотом пишет мне строчки
	    
    }
    cerr << "LOS не сошлось за " << maxIter << " итераций. Текущая невязка = " << relRes << endl;
    return false;
}
void SparseMatrix::ComputeDiagAtA(double* diagAtA) const
{
    for (int i = 0; i < n; ++i) diagAtA[i] = 0.0;

    int kgl = 0;
    int kgu = 0;

    for (int i = 0; i < n; i++)
    {
        double d = di[i];
        diagAtA[i] = d * d;

		int i0 = ig[i], i1 = ig[i + 1];

        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            {
                if (kgl < sizeGGL)
                {
                    double v = ggl[kgl++];
                    diagAtA[j] += v * v;
                }
            }
            else if (j > i)
            {
                if (ggu && kgu < sizeGGU)
                {
                    double v = ggu[kgu++];
                    diagAtA[j] += v * v;
                }
            }
		}
    }
}

void SparseMatrix::ApplyDiagPreconditioner(const double* r, double* z) const
{
    for (int i = 0; i < n; ++i) 
    {
        if (fabs(di[i]) > 1e-15)
        {
            z[i] = r[i] / di[i];
        }
        else
        {
            z[i] = r[i];
        }
    }
}

void SparseMatrix::ApplyILUPreconditioner(const double* r, double* z) const
{
    SolveL(r, z);  // z = L^-1 r

    SolveU(z, z);  // z = U^-1 (L^-1 r)
}
void SparseMatrix::SolveL(const double* b, double* y) const
{
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;

        int i0 = ig[i];
        int i1 = ig[i + 1];

        int idx_l = inxL[i];

        for (int k = i0; k < i1; ++k)
        {
            int j = jg[k];
            if (j < i)
            {
                sum += ggl[idx_l] * y[j];
                idx_l++;
            }
            else break;
        }

        y[i] = (b[i] - sum) / di[i];
    }
}
void SparseMatrix::SolveU(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0.0;

        int i0 = ig[i];
        int i1 = ig[i + 1];

        int idx_u = inxU[i];

        int k = i0;
        while (k < i1 && jg[k] < i) k++;

        for (; k < i1; ++k)
        {
            int j = jg[k];
            if (j > i)
            {
                sum += ggu[idx_u] * x[j];
                idx_u++;
            }
        }

        x[i] -= sum; 
    }
}
void SparseMatrix::SolveLTranspose(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = n - 1; i >= 0; --i)
    {
        x[i] /= di[i];

        int i0 = ig[i];
        int i1 = ig[i + 1];
        int idx_l = inxL[i];

        for (int k = i0; k < i1; ++k)
        {
            int j = jg[k];
            if (j < i)
            {
                x[j] -= ggl[idx_l] * x[i];
                idx_l++;
            }
            else break;
        }
    }
}
void SparseMatrix::SolveUTranspose(const double* y, double* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = 0; i < n; ++i)
    {
        int i0 = ig[i];
        int i1 = ig[i + 1];
        int idx_u = inxU[i];

        int k = i0;
        while (k < i1 && jg[k] < i) k++;

        for (; k < i1; ++k)
        {
            int j = jg[k];
            if (j > i)
            {
                x[j] -= ggu[idx_u] * x[i];
                idx_u++;
            }
        }
    }
}
void SparseMatrix::PreparePreconditioner()
{
    int cnt_l = 0;
    int cnt_u = 0;
    for (int i = 0; i < n; ++i)
    {
        inxL[i] = cnt_l;
        inxU[i] = cnt_u;

        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int k = i0; k < i1; ++k)
        {
            int j = jg[k];
            if (j < i) cnt_l++;
            else if (j > i) cnt_u++;
        }
    }
}