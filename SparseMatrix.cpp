// SparseMatrix.cpp

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

    Atb = new real[n];
    r = new real[n];
    z = new real[n];
    p = new real[n];
    Ap = new real[n];
    tmp = new real[n];
    invDiag = new real[n];
    diagAtA = new real[n];
    inxL = new real[n];
    inxU = new real[n];
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
bool SparseMatrix::ReadrealArray(const string& file, real*& masData, int size)
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
    masData = new real[size];
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
bool SparseMatrix::ReadrealArray(const string& file, real*& masData, int size, int& GetRealSize)
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
    masData = new real[size];
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
    if (!ReadrealArray("ggl.txt", ggl, jgSize, sizeGGL)) return false;
    if (!ReadrealArray("ggu.txt", ggu, jgSize, sizeGGU)) return false;
    if (!ReadrealArray("di.txt", di, n)) return false;
    if (!ReadrealArray("pr.txt", b, n)) return false;
    AllocateWorkVectors();
    return true;
}

void SparseMatrix::GenerateHilbertMatrix(int Size)
{
    this->n = Size;

    delete[] ig; delete[] jg; delete[] ggl; delete[] ggu;
    delete[] di; delete[] b;

    long long nz = (long long)n * (n - 1);

    sizeGGL = nz / 2;
    sizeGGU = nz / 2;

    ig = new int[n + 1];
    jg = new int[nz];
    ggl = new real[sizeGGL];
    ggu = new real[sizeGGU];
    di = new real[n];
    b = new real[n];

    ig[0] = 0;
    int k = 0;      // индекс для jg
    int k_l = 0;    // индекс для ggl
    int k_u = 0;    // индекс для ggu

    for (int i = 0; i < n; i++)
    {
        di[i] = 1.0 / (real)(i + i + 1);

        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;

            jg[k] = j;

            real val = 1.0 / (real)(i + j + 1);

            if (j < i)
            {
                ggl[k_l++] = val;
            }
            else
            {
                ggu[k_u++] = val;
            }
            k++;
        }
        ig[i + 1] = k;
    }

    AllocateWorkVectors();

    real* xTrueTemp = new real[n];
    for (int i = 0; i < n; i++) xTrueTemp[i] = i + 1.0;

    Multiply(xTrueTemp, b);

    delete[] xTrueTemp;
}

void SparseMatrix::Multiply(const real* x, real* y) const
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
void SparseMatrix::MultiplyTranspose(const real* v, real* y) const
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

    real** dense = new real* [n];

    for (int i = 0; i < n; i++)
    {
        dense[i] = new real[n];
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
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(8) << setprecision(3) << dense[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << "b = ";
    for (int i = 0; i < n; i++)
    {
        cout << setw(8) << setprecision(3) << b[i] << " ";
    }
    cout << endl;
}
void SparseMatrix::PrintResults(int iter, real normR, const real* x, const real* xTrue, int n) const
{
    cout << "Iter = " << iter << ", normR = " << scientific << setprecision(6) << normR
        << ", RelRes: " << GetRelativeResidual(x) << endl;

    cout << "\n------------------------------" << endl;
    cout << "|          x_calc" << endl;
    cout << "------------------------------" << endl;

    for (int i = 0; i < n; i++)
    {
        cout << scientific << setprecision(15) << setw(25) << x[i] << endl;
    }
    cout << "\n------------------------------" << endl;
    cout << "|          x error" << endl;
    cout << "------------------------------" << endl;
    for (int i = 0; i < n; i++)
    {
        real diff = xTrue[i] - x[i];
        if (fabs(diff) < 1e-16)
        {
            cout << setw(25) << "0.0";
        }
        else
        {
            double order = floor(log10(fabs(diff)));

            int magnitude = (int)fabs(order);

            int prec = 15 - magnitude;

            if (prec < 1) prec = 1;
            if (prec > 15) prec = 15;

            cout << scientific << setprecision(prec) << setw(25) << diff;
        }

        cout << endl;
    }

    cout << "--------------------------------------------------------\n" << endl;

    cout << defaultfloat << setprecision(6);
}

real SparseMatrix::Dot(const real* a, const real* b) const
{
    real result = 0.0;
    for (int i = 0; i < n; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}
real SparseMatrix::Norm2(const real* a) const
{
    return sqrt(Dot(a, a));
}
bool SparseMatrix::InitialSolution(real* x) const
{
    Multiply(x, tmp);
    for (int i = 0; i < n; ++i) tmp[i] = b[i] - tmp[i];
    MultiplyTranspose(tmp, r);

    real norm0 = Norm2(r);
    if (norm0 < eps)
    {
        cout << "MSG: начальное приближение является решением, 0 итераций.\n";
        return true;
    }
    return false;
}

bool SparseMatrix::SolveMSG(real* x,real*xTrue, const string& prec) const
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
            real v = diagAtA[i];
            if (fabs(v) < EPS) throw runtime_error("\nОШИБКА!!! Диагональный элемент равен нулю. Метод диальногальной предобусловленности не дал результат\n");
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

    real rz_old = Dot(r, z);
    real normAtb = Norm2(Atb);
    //if (normAtb < 1e-16) normAtb = 1.0;
    
    for (int iter = 0; iter < maxIter; ++iter)
    {
        Multiply(p, tmp);
        MultiplyTranspose(tmp, Ap);

        real pAp = Dot(p, Ap);
        if (fabs(pAp) < EPS)
        {
            cout << "Решение сошлось (pAp ~ 0) на итерации " << iter;
            return true;
        }

        real alpha = rz_old / pAp;

        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        real relRes = Norm2(r) / normAtb;
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

        real rz_new = Dot(r, z);
        real beta = rz_new / rz_old;
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
bool SparseMatrix::SolveLOS(real* x, real* xTrue, const string& prec) const
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

    real normB = Norm2(b);
    real normR = Norm2(r);
    real relRes = normR / normB;

    if (relRes < eps)
    {
        cout << "LOS сошёлся на 0 итерации" << endl;
        return true;
    }

    for (int iter = 0; iter < maxIter; ++iter)
    {
        real pDotR = Dot(p, r);
        real pDotP = Dot(p, p);

        if (fabs(pDotP) < EPS)
        {
            cout << "LOS: (p,p) ~ 0 на итерации " << iter << endl;
            return true;
        }

		real alpha = pDotR / pDotP;

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
        real pDotAp = Dot(p, Ap);
		real beta = -pDotAp / pDotP;

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
void SparseMatrix::ComputeDiagAtA(real* diagAtA) const
{
    for (int i = 0; i < n; ++i) diagAtA[i] = 0.0;

    int kgl = 0;
    int kgu = 0;

    for (int i = 0; i < n; i++)
    {
        real d = di[i];
        diagAtA[i] = d * d;

		int i0 = ig[i], i1 = ig[i + 1];

        for (int p = i0; p < i1; p++)
        {
            int j = jg[p];
            if (j < i)
            {
                if (kgl < sizeGGL)
                {
                    real v = ggl[kgl++];
                    diagAtA[j] += v * v;
                }
            }
            else if (j > i)
            {
                if (ggu && kgu < sizeGGU)
                {
                    real v = ggu[kgu++];
                    diagAtA[j] += v * v;
                }
            }
		}
    }
}

void SparseMatrix::ApplyDiagPreconditioner(const real* r, real* z) const
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
void SparseMatrix::ApplyILUPreconditioner(const real* r, real* z) const
{
    SolveL(r, z);  // z = L^-1 r

    SolveU(z, z);  // z = U^-1 (L^-1 r)
}
void SparseMatrix::SolveL(const real* b, real* y) const
{
    for (int i = 0; i < n; ++i)
    {
        real sum = 0.0;

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
void SparseMatrix::SolveU(const real* y, real* x) const
{
    for (int i = 0; i < n; ++i) x[i] = y[i];

    for (int i = n - 1; i >= 0; --i)
    {
        real sum = 0.0;

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
void SparseMatrix::SolveLTranspose(const real* y, real* x) const
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
void SparseMatrix::SolveUTranspose(const real* y, real* x) const
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

void SparseMatrix::WriteSolution(const string& filename, const real* x) const
{
    ofstream fout(filename);
    if (!fout.is_open())
    {
        cerr << "Ошибка открытия файла для записи: " << filename << endl;
        return;
    }

    fout << scientific << setprecision(15);
    for (int i = 0; i < n; i++)
    {
        fout << x[i] << endl;
    }
    fout.close();
    cout << "Решение записано в файл " << filename << endl;
}
real SparseMatrix::GetRelativeResidual(const real* x) const
{
    Multiply(x, r); 
    real normNumerator = 0.0; 
    real normDenominator = 0.0; 

    for (int i = 0; i < n; i++)
    {
        real residue = b[i] - r[i];
        normNumerator += residue * residue;
        normDenominator += b[i] * b[i];
    }

    if (fabs(normDenominator) < EPS) return 0.0; 

    return sqrt(normNumerator) / sqrt(normDenominator);
}