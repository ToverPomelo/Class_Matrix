#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
using namespace std;
#include "matrix.hpp"
#define ACCURACY 0.00000000001

void error(int mark)
{
    switch(mark)
    {
        case 1:std::cout<<"Row or column should be more than zero!"<<std::endl;
            break;
        case 2:std::cout<<"Format error!"<<std::endl;
            break;
        case 3:std::cout<<"Inner bracket error!"<<std::endl;
            break;
        case 4:std::cout<<"Outer bracket error!"<<std::endl;
            break;
        case 5:std::cout<<"Size not match!"<<std::endl;
            break;
        case 6:std::cout<<"Warning: row not exit!"<<std::endl;
            break;
        case 7:std::cout<<"Not invertable!"<<std::endl;
            break;
        case 8:std::cout<<"Not projectable!"<<std::endl;
            break;
        default:break;
    }
    exit(0);
}

Matrix::Matrix()
{
    m_row = 0;
    m_col = 0;
    matrix_p = NULL;
}

Matrix::Matrix(int row,int col)
{
    if(row>0 && col >0)
    {
        m_row = row;
        m_col = col;
        matrix_p = new double*[m_row];

        for(int r=0;r<m_row;r++)
        {
            matrix_p[r] = new double[m_col];

            for(int c=0;c<m_col;c++)
            {
                matrix_p[r][c] = 0;   //Set zero
            }
        }
    }
    else
    {
        error(1);
    }
}

Matrix::Matrix(int row,int col,double *input_matrix)
{
    std::cout<<sizeof(&input_matrix)/sizeof(double)<<std::endl;

    if(row>0 && col >0)
    {
        m_row = row;
        m_col = col;
        matrix_p = new double*[m_row];

        for(int r=0;r<m_row;r++)
        {
            matrix_p[r] = new double[m_col];

            for(int c=0;c<m_col;c++)
            {
                matrix_p[r][c] = input_matrix[r*m_col+c];
            }
        }
    }
    else
    {
        error(1);
    }
}

//"{{*,...,*},...,{*,...,*}}"
Matrix::Matrix(string psd_m)
{
    m_row = 0;
    m_col = 0;
    int col_iter = 0;

    std::vector<double> iter;
    int digit_left = 0;
    int digit_right = 0;
    double num = 0;
    bool in_bracket = 0;
    bool in_bracket_out = 0;
    bool is_double = 0;

    for(int i=0;i<psd_m.length();i++)
    {
        if(psd_m[i] == ' ')
            continue;
        else if(psd_m[i]=='{' && !in_bracket_out && !in_bracket)
        {
            in_bracket_out = 1;
            continue;
        }
        else if(psd_m[i]=='}' && in_bracket_out && !in_bracket)
        {
            m_row++;

            in_bracket_out = 0;
            continue;
        }
        else if(psd_m[i]==',' && !in_bracket && in_bracket_out)
        {
            m_row++;
            continue;
        }

        if(in_bracket_out)
        {
            if(psd_m[i]=='{' && !in_bracket)
            {
                in_bracket = 1;
                continue;
            }
            else if(psd_m[i]=='}' && in_bracket)
            {
                col_iter++;
                if(!m_col) m_col = col_iter;
                if(m_col != col_iter)
                {
                    error(5);
                    m_col = 0;
                    m_row = 0;
                    break;
                }
                col_iter = 0;

                iter.push_back(num);
                digit_left = 0;
                digit_right = 0;
                is_double = 0;
                num = 0;
                in_bracket = 0;
                continue;
            }
            else if(psd_m[i]==',' && in_bracket)
            {
                col_iter++;

                iter.push_back(num);
                digit_left = 0;
                digit_right = 0;
                is_double = 0;
                num = 0;
            }
            else if(psd_m[i]=='.' && in_bracket)
            {
                is_double = 1;
                continue;
            }
            else if(psd_m[i]>='0' && psd_m[i]<='9' && in_bracket && !is_double)
            {
                int iter_num = (int)psd_m[i] - (int)'0';
                num = num*pow(10,digit_left) + iter_num;
                digit_left++;
            }
            else if(psd_m[i]>='0' && psd_m[i]<='9' && in_bracket && is_double)
            {
                digit_right++;
                int iter_num = (int)psd_m[i] - (int)'0';
                num = num + iter_num*pow(10,(-1)*digit_right);
            }
            else if(in_bracket)
            {
                error(3);
                break;
            }
            else
            {
                error(4);
                break;
            }
        }
    }
    //cout<<m_col<<" "<<m_row;    //check
    matrix_p = new double*[m_row];
    for(int r=0;r<m_row;r++)
    {
        matrix_p[r] = new double[m_col];
        for(int c=0;c<m_col;c++)
        {
            matrix_p[r][c] = iter[r*m_col+c];
        }
    }
}

Matrix::Matrix(string mark_s,int mark_i)
{
    if(mark_s=="I")
    {
        m_row = mark_i;
        m_col = mark_i;

        matrix_p = new double*[m_row];
        for(int r=0;r<m_row;r++)
        {
            matrix_p[r] = new double[m_col];
            for(int c=0;c<m_col;c++)
            {
                if(c == r)
                {
                    matrix_p[r][c] = 1;
                }
                else
                {
                    matrix_p[r][c] = 0;
                }
            }
        }
    }
}

Matrix::Matrix(const Matrix &other)
{
    m_row = other.m_row;
    m_col = other.m_col;

    matrix_p = new double*[m_row];
    for(int r=0;r<m_row;r++)
    {
        matrix_p[r] = new double[m_col];
        for(int c=0;c<m_col;c++)
        {
            matrix_p[r][c] = other.matrix_p[r][c];
        }
    }
}

void Matrix::free_matrix()
{
    for(int r=0;r<m_row;r++)
    {
        delete matrix_p[r];
    }
    delete matrix_p;
    m_row = 0;
    m_col = 0;
    //cout<<"Emm~~"<<endl;
}

Matrix::~Matrix()
{
    //cout<<"Wow~~"<<endl;
    //print();
    free_matrix();
    //cout<<"Wooooo~"<<endl;
}

void Matrix::scan()
{
    std::cout<<"Please input your matrix:"<<std::endl;
    for(int r=0;r<m_row;r++)
    {
        for(int c=0;c<m_col;c++)
            std::cin>>matrix_p[r][c];
    }
}

void Matrix::print() const
{
    for(int r=0;r<m_row;r++)
    {
        for(int c=0;c<m_col;c++)
            std::cout<<matrix_p[r][c]<<"\t";
        std::cout<<std::endl;
    }
}

int Matrix::row() const
{
    return m_row;
}

int Matrix::column() const
{
    return m_col;
}

int Matrix::rank() const
{   //similar to to_u()
    Matrix newMatrix(*this);
    int pivot = 0;
    for(int c=0;c<m_col;c++)
    {
        int mark_r = -1;
        double mark_num = (numeric_limits<double>::max)();
        double quotient = 1;
        bool isn_zero[m_row] = {0};
        for(int r=pivot;r<m_row;r++)
        {
            if(newMatrix[r][c]!=0)
            {
                isn_zero[r] = 1;
            }
            if(newMatrix[r][c]<mark_num && newMatrix[r][c]!=0)
            {
                mark_num = newMatrix[r][c];
                mark_r = r;
            }
        }
        if(mark_r != -1) pivot++;
        if(pivot>=m_row || pivot>=m_col) break;
        if(mark_r!=0 && mark_r!=-1) //change to the first line
        {
            for(int i=pivot-1;i<m_col;i++)
            {
                newMatrix[pivot-1][i] = newMatrix[pivot-1][i] + newMatrix[mark_r][i];
                newMatrix[mark_r][i] = newMatrix[pivot-1][i] - newMatrix[mark_r][i];
                newMatrix[pivot-1][i] = newMatrix[pivot-1][i] - newMatrix[mark_r][i];
            }
        }

        if(mark_r != -1)           //multiply and add
        {
            for(int r=pivot;r<m_row;r++)
            {
                if(!isn_zero[r]) continue;
                quotient = newMatrix[r][c]/newMatrix[pivot-1][c];
                for(int i=pivot-1;i<m_col;i++)
                {
                    newMatrix[r][i] -= newMatrix[pivot-1][i]*quotient;
                    if(abs(newMatrix[r][i])<ACCURACY)
                        newMatrix[r][i]=0;
                }
            }
        }
    }
    return pivot;
}

/*****************************Basis_get**************************/

Matrix Matrix::get_transpose() const
{
    Matrix newMatrix(m_col,m_row);
    for(int r=0;r<m_row;r++)
    {
        for(int c=0;c<m_col;c++)
        {
            newMatrix[c][r] = matrix_p[r][c];
        }
    }
    return newMatrix;
}

Matrix Matrix::get_rref() const
{
    Matrix newMatrix(*this);

    double quotient = 1;
    for(int c=0;c<m_col;c++)
    {
        if(newMatrix[c][c] == 0)
        {
            int r2;
            for(r2=c+1;r2<m_row;r2++)
            {
                if(newMatrix[r2][c] != 0)
                    break;
            }
            if(r2 != c)
            {
                for(int i=c;i<m_col;i++)   //because all before c is 0
                {
                    newMatrix[c][i] += newMatrix[r2][i];
                    if(abs(newMatrix[c][i])<ACCURACY)
                        newMatrix[c][i]=0;
                }
            }
        }
//
        for(int r=c+1;r<m_row;r++)
        {
            if(newMatrix[r][c] == 0) continue;
            quotient = newMatrix[r][c]/newMatrix[c][c];
            for(int i=c;i<m_col;i++)
            {
                newMatrix[r][i] -= newMatrix[c][i]*quotient;
                if(abs(newMatrix[r][i])<ACCURACY)
                    newMatrix[r][i]=0;
            }
        }

        if(c+1>=m_row)
            break;
    }
//rref:
    for(int r=0;r<m_row;r++)
    {
        for(int c=r;c<m_col;c++)
        {
            if(newMatrix[r][c] == 0) continue;

            for(int r2=0;r2<r;r2++)
            {
                quotient = newMatrix[r2][c]/newMatrix[r][c];
                for(int i=c;i<m_col;i++)
                {
                    newMatrix[r2][i] -= newMatrix[r][i]*quotient;
                    if(abs(newMatrix[r2][i])<ACCURACY)
                        newMatrix[r2][i]=0;
                }
            }
            break;
        }
    }
    return newMatrix;
}

Matrix Matrix::get_augmented(const Matrix &other) const
{
    if(m_row == other.m_row)
    {
        Matrix newMatrix(m_row,m_col+other.m_col);
        for(int r=0;r<m_row;r++)
        {
            for(int c=0;c<m_col+other.m_col;c++)
            {
                if(c<m_col)
                {
                    newMatrix[r][c] = matrix_p[r][c];
                }
                else
                {
                    newMatrix[r][c] = other.matrix_p[r][c-m_col];
                }
            }
        }
        return newMatrix;
    }
    else
    {
        error(5);
    }
}

Matrix Matrix::get_invert() const
{
    if(m_row==m_col && rank()==m_row)     //inefficient
    {
        if(m_row == 1)
        {
            Matrix newMatrix(1,1);
            newMatrix[0][0] = 1/matrix_p[0][0];
            return newMatrix;
        }
        Matrix i("I",m_row);
        Matrix arg = get_augmented(i);
        arg.to_rref();

        Matrix newMatrix(m_row,m_col);
        for(int r=0;r<m_row;r++)
        {
            for(int c=0;c<m_col;c++)
            {
                newMatrix[r][c] = arg.matrix_p[r][c+m_col];
            }
        }
        return newMatrix;
    }
    else
    {
        error(7);
    }
}

/************************Basis_to******************************/

void Matrix::to_transpose()
{
    *this = (*this).get_transpose();
}

void Matrix::to_u()
{
    for(int c=0;c<m_col;c++)
    {
        double quotient = 1;
        if(matrix_p[c][c] == 0)
        {
            int r2;
            for(r2=c+1;r2<m_row;r2++)
            {
                if(matrix_p[r2][c] != 0)
                    break;
            }
            if(r2 != c)
            {
                for(int i=c;i<m_col;i++)   //because all before c is 0
                {
                    matrix_p[c][i] += matrix_p[r2][i];
                    if(abs(matrix_p[c][i])<ACCURACY)
                        matrix_p[c][i]=0;
                }
            }
        }
//
        for(int r=c+1;r<m_row;r++)
        {
            if(matrix_p[r][c] == 0) continue;
            quotient = matrix_p[r][c]/matrix_p[c][c];
            for(int i=c;i<m_col;i++)
            {
                matrix_p[r][i] -= matrix_p[c][i]*quotient;
                if(abs(matrix_p[r][i])<ACCURACY)
                    matrix_p    [r][i]=0;
            }
        }
        if(c+1>=m_row)
            break;
    }
}

void Matrix::to_rref()
{
    double quotient = 1;
    for(int c=0;c<m_col;c++)
    {
        if(matrix_p[c][c] == 0)
        {
            int r2;
            for(r2=c+1;r2<m_row;r2++)
            {
                if(matrix_p[r2][c] != 0)
                    break;
            }
            if(r2 != c)
            {
                for(int i=c;i<m_col;i++)   //because all before c is 0
                {
                    matrix_p[c][i] += matrix_p[r2][i];
                    if(abs(matrix_p[c][i])<ACCURACY)
                        matrix_p[c][i]=0;
                }
            }
        }
//
        for(int r=c+1;r<m_row;r++)
        {
            if(matrix_p[r][c] == 0) continue;
            quotient = matrix_p[r][c]/matrix_p[c][c];
            for(int i=c;i<m_col;i++)
            {
                matrix_p[r][i] -= matrix_p[c][i]*quotient;
                if(abs(matrix_p[r][i])<ACCURACY)
                    matrix_p[r][i]=0;
            }
        }

        if(c+1>=m_row)
            break;
    }
//rref:
    for(int r=0;r<m_row;r++)
    {
        for(int c=r;c<m_col;c++)
        {
            if(matrix_p[r][c] == 0) continue;

            for(int r2=0;r2<r;r2++)
            {
                quotient = matrix_p[r2][c]/matrix_p[r][c];
                for(int i=c;i<m_col;i++)
                {
                    matrix_p[r2][i] -= matrix_p[r][i]*quotient;
                    if(abs(matrix_p[r2][i])<ACCURACY)
                        matrix_p[r2][i]=0;
                }
            }
            break;
        }
    }
}

void Matrix::to_augmented(const Matrix &other)
{
    if(m_row == other.m_row)
    {
        Matrix newMatrix(m_row,m_col+other.m_col);
        for(int r=0;r<m_row;r++)
        {
            for(int c=0;c<m_col+other.m_col;c++)
            {
                if(c<m_col)
                {
                    newMatrix[r][c] = matrix_p[r][c];
                }
                else
                {
                    newMatrix[r][c] = other.matrix_p[r][c-m_col];
                }
            }
        }
        free_matrix();
        *this = newMatrix;
    }
    else
    {
        error(5);
    }
}

void Matrix::to_invert()
{
    if(m_row==m_col && rank()==m_row)     //inefficient
    {
        if(m_row == 1)
        {
            matrix_p[0][0] = 1/matrix_p[0][0];
            return;
        }
        Matrix i("I",m_row);
        Matrix arg = get_augmented(i);
        arg.to_rref();

        for(int r=0;r<m_row;r++)
        {
            for(int c=0;c<m_col;c++)
            {
                matrix_p[r][c] = arg.matrix_p[r][c+m_col];
            }
        }
    }
    else
    {
        error(7);
    }
}

/**************************************/

Matrix Matrix::get_u() const
{
    Matrix newMatrix(*this);

    for(int c=0;c<m_col;c++)
    {
        double quotient = 1;
        if(newMatrix[c][c] == 0)
        {
            int r2;
            for(r2=c+1;r2<m_row;r2++)
            {
                if(newMatrix[r2][c] != 0)
                    break;
            }
            if(r2 != c)
            {
                for(int i=c;i<m_col;i++)   //because all before c is 0
                {
                    newMatrix[c][i] += newMatrix[r2][i];
                    if(abs(newMatrix[c][i])<ACCURACY)
                        newMatrix[c][i]=0;
                }
            }
        }
//
        for(int r=c+1;r<m_row;r++)
        {
            if(newMatrix[r][c] == 0) continue;
            quotient = newMatrix[r][c]/newMatrix[c][c];
            for(int i=c;i<m_col;i++)
            {
                newMatrix[r][i] -= newMatrix[c][i]*quotient;
                if(abs(newMatrix[r][i])<ACCURACY)
                    newMatrix[r][i]=0;
            }
        }

        if(c+1>=m_row)
            break;
    }
    return newMatrix;
}

Matrix Matrix::get_e() const
{
    Matrix i("I",m_row);
    Matrix arg = get_augmented(i);
    arg = arg.get_u();

    Matrix newMatrix(m_row,m_row);
    for(int r=0;r<m_row;r++)
    {
        for(int c=0;c<m_row;c++)
        {
            newMatrix[r][c] = arg.matrix_p[r][c+m_col];
        }
    }
    return newMatrix;
}

Matrix Matrix::get_l() const
{
    return get_e().get_invert();
}

/************************Determinant************************/

double Matrix::get_det() const
{
    if(m_row == m_col)
    {
        Matrix iter(*this);
        iter.to_u();
        double det = 1;
        for(int i=0;i<m_col;i++)
        {
            det *= iter[i][i];
        }
        return det;
    }
    else
    {
        return 0;
    }

}

/**********************Projection**********************/

Matrix Matrix::get_p() const
{
    if(rank() == m_col)
    {
        return (*this) * ((*this).get_transpose() * (*this)).get_invert() * (*this).get_transpose();
    }
    else
    {
        error(8);
    }
}

Matrix get_projection(const Matrix &a,const Matrix &b)
{
    if(b.m_row == a.m_row)
    {
        return b.get_p()*a;
    }
    else
    {
        error(8);
    }
}

/*********************Eigenvalues************************/

double Matrix::trace() const
{
    if(m_row == m_col)
    {
        double trace = 0;
        for(int i=0;i<m_row;i++)
        {
            trace += matrix_p[i][i];
        }
        return trace;
    }
    else
    {
        error(5);
    }
}

/*********************Operators***********************/

double* &Matrix::operator[](int row)
{
    //cout<<".<"<<row<<">."<<endl;
    if(row>=0 && row<m_row)
        return matrix_p[row];
    else
    {
        error(6);
    }
}

Matrix &Matrix::operator=(const Matrix &other)
{
    free_matrix();
    m_row = other.m_row;
    m_col = other.m_col;

    matrix_p = new double*[m_row];
    for(int r=0;r<m_row;r++)
    {
        matrix_p[r] = new double[m_col];
        for(int c=0;c<m_col;c++)
        {
            matrix_p[r][c] = other.matrix_p[r][c];
        }
    }
}

Matrix &Matrix::operator+=(const Matrix &other)
{
    *this = (*this)+other;
}

Matrix &Matrix::operator-=(const Matrix &other)
{
    *this = (*this)-other;
}

Matrix &Matrix::operator*=(const Matrix &other)
{
    *this = (*this)*other;
}

Matrix &Matrix::operator*=(const double &num)
{
    *this = (*this)*num;
}

Matrix &Matrix::operator/=(const double &num)
{
    *this = (*this)/num;
}

Matrix operator+(const Matrix &m1,const Matrix &m2)
{
    if(m1.m_row==m2.m_row && m1.m_col==m2.m_col)
    {
        Matrix newMatrix;
        newMatrix.m_row = m1.m_row;
        newMatrix.m_col = m1.m_col;

        newMatrix.matrix_p = new double*[m1.m_row];
        for(int r=0;r<m1.m_row;r++)
        {
            newMatrix.matrix_p[r] = new double[m1.m_col];
            for(int c=0;c<m1.m_col;c++)
            {
                newMatrix.matrix_p[r][c] = m1.matrix_p[r][c] + m2.matrix_p[r][c];
            }
        }
        return newMatrix;
    }
    else
    {
        error(5);
    }
}

Matrix operator-(const Matrix &m1,const Matrix &m2)
{
    if(m1.m_row==m2.m_row && m1.m_col==m2.m_col)
    {
        Matrix newMatrix;
        newMatrix.m_row = m1.m_row;
        newMatrix.m_col = m1.m_col;

        newMatrix.matrix_p = new double*[m1.m_row];
        for(int r=0;r<m1.m_row;r++)
        {
            newMatrix.matrix_p[r] = new double[m1.m_col];
            for(int c=0;c<m1.m_col;c++)
            {
                newMatrix.matrix_p[r][c] = m1.matrix_p[r][c] - m2.matrix_p[r][c];
            }
        }
        return newMatrix;
    }
    else
    {
        error(5);
    }
}

Matrix operator*(double num,const Matrix &m)
{
    Matrix newMatrix;
    newMatrix.m_row = m.m_row;
    newMatrix.m_col = m.m_col;

    newMatrix.matrix_p = new double*[m.m_row];
    for(int r=0;r<m.m_row;r++)
    {
        newMatrix.matrix_p[r] = new double[m.m_col];
        for(int c=0;c<m.m_col;c++)
        {
            newMatrix.matrix_p[r][c] = num*m.matrix_p[r][c];
        }
    }
    return newMatrix;
}

Matrix operator*(const Matrix &m,double num)
{
    Matrix newMatrix;
    newMatrix.m_row = m.m_row;
    newMatrix.m_col = m.m_col;

    newMatrix.matrix_p = new double*[m.m_row];
    for(int r=0;r<m.m_row;r++)
    {
        newMatrix.matrix_p[r] = new double[m.m_col];
        for(int c=0;c<m.m_col;c++)
        {
            newMatrix.matrix_p[r][c] = num*m.matrix_p[r][c];
        }
    }
    return newMatrix;
}

Matrix operator*(const Matrix &m1,const Matrix &m2)
{
    if(m1.m_col == m2.m_row)
    {
        Matrix newMatrix(m1.m_row,m2.m_col);
        for(int r=0;r<m1.m_row;r++)
        {
            for(int c=0;c<m2.m_col;c++)
            {
                for(int i=0;i<m1.m_col;i++)
                {
                    newMatrix[r][c] += m1.matrix_p[r][i]*m2.matrix_p[i][c];
                }
            }
        }
        return newMatrix;
    }
    else
    {
        error(5);
    }
}

Matrix operator/(const Matrix &m,double num)
{
    Matrix newMatrix;
    newMatrix.m_row = m.m_row;
    newMatrix.m_col = m.m_col;

    newMatrix.matrix_p = new double*[m.m_row];
    for(int r=0;r<m.m_row;r++)
    {
        newMatrix.matrix_p[r] = new double[m.m_col];
        for(int c=0;c<m.m_col;c++)
        {
            newMatrix.matrix_p[r][c] = m.matrix_p[r][c]/num;
        }
    }
    return newMatrix;
}

int main()
{
    //double e[4] = {1,2,3,4};
    string e = "{{1,2}}";
    Matrix a(e);
    cout<<"a:"<<a.rank()<<endl;
    a = a.get_transpose();
    a.print();
    Matrix p = a.get_p();
    cout<<"p"<<endl;
    p.print();
    cout<<"a_p"<<endl;
    a = get_projection(a,a);
    a.print();

    string r = "{{2,7},{6,5}}";
    Matrix b(r);
    b.to_transpose();
    cout<<"b:"<<b.trace()<<endl;
    b.print();

    double d[4] = {1,2,3,4};
    Matrix c(2,2,d);
    c.print();







    return 0;
}
