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
        default:break;
    }
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

Matrix Matrix::transpose() const
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

Matrix Matrix::to_u() const
{
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
    return newMatrix;
}

Matrix Matrix::to_rref() const
{
    Matrix newMatrix(*this);
    int pivot = 0;
//to_u
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
//to_rref
    int count_p = 0;
    for(int r=0;r<m_row;r++)
    {
        int mark_pivot;
        double quotient = 1;
        bool zero_row = 1;
        for(int c=count_p;c<m_col;c++)     //find pivot line
        {
            if(newMatrix[r][c] != 0)
            {
                mark_pivot = c;
                count_p++;
                zero_row = 0;
                break;
            }
        }
        if(zero_row)   break;            //meet zero row and break;
        else
        {
            quotient = newMatrix[r][mark_pivot];       //set pivot to 1
            for(int c=mark_pivot;c<m_col;c++)
            {
                newMatrix[r][c] /= quotient;
                if(abs(newMatrix[r][c])<ACCURACY)
                    newMatrix[r][c] = 0;
            }
        }
        for(int _r=0;_r<r;_r++)
        {
            if(newMatrix[_r][mark_pivot]==0)  continue;
            quotient = newMatrix[_r][mark_pivot]/newMatrix[r][mark_pivot];
            for(int _c=mark_pivot;_c<m_col;_c++)
            {
                newMatrix[_r][_c] -= newMatrix[r][_c]*quotient;
                if(abs(newMatrix[_r][_c])<ACCURACY)
                    newMatrix[_r][_c] = 0;

            }
        }
        if(count_p>=pivot)   break;
    }

    return newMatrix;
}

double* &Matrix::operator[](int row)
{
    if(row>=0 && row<m_row)
        return matrix_p[row];
    else
    {
        error(6);
    }
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
    string e = "{{1,2,3},{4,5,6}}";

    Matrix a(e);
    cout<<"a:"<<endl;
    a.print();

    e = "{{0},{1},{1}}";
    Matrix b(e);
    cout<<"b:"<<endl;
    b.print();

    Matrix c = a.to_rref();
    cout<<"c:"<<endl;
    c.print();



    return 0;
}
