#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#include "matrix.hpp"

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
        default:break;
    }
}

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
/*
Matrix Matrix::plus(Matrix a,Matrix b) const
{
    if(a.m_row==b.m_row && a.col==b.col)
    {
        Matrix iter(a.row,a.col);
        for(int r=0;r<a.row;r++)
        {
            for(int c=0;c<a.col;c++)
            {
                iter.matrix[r][c] = a.matrix[r][c] + b.matrix[r][c];
                cout<<iter.matrix[r][c]<<endl;
            }
        }
        return iter;
    }
    else
    {
        error(1);
    }
}

Matrix &Matrix::operator=(const Matrix &newmat)
{
    row = newmat.row;
    col = newmat.col;
    for(int r=0;r<newmat.row;r++)
    {
        for(int c=0;c<newmat.col;c++)
            matrix[r][c] = newmat.matrix[r][c];
    }
    return *this;
}
*/
int main()
{
    //double e[4] = {1,2,3,4};
    string e = "{{2,4,5},{125,4.6,7.5}}";

    Matrix a(e);
    cout<<"a:"<<endl;
//    a.scan();
    a.print();



    return 0;
}
