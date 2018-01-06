/***********************************************/
/*          It's  full  of  water!             */
/*              都 是 水 的！                  */
/***********************************************/


#include <iostream>
#include <limits>
#include <cmath>
#define MAX 100
#define ACCURACY 0.00000000001
using namespace std;

int row,col;
double matrix[MAX][MAX];
int pivot = 0; 

void swap(int &a,int &b)
{
    a = a + b;
    b = a - b;
    a = a - b;
}


int main()
{
//input
    cout<<"Row==";
    cin>>row;
    cout<<"Column==";
    cin>>col;
    cout<<"Go on:"<<endl;
    for(int r=0;r<row;r++)
    {
        for(int c=0;c<col;c++)
            cin>>matrix[r][c];
    }
//to U
    for(int c=0;c<col;c++)
    {
        int mark_r = -1;
        double mark_num = (numeric_limits<double>::max)();
        double quotient = 1;
        bool isn_zero[row] = {0};
        for(int r=pivot;r<row;r++)
        {
            if(matrix[r][c]!=0)
            {
                isn_zero[r] = 1;
            }
            if(matrix[r][c]<mark_num && matrix[r][c]!=0)
            {
                mark_num = matrix[r][c];
                mark_r = r;
            }
        }
        if(mark_r != -1) pivot++;
        if(pivot>=row || pivot>=col) break;       
        if(mark_r!=0 && mark_r!=-1) //change to the first line
        {
            for(int i=pivot-1;i<col;i++)
            {
                swap(matrix[pivot-1][i],matrix[mark_r][i]);
            }
        }

        if(mark_r != -1)           //multiply and add
        {
            for(int r=pivot;r<row;r++)
            {
                if(!isn_zero[r]) continue;
                quotient = matrix[r][c]/matrix[pivot-1][c];
                for(int i=pivot-1;i<col;i++)
                {
                    matrix[r][i] -= matrix[pivot-1][i]*quotient;
                    if(abs(matrix[r][i])<ACCURACY)
                        matrix[r][i]=0;
                }
            }
        }  
    }
    cout<<"Rank=="<<pivot<<endl;
    cout<<"U:"<<endl;
    for(int r=0;r<row;r++)
    {
        for(int c=0;c<col;c++)
            cout<<matrix[r][c]<<"\t";
        cout<<endl;
    }
//to rref

    int count_p = 0;
    for(int r=0;r<row;r++)
    {
        int mark_pivot;
        double quotient = 1;
        bool zero_row = 1;
        for(int c=count_p;c<col;c++)     //find pivot line
        {
            if(matrix[r][c] != 0)
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
            quotient = matrix[r][mark_pivot];       //set pivot to 1
            for(int c=mark_pivot;c<col;c++)
            {
                matrix[r][c] /= quotient;
                if(abs(matrix[r][c])<ACCURACY)
                    matrix[r][c] = 0;
            }
        }
        for(int _r=0;_r<r;_r++)
        {
            if(matrix[_r][mark_pivot]==0)  continue;
            quotient = matrix[_r][mark_pivot]/matrix[r][mark_pivot];
            for(int _c=mark_pivot;_c<col;_c++)
            {
                matrix[_r][_c] -= matrix[r][_c]*quotient;
                if(abs(matrix[_r][_c])<ACCURACY)
                    matrix[_r][_c] = 0;

            }
        }
        if(count_p>=pivot)   break; 
    }
    
    cout<<"rref:"<<endl;
    for(int r=0;r<row;r++)
    {
        for(int c=0;c<col;c++)
            cout<<matrix[r][c]<<"\t";
        cout<<endl;
    }

    return 0;
}
