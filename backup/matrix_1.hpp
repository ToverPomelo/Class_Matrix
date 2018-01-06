#define MAXROW 100
#define MAXCOL 100

class Matrix
{
    private:
        int m_row;
        int m_col;
        double **matrix_p;
        //Matrix s_to_m(string psd_m);

    public:
        Matrix();
        Matrix(int row,int col);
        Matrix(int row,int col,double *input_matrix); //not safe
        Matrix(string psd_matrix);
        void scan();
        void print() const;
        Matrix plus(Matrix a,Matrix b) const;

        Matrix &operator=(const Matrix &newmat);



      //  ~Matrix();
};
