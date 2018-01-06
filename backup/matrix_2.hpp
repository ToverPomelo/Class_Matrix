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
        Matrix(const Matrix &other);

        void scan();
        void print() const;
        int row() const;
        int column() const;
        Matrix transpose() const;
        Matrix to_u() const;
        Matrix to_rref() const;

        double* &operator[](int row);

        friend Matrix operator+(const Matrix &m1,const Matrix &m2);
        friend Matrix operator-(const Matrix &m1,const Matrix &m2);
        friend Matrix operator*(double num,const Matrix &m);
        friend Matrix operator*(const Matrix &m,double num);
        friend Matrix operator*(const Matrix &m1,const Matrix &m2);
        friend Matrix operator/(const Matrix &m,double num);



      //  ~Matrix();
};
