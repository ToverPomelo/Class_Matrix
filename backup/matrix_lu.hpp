#define MAXROW 100
#define MAXCOL 100

class Matrix
{
    private:
        int m_row;
        int m_col;
        double **matrix_p;
        void free_matrix();
        //Matrix s_to_m(string psd_m);

    public:
        Matrix();
        Matrix(int row,int col);
        Matrix(int row,int col,double *input_matrix); //not safe
        Matrix(string psd_matrix);
        Matrix(string mark_s,int mark_i);
        Matrix(const Matrix &other);
        ~Matrix();

        void scan();
        void print() const;
        int row() const;
        int column() const;
        int rank() const;
//basis
        Matrix get_transpose() const;
        Matrix eliminate_get_u() const;
        Matrix get_rref() const;
        Matrix get_augmented(const Matrix &other) const;
        Matrix get_invert() const;
        void to_transpose();
        void eliminate_to_u();    //change lines?
        void to_rref();
        void to_augmented(const Matrix &other);
        void to_invert();
//lu factorization
        Matrix get_u() const;
        Matrix get_e() const;
        Matrix get_l() const;

        double* &operator[](int row);
        Matrix &operator=(const Matrix &other);

        friend Matrix operator+(const Matrix &m1,const Matrix &m2);
        friend Matrix operator-(const Matrix &m1,const Matrix &m2);
        friend Matrix operator*(double num,const Matrix &m);
        friend Matrix operator*(const Matrix &m,double num);
        friend Matrix operator*(const Matrix &m1,const Matrix &m2);
        friend Matrix operator/(const Matrix &m,double num);

};
