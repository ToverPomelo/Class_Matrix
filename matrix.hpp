#define MAXROW 100
#define MAXCOL 100

class Matrix
{
    private: 
        int m_row;                         //行数
        int m_col;                         //列数
        double **matrix_p;                 //二维数组，用来放矩阵
        void free_matrix();                //删除矩阵
        //Matrix s_to_m(string psd_m);     //这个。。忘了是用来干什么的。。。

    public:
        Matrix();      
        Matrix(int row,int col);            //指定行数和列数时初始化为0矩阵
        Matrix(int row,int col,double *input_matrix); //not safe?
        Matrix(string psd_matrix);          //以{{*,...,*},...,{*,...,*}}格式的字符串传递矩阵
        Matrix(string mark_s,int mark_i);   //这里构造一些特殊矩阵，如单位矩阵
        Matrix(const Matrix &other);        //复制
        ~Matrix();

        void scan();                        //通过输入给矩阵赋值
        void print() const;                 //打印矩阵
        int row() const;                    //返回矩阵的行数
        int column() const;                 //返回矩阵的列数
        int rank() const;                   //返回矩阵的秩
                                            //（rank这里偷懒直接用了以前写过消元的代码，应该有bug）
//basis
        Matrix get_transpose() const;       //返回矩阵的转置矩阵
        Matrix get_rref() const;            //返回矩阵的rref（简化阶梯形）形式
        Matrix get_augmented(const Matrix &other) const;   //返回该矩阵与矩阵other增广矩阵
        Matrix get_invert() const;          //返回矩阵的逆矩阵
        void to_transpose();                //变为转置矩阵
        void to_u();    //change lines?     //变为上三角（阶梯形）
        void to_rref();                     //变rref
        void to_augmented(const Matrix &other);    //变增广
        void to_invert();                   //变逆矩阵
//lu factorization（LU分解）
        Matrix get_u() const;               //得到U               
        Matrix get_e() const;               //得到E（L的逆）
        Matrix get_l() const;               //得到L
//determinant（行列式）
        double get_det() const;             //返回矩阵行列式
//projection（投影）
        Matrix get_p() const;               //返回P（投影矩阵）
        friend Matrix get_projection(const Matrix &a,const Matrix &b);  //from a to b
									//从a到b的投影
//Eigen~（特征值与特征向量，未完成，据说可以用QR分解求）
        double trace() const;               //返回矩阵的迹

        double* &operator[](int row);                //重载[]，让矩阵可以像二维数组那样用
        Matrix &operator=(const Matrix &other);      //重载复制
        Matrix &operator+=(const Matrix &other);     
        Matrix &operator-=(const Matrix &other);
        Matrix &operator*=(const Matrix &other);
        Matrix &operator*=(const double &num);       //这几个偷懒才这样写的啊
        Matrix &operator/=(const double &num);       //如果矩阵很大的话效率好像会跟不上

        friend Matrix operator+(const Matrix &m1,const Matrix &m2);   //矩阵相加
        friend Matrix operator-(const Matrix &m1,const Matrix &m2);   //相减
        friend Matrix operator*(double num,const Matrix &m);          //数乘
        friend Matrix operator*(const Matrix &m,double num);          //数乘
        friend Matrix operator*(const Matrix &m1,const Matrix &m2);   //矩阵相乘
        friend Matrix operator/(const Matrix &m,double num);          //相除

};
