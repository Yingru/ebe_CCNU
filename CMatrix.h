#ifndef CMATRIXH
#define CMATRIXH
#include<iostream>
#include<cassert>
using namespace std;
/**********Matrix_1D****************************/
template <class mytype>
class CMatrix_1D{
    private:
	mytype *m_data;   //storage data
	void copy(const CMatrix_1D <mytype> & mat1);
    public:
	int ncols;  //total rows and columns
	int col_l, col_h;  //column index

	void init(int,int);
	CMatrix_1D();      //construct a 1*1 matrix
	~CMatrix_1D();
	CMatrix_1D(const CMatrix_1D <mytype> &mat1);
	CMatrix_1D(int n); //construct a 1*n matrix
	CMatrix_1D(int nl, int nh);//construct a (nh-nl+1) matrix
	mytype& operator()(int); 
//	CMatrix_1D <mytype> & operator=(const CMatrix_1D <mytype> &mat1);
};


template < class mytype >
CMatrix_1D<mytype>::~CMatrix_1D(){
    delete [] m_data;
}

template < class mytype >
CMatrix_1D<mytype>::CMatrix_1D(const CMatrix_1D <mytype> &mat1){
this->copy(mat1);
}

template < class mytype >
void CMatrix_1D<mytype>::init(int nl, int nh){
    col_l = nl;
    col_h = nh;
    ncols = nh-nl+1 ;

    m_data = new mytype[ncols];
    assert(m_data != 0);//check allocate 
    for(int i=0; i!=ncols; i++)
	m_data[i]=0; 
}

template < class mytype >
CMatrix_1D <mytype>::CMatrix_1D(){
    init(0,0);
}

template < class mytype >
CMatrix_1D <mytype>::CMatrix_1D(int nl, int nh){
    init(nl,nh);
}



template < class mytype >
mytype& CMatrix_1D<mytype>::operator ()(int j){
    assert(j>=col_l && j<=col_h);
    return m_data[(j-col_l)];
}

template < class mytype >
void CMatrix_1D<mytype>::copy(const CMatrix_1D <mytype> &mat1){
    col_l=mat1.col_l;
    col_h=mat1.col_h;
    ncols=col_h-col_l;
    int idata=col_h-col_l+1;
    m_data = new mytype [idata];
    for(int i=0; i!=idata; i++)
	m_data[i]=mat1.m_data[i];
}

//template < class mytype >
//CMatrix_1D <mytype> & CMatrix_1D::operator=(const CMatrix_1D <mytype> & mat1){
//    if (this == &mat1) return *this;
//    delete [] m_data;
//    this->copy(mat1);
//    return *this;
//}
//

/**********Matrix_1D****************************/


/**********Matrix_2D****************************/
template <class mytype>
class CMatrix_2D{
    private:
	mytype *m_data;   //storage data
    public:

	int nrows, ncols;  //total rows and columns
	int row_l, row_h;  //row index low-high
	int col_l, col_h;  //column index


	void init(int,int,int,int);
	CMatrix_2D();      //construct a 1*1 matrix
	~CMatrix_2D();      //destruct a 1*1 matrix
	CMatrix_2D(int m, int n); //construct a m*n matrix
	CMatrix_2D(int ml, int mh, int nl, int nh);//construct a (mh-ml+1)*(nh-nl+1) matrix
	mytype& operator()(int, int); 
};

template < class mytype >
CMatrix_2D<mytype>::~CMatrix_2D(){
    delete [] m_data;
}


template < class mytype >
void CMatrix_2D<mytype>::init(int ml, int mh, int nl, int nh){
    row_l = ml;
    row_h = mh;
    col_l = nl;
    col_h = nh;
    nrows = mh-ml+1 ;
    ncols = nh-nl+1 ;

    m_data = new mytype[nrows*ncols];
    for(int i=0; i!=nrows*ncols; i++)
	m_data[i]=0; //Here is not safe
}

template < class mytype >
CMatrix_2D <mytype>::CMatrix_2D(){
    init(0,0,0,0);
}

template < class mytype >
CMatrix_2D <mytype>::CMatrix_2D(int m, int n){
    init(0,m,0,n);
}

template < class mytype >
CMatrix_2D <mytype>::CMatrix_2D(int ml, int mh, int nl, int nh){
    init(ml,mh,nl,nh);
}

template < class mytype >
mytype& CMatrix_2D<mytype>::operator ()(int i, int j){
    if( (i<row_l)or(i>row_h)or(j<col_l)or(j>col_h) ){
	cout<<"Out of matrix range!"<<endl;
    }
    else {
	return m_data[(i-row_l)*ncols+(j-col_l)];
    }
}
/**********Matrix_2D****************************/


/**********Matrix_3D****************************/
template <class mytype>
class CMatrix_3D{
    private:
	mytype *m_data;   //storage data

    public:

	int npags, nrows, ncols;  //total rows and columns
	int pag_l, pag_h;
	int row_l, row_h;  //row index low-high
	int col_l, col_h;  //column index
	void init(int,int,int,int,int,int);

	CMatrix_3D();      //construct a 1*1 matrix
	~CMatrix_3D();      //destruct a 1*1 matrix
	CMatrix_3D(int p, int m, int n); //construct a m*n matrix
	CMatrix_3D(int pl, int ph, int ml, int mh, int nl, int nh);//construct a (mh-ml+1)*(nh-nl+1) matrix
	mytype& operator()(int, int, int); 
};

template < class mytype >
CMatrix_3D <mytype>::~CMatrix_3D(){
    delete [] m_data;
}


template < class mytype >
void CMatrix_3D<mytype>::init(int pl, int ph, int ml, int mh, int nl, int nh){
    pag_l = pl;
    pag_h = ph;
    row_l = ml;
    row_h = mh;
    col_l = nl;
    col_h = nh;
    npags = ph-pl+1 ;
    nrows = mh-ml+1 ;
    ncols = nh-nl+1 ;

    m_data = new mytype[npags*nrows*ncols];
    for(int i=0; i!=npags*nrows*ncols; i++)
	m_data[i]=0; //Here is not safe
}

template < class mytype >
CMatrix_3D <mytype>::CMatrix_3D(){
    init(0,0,0,0,0,0);
}

template < class mytype >
CMatrix_3D <mytype>::CMatrix_3D(int p, int m, int n){
    init(0,p,0,m,0,n);
}

template < class mytype >
CMatrix_3D <mytype>::CMatrix_3D(int pl, int ph, int ml, int mh, int nl, int nh){
    init(pl,ph,ml,mh,nl,nh);
}

template < class mytype >
mytype& CMatrix_3D<mytype>::operator ()(int k, int i, int j){
    //if( (k<pag_l)or(k>pag_h)or(i<row_l)or(i>row_h)or(j<col_l)or(j>col_h) ){
	//cout<<"Out of matrix range!"<<endl;
	//assert( (k<pag_l)or(k>pag_h)or(i<row_l)or(i>row_h)or(j<col_l)or(j>col_h) );
    //}
	return m_data[(k-pag_l)*nrows*ncols+(i-row_l)*ncols+(j-col_l)];
}
/**********Matrix_3D****************************/
#endif
