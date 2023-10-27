#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo {

// forward declarations

// accumulator
template <typename T=real_t>
class gsCoefAccumulator;
template <typename T=real_t>
class gsRowAccumulator;
template <typename T=real_t>
class gsColAccumulator ;
template <typename T=real_t>
class gsRecipeAccumulator;

// per point data
template <typename T=real_t>
class gsPointFuncData;
template <typename T=real_t>
class gsPointMapData;

// bilinear ops
template <typename T=real_t>
class gsBilinearOp;
    // for volumes
template <typename T=real_t>
class gsL2ScalarOp;
template <typename T=real_t>
class gsGradGradOp;
template <typename T=real_t>
class gsDivergenceOp;
template <typename T=real_t>
class gsGradientOp;
    // for boundaries
template <typename T=real_t>
class gsBoundaryL2ScalarOp;
template <typename T=real_t>
class gsBoundaryNormalDerValueOp;
template <typename T=real_t>
class gsBoundaryNormalDerNormalDerOp;

    // test -> DEPRECATED
template <typename T=real_t>
class gsL2TestOp;
template <typename T=real_t>
class gsL2TestVecOp;
template <typename T=real_t>
class gsBoundaryL2ScalarTestVecOp;
template <typename T=real_t>
class gsBoundaryL2ScalarTestOp;
template <typename T=real_t>
class gsBoundaryL2TestOp;
template <typename T=real_t>
class gsBoundaryL2TestVecOp;


// real code
template <typename T>
class gsRecipeAccumulator
{
protected:
    gsMatrix<T> &m_data;
    T            m_fact;
public:
    gsRecipeAccumulator(gsMatrix<T> &d, T f=1)
        : m_data(d),m_fact(f)
    {}
    template <typename exprT>
    void operator+= (const exprT ex) {m_data+=ex*m_fact;}
    template <typename exprT>
    void operator-= (const exprT ex) {m_data-=ex*m_fact;}

    const gsMatrix<T> & data() {return m_data;}
    
    gsCoefAccumulator<T> operator() (index_t r,index_t c) {return gsCoefAccumulator<T>(*this,r,c);}
    gsRowAccumulator<T> row (index_t r) {return gsRowAccumulator<T>(*this,r);}
    gsColAccumulator<T> col (index_t c) {return gsColAccumulator<T>(*this,c);}
};


template <typename T>
class gsColAccumulator : public gsRecipeAccumulator<T>
{
private:
    using gsRecipeAccumulator<T>::m_fact;
    using gsRecipeAccumulator<T>::m_data;

    index_t m_col;
public:
    gsColAccumulator(gsRecipeAccumulator<T> o, index_t c)
        : gsRecipeAccumulator<T>(o), m_col(c)
    {}
    template <typename exprT>
    void operator+=(const exprT ex) {m_data.col(m_col)+=ex*m_fact;}
    template <typename exprT>
    void operator-=(const exprT ex) {m_data.col(m_col)-=ex*m_fact;}
};

template <typename T>
class gsRowAccumulator : public gsRecipeAccumulator<T>
{
private:
    using gsRecipeAccumulator<T>::m_fact;
    using gsRecipeAccumulator<T>::m_data;

    index_t m_row;
public:
    gsRowAccumulator(gsRecipeAccumulator<T> o, index_t r)
        : gsRecipeAccumulator<T>(o), m_row(r)
    {}
    template <typename exprT>
    void operator+=(const exprT ex) {m_data.row(m_row)+=ex*m_fact;}
    template <typename exprT>
    void operator-=(const exprT ex) {m_data.row(m_row)-=ex*m_fact;}
};

template <typename T>
class gsCoefAccumulator : public gsRecipeAccumulator<T>
{
private:
    using gsRecipeAccumulator<T>::m_fact;
    using gsRecipeAccumulator<T>::m_data;

    index_t m_row;
    index_t m_col;
public:
    gsCoefAccumulator(gsRecipeAccumulator<T> o, index_t r, index_t c)
        : gsRecipeAccumulator<T>(o), m_row(r), m_col(c)
    {}
    template <typename exprT>
    void operator+=(const exprT ex) {m_data(m_row,m_col)+=ex*m_fact;}
    template <typename exprT>
    void operator-=(const exprT ex) {m_data(m_row,m_col)-=ex*m_fact;}
};



template <typename T>
class gsPointFuncData
{
public:
    typedef typename gsFuncData<T>::matrixView matrixView;
    typedef typename gsFuncData<T>::matrixTransposeView matrixTransposeView;
private:
    const gsFuncData<T>       *m_data;
    index_t              m_k;
public:
    gsPointFuncData(const gsFuncData<T> &data, index_t point_idx=0)
        : m_data(&data), m_k(point_idx)
    {}

    static gsPointFuncData<T> make(const gsFuncData<T> &data)
    {return gsPointFuncData<T>(data);}

    inline void setPoint (index_t p) {m_k=p;}

    const gsFuncData<T> & data() const {return *m_data;}
    const std::pair<short_t, short_t>    &info      () const {return m_data->dim;}
    inline unsigned                   flags     () const {return m_data->flags;}

    inline const gsMatrix<index_t>   &actives   () const {return m_data->actives;}
    inline const matrixView           value     () const {return m_data->eval(m_k);}
    inline const matrixView           deriv     () const {return m_data->deriv(m_k);}
    inline const matrixView           deriv2    () const {return m_data->deriv2(m_k);}
    inline const matrixView           div       () const {return m_data->div(m_k);}
    inline const matrixView           curl      () const {return m_data->curl(m_k);}
    inline const matrixView           laplacian () const {return m_data->laplacian(m_k);}
    inline const matrixTransposeView  jacobians () const {return m_data->values[1].reshapeCol(m_k, m_data->dim.first, m_data->values[1].rows()/m_data->dim.first).transpose(); }
    inline const matrixTransposeView  jacobian  (index_t fId) const {return m_data->jacobian(m_k,fId);}
};

/// This is needed because we decided to split gsMapData from gsFunData.
template <typename T>
class gsPointMapData
{
public:
    typedef typename gsMapData<T>::matrixView matrixView;
    typedef typename gsMapData<T>::matrixTransposeView matrixTransposeView;
    typedef typename gsMapData<T>::constColumn constColumn;
private:
    const gsMapData<T>        *m_data;
    index_t              m_k;
public:
    gsPointMapData(const gsMapData<T> &data, index_t point_idx=0)
        : m_data(&data), m_k(point_idx)
    {}

    static gsPointMapData<T> make(const gsFuncData<T> &data)
    {return gsPointMapData<T>(data);}

    inline void setPoint (index_t p) {m_k=p;}

    const std::pair<short_t, short_t>    &info      () const {return m_data->dim;}
    inline unsigned                   flags     () const {return m_data->flags;}

    inline const gsMatrix<index_t>   &actives   () const {return m_data->actives;}
    inline const matrixView           value     () const {return m_data->eval(m_k);}
    inline const matrixView           deriv     () const {return m_data->deriv(m_k);}
    inline const matrixView           deriv2    () const {return m_data->deriv2(m_k);}
    inline const matrixView           div       () const {return m_data->div(m_k);}
    inline const matrixView           curl      () const {return m_data->curl(m_k);}
    inline const matrixView           laplacian () const {return m_data->laplacian(m_k);}
    inline const matrixTransposeView  jacobian  () const {return m_data->jacobian(m_k,0);}

    inline constColumn          point     () const { return m_data->point(m_k);}
    inline T                    measure   () const { return m_data->measure(m_k);}
    inline matrixView           fundForm  () const { return m_data->fundForm(m_k);}
    inline constColumn          normal    () const { return m_data->normal(m_k);}
    inline constColumn          outNormal () const { return m_data->outNormal(m_k);}
};


/**
  The bilinear operator class is an abstract interface for operators
  in the recipe assembler framework.
  Each operator must implement the full interface.
**/
template <typename T>
class gsBilinearOp
{
public:
    virtual ~gsBilinearOp()
    {}
    /**
      Evaluate the operator between the spaces
    **/
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>     result
            ) const = 0;
    virtual unsigned    testSpaceNeeds()    const  = 0;
    virtual unsigned    unknownSpaceNeeds() const  = 0;
    virtual unsigned    geometryNeeds()     const  {return 0;}
    virtual void        outputSize (unsigned &r, unsigned &c) const { GISMO_UNUSED(r); GISMO_UNUSED(c); }
};


//template <typename T>
//class gsInterfaceOperator : public gsBilinearOp<T>
//{
//public:
//    virtual gsMatrix<T>  pointEval (
//            const gsPointestSpaceuator<T>     &testSpace,
//            const gsPointestSpaceuator<T>     &unknownSpace,
//            const gsPointGeometryEvaluator<T>  &geoEval
//            ) const {GISMO_NO_IMPLEMENTATION}
//    virtual void pointEval (
//            const gsPointFuncData<T>  &testSpace,
//            const gsPointFuncData<T>  &unknownSpace,
//            const gsPointMapData <T>  &geoEval,
//            gsRecipeAccumulator<T>  result
//            ) const {GISMO_NO_IMPLEMENTATION}
//    virtual gsMatrix<T>  pointEval (
//            const gsPointFuncData<T>  &tSpaceMaster,
//            const gsPointFuncData<T>  &uSpaceMaster,
//            const gsPointFuncData<T>  &tSpaceSlave,
//            const gsPointFuncData<T>  &uSpaceSlave,
//            const gsPointMapData <T>  &geoEval,
//            gsRecipeAccumulator<T>  result
//            ) const = 0;

//    virtual unsigned    testSpaceNeeds() const  = 0;
//    virtual unsigned    unknownSpaceNeeds() const  = 0;
//    virtual unsigned    geometryNeeds()    const  = 0;
//    virtual void        outputSize (unsigned &r, unsigned &c) const {}
//};


template <typename T>
class gsLaplaceLaplaceOp : public gsBilinearOp<T>
{
public:
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_LAPLACIAN;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_LAPLACIAN ;}
};

template<typename T>
class gsLaplaceOp : public gsBilinearOp<T>
{
public:
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    //BUG! I Need to pass flags here. This should not be necessary for the test space
    virtual unsigned  testSpaceNeeds()    const {return NEED_VALUE;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_LAPLACIAN;}
};


template <typename T>
class gsGradGradOp : public gsBilinearOp<T>
{
public:
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_GRAD;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_GRAD;}
};

template <typename T>
class gsGenericSecondOrderOp : public gsBilinearOp<T>
{
protected:
    gsFunction<T> *m_A;
    gsFunction<T> *m_b;
    gsFunction<T> *m_c;

public:
    gsGenericSecondOrderOp (gsFunction<T> *A,gsFunction<T> *b=NULL,gsFunction<T> *c=NULL)
        : m_A(A),m_b(b),m_c(c)
    {}
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_VALUE |NEED_GRAD;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE |NEED_GRAD;}
    virtual unsigned  geometryNeeds()     const {return NEED_VALUE ;}
};

template <typename T>
class gsL2ScalarOp : public gsBilinearOp<T>
{
public:
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_VALUE;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE;}
};


template <typename T>
class gsL2TestOp : public gsBilinearOp<T>
{
private:
    const gsFunction<T>    *testFunc;
public:
    gsL2TestOp(const gsFunction<T> &test)
        : testFunc(&test)
    {}
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds() const    {return NEED_VALUE;}
    virtual unsigned  unknownSpaceNeeds() const {return 0;}
    virtual unsigned  geometryNeeds()    const  {return NEED_VALUE;}
    virtual void      outputSize (unsigned & /*r*/,
                                  unsigned & c) const
    { c=testFunc->targetDim(); }
};


template <typename T>
class gsL2TestVecOp : public gsBilinearOp<T>
{
private:
    const gsFunction<T>    *testFunc;
public:
    gsL2TestVecOp(const gsFunction<T> &test)
        : testFunc(&test)
    {}
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_VALUE;}
    virtual unsigned  unknownSpaceNeeds() const {return 0;}
    virtual unsigned  geometryNeeds()     const {return NEED_VALUE;}
    virtual void      outputSize (unsigned & /*r*/,
                                  unsigned & c) const
    { c=1; }
};


template <typename T>
class gsDivergenceOp : public gsBilinearOp<T>
{
public:
    gsDivergenceOp() {}
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds()    const {return NEED_VALUE ;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_DIV;}
};

//JS2: (grad(u), v) not tested
template <typename T>
class gsGradientOp : public gsBilinearOp<T>
{
public:
    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds()    const {return NEED_VALUE;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
};


template <typename T>
class gsLinElastOp : public gsBilinearOp<T>
{
public:

    gsLinElastOp( T YoungsModulus, T PoissonsRatio )
    {
        // maybe superfluous to keep all parameters,
        // but these two extra variables won't hurt.
        m_YoungsModulus=YoungsModulus;
        m_PoissonsRatio=PoissonsRatio;

        m_lambda = m_YoungsModulus * m_PoissonsRatio / ( (1.+m_PoissonsRatio)*(1.-2.*m_PoissonsRatio)) ;
        m_mu     = m_YoungsModulus / (2.*(1.+m_PoissonsRatio)) ;
    }
    virtual void pointEval (

            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned  testSpaceNeeds()    const {return NEED_VALUE | NEED_GRAD | NEED_DIV;}
    virtual unsigned  unknownSpaceNeeds() const {return NEED_VALUE | NEED_GRAD | NEED_DIV;}

    T m_YoungsModulus;
    T m_PoissonsRatio;
    T m_lambda;
    T m_mu;
};


// ///////////////////// boundary

template <typename T>
class gsBoundaryL2TestOp :
        public gsBilinearOp<T>
{
private:
    gsFunction<T>    *testFunc;
    mutable gsMatrix<T>       temp;
public:
    gsBoundaryL2TestOp(gsFunction<T>    &f)
        : testFunc(&f)
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const    {return NEED_VALUE;}
    virtual unsigned    unknownSpaceNeeds() const {return 0;}
    virtual unsigned    geometryNeeds()    const  {return NEED_VALUE;}
    virtual void        outputSize (unsigned & /*r*/,
                                    unsigned & c) const
    { c=testFunc->targetDim(); }
};


template <typename T>
class gsBoundaryL2TestVecOp :
        public gsBilinearOp<T>
{
private:
    gsFunction<T>    *testFunc;
    mutable gsMatrix<T>       temp;

public:
    gsBoundaryL2TestVecOp(gsFunction<T>    &f)
        : testFunc(&f)
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const {return NEED_VALUE;}
    virtual unsigned    unknownSpaceNeeds() const {return 0;}
    virtual unsigned    geometryNeeds()    const  {return NEED_VALUE ;}
    virtual void        outputSize (unsigned & /*r*/,
                                    unsigned & c) const
    { c=1; }
};

template <typename T>
class gsBoundaryL2ScalarOp :
        public gsBilinearOp<T>
{
private:

public:
    gsBoundaryL2ScalarOp()
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const {return NEED_VALUE;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_VALUE;}
};

template <typename T>
class gsBoundaryNormalDerValueOp :
        public gsBilinearOp<T>
{
private:

public:
    gsBoundaryNormalDerValueOp()
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const    {return NEED_VALUE;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    geometryNeeds()      const {return NEED_OUTER_NORMAL;}
};

template <typename T>
class gsBoundaryNormalDerNormalDerOp :
        public gsBilinearOp<T>
{
public:
    gsBoundaryNormalDerNormalDerOp()
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const    {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return NEED_GRAD;}
    virtual unsigned    geometryNeeds()     const {return NEED_OUTER_NORMAL;}
};

template <typename T>
class gsBoundaryNormalDerTestOp :
        public gsBilinearOp<T>
{
private:
    gsFunction<T>    *testFunc;
    mutable gsMatrix<T>       temp;
public:
    gsBoundaryNormalDerTestOp(gsFunction<T>    &f)
        : testFunc(&f)
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds()    const {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return 0;}
    virtual unsigned    geometryNeeds()     const {return NEED_VALUE | NEED_OUTER_NORMAL;}
    virtual void        outputSize (unsigned & /*r*/,
                                    unsigned & c) const
    { c=1; }
};

template <typename T>
class gsBoundaryNormalDerTestNormalOp :
        public gsBilinearOp<T>
{
private:
    gsFunction<T>    *testFunc;
    mutable gsMatrix<T>       temp;
public:
    gsBoundaryNormalDerTestNormalOp(gsFunction<T>    &f)
        : testFunc(&f)
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds()    const {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return 0;}
    virtual unsigned    geometryNeeds()     const {return NEED_VALUE | NEED_OUTER_NORMAL;}
    virtual void        outputSize (unsigned & /*r*/,
                                    unsigned & c) const
    {  c=1; }
};


template <typename T>
class gsBoundaryNormalDerTestVecOp :
        public gsBilinearOp<T>
{
private:
    gsFunction<T>    *testFunc;
    mutable gsMatrix<T>       temp;
public:
    gsBoundaryNormalDerTestVecOp(gsFunction<T>    &f)
        : testFunc(&f)
    {}

    virtual void pointEval (
            const gsPointFuncData<T>  &testSpace,
            const gsPointFuncData<T>  &unknownSpace,
            const gsPointMapData <T>  &geoEval,
            gsRecipeAccumulator<T>  result
            ) const;
    virtual unsigned    testSpaceNeeds() const    {return NEED_GRAD;}
    virtual unsigned    unknownSpaceNeeds() const {return 0;}
    virtual unsigned    geometryNeeds()    const  {return NEED_VALUE | NEED_OUTER_NORMAL;}
    virtual void        outputSize (unsigned & /*r*/,
                                    unsigned & c) const
    { c=1; }
};





} // namespace gismo
