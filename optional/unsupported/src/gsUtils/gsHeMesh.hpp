
#pragma once

namespace gismo
{


//////////////////////////////////////////////////////
// Source
//////////////////////////////////////////////////////

template<class T>
gsHeMesh<T>::~gsHeMesh()
{
    //apaga os vertices
    typename std::vector<gsHalfVertexHandle>::iterator vIter;
    for(vIter = vertex.begin(); vIter != vertex.end(); vIter++) {
        delete *vIter;
    }
    vertex.clear();

    //apaga as arestas
    typename std::vector<gsHalfEdgeHandle>::iterator eIter;
    for(eIter = edge.begin(); eIter != edge.end(); eIter++) {
        delete *eIter;
    }
    edge.clear();

    //apaga as faces
    typename std::vector<gsHalfFaceHandle>::iterator fIter;
    for(fIter = face.begin(); fIter != face.end(); fIter++) {
        delete *fIter;
    }
    face.clear();
}

template<class T>
void gsHeMesh<T>::addVertex(gsHalfVertexHandle v)
{
    vertex.push_back(v);
    numVertices++;
}

template<class T>
void gsHeMesh<T>::initialize()
{  
    setHeMate();
    setValence();
    setFaceClassifications();
    setRegularity(6);
    setBasisFactor();
    setBoxsplineTable();
}

template<class T>
void gsHeMesh<T>::solveEquation(T x0, T x1, T x2, T x3, T y0, T y1, T y2, T y3, T x1y1, T x2y1, T x3y1, T x1y2, T x2y2, T x3y2, T x1y3, T x2y3, T x3y3)
{
    gsHeMesh<T> refmesh;
    gsMatrix<T> loadvectorx0(10, 1);
    gsMatrix<T> loadvectorx1(10, 1);
    gsMatrix<T> loadvectorx2(10, 1);
    gsMatrix<T> loadvectorx3(10, 1);
    gsMatrix<T> loadvectory0(10, 1);
    gsMatrix<T> loadvectory1(10, 1);
    gsMatrix<T> loadvectory2(10, 1);
    gsMatrix<T> loadvectory3(10, 1);
    gsMatrix<T> loadvectorx1y1(10, 1);
    gsMatrix<T> loadvectorx2y1(10, 1);
    gsMatrix<T> loadvectorx3y1(10, 1);
    gsMatrix<T> loadvectorx1y2(10, 1);
    gsMatrix<T> loadvectorx2y2(10, 1);
    gsMatrix<T> loadvectorx3y2(10, 1);
    gsMatrix<T> loadvectorx1y3(10, 1);
    gsMatrix<T> loadvectorx2y3(10, 1);
    gsMatrix<T> loadvectorx3y3(10, 1);
    gsMatrix<T> massmatrix(10, 10);
    refmesh.buildReferenceConfiguration();
    std::vector<gsHalfFaceHandle> refmeshfaces = refmesh.getHalfFaces();
    std::vector<gsHalfVertexHandle> refmeshvertices = refmesh.getVertices();
    massmatrix = getMassMatrix(refmeshfaces, refmeshvertices);
    loadvectorx0 = getLoadVector(refmeshfaces, refmeshvertices, x0);
    loadvectorx1 = getLoadVector(refmeshfaces, refmeshvertices, x1);
    loadvectorx2 = getLoadVector(refmeshfaces, refmeshvertices, x2);
    loadvectorx3 = getLoadVector(refmeshfaces, refmeshvertices, x3);
    loadvectory0 = getLoadVector(refmeshfaces, refmeshvertices, y0);
    loadvectory1 = getLoadVector(refmeshfaces, refmeshvertices, y1);
    loadvectory2 = getLoadVector(refmeshfaces, refmeshvertices, y2);
    loadvectory3 = getLoadVector(refmeshfaces, refmeshvertices, y3);
    loadvectorx1y1 = getLoadVector(refmeshfaces, refmeshvertices, x1y1);
    loadvectorx2y1 = getLoadVector(refmeshfaces, refmeshvertices, x2y1);
    loadvectorx3y1 = getLoadVector(refmeshfaces, refmeshvertices, x3y1);
    loadvectorx1y2 = getLoadVector(refmeshfaces, refmeshvertices, x1y2);
    loadvectorx2y2 = getLoadVector(refmeshfaces, refmeshvertices, x2y2);
    loadvectorx3y2 = getLoadVector(refmeshfaces, refmeshvertices, x3y2);
    loadvectorx1y3 = getLoadVector(refmeshfaces, refmeshvertices, x1y3);
    loadvectorx2y3 = getLoadVector(refmeshfaces, refmeshvertices, x2y3);
    loadvectorx3y3 = getLoadVector(refmeshfaces, refmeshvertices, x3y3);

    gsMatrix<T> solutionx0 = LUDecomposition(massmatrix,loadvectorx0);
    gsMatrix<T> solutionx1 = LUDecomposition(massmatrix,loadvectorx1);
    gsMatrix<T> solutionx2 = LUDecomposition(massmatrix,loadvectorx2);
    gsMatrix<T> solutionx3 = LUDecomposition(massmatrix,loadvectorx3);
    gsMatrix<T> solutiony0 = LUDecomposition(massmatrix,loadvectory0);
    gsMatrix<T> solutiony1 = LUDecomposition(massmatrix,loadvectory1);
    gsMatrix<T> solutiony2 = LUDecomposition(massmatrix,loadvectory2);
    gsMatrix<T> solutiony3 = LUDecomposition(massmatrix,loadvectory3);
    gsMatrix<T> solutionx1y1 = LUDecomposition(massmatrix,loadvectorx1y1);
    gsMatrix<T> solutionx2y1 = LUDecomposition(massmatrix,loadvectorx2y1);
    gsMatrix<T> solutionx3y1 = LUDecomposition(massmatrix,loadvectorx3y1);
    gsMatrix<T> solutionx1y2 = LUDecomposition(massmatrix,loadvectorx1y2);
    gsMatrix<T> solutionx2y2 = LUDecomposition(massmatrix,loadvectorx2y2);
    gsMatrix<T> solutionx3y2 = LUDecomposition(massmatrix,loadvectorx3y2);
    gsMatrix<T> solutionx1y3 = LUDecomposition(massmatrix,loadvectorx1y3);
    gsMatrix<T> solutionx2y3 = LUDecomposition(massmatrix,loadvectorx2y3);
    gsMatrix<T> solutionx3y3 = LUDecomposition(massmatrix,loadvectorx3y3);
    std::cout << "\n" << "\n" <<"Solution of the equation in parts: " << "\n";
    for(index_t i = 0; i < 10; i++)
    {
        std::cout << solutionx0(i,0) << "+" << solutionx1(i,0) << "+" << solutionx2(i,0) << "+" << solutionx3(i,0) << "\n";
       // std::cout << solutiony0(i,0) << "+" << solutiony1(i,0) << "+" << solutiony2(i,0) << "+" << solutiony3(i,0) << "\n";
    }
    T resultx0 = 0;
    T resultx1 = 0;
    T resultx2 = 0;
    T resultx3 = 0;
    T resulty0 = 0;
    T resulty1 = 0;
    T resulty2 = 0;
    T resulty3 = 0;
    T resultx1y1 = 0;
    T resultx2y1 = 0;
    T resultx3y1 = 0;
    T resultx1y2 = 0;
    T resultx2y2 = 0;
    T resultx3y2 = 0;
    T resultx1y3 = 0;
    T resultx2y3 = 0;
    T resultx3y3 = 0;

    int domaincounter = 0;
    gsMatrix<T> basis;
    for ( typename std::vector<gsHalfFaceHandle>::iterator
              it = face.begin(); it!= face.end(); ++it)
    {
        gsHalfFaceHandle fa = *it;
        if(fa->getClassification() == 0)
        {
            basis = basisfunctions(0.5,0.5,0.0,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
            for(index_t i = 0; i < 10; i++)
            {
                resultx0 = resultx0 + solutionx0(i,0)*basis(0,i);
                resultx1 = resultx1 + solutionx1(i,0)*basis(0,i);
                resultx2 = resultx2 + solutionx2(i,0)*basis(0,i);
                resultx3 = resultx3 + solutionx3(i,0)*basis(0,i);
                resulty0 = resulty0 + solutiony0(i,0)*basis(0,i);
                resulty1 = resulty1 + solutiony1(i,0)*basis(0,i);
                resulty2 = resulty2 + solutiony2(i,0)*basis(0,i);
                resulty3 = resulty3 + solutiony3(i,0)*basis(0,i);
                resultx1y1 = resultx1y1 + solutionx1y1(i,0)*basis(0,i);
                resultx2y1 = resultx2y1 + solutionx2y1(i,0)*basis(0,i);
                resultx3y1 = resultx3y1 + solutionx3y1(i,0)*basis(0,i);
                resultx1y2 = resultx1y2 + solutionx1y2(i,0)*basis(0,i);
                resultx2y2 = resultx2y2 + solutionx2y2(i,0)*basis(0,i);
                resultx3y2 = resultx3y2 + solutionx3y2(i,0)*basis(0,i);
                resultx1y3 = resultx1y3 + solutionx1y3(i,0)*basis(0,i);
                resultx2y3 = resultx2y3 + solutionx2y3(i,0)*basis(0,i);
                resultx3y3 = resultx3y3 + solutionx3y3(i,0)*basis(0,i);
            }
            domaincounter++;
        }
    }
    resultx0 = resultx0/domaincounter;
    resultx1 = resultx1/domaincounter;
    resultx2 = resultx2/domaincounter;
    resultx3 = resultx3/domaincounter;
    resulty0 = resulty0/domaincounter;
    resulty1 = resulty1/domaincounter;
    resulty2 = resulty2/domaincounter;
    resulty3 = resulty3/domaincounter;
    resultx1y1 = resultx1y1/domaincounter;
    resultx2y1 = resultx2y1/domaincounter;
    resultx3y1 = resultx3y1/domaincounter;
    resultx1y2 = resultx1y2/domaincounter;
    resultx2y2 = resultx2y2/domaincounter;
    resultx3y2 = resultx3y2/domaincounter;
    resultx1y3 = resultx1y3/domaincounter;
    resultx2y3 = resultx2y3/domaincounter;
    resultx3y3 = resultx3y3/domaincounter;
    std::cout << "Result X: " << resultx0 << " + " << resultx1 << "x + " << resultx2 << "x^2 + " << resultx3 << "x^3" << "\n";
    std::cout << "Result Y: " << resulty0 << " + " << resulty1 << "y + " << resulty2 << "y^2 + " << resulty3 << "y^3" << "\n";
}

template<class T>
gsMatrix<T> gsHeMesh<T>::interpolate(std::vector<T> y, int option)
{
    //option 1=linear, 2=quadratic, 3=cubic
    if(option == 1)
    {
        int n = y.size();
        gsMatrix<T> result(1, 2);
        std::vector<T> x;
        T sumx = 0;
        T sumy = 0;
        T alpha = 0;
        T beta = 0;
        T upper = 0;
        T lower = 0;
        for(int i = 0; i <= n-1; i++)
        {
            x.push_back((double)i/(n-1));
            sumx = sumx + x[i];
            sumy = sumy + y[i];
        }
        for(int i = 0; i <= n-1; i++)
        {
            upper = upper + x[i]*y[i];
            lower = lower + x[i]*x[i];
        }
        upper = upper / n;
        lower = lower / n;
        upper = upper - sumx*sumy;
        lower = lower - sumx*sumx;
        beta = upper/lower;
        alpha = sumy - beta*sumx;
        result(0,0) = alpha;
        result(0,1) = beta;
        std::cout << result(0,1) << " x + " << result(0,0) << "\n";
        return result;
    }
    if(option == 2)
    {
        int n = y.size();
        gsMatrix<T> result(1, 3);
        std::vector<T> x;
        T s40 = 0;
        T s30 = 0;
        T s20 = 0;
        T s10 = 0;
        T s00 = 0;
        T s21 = 0;
        T s11 = 0;
        T s01 = 0;
        T a = 0;
        T b = 0;
        T c = 0;
        T det = 0;
        for(int i = 0; i <= n-1; i++)
        {
            x.push_back((double)i/(n-1));
        }
        for(int i = 0; i <= n-1; i++)
        {
            s40 += pow(x[i],4);
            s30 += pow(x[i],3);
            s20 += pow(x[i],2);
            s10 += x[i];
            s00 += 1;
            s21 += pow(x[i],2)*y[i];
            s11 += x[i]*y[i];
            s01 += y[i];
        }
        //Cramer's law
        det = s40*(s20*s00 - s10*s10) - s30*(s30*s00 - s10*s20) + s20*(s30*s10 - s20*s20);
        a = (s21*(s20*s00 - s10*s10) - s11*(s30*s00 - s10*s20) + s01*(s30*s10 - s20*s20))/det;
        b = (s40*(s11*s00 - s01*s10) - s30*(s21*s00 - s01*s20) + s20*(s21*s10 - s11*s20))/det;
        c = (s40*(s20*s01 - s10*s11) - s30*(s30*s01 - s10*s21) + s20*(s30*s11 - s20*s21))/det;
        result(0,2) = a;
        result(0,1) = b;
        result(0,0) = c;
        std::cout << result(0,2) << " x^2 + " << result(0,1) << " x + " << result(0,0) << "\n";
        return result;
    }
    if(option == 3)
    {
        //Zheng Yan algorithm
        int n = y.size();
        std::vector<T> x;
        gsMatrix<T> result(n, 4);
        std::vector<T> d;
        std::vector<T> s;
        std::vector<T> delta;
        std::vector<T> deltax;
        for(int i = 0; i <= n-1; i++)
        {
            x.push_back((double)i/(n-1));
        }
        for(int i = 0; i <= n-2; i++)
        {
            delta.push_back((y[i+1]-y[i])/(x[i+1]-x[i]));
            deltax.push_back(x[i+1]-x[i]);
        }
    //    delta.push_back(y[n-1]/x[n-1]);
    //    std::cout << delta[n-1] << "\n";
    //    deltax.push_back(x[n-1]);
        d.push_back(0); //??
        s.push_back(d[0]);
        for(int i = 1; i <= n-3; i++)
        {
            d.push_back(delta[i-1]+(delta[i]-delta[i-1])/(deltax[i]+deltax[i-1])*deltax[i-1]
                    +((delta[i]-delta[i-1])/(deltax[i]+deltax[i-1])-(delta[i+1]-delta[i])/(deltax[i+1]+deltax[i]))/(deltax[i+1]+deltax[i]+deltax[i-1])*deltax[i-1]*deltax[i]);

            if(d[i]*delta[i]<=0)
                s.push_back(0);
            else
                s.push_back(d[i]);
            std::cout<<s[i]<<"\n";

        }
     //   d.push_back(d[n-2]+(d[n-2]-d[n-3]));
     //   s.push_back(d[n-2]+(d[n-2]-d[n-3]));
        T resultx3 = 0;
        T resultx2 = 0;
        T resultx1 = 0;
        T resultx0 = 0;
        std::cout << "{";
        for(int i = 0; i <= n-4; i++)
        {
            result(i,0) = ((s[i]+s[i+1]-2*delta[i])/pow(deltax[i],2))*(pow(-x[i],3)) + ((-2*s[i]-s[i+1]+3*delta[i])/deltax[i]) * pow(x[i],2) + s[i]*(-x[i]) + y[i];
            result(i,1) = ((s[i]+s[i+1]-2*delta[i])/pow(deltax[i],2))*3*(pow(x[i],2)) + ((-2*s[i]-s[i+1]+3*delta[i])/deltax[i]) * 2 *(-x[i]) + s[i];
            result(i,2) = ((s[i]+s[i+1]-2*delta[i])/pow(deltax[i],2))*3*(-x[i]) + ((-2*s[i]-s[i+1]+3*delta[i])/deltax[i]);
            result(i,3) = ((s[i]+s[i+1]-2*delta[i])/pow(deltax[i],2));
            resultx3 = resultx3 + result(i,3);
            resultx2 = resultx2 + result(i,2);
            resultx1 = resultx1 + result(i,1);
            resultx0 = resultx0 + result(i,0);
            std::cout << "{" << result(i,3) << " x^3 + " << result(i,2) << " x^2 + " << result(i,1) << " x + " << result(i,0) << "," << x[i] << "," << x[i+1] << "}," << "\n";
        }
        resultx3 = resultx3/(n-3);
        resultx2 = resultx2/(n-3);
        resultx1 = resultx1/(n-3);
        resultx0 = resultx0/(n-3);
        std::cout << "\n";
        std::cout << resultx3 << " x^3 + " << resultx2 << " x^2 + " << resultx1 << " x + " << resultx0 << "\n";
        return result;
    }
}

template<class T>
gsMatrix<T> gsHeMesh<T>::LUDecomposition(gsMatrix<T> A, gsMatrix<T> b)
{
    int msize = 10;
    gsMatrix<T> L(msize,msize);
    gsMatrix<T> U(msize,msize);
    U = A;

    for(index_t k = 0; k < msize; k++)
    {
        L(k,k) = 1;
        for(index_t j = k+1; j < msize; j++)
        {
            L(j,k)=U(j,k)/U(k,k);
            for(index_t i = k; i < msize; i++)
            {
                U(j,i)=U(j,i)-L(j,k)*U(k,i);
            }
        }
    }

    T temp = b(0,0);
    gsMatrix<T> y(msize,1);
    y(0,0) = temp;
    for(index_t i = 1; i < msize; i++)
    {
        temp = 0;
        for(index_t j = i-1; j >= 0; j--)
        {
            temp = temp + L(i,j)*y(j,0);
        }
        y(i,0) = b(i,0) - temp;
    }

    gsMatrix<T> x(msize,1);
    x(msize-1,0) = y(msize-1,0)/U(msize-1,msize-1);
    for(index_t i = msize-2; i >= 0; i--)
    {
        temp = 0;
        for(index_t j = msize-1; j > i; j--)
        {
            temp = temp + U(i,j) * x(j,0);
        }
        x(i,0) = (y(i,0)-temp)/U(i,i);
    }
    return(x);
}

template<class T>
void gsHeMesh<T>::setHeMate()
{
    // Method in the meantime: consider each pair of halfedges and see if they are mates
    // Todo: consider each pair of faces
    unsigned int noMate(0); // number of mates
    bool outer = true;
    std::vector<gsHalfEdgeHandle> boundarycandidate;
    std::vector<gsHalfVertexHandle> ver;
    std::vector<gsHalfEdgeHandle> edgecopy;
    edgecopy = edge;

    gsHalfFaceHandle f = new gsHalfFace<T>();
    f->setClassification(2);
    for (typename std::vector<gsHalfEdgeHandle>::iterator it1=edge.begin(); it1!=edge.end(); ++it1)
    {
        for (typename std::vector<gsHalfEdgeHandle>::iterator it2=edge.begin(); it2!=edge.end(); ++it2)
        {
             // a pair of half edges are mates iff the source of one of them is the target of the orther
            gsHalfVertexHandle source1, source2, target1, target2;
            source1 = (*it1)->getSource();
            source2 = (*it2)->getSource();
            target1 = (*it1)->getNext()->getSource();
            target2 = (*it2)->getNext()->getSource();
            if ((source1 == target2
                && source2 == target1) &&(source1 != source2 && target1 != target2))
            {
                noMate++;
                outer = false;
                (*it1)->setMate(*it2);
                (*it2)->setMate(*it1);
            }
        }
        if(outer)
        {
            gsHalfEdgeHandle heh = this->makeHalfEdge((*it1)->getNext()->getSource(),f);
            ver.push_back((*it1)->getSource());
            noMate++;
            (*it1)->setMate(heh);
            heh->setMate(*it1);
            heh->setFace(f);
            heh->setId(numHalfEdges++);
            boundarycandidate.push_back(heh);
            f->setVertex(heh->getSource());
        }
        outer = true;
    }
    int i = 0;
    for (typename std::vector<gsHalfEdgeHandle>::iterator it=boundarycandidate.begin(); it!=boundarycandidate.end(); ++it)
    {
        gsHalfEdgeHandle heh = *it;
        for (typename std::vector<gsHalfEdgeHandle>::iterator it2=boundarycandidate.begin(); it2!=boundarycandidate.end(); ++it2)
        {
            if((*it2)->getSource() == ver[i])
                heh->setNext((*it2));
        }
        edge.push_back(heh);
        i++;
    }
    face.push_back(f);
    f->setId(numHalfFaces++);
}

template<class T>
typename gsHeMesh<T>::gsHalfVertexHandle gsHeMesh<T>::addVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z)
{
    gsHalfVertexHandle v = this->makeHeVertex(x,y,z);
    vertex.push_back(v );
    v->setId(numVertices++);
    return v;
}

template<class T>
void gsHeMesh<T>::addHalfEdge(gsHalfEdgeHandle he)
{
    edge.push_back(he);
    he->setId(numHalfEdges++);
}

template<class T>
typename gsHeMesh<T>::gsHalfFaceHandle gsHeMesh<T>::addHalfFace(std::vector<gsHalfVertexHandle> V)
{
    gsHalfFaceHandle f = new gsHalfFace<T>();
    std::vector<gsHalfEdgeHandle> boundarycandidate;

    std::vector<gsHalfEdgeHandle> H;
    std::vector<gsHalfEdgeHandle> tempHeNext;
    for ( typename std::vector<gsHalfVertexHandle>::iterator
          it = V.begin(); it!= V.end(); ++it)
    {
        f->setVertex(*it);
        // initialize the half-edge associated with the vertex *it and the face f
        gsHalfEdgeHandle tempHe = this->makeHalfEdge(*it,f);
        (*it)->setHalfEdge(tempHe);
        tempHeNext.push_back(tempHe);
    }
    gsHalfEdgeHandle next = 0;
    for ( typename std::vector<gsHalfEdgeHandle>::iterator
          it = tempHeNext.begin(); it!= tempHeNext.end(); ++it)
    {
        gsHalfEdgeHandle tempHe = *it;
        it++;
        if(it==tempHeNext.end())
            next = *(tempHeNext.begin());
        else
            next = *it;
        it--;
        tempHe->setNext(next);
        tempHe->setId(numHalfEdges++);
        boundarycandidate.push_back(tempHe);
        edge.push_back(tempHe);
    }
    f->setBoundary(boundarycandidate);
    face.push_back(f);
    f->setId(numHalfFaces++);
    return f;
}

template<class T>
int gsHeMesh<T>::getValence(gsHalfVertexHandle s)
{
    int val = 0;
    gsHalfEdgeHandle control = s->getHalfEdge();
    gsHalfEdgeHandle heh = control;
    do
    {
        heh = heh->getMate()->getNext();
        val++;
    }
    while(control != heh);
    return val;
}

template<class T>
int gsHeMesh<T>::getValence(gsHalfVertexHandle s, std::list<int> indices)
{
    int val = 0;
    bool redundant = false;
    gsHalfEdgeHandle control = s->getHalfEdge();
    gsHalfEdgeHandle heh = control;
    do
    {
        heh = heh->getMate()->getNext();
        for(typename std::list<int>::iterator
            it = indices.begin(); it!= indices.end(); ++it)
        {
            if(heh->getMate()->getSource()->getId() == (*it))
            {
                redundant = true;
                break;
            }
        }
        if(redundant)
            val++;
        redundant = false;
    }
    while(control != heh);
    return val;

}

template<class T>
void gsHeMesh<T>::setValence()
{
    for(typename std::vector<gsHalfVertexHandle>::iterator
        it = vertex.begin(); it!= vertex.end(); ++it)
    {
        gsHalfVertexHandle tempvertex = *it;
        tempvertex->setValence(getValence(tempvertex));
    }
}


template<class T>
std::list<int> gsHeMesh<T>::get1Ring(gsHalfFaceHandle fa)
{
    std::list<int> indices;
    std::list<int> tempindices;
    std::vector<gsHalfVertexHandle> v = fa->getVertices();
    for ( typename std::vector<gsHalfVertexHandle>::iterator
              it = v.begin(); it!= v.end(); ++it)
    {
        gsHalfVertexHandle hvh = *it;
        tempindices = get1Ring(hvh);
        for ( typename std::list<int>::iterator
                  it2 = tempindices.begin(); it2!= tempindices.end(); ++it2)
        {
             indices.push_back(*it2);
        }
    }
    indices.sort(); //sort do detect multiple equal elements
    indices.unique(); //delete equal elements
    return indices;
}

template<class T>
std::list<int> gsHeMesh<T>::get1Ring(gsHalfVertexHandle ve)
{
    std::list<int> indices;
    gsHalfEdgeHandle control = ve->getHalfEdge();
    gsHalfEdgeHandle heh = control;
    do
    {
        heh = heh->getMate();
        indices.push_back(heh->getSource()->getId());
        heh = heh->getNext();
    }
    while(control != heh);
    return indices;
}

template<class T>
std::list<int> gsHeMesh<T>::get2Ring(gsHalfFaceHandle fa)
{
    std::list<int> indices;
    std::list<int> tempindices1;
    std::list<int> tempindices2;
    std::vector<gsHalfVertexHandle> v = fa->getVertices();
    for ( typename std::vector<gsHalfVertexHandle>::iterator
              it = v.begin(); it!= v.end(); ++it)
    {
        gsHalfVertexHandle hvh = *it;
        tempindices1 = get1Ring(hvh);
        for ( typename std::list<int>::iterator
                  it2 = tempindices1.begin(); it2!= tempindices1.end(); ++it2)
        {
            tempindices2 = get1Ring(getVertexAtId(*it2));
            for ( typename std::list<int>::iterator
                      it3 = tempindices2.begin(); it3!= tempindices2.end(); ++it3)
            {
                 indices.push_back(*it3);
            }
        }
    }
    indices.sort(); //sort do detect multiple equal elements
    indices.unique(); //delete equal elements
    return indices;
}

template<class T>
void gsHeMesh<T>::setFaceClassifications()
{
    bool domain = true;
    for ( typename std::vector<gsHalfFaceHandle>::iterator
          it = face.begin(); it!= face.end(); ++it)
    {
        gsHalfFaceHandle hfh = *it;
        if(hfh->getClassification() == -1)
        {
            std::vector<gsHalfEdgeHandle> boundary = hfh->getBoundary();
            for ( typename std::vector<gsHalfEdgeHandle>::iterator
                  it2 = boundary.begin(); it2!= boundary.end(); ++it2)
            {
                gsHalfEdgeHandle tempHe = *it2;
                if(tempHe->getSource()->getValence()!=6)
                    domain = false;
            }
            if(domain)
                hfh->setClassification(0);
            else
                hfh->setClassification(1);
            domain = true;
        }
    }
}

template<class T>
std::vector<typename gsHeMesh<T>::gsHalfFaceHandle> gsHeMesh<T>::getFacesByClassification(int classification)
{
    std::vector<gsHalfFaceHandle> fa;
    for ( typename std::vector<gsHalfFaceHandle>::iterator
          it = face.begin(); it!= face.end(); ++it)
    {
        gsHalfFaceHandle hfh = *it;
        if(hfh->getClassification() == classification)
        {
            fa.push_back(hfh);
        }
    }
    if(fa.size() > 0)
        return fa;
    else
    {
        std::cout << "No Faces found! Return all faces.";
        return face;
    }
}

template<class T>
void gsHeMesh<T>::setRegularity(int check)
{
    bool regtest = true;
    std::vector<gsHalfFaceHandle > fa;
    for ( typename std::vector<gsHalfFaceHandle>::iterator
          it = face.begin(); it!= face.end(); ++it)
    {
        std::vector<gsHalfVertexHandle> v = (*it)->getVertices();
        for ( typename std::vector<gsHalfVertexHandle>::iterator
              it2 = v.begin(); it2!= v.end(); ++it2)
        {
            if((*it2)->getValence() != check)
                regtest = false;
        }
        if(regtest)
        {
            (*it)->setRegular(true);
            fa.push_back((*it));
        }
        regtest = true;
    }
    setRegularFaces(fa);
}

template<class T>
void gsHeMesh<T>::setBasisFactor()
{
    int l = 0;
    for ( typename std::vector<gsHalfFaceHandle>::iterator
          it = regface.begin(); it!= regface.end(); ++it)
    {
        std::list<int> ve = get1Ring((*it));
        for ( typename std::list<int>::iterator
              it2 = ve.begin(); it2!= ve.end(); ++it2)
        {
            gsHalfVertexHandle v = getVertexAtId(*it2);
            if(v->getClassification() == -1)
            {
                v->setClassification(l);
                l++; //? maybe after for this increment
            }
        }
    }
}
template<class T>
int gsHeMesh<T>::getVertexType(typename gsHeMesh<T>::gsHalfVertexHandle v, std::list<int> indices)
{
    int type = -1;

    if(getValence(v,indices) == 3)
    {
        gsHalfEdgeHandle tempedge = v->getHalfEdge();
        while(getValence(tempedge->getNext()->getSource(),indices) != 6)
        {
            tempedge = tempedge->getMate()->getNext();
        }

        if(getValence(tempedge->getNext()->getNext()->getSource(),indices) == 3)
            type = 3;
        if(getValence(tempedge->getNext()->getNext()->getSource(),indices) == 4)
            type = 1;
    }
    if(getValence(v,indices) == 4)
        type = 2;
    if(getValence(v,indices) == 6)
        type = 0;
    return type;
}

template<class T>
T gsHeMesh<T>::VertexDistance(gsHalfVertexHandle v1, gsHalfVertexHandle v2)
{
    return ((v1->x()-v2->x)^2 + (v1->y()-v2->y)^2 + (v1->z()-v2->z)^2)^(1/2);
}

template<class T>
gsMatrix<T> gsHeMesh<T>::evaluate1Ring(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfFaceHandle hf, std::vector<gsHalfFaceHandle> refmeshfaces)
{
     std::vector<std::vector<T> > refvalues;
     std::list<int> indices;
     T x = 0;
     T y = 0;
     T testsum = 0;

     indices = get1Ring(hf);
     for ( typename std::list<int>::iterator
               it = indices.begin(); it!= indices.end(); ++it)
     {
         refvalues.push_back(transformPoint(bary1, bary2, bary3, v1, v2, v3, getVertexAtId(*it),refmeshfaces,indices));
     }

     gsBoxSplineBasis<2,real_t> boxspline(2,2,2);
     gsMatrix<T> u(2,refvalues.size());
     for(index_t i = 0; i < refvalues.size(); i++)
     {
        // std::cout << (i+1) << " X: " << (refvalues[i])[0] << " Y: " << (refvalues[i])[1] << "\n";
         u(0,i) = (refvalues[i])[0];
         u(1,i) = (refvalues[i])[1];
     }

     gsMatrix<T> result = boxspline.evalresult(u);
/*     gsVector<int> r(2);
     r << 0,1;
     gsMatrix<> result = boxspline.derivativeresult(u,r); */
     for(index_t i = 0; i < result.cols(); i++)
     {
         std::cout << result(0,i) << "\n";
         testsum = testsum + result(0,i);
     }
     if(testsum < 1-0.1 || testsum > 1+0.1)
        std::cout <<"Warning! Testsum is unequal to expected value 1 - Testsum has value " << testsum << "\n";
     int j = 0;
     for ( typename std::list<int>::iterator
               it = indices.begin(); it!= indices.end(); ++it)
     {
         std::cout << "PHI: " << result(0,j) << "\n";

         x = x + getVertexAtId(*it)->x()*result(0,j);
         y = y + getVertexAtId(*it)->y()*result(0,j);
 //                  std::cout << (j+1) << " Y: " << y <<"\n";
         j++;
     }
    // gsHalfVertexHandle evalpoint = this->makeHeVertex(x,y,0);
  //   std::cout << std::setprecision(15) << "Evaluationpoint: X: " << evalpoint->x() << " Y: " << evalpoint->y() << "\n";
  //   std::cout << std::setprecision(15) << "Evaluationpoint: X: " << x << " Y: " << y << "\n";
       std::cout << std::setprecision(15) << "{" << x << "," << y << "}" << "\n";
     return result;
}

template<class T>
gsMatrix<T> gsHeMesh<T>::basisfunctions(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfFaceHandle hf, std::vector<gsHalfFaceHandle> refmeshfaces)
{
     std::vector<std::vector<T> > refvalues;
     std::list<int> indices = get1Ring(hf);
     for ( typename std::list<int>::iterator
               it = indices.begin(); it!= indices.end(); ++it)
     {
         refvalues.push_back(transformPoint(bary1, bary2, bary3, v1, v2, v3, getVertexAtId(*it),refmeshfaces,indices));

     }

     gsBoxSplineBasis<2,real_t> boxspline(2,2,2);
     gsMatrix<T> u(2,refvalues.size());
     for(index_t i = 0; i < u.cols(); i++)
     {
         u(0,i) = refvalues[i][0];
         u(1,i) = refvalues[i][1];
     }
     gsMatrix<T> tempresult = boxspline.evalresult(u);
     gsMatrix<T> result(1,refvalues.size()-2);
     int ind = 0;
     for(index_t j = 0; j < tempresult.cols(); j++)
     {
         if(tempresult(0,j) > 0.00001)
         {
             result(0,ind) = tempresult(0,j);
             ind++;
         }
     }
 /*    gsVector<int> r(2);
          r << 1,1;
     gsMatrix<T> tempresult = boxspline.derivativeresult(u,r);
     gsMatrix<T> result(1,refvalues.size()-2);
     int ind = 0;
     std::cout<< "Derivative: " << "\n";
     for(index_t j = 0; j < tempresult.cols(); j++)
     {
         std::cout << tempresult(0,j) << " next: ";
         if(tempresult(0,j) > 0.00001)
         {
             result(0,ind) = tempresult(0,j);
             ind++;
         }
     } */
     return result;
}

template<class T>
gsMatrix<T> gsHeMesh<T>::getMassMatrix(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices)
{
    gsMatrix<T> massmatrix(10, 10); //(getnumVertices(),getnumVertices())
    gsMatrix<T> basis;
    std::list<int> indices;
    T Gx1=0;
    T Gx2=0;
    T Gy1=0;
    T Gy2=0;
    std::cout << "massmatrix: \n";
    for(index_t i = 0; i < massmatrix.rows(); i++)
    {
        for(index_t j = 0; j < massmatrix.cols(); j++)
        {
            massmatrix(i,j) = 0;
//            for ( typename std::vector<gsHalfFaceHandle>::iterator
//                      it = refmeshfaces.begin(); it!= refmeshfaces.end(); ++it)
            for ( typename std::vector<gsHalfFaceHandle>::iterator
                      it = face.begin(); it!= face.end(); ++it)
            {
                gsHalfFaceHandle fa = *it;
                if(fa->getClassification() == 0)
                {
                    Gx1=0;
                    Gx2=0;
                    Gy1=0;
                    Gy2=0;
                    indices = get1Ring(fa);
                    index_t index = 0;
                    index_t ind = 0;
                    for ( typename std::list<int>::iterator
                              it2 = indices.begin(); it2!= indices.end(); ++it2)
                    {
                        index = *it2;
                        Gx1=Gx1+refmeshvertices[index]->x()*(getPhi1())[ind];
                        Gx2=Gx2+refmeshvertices[index]->x()*(getPhi2())[ind];
                        Gy1=Gy1+refmeshvertices[index]->y()*(getPhi1())[ind];
                        Gy2=Gy2+refmeshvertices[index]->y()*(getPhi2())[ind];
                        ind++;
                    }
                    basis = basisfunctions(0.5,0.5,0,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    massmatrix(i,j) = massmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt( math::abs(Gx1*Gy2-Gx2*Gy1));
                    basis = basisfunctions(0.0,0.5,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    massmatrix(i,j) = massmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1));
                    basis = basisfunctions(0.5,0,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    massmatrix(i,j) = massmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1)); }
            }
            std::cout << massmatrix(i,j) << " ";

        }
        std::cout << "\n";
    }
    return massmatrix;
}

template<class T>
gsMatrix<T> gsHeMesh<T>::getStiffnessMatrix(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices)
{
    //not ready
    gsMatrix<T> stiffnessmatrix(10, 10);

    gsMatrix<T> basis;
    std::list<int> indices;
    T Gx1=0;
    T Gx2=0;
    T Gy1=0;
    T Gy2=0;
    for(index_t i = 0; i < stiffnessmatrix.rows(); i++)
    {
        for(index_t j = 0; j < stiffnessmatrix.cols(); j++)
        {
            stiffnessmatrix(i,j) = 0;
//            for ( typename std::vector<gsHalfFaceHandle>::iterator
//                      it = refmeshfaces.begin(); it!= refmeshfaces.end(); ++it)
            for ( typename std::vector<gsHalfFaceHandle>::iterator
                      it = face.begin(); it!= face.end(); ++it)
            {
                gsHalfFaceHandle fa = *it;
                if(fa->getClassification() == 0)
                {
                    Gx1=0;
                    Gx2=0;
                    Gy1=0;
                    Gy2=0;
                    indices = get1Ring(fa);
                    index_t index = 0;
                    index_t ind = 0;
                    for ( typename std::list<int>::iterator
                              it2 = indices.begin(); it2!= indices.end(); ++it2)
                    {
                        index = *it2;
                        Gx1=Gx1+refmeshvertices[index]->x()*(getPhi1())[ind];
                        Gx2=Gx2+refmeshvertices[index]->x()*(getPhi2())[ind];
                        Gy1=Gy1+refmeshvertices[index]->y()*(getPhi1())[ind];
                        Gy2=Gy2+refmeshvertices[index]->y()*(getPhi2())[ind];
                        ind++;
                    }
                    basis = basisfunctions(0.5,0.5,0,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    stiffnessmatrix(i,j) = stiffnessmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt( math::abs(Gx1*Gy2-Gx2*Gy1));
                    basis = basisfunctions(0.0,0.5,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    stiffnessmatrix(i,j) = stiffnessmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1));
                    basis = basisfunctions(0.5,0,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                    stiffnessmatrix(i,j) = stiffnessmatrix(i,j) + 1.0/6.0*basis(0,i)*basis(0,j)*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1)); }
            }
            std::cout << stiffnessmatrix(i,j) << " ";

        }
        std::cout << "\n";
    }
    return stiffnessmatrix;
}

template<class T>
gsMatrix<T> gsHeMesh<T>::getLoadVector(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices, T functionvalue)
{
    gsMatrix<T> loadvector(10, 1);
    gsMatrix<> basis;
    std::list<int> indices;
    T Gx1=0;
    T Gx2=0;
    T Gy1=0;
    T Gy2=0;
    std::cout << "loadvector: \n";
    for(index_t i = 0; i < loadvector.rows(); i++)
    {
        loadvector(i,0) = 0;
        for ( typename std::vector<gsHalfFaceHandle>::iterator
                  it = face.begin(); it!= face.end(); ++it)
        {
            gsHalfFaceHandle fa = *it;
            if(fa->getClassification() == 0)
            {
                Gx1=0;
                Gx2=0;
                Gy1=0;
                Gy2=0;
                indices = get1Ring(fa);
                int index = 0;
                int ind = 0;
                for ( typename std::list<int>::iterator
                          it2 = indices.begin(); it2!= indices.end(); ++it2)
                {
                    index = *it2;
                    Gx1=Gx1+refmeshvertices[index]->x()*(getPhi1())[ind];
                    Gx2=Gx2+refmeshvertices[index]->x()*(getPhi2())[ind];
                    Gy1=Gy1+refmeshvertices[index]->y()*(getPhi1())[ind];
                    Gy2=Gy2+refmeshvertices[index]->y()*(getPhi2())[ind];
                    ind++;
                }


                basis = basisfunctions(0.5,0.5,0,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                loadvector(i,0) = loadvector(i,0) + 1.0/6.0*basis(0,i)*functionvalue*math::sqrt( math::abs(Gx1*Gy2-Gx2*Gy1));
                basis = basisfunctions(0.0,0.5,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                loadvector(i,0) = loadvector(i,0) + 1.0/6.0*basis(0,i)*functionvalue*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1));
                basis = basisfunctions(0.5,0,0.5,(fa->getBoundary()[0])->getSource(),(fa->getBoundary()[1])->getSource(),(fa->getBoundary()[2])->getSource(),fa,refmeshfaces);
                loadvector(i,0) = loadvector(i,0) + 1.0/6.0*basis(0,i)*functionvalue*math::sqrt(math::abs(Gx1*Gy2-Gx2*Gy1));
            }
        }
        std::cout << "\n" << loadvector(i,0);

    }
    return loadvector;
}

template<class T>
void gsHeMesh<T>::setBoxsplineTable()
{
    phi.push_back(0.328125);
    phi.push_back(0.328125);
    phi.push_back(0.13542);
    phi.push_back(0.015625);
    phi.push_back(0.00521);
    phi.push_back(0.015625);
    phi.push_back(0.13542);
    phi.push_back(0.015625);
    phi.push_back(0.00521);
    phi.push_back(0.015625);
    phi.push_back(0);
    phi.push_back(0);

    phi1.push_back(-0.54167);
    phi1.push_back(0.54167);
    phi1.push_back(0);
    phi1.push_back(-0.08333);
    phi1.push_back(-0.04167);
    phi1.push_back(-0.08333);
    phi1.push_back(0);
    phi1.push_back(0.08333);
    phi1.push_back(0.04167);
    phi1.push_back(0.08333);
    phi1.push_back(0);
    phi1.push_back(0);

    phi2.push_back(-0.27083);
    phi2.push_back(0.27083);
    phi2.push_back(0.375);
    phi2.push_back(0.02083);
    phi2.push_back(-0.02083);
    phi2.push_back(-0.10416);
    phi2.push_back(-0.375);
    phi2.push_back(-0.02083);
    phi2.push_back(0.02083);
    phi2.push_back(0.10416);
    phi2.push_back(0);
    phi2.push_back(0);
}

template<class T>
std::vector<T> gsHeMesh<T>::transformPoint(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfVertexHandle veval, std::vector<gsHalfFaceHandle> f, std::list<int> indices)
{
  //  gsHalfVertexHandle result = new gsHalfVertex<T>();
    std::vector<T> resultcoordinates;
    gsHalfVertexHandle tempcoordinate1 = NULL;
    gsHalfVertexHandle tempcoordinate2 = NULL;
    gsHalfVertexHandle tempcoordinate3 = NULL;
    gsHalfEdgeHandle controledge = veval->getHalfEdge();
    int type = getVertexType(veval,indices);
    int counter = 0;

    if(bary1+bary2+bary3 < 1.1 && bary1+bary2+bary3 > 0.9 ) //property for being a barycentric coordinate with treat of rounding errors (very high tolerance)
    {
        while((!controledge->getNext()->getSource()->equal(v1)) && (!controledge->getNext()->getSource()->equal(v2)) && (!controledge->getNext()->getSource()->equal(v3)))
        {
            controledge = controledge->getMate()->getNext();
            if(counter<6)
                counter++;
            else
            {
                if(counter==6)
                {
                    controledge = veval->getHalfEdge()->getNext();
                    counter++;
                }
                if(counter>6 && counter <12)
                    counter++;
                else
                {
                    counter = 0;
                    controledge = veval->getHalfEdge()->getNext()->getNext();
                }
            }
        }
        if(type == 0)
        {
      //      gsHalfFaceHandle fa = tempmesh.getHalfFaceAtIndex(10); //maybe not necessary
            gsHalfFaceHandle fa = f[10];
            std::vector<gsHalfVertexHandle> tempvertices = fa->getVertices();
            for ( typename std::vector<gsHalfVertexHandle>::iterator
                      it = tempvertices.begin(); it!= tempvertices.end(); ++it)
            {
                gsHalfVertexHandle tempv = *it;
                if(!veval->equal(v1)&& !veval->equal(v2) && !veval->equal(v3))
                {
                    if(controledge->getNext()->getSource()->equal(v1)) //no rotation
                    {
                        if(tempv->getId() == 1) //(2,2) old Ids: (9,5,10); new Ids: (1,2,3)
                        {
                            tempcoordinate1 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate2 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate3 = tempv;
                        }
                     //   std::cout << "Type0 - Number 1" << "\n";
                    }
                    if(controledge->getNext()->getSource()->equal(v2)) //1 left rotation
                    {
                        if(tempv->getId() == 1) //(2,2)
                        {
                            tempcoordinate3 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate1 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate2 = tempv;
                        }
                     //   std::cout << "Type0 - Number 2" << "\n";
                    }
                    if(controledge->getNext()->getSource()->equal(v3)) // 2 left rotations
                    {
                        if(tempv->getId() == 1) //(2,2)
                        {
                            tempcoordinate2 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate3 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate1 = tempv;
                        }
                     //   std::cout << "Type0 - Number 3" << "\n";
                    }
                }
                else
                {
                    if(veval->equal(v1))
                    {
                        if(tempv->getId() == 1) //(2,2)
                        {
                            tempcoordinate1 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate2 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate3 = tempv;
                        }
                   //     std::cout << "Type0 - Number 1" << "\n";
                    }
                    if(veval->equal(v2))
                    {
                        if(tempv->getId() == 1) //(2,2)
                        {
                            tempcoordinate3 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate1 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate2 = tempv;
                        }
                    //    std::cout << "Type0 - Number 2" << "\n";
                    }
                    if(veval->equal(v3))
                    {
                        if(tempv->getId() == 1) //(2,2)
                        {
                            tempcoordinate2 = tempv;
                        }
                        if(tempv->getId() == 2) //(2,1)
                        {
                            tempcoordinate3 = tempv;
                        }
                        if(tempv->getId() == 3) //(3,2)
                        {
                            tempcoordinate1 = tempv;
                        }
                    //    std::cout << "Type0 - Number 3" << "\n";
                    }
                }
            }
            resultcoordinates.push_back(bary1*tempcoordinate1->x()+bary2*tempcoordinate2->x()+bary3*tempcoordinate3->x());
            resultcoordinates.push_back(bary1*tempcoordinate1->y()+bary2*tempcoordinate2->y()+bary3*tempcoordinate3->y());
            resultcoordinates.push_back(bary1*tempcoordinate1->z()+bary2*tempcoordinate2->z()+bary3*tempcoordinate3->z());
            return(resultcoordinates);
        //    std::cout << "Type 0: X: " << result->x() << " Y: " << result->y() << "\n";
        }

        if(type == 1)
        {
            //gsHalfFaceHandle fa = tempmesh.getHalfFaceAtIndex(4);
            gsHalfFaceHandle fa = f[4];
            std::vector<gsHalfVertexHandle> tempvertices = fa->getVertices();
            for ( typename std::vector<gsHalfVertexHandle>::iterator
                      it = tempvertices.begin(); it!= tempvertices.end(); ++it)
            {
                gsHalfVertexHandle tempv = *it;
                if(controledge->getNext()->getSource()->equal(v1))
                {
                    if(tempv->getId() == 2) //(2,1) old Ids: (5,4,6); new Ids: (2,10,11)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 11) //(3,1)
                    {
                        tempcoordinate3 = tempv;
                    }
                   // std::cout << "Type1 - Number 1" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v3))
                {
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 11) //(3,1)
                    {
                        tempcoordinate1 = tempv;
                    }
                 //   std::cout << "Type1 - Number 2" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v2))
                {
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 11) //(3,1)
                    {
                        tempcoordinate2 = tempv;
                    }
                 //   std::cout << "Type1 - Number 3" << "\n";
                }
            }
            resultcoordinates.push_back(bary1*tempcoordinate1->x()+bary2*tempcoordinate2->x()+bary3*tempcoordinate3->x());
            resultcoordinates.push_back(bary1*tempcoordinate1->y()+bary2*tempcoordinate2->y()+bary3*tempcoordinate3->y());
            resultcoordinates.push_back(bary1*tempcoordinate1->z()+bary2*tempcoordinate2->z()+bary3*tempcoordinate3->z());
            return(resultcoordinates);
        }

        if(type == 2)
        {
            //gsHalfFaceHandle fa = tempmesh.getHalfFaceAtIndex(9);
            gsHalfFaceHandle fa = f[9];
            std::vector<gsHalfVertexHandle> tempvertices = fa->getVertices();
            for ( typename std::vector<gsHalfVertexHandle>::iterator
                      it = tempvertices.begin(); it!= tempvertices.end(); ++it)
            {
                gsHalfVertexHandle tempv = *it;
                if(controledge->getNext()->getSource()->equal(v1))
                {
                    if(tempv->getId() == 11) //(3,1) old Ids: (6,10,5); new Ids: (11,3,2)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 3) //(3,2)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate1 = tempv;
                    }
                 //   std::cout << "Type2 - Number 1" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v3))
                {
                    if(tempv->getId() == 11) //(3,1)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 3) //(3,2)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate2 = tempv;
                    }
                 //   std::cout << "Type2 - Number 2" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v2))
                {
                    if(tempv->getId() == 11) //(3,1)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 3) //(3,2)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate3 = tempv;
                    }
                 //   std::cout << "Type2 - Number 3" << "\n";
                }
            }
            resultcoordinates.push_back(bary1*tempcoordinate1->x()+bary2*tempcoordinate2->x()+bary3*tempcoordinate3->x());
            resultcoordinates.push_back(bary1*tempcoordinate1->y()+bary2*tempcoordinate2->y()+bary3*tempcoordinate3->y());
            resultcoordinates.push_back(bary1*tempcoordinate1->z()+bary2*tempcoordinate2->z()+bary3*tempcoordinate3->z());
            return(resultcoordinates);
        }

        else //=type 3
        {
            //gsHalfFaceHandle fa = tempmesh.getHalfFaceAtIndex(2);
            gsHalfFaceHandle fa = f[2];
            std::vector<gsHalfVertexHandle> tempvertices = fa->getVertices();
            for ( typename std::vector<gsHalfVertexHandle>::iterator
                      it = tempvertices.begin(); it!= tempvertices.end(); ++it)
            {
                gsHalfVertexHandle tempv = *it;
                if(controledge->getNext()->getSource()->equal(v1))
                {
                    if(tempv->getId() == 9) //(1,0) old Ids: (1,4,5); new Ids: (9,10,2)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate1 = tempv;
                    }
                 //   std::cout << "Type3 - Number 1" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v3))
                {
                    if(tempv->getId() == 9) //(1,0)
                    {
                        tempcoordinate3 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate2 = tempv;
                    }
                 //   std::cout << "Type3 - Number 2" << "\n";
                }
                if(controledge->getNext()->getSource()->equal(v2))
                {
                    if(tempv->getId() == 9) //(1,0)
                    {
                        tempcoordinate1 = tempv;
                    }
                    if(tempv->getId() == 10) //(2,0)
                    {
                        tempcoordinate2 = tempv;
                    }
                    if(tempv->getId() == 2) //(2,1)
                    {
                        tempcoordinate3 = tempv;
                    }
                 //   std::cout << "Type3 - Number 3" << "\n";
                }
            }
            resultcoordinates.push_back(bary1*tempcoordinate1->x()+bary2*tempcoordinate2->x()+bary3*tempcoordinate3->x());
            resultcoordinates.push_back(bary1*tempcoordinate1->y()+bary2*tempcoordinate2->y()+bary3*tempcoordinate3->y());
            resultcoordinates.push_back(bary1*tempcoordinate1->z()+bary2*tempcoordinate2->z()+bary3*tempcoordinate3->z());
            return(resultcoordinates);
        }
    }
    GISMO_ERROR("No barycentric coordinates as input");
}

template<class T>
void gsHeMesh<T>::buildReferenceConfiguration()
{
    //typedef gsMeshElement<>::gsHeVertexHandle gsVertexHandle;
    //typedef gsMeshElement<>::gsHalfFaceHandle gsHalfFaceHandle;
    std::vector<gsHalfVertexHandle> vh;
    (*this).addVertex(0,0); //0
    (*this).addVertex(1,0); //1
    (*this).addVertex(1,1); //2
    (*this).addVertex(0,1); //3
    (*this).addVertex(2,0); //4
    (*this).addVertex(2,1); //5
    (*this).addVertex(3,1); //6
    (*this).addVertex(0,2); //7
    (*this).addVertex(1,2); //8
    (*this).addVertex(2,2); //9
    (*this).addVertex(3,2); //10
    (*this).addVertex(4,2); //11
    (*this).addVertex(1,3); //12
    (*this).addVertex(2,3); //13
    (*this).addVertex(3,3); //14
    (*this).addVertex(4,3); //15
    (*this).addVertex(2,4); //16
    (*this).addVertex(3,4); //17
    (*this).addVertex(4,4); //18

    vh.push_back((*this).getVertexAtIndex(0));vh.push_back((*this).getVertexAtIndex(1));vh.push_back((*this).getVertexAtIndex(2));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(2));vh.push_back((*this).getVertexAtIndex(3));vh.push_back((*this).getVertexAtIndex(0));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(1));vh.push_back((*this).getVertexAtIndex(4));vh.push_back((*this).getVertexAtIndex(5));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(5));vh.push_back((*this).getVertexAtIndex(2));vh.push_back((*this).getVertexAtIndex(1));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(4));vh.push_back((*this).getVertexAtIndex(6));vh.push_back((*this).getVertexAtIndex(5));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(3));vh.push_back((*this).getVertexAtIndex(2));vh.push_back((*this).getVertexAtIndex(8));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(8));vh.push_back((*this).getVertexAtIndex(7));vh.push_back((*this).getVertexAtIndex(3));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(2));vh.push_back((*this).getVertexAtIndex(5));vh.push_back((*this).getVertexAtIndex(9));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(9));vh.push_back((*this).getVertexAtIndex(8));vh.push_back((*this).getVertexAtIndex(2));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(5));vh.push_back((*this).getVertexAtIndex(6));vh.push_back((*this).getVertexAtIndex(10));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(10));vh.push_back((*this).getVertexAtIndex(9));vh.push_back((*this).getVertexAtIndex(5));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(6));vh.push_back((*this).getVertexAtIndex(11));vh.push_back((*this).getVertexAtIndex(10));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(7));vh.push_back((*this).getVertexAtIndex(8));vh.push_back((*this).getVertexAtIndex(12));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(8));vh.push_back((*this).getVertexAtIndex(9));vh.push_back((*this).getVertexAtIndex(13));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(13));vh.push_back((*this).getVertexAtIndex(12));vh.push_back((*this).getVertexAtIndex(8));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(9));vh.push_back((*this).getVertexAtIndex(10));vh.push_back((*this).getVertexAtIndex(14));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(14));vh.push_back((*this).getVertexAtIndex(13));vh.push_back((*this).getVertexAtIndex(9));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(10));vh.push_back((*this).getVertexAtIndex(11));vh.push_back((*this).getVertexAtIndex(15));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(15));vh.push_back((*this).getVertexAtIndex(14));vh.push_back((*this).getVertexAtIndex(10));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(12));vh.push_back((*this).getVertexAtIndex(13));vh.push_back((*this).getVertexAtIndex(16));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(13));vh.push_back((*this).getVertexAtIndex(14));vh.push_back((*this).getVertexAtIndex(17));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(17));vh.push_back((*this).getVertexAtIndex(16));vh.push_back((*this).getVertexAtIndex(13));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(14));vh.push_back((*this).getVertexAtIndex(15));vh.push_back((*this).getVertexAtIndex(18));
    (*this).addHalfFace(vh);vh.clear();

    vh.push_back((*this).getVertexAtIndex(18));vh.push_back((*this).getVertexAtIndex(17));vh.push_back((*this).getVertexAtIndex(14));
    (*this).addHalfFace(vh);vh.clear();

    (*this).initialize();
    std::vector<gsHalfVertexHandle > ve1 = (*this).getVertices();
    std::vector<gsHalfVertexHandle > ve2;
    for ( typename std::vector<gsHalfVertexHandle>::iterator
          it = ve1.begin(); it!= ve1.end(); ++it)
    {
        gsHalfVertexHandle ve = *it;
        int veId = ve->getId();
        switch(veId)
        {
            case 0: ve->setId(0); break;
            case 1: ve->setId(9); break;
            case 2: ve->setId(8); break;
            case 3: ve->setId(18); break;
            case 4: ve->setId(10); break;
            case 5: ve->setId(2); break;
            case 6: ve->setId(11); break;
            case 7: ve->setId(17); break;
            case 8: ve->setId(7); break;
            case 9: ve->setId(1); break;
            case 10: ve->setId(3); break;
            case 11: ve->setId(12); break;
            case 12: ve->setId(16); break;
            case 13: ve->setId(6); break;
            case 14: ve->setId(5); break;
            case 15: ve->setId(4); break;
            case 16: ve->setId(15); break;
            case 17: ve->setId(14); break;
            case 18: ve->setId(13); break;
        }
        ve2.push_back(ve);
    }
    (*this).setVertices(ve2);
}

template<class T>
void gsHeMesh<T>::getAllEdgesPrinted()
{
    int i = 1;
    for ( typename std::vector<gsHalfEdgeHandle>::iterator
          it = edge.begin(); it!= edge.end(); ++it)
    {
        gsHalfEdgeHandle heh = *it;
        std::cout << i << ") From x: " << heh->getSource()->x() << " y: " << heh->getSource()->y()  << "\n";
        std::cout << "to x: " << heh->getNext()->getSource()->x() << " y: " << heh->getNext()->getSource()->y()  << "\n";
        i++;
    }
}

template<class T>
void gsHeMesh<T>::Testprint()
{
    int i = 1;

        gsHalfEdgeHandle heh = edge[6];
        gsHalfEdgeHandle control = heh;
        do
        {
            try
            {
            heh = heh->getNext();
            std::cout << i << ") From x: " << heh->getSource()->x() << " y: " << heh->getSource()->y()  << "\n";
            std::cout << "to x: " << heh->getNext()->getSource()->x() << " y: " << heh->getNext()->getSource()->y()  << "\n";
            i++;
            if(i>100)
                break;
            }
            catch(int e)
            {
                std::cout << "Error at " << control->getId() << " Errortyp:" << e;
                break;
            }
        }
        while(control != heh);

}

/*
template<class T>
typename gsHeMesh<T>::gsHalfEdgeHandle gsHeMesh<T>::changeOrientation(typename gsHeMesh<T>::gsHalfEdgeHandle heh)
{
    gsHalfEdgeHandle tempHe = *heh;
    heh->setSource(heh->getNext()->getSource());
    while(tempHe->getNext())
    {
        tempHe = tempHe->getNext();
    }
    heh->setNext(tempHe);
    return heh;
} */

/*
template<class T>
std::vector<typename gsHeMesh<T>::gsHalfEdgeHandle> gsHeMesh<T>::changeOrientation(std::vector<gsHalfEdgeHandle> heh)
{
    std::vector<gsHalfEdgeHandle> handles;
    for ( typename std::vector<gsHalfEdgeHandle>::iterator
          it = heh.begin(); it!= heh.end(); ++it)
    {
        gsHalfEdgeHandle tempHe = *it;
        gsHalfVertexHandle vh = tempHe->getSource();
        tempHe->setSource(tempHe->getNext()->getSource());
        tempHe->getNext()->setSource(vh);
        handles.push_back(tempHe);
    }
    return handles;
}

template<class T>
bool gsHeMesh<T>::checkRightBoundaryOrientation()
{
    gsHalfEdgeHandle heh = globalboundary[0];
    for ( typename std::vector<gsHalfEdgeHandle>::iterator
          it = edge.begin(); it!= edge.end(); ++it)
    {
        gsHalfEdgeHandle tempHe = *it;
        if(heh->getSource() == tempHe->getNext()->getSource()&& heh->getNext()->getSource() == tempHe->getSource())
            return true;
    }
    return false;
 }
*/
}
