
#include <gsIO/gsLsdynaUtils.h>


namespace gismo
{

namespace lsdyna    
{

template <typename T>    
class gsLsdynaGeometry
{
public:   
    gsLsdynaGeometry(const gsGeometry<T>& geom)
        :
        mGeom(geom),
        materialId(1),
        massDensity(1),
        youngsModulus(1e+7),
        poissonsRatio(3e-1),
        shearCorrection(1),
        numIntPoints(2),
        shellThickness(5e-2)
    {
        
    }

    ~gsLsdynaGeometry()
    {
        
    }


    // --------------------------------------------------
    // getters
    // --------------------------------------------------

    const gsGeometry<T>& getGeometry() const
    {
        return mGeom;
    }

    lsInt getMaterialId() const
    {
        return materialId;
    }

    lsFloat getMassDensity() const
    {
        return massDensity;
    }

    lsFloat getYoungsModulus() const
    {
        return youngsModulus;
    }

    lsFloat getPoissonsRatio() const
    {
        return poissonsRatio;
    }

    lsFloat getShearCorrection() const
    {
        return shearCorrection;
    }

    lsInt getNumIntegrationPoints() const
    {
        return numIntPoints;
    }

    lsFloat getShellThickness() const
    {
        return shellThickness;
    }
    

    // --------------------------------------------------
    // setters
    // --------------------------------------------------
    
    void setMaterialId(const lsInt id)
    {
        materialId = id;
    }

    void setMassDensity(const lsFloat density)
    {
        massDensity = density;
    }

    void setYoungsModulus(const lsFloat modulus)
    {
        youngsModulus = modulus;
    }

    void setPoissonsRatio(const lsFloat ratio)
    {
        poissonsRatio = ratio;
    }

    void setShearCorrection(const lsFloat correction)
    {
        shearCorrection = correction;
    }

    void setNumIntegrationPoints(const lsInt points)
    {
        numIntPoints = points;
    }

    void setShellThickness(const lsFloat thickness)
    {
        shellThickness = thickness;
    }


    // --------------------------------------------------
    // data memebers
    // --------------------------------------------------
private:

    // Input geometry
    const gsGeometry<T>& mGeom;

    // material id
    lsInt materialId;

    // mass density
    lsFloat massDensity;

    // youngs modulus
    lsFloat youngsModulus;

    // poissons ratio
    lsFloat poissonsRatio;

    // shear correction
    lsFloat shearCorrection;

    // number of integration points
    lsInt numIntPoints;

    // shell thickness
    lsFloat shellThickness;
};

    
template <typename T>
class gsLsdynaFEMShell : public gsLsdynaGeometry<T>
{
public:

    // --------------------------------------------------
    // constructors - destructors
    // --------------------------------------------------

    gsLsdynaFEMShell(const gsGeometry<T>& geom)
        :
        gsLsdynaGeometry<T>(geom),
        partId(1001),
        elementForm(16)
    {
        
    }

    ~gsLsdynaFEMShell()
    {

    }

public:

    void write(std::ostream& out,
               const gsVector<index_t, 2>& gridDimension) const
    {
        setGlobalState(out);

        write::DEFINE_CURVE(out, 123, 0, 5e-3, 0, 1);

        write::MAT_ELASTIC(out, this->getMaterialId(), this->getMassDensity(),
                           this->getYoungsModulus(), this->getPoissonsRatio());

        write::PART(out, partId, partId, this->getMaterialId());

        write::SECTION_SHELL(out, partId, elementForm,
                             this->getShearCorrection(),
                             this->getNumIntegrationPoints(),
                             this->getShellThickness());

        writeNodes(out, gridDimension);

        writeElementShell(out, gridDimension);
    }

    void writeNodes(std::ostream& out,
                    const gsVector<index_t, 2>& gridDimension) const
    {
        index_t numGridPts = gridDimension.prod();

        gsMatrix<T> points(3, numGridPts);
        gsVector<bool> onBoundary(numGridPts);
        
        gsMatrix<T> pt(3, 1); // temporary variable
        gsMatrix<T> param(2, 1); // temporary variable

        // last index in given direction
        const index_t lastIndex0 = gridDimension(0) - 1;
        const index_t lastIndex1 = gridDimension(1) - 1;

        gsGridIterator<index_t,CUBE,2> it(gridDimension);
        index_t flatIndex = 0;
        for (; it; ++it, ++flatIndex)
        {
            const gsVector<index_t, 2>& index = *it;

            // check if it is on boundary
            if (index(0) == 0 || index(0) == lastIndex0 ||
                index(1) == 0 || index(1) == lastIndex1)
            {
                onBoundary(flatIndex) = true;
            }
            else
            {
                onBoundary(flatIndex) = false;
            }

            param(0, 0) = (index(0) * 1.0) / lastIndex0;
            param(1, 0) = (index(1) * 1.0) / lastIndex1;
            
            this->getGeometry().eval_into(param, pt);

            points.col(flatIndex) = pt.col(0);
        }

        write::NODES(out, points, onBoundary, 1000001);
    }

    void writeElementShell(std::ostream& out,
                           const gsVector<index_t, 2>& gridDim) const
    {
        gsVector<index_t, 2> elementDim;
        elementDim(0) = gridDim(0) - 1;
        elementDim(1) = gridDim(1) - 1;

        gsVector<index_t> indices(4);
        gsGridIterator<index_t,CUBE,2> it(elementDim);

        index_t flatIndex = 0;
        for (; it; ++it, ++flatIndex)
        {
            gsVector<index_t, 2> index = *it;

            indices(0) = getFlatIndex(index, gridDim);
            index(0) += 1;

            indices(1) = getFlatIndex(index, gridDim);
            index(1) += 1;

            indices(2) = getFlatIndex(index, gridDim);
            index(0) -= 1;

            indices(3) = getFlatIndex(index, gridDim);

            write::ELEMENT_SHELL(out, 200001 + flatIndex, partId,
                                 indices, 1000001);
        }
    }

    static
    index_t getFlatIndex(const gsVector<index_t, 2>& index,
                         const gsVector<index_t, 2>& gridDim) 
    {
        return gridDim(0) * index(1) + index(0);
    }
    
    // --------------------------------------------------
    // setters
    // --------------------------------------------------

    void setElementForm(const lsInt form)
    {
        elementForm = form;
    }

private:

    void setGlobalState(std::ostream& out) const
    {
        out << std::uppercase;
        out << std::setprecision(16);
        out << std::scientific;
    }

    // disable copy constructor and assignment operator
    gsLsdynaFEMShell(const gsLsdynaFEMShell& other);
    const gsLsdynaFEMShell operator=(const gsLsdynaFEMShell& other);

    // --------------------------------------------------
    // data memebers
    // --------------------------------------------------
private:

    // part id
    lsInt partId;

    // element formulation
    lsInt elementForm;
};


template <typename T>
class gsLsdynaIGAShell: public gsLsdynaGeometry<T>
{
public:
    // --------------------------------------------------
    // constructors - destructors
    // --------------------------------------------------

    gsLsdynaIGAShell(const gsGeometry<T>& geom)
        :
        gsLsdynaGeometry<T>(geom),
        numElements(),
        startPartId(1001),
        startInterpolPartId(8001),
        startNodeId(1000001),
        lumpingMassMatOption(0),
        shellFormulation(1)
    {
        numElements = geom.basis().numElements();
    }

    ~gsLsdynaIGAShell()
    {

    }

    // --------------------------------------------------
    // data memebers
    // --------------------------------------------------
public:
    void write(std::ostream& out) const
    {
        setGlobalState(out);

        write::DEFINE_CURVE(out, 123, 0, 5e-3, 0, 1);

        writePartAndSection(out);

        writeInterpolationPartAndSection(out);

        writeNodes(out);

        writeGeneralizedShell(out);
    }

private:
    void setGlobalState(std::ostream& out) const
    {
        out << std::uppercase;
        out << std::setprecision(16);
        out << std::scientific;
    }

    void writePartAndSection(std::ostream& out) const
    {
        write::MAT_ELASTIC(out, this->getMaterialId(), this->getMassDensity(),
                           this->getYoungsModulus(), this->getPoissonsRatio());

        for (index_t el = 0; el != numElements; el++)
        {
            lsInt partId = startPartId + el;

            write::PART(out, partId, partId, this->getMaterialId());

            write::SECTION_SHELL(out, partId, partId,
                                 this->getShearCorrection(),
                                 this->getNumIntegrationPoints(),
                                 this->getShellThickness());
        }
    }

    void writeInterpolationPartAndSection(std::ostream& out) const
    {
        write::MAT_ELASTIC(out, 2, 1e-10, 1e-10, 1e-10);

        write::PART(out, startInterpolPartId, startInterpolPartId, 2);

        // 98 means interpolation shell
        write::SECTION_SHELL(out, startInterpolPartId, 98,
                             this->getShearCorrection(),
                             this->getNumIntegrationPoints(),
                             this->getShellThickness());
    }

    void writeNodes(std::ostream& out) const
    {
        gsMatrix<T> coefs = this->getGeometry().coefs();
        coefs.transposeInPlace();
        gsMatrix<index_t> pBoundary = this->getGeometry().basis().allBoundary();

        gsVector<bool> onBoundary(coefs.cols());
        onBoundary.setZero();
        
        for (index_t row = 0; row != pBoundary.rows(); row++)
        {
            onBoundary((pBoundary)(row)) = true;
        }

        write::NODES(out, coefs, onBoundary, startNodeId);
    }

    void writeGeneralizedShell(std::ostream& out) const
    {
        const gsBasis<T>& basis = this->getGeometry().basis();
        gsDomainIterator<T>* pDomIt = basis.makeDomainIterator().release();

        gsVector<index_t> quadDimension(2); // dimension of quadrature points
        quadDimension[0] = basis.degree(0) + 1;
        quadDimension[1] = basis.degree(1) + 1;

        gsGaussRule<T> quRule( quadDimension );
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;

        gsFuncData<T> bdata(NEED_DERIV|NEED_ACTIVE);

        for (lsInt elementId = 0; pDomIt->good(); pDomIt->next(), elementId++)
        {
            quRule.mapTo(pDomIt->lowerCorner(), pDomIt->upperCorner(), quNodes, quWeights);
            basis.compute(quNodes, bdata);

            writeElementGeneralizedShell(out, elementId, bdata, quWeights);
        }

        delete pDomIt;
    }

    void writeElementGeneralizedShell(std::ostream& out,
                                      const lsInt elementId,
                                      const gsFuncData<T> & bdata,
                                      const gsVector<T> & quWeights) const
    {
        write::DEFINE_ELEMENT_GENERALIZED_SHELL_header(
            out, startPartId + elementId, bdata.allValues().cols(),
            bdata.actives.rows(), lumpingMassMatOption, shellFormulation);


        write::DEFINE_ELEMENT_GENERALIZED_SHELL_values<T>(
            out, quWeights, bdata.values[0], bdata.values[1]);

        checkValuesAndDerivatives(bdata.values[0], bdata.values[1]);
    }

    

    void checkValuesAndDerivatives(const gsMatrix<T>& values,
                                   const gsMatrix<T>& derivs) const
    {
        for (index_t k = 0; k != values.cols(); k++)
        {
            if (1e-10 < math::abs(values.col(k).sum() - 1))
            {
                gsWarn << "& ** Warning **\n"
                       << "& at gauss point: " << k + 1 << "\n"
                       << "& sum of values of basis functions in not equal to 1\n"
                       << "& sum is: " << values.col(k).sum() << "\n";
            }

            gsMatrix<T> der = derivs.col(k);
            der.resize(2, values.rows());

            for (index_t u = 0; u != der.rows(); ++u)
            {
                if (1e-10 < math::abs(der.row(u).sum()))
                {
                    std::string dir = (u == 0) ? "x" : "y";

                    gsWarn << "& ** Warning **\n"
                           << "& at gauss point: " << k + 1 << "\n"
                           << "& sum of derivative is not equal to 0 in " << dir
                           << " direction\n"
                           << "& sum is: " << der.row(u).sum() << "\n";
                }
            }


        }
    }
    
    // --------------------------------------------------
    // data memebers
    // --------------------------------------------------
private:

    // number of elements
    index_t numElements;

    // beggining of the part id
    lsInt startPartId;

    // beggining of the interpolation part id
    lsInt startInterpolPartId;

    // beggining of the node id
    lsInt startNodeId;

    // Option for lumping of mass matrix:
    lsInt lumpingMassMatOption;

    // Shell formulation to be used
    lsInt shellFormulation;
};
    

} // namespace lsdyna

} // namespace gismo




