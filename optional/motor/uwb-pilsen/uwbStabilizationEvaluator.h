/** @file uwbStabilizationEvaluator.h

    Author(s): E. Turnerova

*/

#pragma once

namespace gismo
{

template <class T>
class uwbStabilizationEvaluator
{
public:
    uwbStabilizationEvaluator()
    { }

    ~uwbStabilizationEvaluator()
    { }

    void initialize(index_t nPoints, index_t dim)
    {
        m_nPoints = nPoints;
        m_dim = dim;
        initMembers();
    }

    void initMembers()
    {
        m_diffusionCoeff.setZero(m_nPoints);
        m_advectionCoeff.resize(m_nPoints);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s].setZero(m_dim);
        m_reactionCoeff.setZero(m_nPoints);
        m_deg = 0;
        m_timeStep = 0.;
        m_elemDiam = 0.;
        m_alpha = 1.;
        m_srbavScaleFactorH = 1.;
        m_bTauDeg = false;
    }

    void initAtElement(gsVector<T> advectionCoeff, T diffusionCoeff, bool tauDeg = false)
    {
        m_diffusionCoeff.setConstant(m_nPoints, diffusionCoeff);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionCoeff;
        m_reactionCoeff.setZero(m_nPoints);
        m_bTauDeg = tauDeg;
    }

    void initAtElement(gsVector<T> advectionCoeff, T diffusionCoeff, T reactionCoeff, bool tauDeg = false)
    {
        m_diffusionCoeff.setConstant(m_nPoints, diffusionCoeff);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionCoeff;
        m_reactionCoeff.setConstant(m_nPoints, reactionCoeff);
        m_bTauDeg = tauDeg;
    }

    void initAtElement(gsMatrix<T>& advectionSolVals, bool tauDeg = false)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff.setZero(m_nPoints);
        m_reactionCoeff.setZero(m_nPoints);
        m_bTauDeg = tauDeg;
    }

    void initAtElement(gsMatrix<T>& advectionSolVals, gsMatrix<T>& diffusionSolVals, bool tauDeg = false)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff = diffusionSolVals.row(0);
        m_reactionCoeff.setZero(m_nPoints);
        m_bTauDeg = tauDeg;
    }

public:
    void setElemLength(T elemDiam, T advElemDiam, gsVector<T> diffElemDiam)
    {
        if (elemDiam <= 0.)
            GISMO_ERROR("'elemDiam' equals zero. It is not set in the stabilization evaluator!");
        m_elemDiam = elemDiam;
        m_advElemLength.resize(diffElemDiam.size());
        for(int var = 0; var < diffElemDiam.size(); var++)
            m_advElemLength[var].setConstant(m_nPoints, advElemDiam);
        m_diffElemLength.resize(diffElemDiam.size());
        for(int var = 0; var < diffElemDiam.size(); var++)
            m_diffElemLength[var].setConstant(m_nPoints, diffElemDiam(var));
    }

    void setElemLength(T elemDiam, std::vector<gsVector<T> > advElemDiam, std::vector<gsVector<T> > diffElemDiam)
    {
        if (elemDiam <= 0.)
            GISMO_ERROR("'elemDiam' equals zero. It is not set in the stabilization evaluator!");
        m_elemDiam = elemDiam;
        m_advElemLength = advElemDiam;
        m_diffElemLength = diffElemDiam;
    }

    void setElemLength(T elemDiam, std::vector<gsVector<T> > advElemDiam)
    {
        if (elemDiam <= 0.)
            GISMO_ERROR("'elemDiam' equals zero. It is not set in the stabilization evaluator!");
        m_elemDiam = elemDiam;
        m_advElemLength = advElemDiam;
        m_diffElemLength = advElemDiam;
    }

    void setSUPGvars(index_t tauStabType, index_t deg, T timeStep)
    {
        m_tauStabType = tauStabType;
        m_deg = deg;
        m_timeStep = timeStep;
    }

    void setTauType(index_t tauType) { m_tauStabType = tauType; }

    void setCrosswindVars(index_t crosswindType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual)
    {
        m_crosswindType = crosswindType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
    }

    void setSRBAVvars(index_t srbavType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T alpha = 1., T srbavScaleFactorH = 1.)
    {
        m_srbavType = srbavType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_alpha = alpha;
        m_srbavScaleFactorH = srbavScaleFactorH;
    }

    void setIsoADvars(index_t isoType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T alpha = 1.)
    {
        m_isoType = isoType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_alpha = alpha;
    }

    void setIsoADVars(std::vector< gsVector<T> >& residual, std::vector< gsMatrix<T> >& solGrads)
    {
        m_residual = residual;
        m_solGrads = solGrads;
    }

    void setADVars(index_t tauType, index_t deg, T timeStep)
    {
        m_tauStabType = tauType;
        m_deg = deg;
        m_timeStep = timeStep;
    }

    const gsVector<T> & getDiffusionCoefficient() const { return m_diffusionCoeff; }
    const T getDiffusionCoefficient(const int i) const { return m_diffusionCoeff(i); }
    const std::vector<gsVector<T> > & getAdvectionCoefficient() const { return m_advectionCoeff; }
    const gsVector<T> getAdvectionCoefficient(const int i) const { return m_advectionCoeff[i]; }
    const gsVector<T> & getReactionCoefficient() const { return m_reactionCoeff; }
    const T getReactionCoefficient(const int i) const { return m_reactionCoeff(i); }

    T getTauS(index_t k, T diffCoeff = -1., T reactCoeff = 0., index_t var = 0)
    {
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);

        m_reactionCoeff(k) = math::max(reactCoeff, 0.);

        T tau_s = 0.;
        T normAdvectionCoeff = m_advectionCoeff[k].norm();
        T P = 0.;

        int deg = 1;
        if (m_bTauDeg)
            deg = m_deg;

        switch (m_tauStabType)
        {
        case 0:
            P = m_advElemLength[var](k) * normAdvectionCoeff / (2. * m_diffusionCoeff(k));
            tau_s = (m_advElemLength[var](k) / (2. * deg * normAdvectionCoeff)) * (1. / tanh(P) - 1. / P);

            if (P == 0. || normAdvectionCoeff == 0.)
                tau_s = 0.;

          break;

        case 1:
            tau_s = m_advElemLength[var](k) / (2. * deg * normAdvectionCoeff);

            if (normAdvectionCoeff == 0.)
                tau_s = 0.;
            break;

        case 2:
            tau_s = 1. / (math::sqrt(math::pow(2. * deg * normAdvectionCoeff / m_advElemLength[var](k), 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemLength[var](k), 2), 2)));
            break;

        case 3:
            tau_s = 1. / (math::sqrt(math::pow(m_reactionCoeff(k), 2)
                  + math::pow(2. * deg * normAdvectionCoeff / m_advElemLength[var](k), 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemLength[var](k), 2), 2)));
            break;

        case 4:
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * deg * normAdvectionCoeff / m_advElemLength[var](k), 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemLength[var](k), 2), 2)));
            break;

        case 5:
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep + m_reactionCoeff(k), 2)
                  + math::pow(2. * deg * normAdvectionCoeff / m_advElemLength[var](k), 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemLength[var](k), 2), 2)));
            break;

        case 6: //Codina+deg
            //GISMO_ASSERT(m_deg > 0 && reactCoeff >= 0., "'m_deg equals zero or reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (2.* deg * normAdvectionCoeff / m_advElemLength[var](k) +
                    4. * math::abs(m_diffusionCoeff(k)) / math::pow(m_diffElemLength[var](k), 2) +
                    math::abs(1. / m_timeStep + m_reactionCoeff(k)));
            break;

        case 7: //KLR+deg
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = math::min(math::min(m_advElemLength[var](k) / (2. * deg * normAdvectionCoeff),
                                        1. / math::abs(1. / m_timeStep + m_reactionCoeff(k))),
                              math::pow(m_diffElemLength[var](k), 2) / (math::pow(deg, 4) * math::abs(m_diffusionCoeff(k))));
            break;

        default: //default tau            
            gsInfo << "Wrong or no type of the stabilization parameter set in the TM visitor. Default tau is used!\n";

            tau_s = m_advElemLength[var](k) / (2. * deg * normAdvectionCoeff);

            if (normAdvectionCoeff == 0.)
                tau_s = 0.;
        }

        return tau_s;
    }

    T getCrosswindStabParam(index_t k, index_t variable = 0, T diffCoeff = -1., T reactCoeff = 0.)
    {
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);

        T cwStabParam = 0.;

        T normAdvection = m_advectionCoeff[k].norm();
        T denom;
        //T C;

        switch (m_crosswindType) {
          case 0: // John and Knobloch 2005
            denom = normAdvection * (m_solGrads[k].row(variable)).norm() + math::abs(m_residual[variable](k));
            cwStabParam = getTauS(k, m_diffusionCoeff(k), reactCoeff, variable) * math::pow(normAdvection, 2) * math::abs(m_residual[variable](k)) /
                          denom;
            if (denom == 0.)
                cwStabParam = 0.;

            break;

        case 1:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          cwStabParam = (m_elemDiam / 2.) * math::abs(m_residual[variable](k)) / m_solGrads[k].row(variable).norm();
          if (m_solGrads[k].row(variable).norm() == 0.)
              cwStabParam = 0.;

          break;

        case 2:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          cwStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(m_residual[variable](k)), 2) / m_solGrads[k].row(variable).norm();
          if (m_solGrads[k].row(variable).norm() == 0.)
              cwStabParam = 0.;

          break;

        /*case 4:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          cwStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(m_residual[variable](k)), 2);

          break;

        case 5:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          cwStabParam = math::pow(m_elemDiam / 2., 2) *
                  math::abs(m_residual[variable](k)) / math::abs(m_advectionCoeff[k](variable));
          if (m_advectionCoeff[k](variable) == 0.)
              cwStabParam = 0.;

          break;

        case 6:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          //T elemdiam = m_elemDiam;

          cwStabParam = (m_elemDiam / 2.) *
                  math::pow(math::tanh(m_residual[variable](k)), 2) / normAdvection;
          if (normAdvection == 0.)
              cwStabParam = 0.;

          break;*/

        default:
            GISMO_ERROR("Wrong crosswind type set in ADR evaluator.");
            break;
        }

        return cwStabParam;
    }

    gsMatrix<T> getCrosswindProjection(index_t k)
    {
        gsMatrix<T> mIdentity;
        gsVector<T> advection = m_advectionCoeff[k];
        gsMatrix<T> tensorProduct(m_dim, m_dim);
        if (m_dim == 2)
            tensorProduct << advection(0)*advection(0), advection(0)*advection(1),
                             advection(1)*advection(0), advection(1)*advection(1);
        else if (m_dim == 3)
            tensorProduct << advection(0)*advection(0), advection(0)*advection(1), advection(0)*advection(2),
                             advection(1)*advection(0), advection(1)*advection(1), advection(1)*advection(2),
                             advection(2)*advection(0), advection(2)*advection(1), advection(2)*advection(2);
        else
            GISMO_ERROR("Wrong dimension in ADR evaluator. Computation possible only in 2D or 3D.");

        mIdentity.setIdentity(m_dim, m_dim);
        gsMatrix<T> projMatrix;
        if (advection.norm() == 0)
            projMatrix.setZero(m_dim, m_dim);
        else
            projMatrix = mIdentity - tensorProduct/math::pow(advection.norm(), 2);

        return projMatrix;
    }

    T getSRBAVstabParam(index_t k, index_t variable = 0, T diffCoeff = -1., T reactCoeff = 0.)
    {
        /*if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);
        if (reactCoeff < 0.)
            m_reactionCoeff(k) = math::max(m_reactionCoeff(k), 0.);
        else
            m_reactionCoeff(k) = math::max(reactCoeff, 0.);*/

        T normAdvectionCoeff = m_advectionCoeff[k].norm();

        //==============================
        //T L_step = 0.0127;//8. * 0.0127;
        //T L_step = 1.;
        //T L_circleDiam = 2*0.05;
        //L_step = L_circleDiam;
        //T L_step = 1.;
        //m_srbavScaleFactorH = 1./L_step;
        //==============================

        T srbavStabParam = 0.;
        switch (m_srbavType)
        {
          case 0: //srbav h^alpha (tanh-CSD)
            srbavStabParam = math::pow(m_srbavScaleFactorH * m_advElemLength[variable](k), m_alpha) * getTauS(k, diffCoeff, reactCoeff, variable) * math::pow(math::tanh(m_residual[variable](k)), 2);
            //srbavStabParam = math::pow(m_elemDiam, m_alpha) * getTauS(k, diffCoeff, reactCoeff, variable) * math::pow(math::tanh(m_residual[variable](k)), 2);

          break;

          case 1: //Nazarov
            srbavStabParam = math::pow(m_srbavScaleFactorH * m_advElemLength[variable](k), m_alpha) * math::abs(m_residual[variable](k)) / math::pow(normAdvectionCoeff, 2);
            //srbavStabParam = math::pow(m_elemDiam, m_alpha) * math::abs(m_residual[variable](k)) / math::pow(normAdvectionCoeff, 2);
            if (normAdvectionCoeff == 0.)
                srbavStabParam = 0.;

          break;

          /*case 2: //Nazarov + tanh^2(Res)
            srbavStabParam = math::pow(m_advElemLength[variable](k) / m_srbavScaleFactorH, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2) / math::pow(normAdvectionCoeff, 2);
            //srbavStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2) / math::pow(normAdvectionCoeff, 2);
            if (normAdvectionCoeff == 0.)
                srbavStabParam = 0.;

          break;

          case 3: //Nazarov + tanh^2(Res)
            srbavStabParam = math::pow(m_advElemLength[variable](k) / m_srbavScaleFactorH, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);
            //srbavStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);

          break;

          case 4: //Nazarov + tanh^2(Res)
            srbavStabParam = math::pow(m_advElemLength[variable](k) / m_srbavScaleFactorH, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2) / normAdvectionCoeff;
            //srbavStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2) / normAdvectionCoeff;
            if (normAdvectionCoeff == 0.)
                srbavStabParam = 0.;

          break;*/

          default:
              GISMO_ERROR("Wrong SRBAVtype set in stabilization evaluator.");
              break;
          }

        return srbavStabParam;
    }

    //------- original
    T getIsoADStabParamOriginal(index_t k, index_t variable = 0, T diffCoeff = -1.)
    {
        /*T alpha = 3./2.;
        T nu = 1.999;
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);*/
        // Johnson 1990: alpha and nu from (3/2, 2) such that Johnson suggested nu close to 2
        //T isoStabParam = math::max(0., alpha * math::pow(m_elemDiam, nu) * math::abs(m_residual[variable](k)) - m_diffusionCoeff(k));

        T denom = (m_solGrads[k].row(variable)).norm();
        //T isoStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(m_residual[variable](k)), 2) / denom;
        //T isoStabParam = m_elemDiam * math::abs(m_residual[variable](k)) / denom;
        T isoStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(math::abs(m_residual[variable](k)) / denom), 2);

        if (denom == 0.)
            isoStabParam = 0.;

        return isoStabParam;
    }

    T getRansIsoADStabParam(index_t k, index_t variable = 0)
    {
        T isoStabParam;

        switch (m_isoType)
        {
          case 3: //Nazarov h^alpha |Res|
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::abs(m_residual[variable](k));

            break;
          case 4: //Nazarov h^alpha tanh(Res)
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;

          default:
              GISMO_ERROR("Wrong isoADType set in stabilization evaluator.");
              break;
          }

        return isoStabParam;
    }

    T maxNorm()
    {
        T maxVal = 0.;
        for (int k = 0; k < m_nPoints; ++k)
        {
            maxVal = math::max(m_advectionCoeff[k].norm(), maxVal);

        }
        return maxVal;
    }

    T getIsoADStabParam(index_t k, index_t variable = 0, T diffCoeff = -1., index_t isoADtype = 6, index_t tauStabType = 2, T timeStep = -1.,
                        index_t deg = 0, T reactCoeff = 0.)
    {
        m_tauStabType = tauStabType;
        m_timeStep = timeStep;
        m_deg = deg;

        T normAdvectionCoeff = m_advectionCoeff[k].norm();

        T isoStabParam;

        switch (isoADtype)
        {
          case 1:
            isoStabParam = (1./getTauS(k, diffCoeff, reactCoeff, variable)) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;

          case 2:
            isoStabParam = getTauS(k, diffCoeff, reactCoeff, variable) * math::pow(math::tanh(m_residual[variable](k)), 2) * math::pow(normAdvectionCoeff, 2);

            break;

          case 3: //Nazarov h^alpha |Res|
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::abs(m_residual[variable](k));

            break;
          case 4: //Nazarov h^alpha tanh(Res)
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;
          default:
              GISMO_ERROR("Wrong isoADType set in stabilization evaluator.");
              break;
          }

        return isoStabParam;
    }

protected:
    index_t     m_nPoints;
    index_t     m_dim;
    index_t     m_deg;
    index_t     m_tauStabType;
    index_t     m_crosswindType;
    index_t     m_srbavType;
    index_t     m_isoType;
    T           m_timeStep;
    T           m_elemDiam;
    std::vector<gsVector<T> > m_advElemLength;
    std::vector<gsVector<T> > m_diffElemLength;
    T           m_alpha;
    T           m_srbavScaleFactorH;
    std::vector< gsMatrix<T> > m_solGrads;
    std::vector< gsVector<T> > m_residual;
    gsVector<T> m_diffusionCoeff;
    std::vector<gsVector<T> > m_advectionCoeff;
    gsVector<T> m_reactionCoeff;
    bool        m_bTauDeg;
};

//=====================================================================================================================================

} // namespace gismo

