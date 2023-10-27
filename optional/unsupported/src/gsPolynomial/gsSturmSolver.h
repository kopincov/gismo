

#include <gsPolynomial/gsMonomialPoly.h>


namespace gismo{

template<class T>
class gsSturmSolver
{
    

};

/// entry function
/// needs a polynomial and a boundary eps as an input and returns a
/// sorted matrix, where the intervals of the real roots are stored
	template<typename T>
	gsMatrix<T> FindRootIntervals(gsMonomialPoly<T> & poly, T eps)
    {
		gsMatrix<T> root(1, 3);
		root.setZero();

		root = CalcIntervals(poly, eps, root);
		root.conservativeResize(root.rows() - 1, 3);
		root = SortRoots(root);

		return root;
	}
	

	/// calculates the sturm sequence and the boundaries and calls
	/// then the function RootIsolation
	template<typename T>
	gsMatrix<T> CalcIntervals(const gsMonomialPoly<T> & poly, T eps, gsMatrix<T> rootIntervals)
	{
		gsMonomialPoly<T> copyPoly = poly;
		gsVector<T> coeff = copyPoly.coefs();
		int degree = copyPoly.basis().size() - 1;

		// checks if the input polynomial has 0 as a root
		int divByX = 0;

		while (!copyPoly.isNull() && copyPoly.coefs()(0) == 0)
        {
			rootIntervals(rootIntervals.rows() - 1, 0) = 0;
			rootIntervals(rootIntervals.rows() - 1, 1) = 0;
			rootIntervals(rootIntervals.rows() - 1, 2) += 1;

			gsMonomialBasis<T> monomialBasis(1);
			gsMatrix<T> coefs(2, 1);
			coefs << 0, 1;
			gsMonomialPoly<T> x(monomialBasis, coefs);
			copyPoly = PolyDivision(copyPoly, x, 1);
			divByX = 1;
		}

		if (divByX == 1)
        {
			rootIntervals.conservativeResize(rootIntervals.rows() + 1, 3);
		}

		//if the polynomial is not constant, the sturm sequence is calculated
		if (! copyPoly.isConstant() )
        {
			gsMatrix<T> sturmSeq(degree + 1, degree + 1);
			sturmSeq = SturmSequence(copyPoly);

			//checks if the polynomial has multiple roots
			if (sturmSeq.bottomRows(1).isZero())
			{
				gsVector< T >vec = sturmSeq.row(sturmSeq.rows() - 2);
				gsMonomialPoly<T >poly(copyPoly.basis(), vec);
				rootIntervals = CalcIntervals(poly, eps, rootIntervals);
				rootIntervals = CalcIntervals(PolyDivision(copyPoly, poly, 1), eps, rootIntervals);
			}
			// the polynomial doesn't have multiple roots
			else{
				RootIsolation(T(-CauchyBound(copyPoly)), CauchyBound(copyPoly),
                              eps, sturmSeq, rootIntervals);

			}
		}
		return rootIntervals;
	}

	
	/// Algorithm for the isolation of the roots
	template <typename T>
	void RootIsolation(T leftBound, T rightBound, T eps, gsMatrix<T> & sturmSeq, gsMatrix<T> & rootIntervals){

		int roots = SignChanges(HornerEval(sturmSeq, leftBound)) - SignChanges(HornerEval(sturmSeq, rightBound));

		//no roots in this interval
		if (roots == 0) { return; }

		//one root and the length of the intervals is <= 0.01
		else if (roots == 1 && rightBound - leftBound <= 0.01)
        {
			rootIntervals(rootIntervals.rows() - 1, 0) = leftBound;
			rootIntervals(rootIntervals.rows() - 1, 1) = rightBound;
			rootIntervals(rootIntervals.rows() - 1, 2) = 1;
			rootIntervals.conservativeResize(rootIntervals.rows() + 1, 3);
		}
		//more roots and the length of the interval is smaller than the border eps
		else if (roots > 1 && rightBound - leftBound < eps)
        {
			rootIntervals(rootIntervals.rows() - 1, 0) = leftBound;
			rootIntervals(rootIntervals.rows() - 1, 1) = rightBound;
			rootIntervals(rootIntervals.rows() - 1, 2) = roots;
			rootIntervals.conservativeResize(rootIntervals.rows() + 1, 3);
		}

		//more roots and the length of the interval is greater than eps
		else
        {
			T bound = (rightBound + leftBound) / 2;
			gsVector<T> f = sturmSeq.row(0);
			T root_check = HornerEvalPoly(f, bound);
            
			// the middle of the interval is a root
			if (root_check == 0)
			{
				roots = SignChanges(HornerEval(sturmSeq, T(bound - eps/2)))
                      - SignChanges(HornerEval(sturmSeq, T(bound + eps/2)));
				if (roots > 1)
				{
					rootIntervals(rootIntervals.rows() - 1, 0) = bound - eps / 2;
					rootIntervals(rootIntervals.rows() - 1, 1) = bound + eps / 2;
					rootIntervals(rootIntervals.rows() - 1, 2) = roots;
					rootIntervals.conservativeResize(rootIntervals.rows() + 1, 3);
				}
				else
				{
					rootIntervals(rootIntervals.rows() - 1, 0) = bound;
					rootIntervals(rootIntervals.rows() - 1, 1) = bound;
					rootIntervals(rootIntervals.rows() - 1, 2) = 1;
					rootIntervals.conservativeResize(rootIntervals.rows() + 1, 3);
				}
				RootIsolation(leftBound, T(bound - eps/2), eps, sturmSeq, rootIntervals);
				RootIsolation(T(bound + eps/2), rightBound, eps, sturmSeq, rootIntervals);
			}
			// the middle of the interval is not a root, so a bisection is done
			else
			{
				RootIsolation(leftBound, bound, eps, sturmSeq, rootIntervals);
				RootIsolation(bound, rightBound, eps, sturmSeq, rootIntervals);

			}
		}
	}


	/// returns the Cauchy Bound of a Polynomial, which is 1+max(abs(a_i/a_n))
	template<typename T>
	T CauchyBound(gsMonomialPoly<T> &poly)
	{
		gsVector<T> vec = poly.coefs();

		int lead = vec.size() - 1; // right-most non-zero
		for (int i = vec.size() - 1; i >= 0; i--)
        {
			if (vec(i) != 0)
            {
				lead = i;
				break;
			}
		}
        return 1 + (vec/vec[lead]).array().abs().maxCoeff();
	}

	///returns the sturm sequence in a gsMatrix
	template<typename T>
	gsMatrix<T> SturmSequence(gsMonomialPoly<T> & poly)
	{
		int deg =  poly.basis().size() - 1;
		gsMatrix<T> sturm(deg + 1, deg + 1);
		sturm.setZero();

		//first row = f
		sturm.row(0) = poly.coefs().transpose();

		//second row = f'
		sturm.row(1) = PolyDerivWithSameDegree(poly).coefs().transpose();

		if ( PolyDerivWithSameDegree(poly).isConstant() )
        {
			sturm.conservativeResize(2, sturm.cols());
			return sturm;
		}

		// i-th row = -rem(f_(i-2), f_(i-1)) 
		// until it is constant
		gsVector<T> dividend;
		gsVector<T> divisor;

		for (int i = 2; i <= deg; i++)
		{
			dividend = sturm.row(i - 2);
			divisor = sturm.row(i - 1);
			gsMonomialPoly<T> poly1(poly.basis(), dividend);
			gsMonomialPoly<T> poly2(poly.basis(), divisor);
			gsMonomialPoly<T> poly3 = PolyDivision(poly1, poly2, 2);

			if ( poly3.isNull() )
            {
				sturm.conservativeResize(i + 1, sturm.cols());
				return sturm;
			}


			for (int j = 0; j < poly3.coefs().size(); j++)
            {
				sturm(i, j) -= poly3.coefs()(j);
			}

			if ( poly3.isConstant() )
            {
				sturm.conservativeResize(i + 1, sturm.cols());
				return sturm;
			}

		}

		return sturm;
	}
	
	///Returns the derivation of the input polynomial including a leading 0, s.t. the
	///output polynomial is of the same degree as the input polynomial
	template<typename T>
	gsMonomialPoly<T> PolyDerivWithSameDegree(const gsMonomialPoly<T> & poly)
	{
		gsVector<T> deriv(poly.deg()+1);        
        //AM:
        //deriv.topRows(poly.deg()) = gsVector<int>::LinSpaced(poly.deg(),1,poly.deg()).array()
        //    * poly.coefs().bottomRows(poly.deg()).array();
        
		for (int i = 0; i < deriv.size() - 1; i++)
		{
			deriv(i) = poly.coefs()(i + 1)*(i + 1);
		}
		deriv(deriv.size() - 1) = 0;
        return gsMonomialPoly<T>(poly.basis(), deriv);
	}


	/// Polynomial Division
	/// the input num dedicates if the quotient or the remainder is returned
	/// num = 1 returns the qoutient
	/// num = 2 returns the remainder
	template<typename T>
	gsMonomialPoly<T> PolyDivision(const gsMonomialPoly<T> & dividendPoly,
                                   const gsMonomialPoly<T> & divisorPoly, int num)
	{
		if (num != 1 && num != 2)
        {
			std::cout << "ERROR: The number is not in the corresponding range. The dividend will be returned" << std::endl;
			return dividendPoly;
		}
		gsMonomialPoly<T> dividendCopy = dividendPoly;
		gsMonomialPoly<T> divisorCopy = divisorPoly;

		// storing the coefs of dividend and divisor seperately in a gsVector/gsMatrix
		gsMatrix<T> dividend = dividendCopy.coefs().transpose();

		gsVector<T> divisor(dividend.size());
		divisor.setZero();
		for (int i = 0; i < divisorCopy.coefs().size(); i++)
        {
			divisor(i) = divisorCopy.coefs()(i);
		}


		// getting the real degrees of dividend and divisor
		int deg_dividend = dividend.size() - 1;

		while (dividend(deg_dividend) == 0)
		{
			deg_dividend--;
		}


		int deg_divisor = divisor.size() - 1;

		while (divisor(deg_divisor) == 0)
		{
			deg_divisor--;
		}

		if (deg_divisor > deg_dividend)
        {
			std::cout << "Error, the division is not possible. The degree of the divisor is greater than the degree of the dividend. The dividend will be returned";
			return dividendPoly;
		}


		// Algorithm for polynomial long division
		gsVector<T> quot(dividend.size());
		gsVector<T> rem(dividend.size());

		quot.setZero();
		rem = dividend.transpose();
		int deg_rem = deg_dividend;
		
		int lead = 0;
		T lead_coeff = 0;

        gsVector<T> dummy;
		while (!(rem.maxCoeff() == 0 && rem.minCoeff() == 0) && deg_rem >= deg_divisor)
        {
			lead       = deg_rem - deg_divisor;
			lead_coeff = rem(deg_rem) / divisor(deg_divisor);
			quot(lead) = quot(lead) + lead_coeff;
			dummy      = ShiftRight(divisor, lead);
			rem        = rem - dummy * lead_coeff; 
			deg_rem    = deg_dividend;

			while (rem(deg_rem) == 0 && deg_rem > 0)
			{
				deg_rem--;
			}
		}

		// returning the right output
		if (num == 1){

			gsMonomialBasis<T> monomial_basis(deg_dividend - deg_divisor);
			gsMonomialPoly<T> polyQuot(monomial_basis, quot.head(deg_dividend - deg_divisor + 1));

			return polyQuot;
		}
		else{

			gsMonomialBasis<T> monomial_basis(deg_rem);
			gsMonomialPoly<T> polyRem(monomial_basis, rem.head(deg_rem + 1));

			return polyRem;
		}


	}

	///returns a gsVector where the elements are shifted "num"-spaces to the right
	/// the left side is filled with zeros
	template <typename T>
	gsVector<T> ShiftRight(const gsVector<T> & vec, int num)
	{        
		gsVector<T> shift_vec(vec.size());
        shift_vec.topRows(num).setZero();
        const index_t b = vec.size()-num;
        shift_vec.bottomRows(b) = vec.topRows(b);
		return shift_vec;
	}


	//returns the number of sign changes of a vector
	template <typename T>
	int SignChanges(const gsVector<T> & vec)
	{
        // AM:
        //int result = 0;
        //for (int i = 1; i < vec.size(); ++i)
        //    if (vec[i-1]*vec[i] .. 0) ++result;
            
		gsVector<T> vec_copy = vec;
		for (int i = 0; i < vec.size(); i++)
		{
			vec_copy(i) = vec(i) / math::abs(vec(i));// 0, +1 or -1
		}

		int result = 0;
		for (int i = 0; i < vec_copy.size() - 1; i++)
		{
			if (math::abs(vec_copy(i) + vec_copy(i + 1)) < 2) { result++; }

		}
		return result;
	}


	///in the input gsMatrix the polynomials are stored row by row with their coefficiens
	/// HornerEval returns a vector where the matrix is evaluated at position num
	template <typename T>
	gsVector<T> HornerEval(gsMatrix<T> & mat, T num)
	{
		int col = mat.cols();
		gsVector<T> eval(mat.rows());
		eval = mat.col(col - 1);

		for (int i = 2; i <= mat.cols(); i++)
		{
			eval = eval*num + mat.col(col - i);

		}
		return eval;
	}


	///evaluates a gsVector (which represents the coefficients of a polynomial) at value num
	template <typename T>
	T HornerEvalPoly(const gsVector<T> & vec, T num)
	{
		T eval = vec(vec.size() - 1);

		for (int i = 2; i <= vec.size(); i++)
		{
			eval = eval*num + vec(vec.size() - i);
		}
		return eval;
	}


	/// sorts the Matrix, where the roots are stored by their left bound
	/// and consolidates equal intervals
	template<typename T>
	gsMatrix<T> SortRoots(gsMatrix<T> & unsort)
    {
		gsMatrix<T> sort(unsort.rows(), 3);
		sort.setZero();
		T min = -1;
		int position = -1;

		for (int i = 0; i < unsort.rows(); i++)
		{

			min = unsort.col(0).minCoeff();
			for (int j = 0; j < unsort.rows(); j++)
            {
				if (min == unsort(j, 0)){
					position = j;
					break;
				}
			}

			sort.row(i) = unsort.row(position);
			unsort(position, 0) = INT_MAX;
		}
		
		for (int i = 0; i < sort.rows()-1; i++)
		{

			if ((sort(i, 0) == sort(i + 1, 0) && sort(i, 1) <= sort(i + 1, 1)) ||
				(sort(i, 0) <= sort(i + 1, 0) && sort(i, 1) == sort(i + 1, 1)))
            {
				sort(i, 2) += 1;
				for (int j = i+1; j < sort.rows()-1; j++){
					sort.row(j) = sort.row(j + 1);
				}
				
				sort.conservativeResize(sort.rows()- 1, 3);
				--i;
			}	
		}
		return sort;
	}

}




