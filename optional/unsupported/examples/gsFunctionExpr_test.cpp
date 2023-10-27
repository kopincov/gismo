
#include <iostream>

# include <gismo.h>


using namespace gismo;


/// Define the source function  -exp(x+z).*sin(y)
class source_function : public gismo::gsFunction<> 
{
public:
    /// Shared pointer for source_function
    typedef memory::shared_ptr< source_function > Ptr;

    /// Unique pointer for source_function
    typedef memory::unique_ptr< source_function > uPtr;

    GISMO_CLONE_FUNCTION(source_function)
    
    short_t domainDim() const {return 3;}
    
    /// Computes  -exp(x+z).*sin(y)
    void eval_into(const gsMatrix<>& u, gsMatrix<>& result) const
    { 
        result.resize(1, u.cols()) ;
        
        for( index_t i=0; i< result.cols(); ++i )
            result(i)= -exp( u(0,i)+u(2,i) ) * sin( u(1,i) );
    } 
    
    /// Computes  the gradient
    void deriv_into(const gsMatrix<>& u, gsMatrix<>& result) const
    { 
        result.resize(1, 3*u.cols());
        
        for( index_t i=0; i < u.cols(); ++i )
        {
            result(3*i+0)= -exp( u(0,i)+u(2,i) ) * sin( u(1,i) );
            result(3*i+1)= -exp( u(0,i)+u(2,i) ) * cos( u(1,i) );
            result(3*i+2)= -exp( u(0,i)+u(2,i) ) * sin( u(1,i) );
        }
        result.resize(3, u.cols());
    }
    
  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    { os << "F(x,y,z)= -exp(x+z)*sin(y)"; return os; };
};

/// Define the source function  -exp(x+z).*sin(y)
class source_function2 : public gismo::gsFunction<> 
{
public:
    /// Shared pointer for source_function2
    typedef memory::shared_ptr< source_function2 > Ptr;

    /// Unique pointer for source_function2
    typedef memory::unique_ptr< source_function2 > uPtr;

    GISMO_CLONE_FUNCTION(source_function2)

    short_t domainDim() const {return 2;}

    void eval_into(const gsMatrix<>& u, gsMatrix<>& result) const
    { 
      result.resize(1, u.cols());

      for( index_t i=0; i< result.cols(); ++i )
        result(i)= 
            (8.0-9.0*math::sqrt(u(0,i)*u(0,i)+u(1,i)*u(1,i)))*sin(2*math::atan2(u(1,i),u(0,i))) // atan2(x,y) or atan(y/x)
	  / (u(0,i)*u(0,i)+u(1,i)*u(1,i));
    }

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    { os << "8.0-9.0*sqrt(x^2 + y^2)*sin(2*atan(y/x)) / (x^2+y^2)"; return os; };
};

/// Define the source function  2*pi^2 sin(2*pi*x) * sin(2*pi*y)
class green_function : public gismo::gsFunction<>
{
public:
    /// Shared pointer for green_function
    typedef memory::shared_ptr< green_function > Ptr;

    /// Unique pointer for green_function
    typedef memory::unique_ptr< green_function > uPtr;

    GISMO_CLONE_FUNCTION(green_function)

    short_t domainDim() const {return 2;}

  /// Computes the green function
  // u: first 2 columns= boundary point , last 2 columns= inner point
  void eval_into(const gsMatrix<>& u, gsMatrix<>& result) const
    { 
      // Evaluates G
      result.resize(1,u.cols()) ;

      for( index_t i=0; i< result.cols(); ++i )
          result(i)= math::log( math::sqrt( math::pow(u(0,i)-5,2) + math::pow(u(1,i)-2,2)  ) ) / (2*EIGEN_PI);
    }; 

   // Evaluates partial derivatives of G
    void deriv_into(const gsMatrix<>& u, gsMatrix<>& result) const
    {
      result.resize(1, 2*u.cols());

      for( index_t i=0; i < u.cols(); ++i )
      {
          // Deriv wrt u(0,k)
          result(2*i+0)=  ( u(0,i) - 5 )  / (2*EIGEN_PI * ( math::pow(u(0,i)-5,2) + math::pow(u(1,i)-2,2) ) ) ;
          // Deriv wrt u(1,k)
          result(2*i+1)=  ( u(1,i) - 2 )  / (2*EIGEN_PI * ( math::pow(u(0,i)-5,2) + math::pow(u(1,i)-2,2) ) ) ;
      }
      result.resize(2,u.cols());

    }

  /// Prints the object as a string.
  std::ostream &print(std::ostream &os) const
    { os << "G(x,y)= ln(r)/2*pi"; return os; };
};

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit
    
	gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
	
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Evaluation points
    gsMatrix<> u(3,3);
    u<< 0,0,0, 1,2,3, 1.1,4.2,5.7 ;


    // Source function 1
    source_function f;
    gsFunctionExpr<> ff("-exp(x+z)*sin(y)",3);

    gsInfo<<"Manual Function: "<< f <<".\n" << "\n";
    gsInfo<<"Expression Function: "<< ff <<".\n" << "\n";
    gsInfo<<"nvars: "<< ff.domainDim() <<".\n" << "\n";
    gsInfo<<"eval f: "<< f.eval(u) <<".\n" << "\n";
    gsInfo<<"eval ff: "<< ff.eval(u) <<".\n" << "\n";
    gsInfo<<"grad f: \n"<< f.deriv(u) <<".\n" << "\n";
    gsInfo<<"grad ff: \n"<< ff.deriv(u) <<".\n" << "\n";

    // Source function 2
    gsMatrix<> u2(2,3);
    u2<< 1.10204 , 1.09836  , 1,
        0.0894445, 0.0219756, 2;

    source_function2 f2;
    gsInfo<<"points:\n"<< u2 <<".\n" << "\n";

gsInfo<<"atan:\n"<< atan2(u2(0,0),u2(1,0) ) <<".\n" << "\n";

gsFunctionExpr<> ff2("(8.0-9.0*sqrt(x^2 + y^2))*sin(2*atan(y/x))/(x^2+y^2)", 2);
    gsInfo<<"Manual Function: "<< f2 <<".\n" << "\n";
    gsInfo<<"Expression Function: "<< ff2 <<".\n" << "\n";
    gsInfo<<"nvars: "<< ff2.domainDim() <<".\n" << "\n";
    gsInfo<<"eval f2: "<< f2.eval(u2) <<".\n" << "\n";
    gsInfo<<"eval ff2: "<< ff2.eval(u2) <<".\n" << "\n";
    //u2 = u2.colwise().reverse().eval();
    //gsInfo<<"eval f2: "<< f2.eval(u2) <<".\n" << "\n";
    //gsInfo<<"eval ff2: "<< ff2.eval(u2) <<".\n" << "\n";

    gsMatrix<> v =  u.block(0,0,2,3);
    green_function g;
    gsFunctionExpr<> gg("log( sqrt( (x-5)^2 + (y-2)^2 ) ) / (2*pi)",2);

    gsInfo<<"Manual Function: "<< g <<".\n" << "\n";
    gsInfo<<"Expression Function: "<< gg <<".\n" << "\n";
    gsInfo<<"nvars: "<< gg.domainDim() <<".\n" << "\n";
    gsInfo<<"eval g: "<< g.eval(v) <<".\n" << "\n";
    gsInfo<<"eval gg: "<< gg.eval(v) <<".\n" << "\n";
    gsInfo<<"grad g: \n"<< g.deriv(v) <<".\n" << "\n";
    gsInfo<<"grad gg: \n"<< gg.deriv(v) <<".\n" << "\n";

    gsFunctionExpr<> ggcopy = gg;
    gsInfo<<"COPY Expression Function: "<< ggcopy <<".\n" << "\n";
    gsInfo<<"nvars: "<< ggcopy.domainDim() <<".\n" << "\n";
    gsInfo<<"eval g: "<< ggcopy.eval(v) <<".\n" << "\n";
    gsInfo<<"eval gg: "<< ggcopy.eval(v) <<".\n" << "\n";
    gsInfo<<"grad g: \n"<< ggcopy.deriv(v) <<".\n" << "\n";
    gsInfo<<"grad gg: \n"<< ggcopy.deriv(v) <<".\n" << "\n";

    //gsFunctionExpr<> fcopy( ff );

    gsFunctionExpr<> ss("sin(pi*x)*sin(pi*y)",2) ; 
    u2 << 0.1, 1.5, 1.5,
          0.1, 1.5, 0.5;
    gsInfo<<"Expression Function: "<< ss <<".\n" << "\n";
    gsInfo<<"points:\n"<< u2 <<".\n" << "\n";
    gsInfo<<"eval : "  << ss.eval(u2) <<".\n" << "\n";
    gsInfo<<"grad : \n"<< ss.deriv(u2) <<".\n" << "\n";
    
    //gsFunctionExpr<TT>::gsFunctionExpr(const gsFunctionExpr& other)
    
    gsFunctionExpr<> f_vec("x+y","2*x-y",2) ;
    gsFunctionExpr<> g_vec(f_vec);
    gsInfo<<"Expression Vector Function: "<< f_vec <<".\n" << "\n";
    gsInfo<<"Expression Vector Function: "<< g_vec <<".\n" << "\n";
    gsInfo<<"points:\n"<< f_vec <<".\n" << "\n";
    gsInfo<<"eval : "  << f_vec.eval(u2) <<".\n" << "\n";
    gsInfo<<"eval : "  << g_vec.eval(u2) <<".\n" << "\n";


    #ifdef GISMO_WITH_ADIFF
    {
        gsFunctionExpr<> f("x",2);
        gsMatrix<> result;
        gsMatrix<> correct;
        correct.setZero(3,4);
        gsMatrix<> points(2,4);
        points<<0,1,2,3,4,5,6,7;
        f.deriv2_into(points,result);
        if  ( !(result.array()==correct.array()).all() )
            gsInfo<<"second derivative of x->x :\n"<<result<<".\n" << "\n";
    }
    #endif

  if (plot) 
      {
	// Plot solution in paraview
	gsInfo<<"Plotting in Paraview...\n";
        gsMatrix<> domain(2,2);
        domain << 0, 2,
                  0 ,2 ;
	gsWriteParaview( ss, domain, "func", 400 );
	char cmdParaview[100];
	strcpy(cmdParaview,"paraview func.vts\0");
	strcat(cmdParaview," &");
	return system(cmdParaview);	
      }
  else
    {
      return 0;
    }

};
