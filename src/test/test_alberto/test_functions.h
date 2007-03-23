class TestFunction
{
public:
	virtual std::string description() = 0;
	virtual double    f(double x, double y, double z)	 = 0;	
	virtual double dfdx(double x, double y, double z)	 = 0;
	virtual double dfdy(double x, double y, double z)	 = 0;
	virtual double dfdz(double x, double y, double z)	 = 0;
	virtual double d2fdx2(double x, double y, double z)	 = 0;
	virtual double d2fdy2(double x, double y, double z)	 = 0;
	virtual double d2fdz2(double x, double y, double z)	 = 0;
	virtual double d2fdxdy(double x, double y, double z) = 0;
	virtual double d2fdxdz(double x, double y, double z) = 0;
	virtual double d2fdydz(double x, double y, double z) = 0;	
};

class TestFunction01 : 	public TestFunction	// linear function
{
public:
	std::string description()	{return "Linear function";}
	double    f(double x, double y, double z)			{return 2.*x + 2.*y + 2.*z;}
	double dfdx(double x, double y, double z)			{return 2.;}
	double dfdy(double x, double y, double z)			{return 2.;}
	double dfdz(double x, double y, double z)			{return 2.;}
	double d2fdx2(double x, double y, double z)			{return 0.;}
	double d2fdy2(double x, double y, double z)			{return 0.;}
	double d2fdz2(double x, double y, double z)			{return 0.;}
	double d2fdxdy(double x, double y, double z)		{return 0.;}
	double d2fdxdz(double x, double y, double z)		{return 0.;}
	double d2fdydz(double x, double y, double z)		{return 0.;}	
};

class TestFunction02 : 	public TestFunction	// quadratic function
{
public:
	std::string description()	{return "Quadratic function";}
	double    f(double x, double y, double z)		{return x*x + y*y + z*z; }
	double dfdx(double x, double y, double z)		{return 2.*x;}
	double dfdy(double x, double y, double z)		{return 2.*y;}
	double dfdz(double x, double y, double z)		{return 2.*z;}
	double d2fdx2(double x, double y, double z)		{return 2.;}
	double d2fdy2(double x, double y, double z)		{return 2.;}
	double d2fdz2(double x, double y, double z)		{return 2.;}
	double d2fdxdy(double x, double y, double z)	{return 0.;}
	double d2fdxdz(double x, double y, double z)	{return 0.;}
	double d2fdydz(double x, double y, double z)	{return 0.;}		
};

class TestFunction03 : 	public TestFunction	// quadratic function + constant
{
public:
	std::string description()	{return "Quadratic function + constant";}
	double    f(double x, double y, double z)		{return x*x + y*y + z*z + 1; }
	double dfdx(double x, double y, double z)		{return 2.*x;}
	double dfdy(double x, double y, double z)		{return 2.*y;}
	double dfdz(double x, double y, double z)		{return 2.*z;}
	double d2fdx2(double x, double y, double z)		{return 2.;}
	double d2fdy2(double x, double y, double z)		{return 2.;}
	double d2fdz2(double x, double y, double z)		{return 2.;}
	double d2fdxdy(double x, double y, double z)	{return 0.;}
	double d2fdxdz(double x, double y, double z)	{return 0.;}
	double d2fdydz(double x, double y, double z)	{return 0.;}		
};


class TestFunction04 : 	public TestFunction	// quadratic function with mixed terms
{
public:
	std::string description()	{return "Quadratic function with mixed terms";}
	double    f(double x, double y, double z)		{return 3.*x*y + 3.*x*z + 3.*y*z;}
	double dfdx(double x, double y, double z)		{return 3.*y + 3.*z;}
	double dfdy(double x, double y, double z)		{return 3.*x + 3.*z;}
	double dfdz(double x, double y, double z)		{return 3.*x + 3.*y;}
	double d2fdx2(double x, double y, double z)		{return 0.;}
	double d2fdy2(double x, double y, double z)		{return 0.;}
	double d2fdz2(double x, double y, double z)		{return 0.;}
	double d2fdxdy(double x, double y, double z)	{return 3.;}
	double d2fdxdz(double x, double y, double z)	{return 3.;}
	double d2fdydz(double x, double y, double z)	{return 3.;}		
};

class TestFunction05 : 	public TestFunction	// sinusoidal functions - linear combination
{
public:
	std::string description()	{return "Sinusoidal functions - linear combination";}
	double    f(double x, double y, double z)		{return 6.*sin(x) + 6.*sin(y) + 6.*sin(z);}
	double dfdx(double x, double y, double z)		{return 6.*cos(x);}
	double dfdy(double x, double y, double z)		{return 6.*cos(y);}
	double dfdz(double x, double y, double z)		{return 6.*cos(z);}
	double d2fdx2(double x, double y, double z)		{return -6.*sin(x);}
	double d2fdy2(double x, double y, double z)		{return -6.*sin(y);}
	double d2fdz2(double x, double y, double z)		{return -6.*sin(z);}
	double d2fdxdy(double x, double y, double z)	{return 0.;}
	double d2fdxdz(double x, double y, double z)	{return 0.;}
	double d2fdydz(double x, double y, double z)	{return 0.;}			
};

class TestFunction06 : 	public TestFunction	// sinusoidal functions - linear combination + constant
{
public:
	std::string description()	{return "Sinusoidal functions - linear combination + constant";}
	double    f(double x, double y, double z)		{return 6.*sin(x) + 6.*sin(y) + 6.*sin(z) + 1;}
	double dfdx(double x, double y, double z)		{return 6.*cos(x);}
	double dfdy(double x, double y, double z)		{return 6.*cos(y);}
	double dfdz(double x, double y, double z)		{return 6.*cos(z);}
	double d2fdx2(double x, double y, double z)		{return -6.*sin(x);}
	double d2fdy2(double x, double y, double z)		{return -6.*sin(y);}
	double d2fdz2(double x, double y, double z)		{return -6.*sin(z);}		
	double d2fdxdy(double x, double y, double z)	{return 0.;}
	double d2fdxdz(double x, double y, double z)	{return 0.;}
	double d2fdydz(double x, double y, double z)	{return 0.;}		
};


class TestFunction07 : 	public TestFunction // sinusoidal functions - quadratic combination
{
public:
	std::string description()	{return "Sinusoidal functions - quadratic combination";}
	double    f(double x, double y, double z)		{return  2.*sin(x)*sin(y) + 2.*sin(x)*sin(z) + 2.*sin(y)*sin(z);}
	double dfdx(double x, double y, double z)		{return  2.*cos(x)*(sin(y)+sin(z));}
	double dfdy(double x, double y, double z)		{return  2.*cos(y)*(sin(x)+sin(z));}
	double dfdz(double x, double y, double z)		{return  2.*cos(z)*(sin(x)+sin(y));}
	double d2fdx2(double x, double y, double z)		{return -2.*sin(x)*(sin(y)+sin(z));}
	double d2fdy2(double x, double y, double z)		{return -2.*sin(y)*(sin(x)+sin(z));}
	double d2fdz2(double x, double y, double z)		{return -2.*sin(z)*(sin(x)+sin(y));}		
	double d2fdxdy(double x, double y, double z)	{return  2.*cos(y)*cos(x);}
	double d2fdxdz(double x, double y, double z)	{return  2.*cos(z)*cos(x);}
	double d2fdydz(double x, double y, double z)	{return  2.*cos(y)*cos(z);;}	
};


class TestFunction08 : 	public TestFunction // exponential function - cubic combination
{
public:
	std::string description()	{return "Exponential function - cubic combination";}
	double    f(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y) + 2.*exp(x)*exp(z) + 2.*exp(y)*exp(z);}
	double dfdx(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y) + 2.*exp(x)*exp(z);}
	double dfdy(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y) + 2.*exp(y)*exp(z);}
	double dfdz(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(z) + 2.*exp(y)*exp(z);}
	double d2fdx2(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y) + 2.*exp(x)*exp(z);}
	double d2fdy2(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y) + 2.*exp(y)*exp(z);}
	double d2fdz2(double x, double y, double z)		{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(z) + 2.*exp(y)*exp(z);}	
	double d2fdxdy(double x, double y, double z)	{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(y);}
	double d2fdxdz(double x, double y, double z)	{return  exp(x)*exp(y)*exp(z) + 2.*exp(x)*exp(z);}
	double d2fdydz(double x, double y, double z)	{return  exp(x)*exp(y)*exp(z) + 2.*exp(y)*exp(z);}		
};

class AnalyticalFunction
{
private:
	TestFunction *ptTestFunction;
public:
	void   setup(TestFunction *testFunction)			{ptTestFunction = testFunction;}
	std::string description()						{return ptTestFunction->description();}
	double    f(double x, double y, double z)			{return ptTestFunction->f(x,y,z);}
	double dfdx(double x, double y, double z)		{return ptTestFunction->dfdx(x,y,z);}
	double dfdy(double x, double y, double z)		{return ptTestFunction->dfdy(x,y,z);}
	double dfdz(double x, double y, double z)		{return ptTestFunction->dfdz(x,y,z);}
	double d2fdx2(double x, double y, double z)		{return ptTestFunction->d2fdx2(x,y,z);}
	double d2fdy2(double x, double y, double z)		{return ptTestFunction->d2fdy2(x,y,z);}
	double d2fdz2(double x, double y, double z)		{return ptTestFunction->d2fdz2(x,y,z);}
	double d2fdxdy(double x, double y, double z)	{return ptTestFunction->d2fdxdy(x,y,z);}
	double d2fdxdz(double x, double y, double z)	{return ptTestFunction->d2fdxdz(x,y,z);}
	double d2fdydz(double x, double y, double z)	{return ptTestFunction->d2fdydz(x,y,z);}	
};
