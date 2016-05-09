class Exception: public std::exception
{
public:
	explicit Exception(const char* message):
		msg_(message) {}
	explicit Exception(const std::string& message):
		msg_(message) {}
	virtual ~Exception() throw (){}
	virtual const char* msg() const throw (){
		return msg_.c_str();
	}
protected:
	std::string msg_;
};

inline void exists (const std::string& name)
{
	std::ifstream f(name.c_str());
	if (!f.good())
		throw Exception("File doesn't exists:" + name);
}

inline void check_error(const int info)
{
	if (info != 0)
	{
		std::stringstream ss;
		if (info < 0)
			ss << "Lapack: the " << -info << "-th argument had an illegal value";
		else
			ss << "Lapack: U(" << info << "," << info << ") is exactly zero.";
		throw Exception(ss.str());
	}
}
