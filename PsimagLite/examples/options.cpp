#include <PsimagLite/Options.h>
#include <iostream>

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "USAGE is " << argv[0] << " comma,separated,list,of,options\n";
		return 1;
	}
	PsimagLite::Vector<PsimagLite::String>::Type register_opts;
	register_opts.push_back("fast");
	register_opts.push_back("verbose");
	register_opts.push_back("hasthreads");
	PsimagLite::Options::Writeable opt_writeable(register_opts,
	                                             PsimagLite::Options::Writeable::STRICT);

	PsimagLite::String            myoptions(argv[1]);
	PsimagLite::Options::Readable opts_readable(opt_writeable, myoptions);
	std::cout << "fast=" << opts_readable.isSet("fast") << "\n";
}
