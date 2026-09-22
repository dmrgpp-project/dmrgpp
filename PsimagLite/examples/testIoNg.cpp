#include <PsimagLite/Concurrency.h>
#include <PsimagLite/Io/IoNg.h>

void test1()
{
	PsimagLite::IoNg::Out io_out("hello.hdf5", PsimagLite::IoNg::ACC_TRUNC);

	std::vector<double> v(10, 42.0);
	io_out.write(v, "MyVector");

	io_out.close();

	PsimagLite::IoNg::In io_in("hello.hdf5");

	std::vector<double> w;

	io_in.read(w, "MyVector");

	std::cout << w;
}

int main(int argc, char* argv[])
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	test1();
}
