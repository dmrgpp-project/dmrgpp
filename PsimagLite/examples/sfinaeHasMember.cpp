#include <iostream>
#include <type_traits>

// SFINAE test
// source:
// https://stackoverflow.com/questions/257288/is-it-possible-to-write-a-template-to-check-for-a-functions-existence
template <typename T> class HasHelloworld {
	typedef char One;
	struct Two {
		char x[2];
	};

	template <typename C> static One test(typeof(&C::helloworld));
	template <typename C> static Two test(...);

public:

	enum
	{
		VALUE = sizeof(test<T>(0)) == sizeof(char)
	};
};

class A {
public:

	int helloworld() const { return 42; }
};

class B { };

template <typename T> typename std::enable_if<HasHelloworld<T>::VALUE, int>::type f(const T& t)
{
	return t.helloworld();
}

int main()
{
	A a;
	std::cout << f(a) << "\n";
	// B b;
	// f(b);
}
