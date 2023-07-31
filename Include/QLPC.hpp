#include <array>
#include <complex>

/*Квантование коэффициентов линейного предсказания*/
class QLPC{

public:
	QLPC();

	void msvq(const std::array<float, 10>& koef_lpc,
			const std::array<float, 10>& LSF,
			std::array<float, 4>& out);
private:

	std::array<std::array<float, 15>, 9> d;
	std::array<std::array<float, 15>, 9> e;
	std::array<float, 10> w;
	std::array<float, 10> wreal;

	std::array<float, 10> delta;
	std::array<float, 10> temp;
	std::array<std::complex<float>, 10> a;
	std::array<float, 10> b;

};
