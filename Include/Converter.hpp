#include <array>


class Converter{
public:
	Converter();
	void lpc2lsf(const std::array<float, 10>& koef_lpc,
								std::array<float, 10>& out);
	/*void lsf2lpc(const std::array<float, 10>& lsfs,
								std::array<float, 10>& out);
								*/
	void melp_lsf2lpc(std::array<float, 10> lsfs, std::array<float, 10> out);
private:
	using array_temp3f = std::array<float, 3>;
	using array_temp6f = std::array<float, 6>;

	// to lsf
	std::array<array_temp6f, 2> P;
	std::array<array_temp6f, 6> b;

	std::array<float, 10> f1;
	array_temp6f tmp;
	array_temp6f tmp2;

	// to lpc
	std::array<float, 10> new_lsfs;
	std::array<std::array<float, 5>, 2> w;



	// function
	void clear_lpc2lsf();
	void clear_lsf2lpc();
};
