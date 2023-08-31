#include <array>

template <int StepSize>
struct Hamming{
public:
	Hamming(): n({}){ }

	void use_hamming(int M, std::array<float, StepSize>& output);
private:
	std::array<float, StepSize> n;

};

template<int StepSize>
void Hamming<StepSize>::use_hamming(int M, typename std::array<float, StepSize>& output){

	for(int i = 0; i < StepSize; i++)
		n[i] = 1 - M + 2 * i;

	for(int i = 0; i < StepSize; i++)
		output[i] = 0.54 + 0.46 * std::cos(3.14159265358979323846 * n[i] / (M - 1));

}
