
#include <array>
#include <numeric>
#include "Hamming.hpp"

template<int StepSize, int N>
class Lpc{

public:
	Lpc();

	void lpc(std::array<float, StepSize>& signal,
			 std::array<float, N>& output);

	void convolution(const std::array<float, 11>& b,
			const std::array<float, 360>& signal,
			std::array<float, 360>& output);
private:

	using farray_n= std::array<float, N + 1>;
	using stHamming = Hamming<StepSize>;
	using HammingArray = std::array<float, StepSize>;
	stHamming hamming;

	std::array<float, 10> e;
	std::array<float, N> k;
	HammingArray array_hamming;
	HammingArray temp;
	farray_n r;
	farray_n a_new;
	farray_n a_temp;
	farray_n e_new;

};

template<int StepSize, int N>
Lpc<StepSize, N>::Lpc()
	: e({}), k({}), r({}),
	 a_new({}), a_temp({}), e_new({}){

}

template<int StepSize, int N>
void Lpc<StepSize, N>::lpc(typename std::array<float, StepSize>& signal , typename std::array<float, N>& output){

	hamming.use_hamming(StepSize, array_hamming);

	for(int i = 0; i < StepSize; i++){
		temp[i] = signal[i] * array_hamming[i];
	}

    for(int i = 0; i <= N; i++){
    	r[i] = 0;
    	for(int j = 0; j < StepSize - i; j++)
    		r[i] += temp[j] * temp[j + i];
    }

    a_temp[0] = 1;
    e[0] = r[0];
    for (int i = 1; i <= N; i++) {
        k[i - 1] = r[i];

        for(int j = 1; j < i; j++){
        	k[i - 1] -= a_temp[j] * r[i - j];
        }

        k[i - 1] /= e[i - 1];

        a_new[i] = k[i - 1];

        for (int j = 1; j < i; j++)
              a_new[j] =  a_temp[j] -  k[i - 1] *  a_temp[i - j];

         for (int j = 0; j <= i; j++)
              a_temp[j] =  a_new[j];

        e_new[i] = (1 - k[i - 1] * k[i - 1]) * e[i - 1];

        e[i] = e_new[i];
    }

    for(int i = 1; i<= N; i++){
    	output[i - 1] = -a_new[i];
    }
}

template<int StepSize, int N>
void Lpc<StepSize, N>::convolution(const std::array<float, 11>& b,
		const std::array<float, 360>& input,
		std::array<float, 360>& output) {
    for (int n = 0; n < output.size(); n++) {
        output[n] = 0;
        for (int k = 0; k < b.size(); k++) {
            if (n - k >= 0 && n - k < input.size()) {
            	output[n] += input[n - k] * b[k];
            }
        }
    }
}
