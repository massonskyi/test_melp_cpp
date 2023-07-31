#ifndef TEST_INCLUDE_MYFILTER_HPP_
#define TEST_INCLUDE_MYFILTER_HPP_


#include <array>
#include <iterator>

template<int NumSize, int DenumSize, int StepSize>
class LFilter{
public:

	LFilter(const std::array<float, NumSize>& numerator,
			const std::array<float, DenumSize>& denuminator);

	void step(const std::array<float, 180>::iterator signal,
				std::array<float, 180>::iterator output);


private:

	const std::array<float, NumSize> b;
	const std::array<float, DenumSize> a;

	std::array<float, NumSize - 1> state;

};

template<int NumSize, int DenumSize, int StepSize>
LFilter<NumSize,DenumSize,StepSize>::LFilter(const std::array<float, NumSize>& numerator,
		const std::array<float, DenumSize>& denominator)
: b(numerator),  a(denominator), state({}) { }


template<int NumSize, int DenumSize, int StepSize>
void LFilter<NumSize, DenumSize, StepSize>::step(const std::array<float, 180>::iterator signal, std::array<float, 180>::iterator output){
	int num_iter, denum_iter;
	std::array<float, 180>::iterator it_s = state.begin();

	if(b.size() > 1){
		for(int k = 0; k < StepSize; k++)
		{
			num_iter = 0;
			denum_iter = 0;
			it_s = state.begin();
			*(output + k) = *it_s + b[num_iter] * *(signal + k); /* Calculate first delay (output) */

			num_iter++;
			denum_iter++;

			/* Fill in middle delays */
			for(int n = 0; n < int(b.size() - 2); n++){
				*it_s = it_s[1] + *(signal + k) * b[num_iter] - *(output + k) * a[denum_iter];
				num_iter++;
				denum_iter++;
				it_s++;
			}

				/* Calculate last delay */
				*it_s = *(signal + k) * b[num_iter] - *(output + k) * a[denum_iter];


		}
	}else{
		num_iter = 0;

		for(int k = 0; k < StepSize; k++)
		{
			*(output + k) = *(signal + k) * b[num_iter];

		}
	}
}

#endif /* TEST_INCLUDE_MYFILTER_HPP_ */
