#ifndef TEST_INCLUDE_MYFILTER_HPP_
#define TEST_INCLUDE_MYFILTER_HPP_

#include <array>
#include <iterator>

/**
 * @brief A template class that implements a linear filter.
 *
 * The LFilter class takes in two arrays of floats, numerator and denominator,
 * which represent the coefficients of the filter. It provides a step method
 * that applies the filter to an input signal and produces an output signal.
 * The filter is implemented using a difference equation, and the state of the
 * filter is stored in an array called state.
 *
 * @tparam NumSize The size of the numerator array.
 * @tparam DenumSize The size of the denominator array.
 * @tparam StepSize The size of the input and output signal arrays.
 */
template <int NumSize, int DenumSize, int StepSize>
class LFilter
{
public:
	/**
	 * @brief Constructs an LFilter object with the given numerator and denominator arrays.
	 *
	 * @param numerator The array of floats representing the numerator coefficients of the filter.
	 * @param denominator The array of floats representing the denominator coefficients of the filter.
	 */
	LFilter(const std::array<float, NumSize> &numerator,
			const std::array<float, DenumSize> &denominator) : b(numerator), a(denominator), state({}) {};
	LFilter(){}
	/**
	 * @brief Applies the filter to the input signal and produces the output signal.
	 *
	 * This method updates the state of the filter after each step.
	 *
	 * @param signal An iterator pointing to the input signal array.
	 * @param output An iterator pointing to the output signal array.
	 */
	void set_coeffs(std::array<float, NumSize> &numerator,
			std::array<float, DenumSize> &denominator){
		std::copy(numerator.begin(), numerator.end(), b.begin());
		std::copy(denominator.begin(), denominator.end(), a.begin());
	}
	void step(const std::array<float, StepSize>::iterator signal,
			std::array<float, StepSize>::iterator output) {
		int num_iter, denum_iter;
		typename std::array<float, NumSize-1>::iterator it_s = state.begin();
		if (DenumSize > 2)
		{
			for (int k = 0; k < StepSize; k++)
			{
				num_iter = 0;
				denum_iter = 0;
				it_s = state.begin();
				*(output + k) = *it_s + b[num_iter] * *(signal + k); /* Calculate first delay (output) */

				num_iter++;
				denum_iter++;

				/* Fill in middle delays */
				for (int n = 0; n < int(b.size() - 2); n++)
				{
					*it_s = it_s[1] + *(signal + k) * b[num_iter] - *(output + k) * a[denum_iter];
					num_iter++;
					denum_iter++;
					it_s++;
				}

				/* Calculate last delay */
				*it_s = *(signal + k) * b[num_iter] - *(output + k) * a[denum_iter];
			}
		}
		else
		{
		    for (int n = 0; n < StepSize; n++) {
		    	*(output + n) = 0;
		        for (int k = 0; k < NumSize; k++) {
		            if (n - k >= 0 && n - k < StepSize) {
		            	*(output + n) += *(signal + (n - k)) * b[k];
		            }
		        }
		    }
		}
	}


private:
	std::array<float, NumSize> b;	  /**< The array of numerator coefficients. */
	std::array<float, DenumSize> a; /**< The array of denominator coefficients. */
	std::array<float, NumSize - 1> state; /**< The state of the filter. */
};

#endif /* TEST_INCLUDE_MYFILTER_HPP_ */
