#include "QLPC.hpp"
#include "stage.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>

const int N = 128;
const int M = 10;

QLPC::QLPC()
	: d({{}}), e({{}}), w({}),wreal({}),
	  delta({}), temp({}), a({}), b({}){
	for(int i  = 0; i < 9; i++)
		std::fill(d[i].begin(), d[i].end(), 10000000.0f);
}

void QLPC::msvq(const std::array<float, 10>& koef_lpc,
		const std::array<float, 10>& LSF,
		std::array<float, 4>& out)
{

	float tmp;
	for(int k = 0; k < 10;k ++){
		for(int i = 0; i < 10; i++){
			a[i] = std::exp(std::complex<float>(0.0f, -1.0f) *std::complex<float>(LSF[k]) * std::complex<float>(i + 1));
		}

		float sum;
		for(int j = 0; j < 10; j++){
			sum += a[j].real() * koef_lpc[j];
		}
		w[k] = 1.0f + sum;
	}
    for(int i = 0; i < 10; i++) {
    	w[i] = std::abs(w[i]);
    	w[i] = std::pow(w[i], 0.3);
    }

    w[8] *= 0.64;
    w[9] *= 0.16;
    /*
     * d(m,0)    - оценка
     * d(m,1:11) - разность вектора и кодового слова
     * d(m,11:15)- кодовое слово
     */

    // Определение индекса вектора КК первого уровня

    for(int s = 0; s < 128; s++){
    	for(int i = 0, j = s * 10; i < 10;i++, s++){
    		delta[i] = LSF[i] - st1[j];
    	}
    	tmp = 0;
    	for(int i = 0; i < 10; i++){
    		tmp += w[i] * std::pow(delta[i] , 2);
    	}
    	int m = 0;
    	while (m < 9){
    		if (tmp < d[m][0]){

    			// std::copy(d.begin() + m,d.end(), d.begin() + (m - 1));
                for(int i = 0; i < 15; i++){
                	d[m + 1][i] = d[m][i];
                }

    			d[m][0] = tmp;

                // std::copy(d[m - 1].begin() + 1, d[m - 1].begin()+ 11, delta.begin());
                for(int i = 1; i < 11; i++){
                	d[m][i] = delta[i - 1];
                }
    			d[m][12] = s;

                break;
    		}
    		m++;
    	}
    }

    for(int s = 0; s < 3; s++){

    	for(int i = 0;i < 9; i++){
    		for(int j = 0; j < 15; j++){
    			e[i][j] = d[i][j];
    		}
    	}

    	for(int i = 1; i < 11; i++){
    		d[0][i] = e[0][i] - st2[s][i - 1];
    	}
    	d[0][0] = 0;
    	for(int i = 0; i < 10; i++){
    		d[0][0] += w[i] * std::pow(d[0][i], 2);
    	}

    	// std::copy(d[0].begin() + 11, d[0].begin() + 12 + s, e[0].begin() + 11);/
    	for(int i = 11, j = 11; i < 11 + s + 1; i++, j++){
    		d[0][i] = e[0][j];
    	}
    	d[0][12 + s + 1] = 1;

    	for(int m = 0; m < 8; m++){
    		for(int i = 1; i < 11; i++){
    			delta[i - 1] = e[m][i] - st2[0][i - 1];
    		}
    		tmp = 0;
        	for(int i = 0; i < 10; i++){
        		tmp += w[i] * std::pow(delta[i],2);
        	}
    		for(int num = 0; num < m + 1; num++){

    			if(tmp < d[num][1]){

    				// std::copy(d.begin() + num + 1, d.end(), d.begin() + num);
    				for(int i = 0; i < 15; i++){
    					d[num + 1][i] = d[num][i];
    				}

    				d[num][1] = tmp;

    				// std::copy(delta.begin(), delta.end(), d[num].begin() + 1); // Assign delta to d(num,2:11)
    				for(int i = 1; i < 11; i++){
    					d[num][i] = delta[i - 1];
    				}
    				// std::copy(e[0].begin(), e[0].begin() + s, d[num].begin() + 12); // Assign e(1,12:11+s) to
    				for(int i = 11, j = 11; i < 11 + s + 1; i++, j++){
    					d[num][i] = e[0][j];
    				}
    				d[num][12 + s] = 1;

                    break;
    			}
    		}
    		if (tmp >= d[m - 1][0]){
    			// std::copy(d[m].begin() + 1,d[m].begin()+11,delta.begin());
    			for(int i = 1; i < 11; i++){
    				d[m][i] = delta[i - 1];
    			}
    			d[m][0] = tmp;
    			// std::copy(d[m].begin() + 11,d[m].begin() + 12 + s, e[0].begin() + 11 + s);

    			for(int i = 11, j = 11; i < 11 + s + 1;i++, j++){
    				d[m][i] = e[0][j];
    			}
    			d[m][12 + s] = 1;
    		}
    	}
    	for(int j = 0; j < 8; j++){
    		for(int k = 0; k < 64; k++){
    			for(int i = 0, kj = k * 10; i < 11; i++, kj++){
    				delta[i] = e[j][i] - st2[s][kj];
    			}
    			tmp = 0;
            	for(int i = 0; i < 10; i++){
            		tmp += w[i] * std::pow(delta[i],2);
            	}
    			for(int n = 0; n < 8; n++){
    				if(tmp < d[n][0]){
    					// std::copy(d[n + 1].begin(),d[n + 1].begin() + 9,d[n].begin());
    					for(int i = 0; i < 15; i++){
    						d[n + 1][i] = d[n][i];
    					}
    					d[n][0] = tmp;

    					// std::copy(d[n].begin() + 1, d[n].begin() + 11, delta.begin());
    					for(int i = 1; i < 11; i++){
    						d[n][i] = delta[i - 1];
    					}
    					// std::copy(d[n].begin() + 11, d[n].begin() + 12 + s, e[j].begin() + 11 + s);
    					for(int i = 11; i < 11 + s + 1; i++){
    						d[n][i] = e[j][i];
    					}
    					d[n][12 + s] = k;
    				}
    			}
    		}
    	}
    }
    std::copy(d[0].begin() + 11, d[0].end(),out.begin());
    int x = 43;
}
