#include <cmath>
#include <complex>
#include <array>
#include <cassert>
#include <complex>

class MVSQ{
public:
	MVSQ() : w({}), temp_out({}), d({{}}), delta({}), tdelta({}), e({}){};
	void melp_msvq(const std::array<float,10>& lpcs,
					  const std::array<float,10>& f,
					  const std::array<float, 1280>& stage1,
					  const std::array<std::array<float, 640>, 3>& stage2,
					  std::array<float,4>& out)
	 {
	    reset_mvsq();
	    calculate_w(f, lpcs);
	    calculate_d(f, stage1,stage2);
	    for (int i = 11; i < 15; i++) {
	        out[i - 11] = d[0][i];
	    }
	}
private:
    std::array<float, 10> w;
    std::array<float, 10> temp_out;
    std::array<std::array<float, 15>, 9> d;
    std::array<float, 10> delta;
    std::array<float, 10> tdelta;
    std::array<std::array<float, 15>, 9> e;

	float msvq_powr(float x, float y)
	{
	    return std::pow(x, y);
	}
	void reset_mvsq()
	 {
	    for (int i = 0; i < 10; i++) {
	        w[i] = 0;
	        temp_out[i] = 0;
	    }
	    for (int i = 0; i < 9; i++) {
	        for (int j = 0; j < 15; j++) {
	            d[i][j] = 10000000;
	        }
	    }
	}
	void calculate_w(const std::array<float,10>& f, const std::array<float,10>& lpcs)
	{
	    std::complex<float> a[10] = {0};
	    for (int i = 0; i < 10; i++) {
	        for (int j = 0; j < 10; j++) {
	            a[j] = std::exp(std::complex<float>(0.0f, -1.0f)* std::complex<float>(f[i]) * std::complex<float>(j+1));
	        }
	        for (int j = 0; j < 10; j++) {
	            temp_out[i] += a[j].real() * lpcs[j];
	        }
	        temp_out[i] += 1;
	    }
	    for (int i = 0; i < 10; i++) {
	        w[i] = std::abs(temp_out[i]);
	    }
	    for (int i = 0; i < 10; i++) {
	        w[i] = std::pow(w[i], 0.3);
	    }
	    w[8] *= 0.64;
	    w[9] *= 0.16;
	}
	void calculate_d(const std::array<float,10>& f,
			const std::array<float, 1280>& stage1,
			const std::array<std::array<float, 640>, 3>& stage2)
	{
	    float temp, m;
	    for (int i = 0; i < 128; i++) {
	        for (int j = i * 10, index = 0; j < (i + 1) * 10; j++, index++) {
	            delta[index] = f[index] - stage1[j];
	        }
	        temp = 0;
	        for (int j = 0; j < 10; j++) {
	            temp += w[j] * std::pow(delta[j], 2);
	        }
	        m = 0;
	        while (m < 9) {
	            if (temp < d[m][0]) {
	                int j = 8;
	                while (j > m) {
	                    for (int k = 0; k < 15; k++) {
	                        d[j][k] = d[j-1][k];
	                    }
	                    j--;
	                }
	                d[m][0] = temp;
	                for (int j = 1; j < 11; j++) {
	                    d[m][j] = delta[j - 1];
	                }
	                d[m][11] = i + 1;
	                break;
	            }
	            m += 1;
	        }
	    }
	    for (int s = 0; s < 3; s++) {
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 15; j++) {
					e[i][j] = d[i][j];
				}
			}
			for (int i = 1; i < 11; i++) {
				d[0][i] = e[0][i] - stage2[s][i - 1];
			}
			d[0][0] = 0;
			for (int i = 0; i < 10; i++) {
				d[0][0] += w[i] * std::pow(d[0][i + 1], 2);
			}
			for (int i = 11; i < 12 + s; i++) {
				d[0][i] = e[0][i];
			}
			for (int i = 12; i < 13 + s; i++) {
				d[0][i] = 1;
			}
			for (m = 1; m < 8; m++) {
				for (int k = 0; k < 10; k++) {
					delta[k] = e[m][k + 1] - stage2[0][k];
				}
				temp = 0;
				for (int k = 0; k < 10; k++) {
					temp += w[k] * std::pow(delta[k], 2);
				}
				for (int num = 0; num < m; num++) {
					if (temp < d[num][0]) {
						int j = num + 1, n = num;
						while (n < 8) {
						    for (int i = 0; i < 15; i++) {
						        d[j][i] = d[n][i];
						    }
						    j++;
						    n++;
						}

						d[num][0] = temp;
						for (int i = 1; i < 11; i++) {
							d[num][i] =delta[i - 1];
						}
						for (int i = 11; i < 12 + s; i++) {
							d[num][i] = e[0][i];
						}
						for (int i = 12; i < 12 + s; i++) {
							d[num][i] = 1;
						}
						break;
					}
				}
				if (temp >= d[m - 1][0]) {
					for (int i = 1; i < 11; i++) {
						d[m][i] = delta[i - 1];
					}
					d[m][0] = temp;
					for (int i = 11; i < 12 + s; i++) {
						d[m][i] = e[0][i];
					}
					for (int i = 12; i < 13 + s; i++) {
						d[0][i] = 1;
					}
				}
			}
			for (int j = 0; j < 8; j++) {
				for (int k = 0; k < 64; k++) {
					for (int i = k * 10, index = 0; i < (k + 1) * 10; i++, index++) {
						tdelta[index] = stage2[s][i];
					}
					for (int i = 1; i < 11; i++) {
						delta[i - 1] = e[j][i] - tdelta[i - 1];
					}
					temp = 0;
					for (int u = 0; u < 10; u++) {
						temp += w[u] * std::pow(delta[u], 2);
					}
					for (int n = 0; n < 8; n++) {
						if (temp < d[n][0]) {
							int u = m + 1, nu = m;
							while (nu < 8) {
							    for (int i = 0; i < 15; i++) {
							        d[u][i] = d[nu][i];
							    }
							    u++;
							    nu++;
							}
							d[n][0] = temp;
							for (int i = 1; i < 11; i++) {
								d[n][i] = delta[i - 1];
							}
							for (int i = 11; i < 12 + s; i++) {
								d[n][i] = e[j][i];
							}
							d[n][12 + s] = k + 1;
							break;
						}
					}
				}
			}
		}
	}
};
