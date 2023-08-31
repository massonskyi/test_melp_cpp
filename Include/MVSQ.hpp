#include <cmath>
#include <complex>
#include <array>
#include <cassert>
#include <complex>
struct mvsq {
    std::array<float, 10> w;
    std::array<float, 10> temp_out;
    std::array<std::array<float, 15>, 9> d;
    std::array<float, 10> delta;
    std::array<float, 10> tdelta;
    std::array<std::array<float, 15>, 9> e;

};

float msvq_powr(float x, float y) {
    return std::pow(x, y);
}

void reset_mvsq(mvsq* mvsq) {
    for (int i = 0; i < 10; i++) {
        mvsq->w[i] = 0;
        mvsq->temp_out[i] = 0;
    }
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 15; j++) {
            mvsq->d[i][j] = 10000000;
        }
    }
}
void calculate_w(mvsq* mvsq, const std::array<float,10>& f, const std::array<float,10>& lpcs) {
    std::complex<double> a[10] = {0};
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            a[j] = std::exp(std::complex<float>(0.0f, -1.0f)* std::complex<float>(f[i]) * std::complex<float>(j+1));
        }
        for (int j = 0; j < 10; j++) {
            mvsq->temp_out[i] += a[j].real() * lpcs[j];
        }
        mvsq->temp_out[i] += 1;
    }
    for (int i = 0; i < 10; i++) {
        mvsq->w[i] = std::abs(mvsq->temp_out[i]);
    }
    for (int i = 0; i < 10; i++) {
        mvsq->w[i] = std::pow(mvsq->w[i], 0.3);
    }
    mvsq->w[8] *= 0.64;
    mvsq->w[9] *= 0.16;
}

void calculate_d(mvsq* mvsq,
		const std::array<float,10>& f,
		const std::array<float, 1280>& stage1,
		const std::array<std::array<float, 640>, 3>& stage2)
{
    double temp, m;
    for (int i = 0; i < 128; i++) {
        for (int j = i * 10, index = 0; j < (i + 1) * 10; j++, index++) {
            mvsq->delta[index] = f[index] - stage1[j];
        }
        temp = 0;
        for (int j = 0; j < 10; j++) {
            temp += mvsq->w[j] * std::pow(mvsq->delta[j], 2);
        }
        m = 0;
        while (m < 9) {
            if (temp < mvsq->d[m][0]) {
                int j = 8;
                while (j > m) {
                    for (int k = 0; k < 15; k++) {
                        mvsq->d[j][k] = mvsq->d[j-1][k];
                    }
                    j--;
                }
                mvsq->d[m][0] = temp;
                for (int j = 1; j < 11; j++) {
                    mvsq->d[m][j] = mvsq->delta[j - 1];
                }
                mvsq->d[m][11] = i + 1;
                break;
            }
            m += 1;
        }
    }
    for (int s = 0; s < 3; s++) {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 15; j++) {
				mvsq->e[i][j] = mvsq->d[i][j];
			}
		}
		for (int i = 1; i < 11; i++) {
			mvsq->d[0][i] = mvsq->e[0][i] - stage2[s][i - 1];
		}
		mvsq->d[0][0] = 0;
		for (int i = 0; i < 10; i++) {
			mvsq->d[0][0] += mvsq->w[i] * std::pow(mvsq->d[0][i + 1], 2);
		}
		for (int i = 11; i < 12 + s; i++) {
			mvsq->d[0][i] = mvsq->e[0][i];
		}
		for (int i = 12; i < 13 + s; i++) {
			mvsq->d[0][i] = 1;
		}
		for (m = 1; m < 8; m++) {
			for (int k = 0; k < 10; k++) {
				mvsq->delta[k] = mvsq->e[m][k + 1] - stage2[0][k];
			}
			temp = 0;
			for (int k = 0; k < 10; k++) {
				temp += mvsq->w[k] * std::pow(mvsq->delta[k], 2);
			}
			for (int num = 0; num < m; num++) {
				if (temp < mvsq->d[num][0]) {
					int j = num + 1, n = num;
					while (n < 8) {
					    for (int i = 0; i < 15; i++) {
					        mvsq->d[j][i] = mvsq->d[n][i];
					    }
					    j++;
					    n++;
					}

					mvsq->d[num][0] = temp;
					for (int i = 1; i < 11; i++) {
						mvsq->d[num][i] = mvsq->delta[i - 1];
					}
					for (int i = 11; i < 12 + s; i++) {
						mvsq->d[num][i] = mvsq->e[0][i];
					}
					for (int i = 12; i < 12 + s; i++) {
						mvsq->d[num][i] = 1;
					}
					break;
				}
			}
			if (temp >= mvsq->d[m - 1][0]) {
				for (int i = 1; i < 11; i++) {
					mvsq->d[m][i] = mvsq->delta[i - 1];
				}
				mvsq->d[m][0] = temp;
				for (int i = 11; i < 12 + s; i++) {
					mvsq->d[m][i] = mvsq->e[0][i];
				}
				for (int i = 12; i < 13 + s; i++) {
					mvsq->d[0][i] = 1;
				}
			}
		}
		for (int j = 0; j < 8; j++) { //TODO this block is not currently working after 2 iteration
			for (int k = 0; k < 64; k++) {
				for (int i = k * 10, index = 0; i < (k + 1) * 10; i++, index++) {
					mvsq->tdelta[index] = stage2[s][i];
				}
				for (int i = 1; i < 11; i++) {
					mvsq->delta[i - 1] = mvsq->e[j][i] - mvsq->tdelta[i - 1];
				}
				temp = 0;
				for (int u = 0; u < 10; u++) {
					temp += mvsq->w[u] * std::pow(mvsq->delta[u], 2);
				}
				for (int n = 0; n < 8; n++) {
					if (temp < mvsq->d[n][0]) {
						int u = m + 1, nu = m;
						while (nu < 8) {
						    for (int i = 0; i < 15; i++) {
						        mvsq->d[u][i] = mvsq->d[nu][i];
						    }
						    u++;
						    nu++;
						}
						mvsq->d[n][0] = temp;
						for (int i = 1; i < 11; i++) {
							mvsq->d[n][i] = mvsq->delta[i - 1];
						}
						for (int i = 11; i < 12 + s; i++) {
							mvsq->d[n][i] = mvsq->e[j][i];
						}
						mvsq->d[n][12 + s] = k + 1;
						break;
					}
				}
			}
		}
	}
}

void melp_msvq(mvsq *mvsq,
                  const std::array<float,10>& lpcs,
				  const std::array<float,10>& f,
				  const std::array<float, 1280>& stage1,
				  const std::array<std::array<float, 640>, 3>& stage2,
				  std::array<float,4>& out) {
    reset_mvsq(mvsq);
    calculate_w(mvsq, f, lpcs);
    calculate_d(mvsq,f, stage1,stage2);
    for (int i = 11; i < 15; i++) {
        out[i - 11] = mvsq->d[0][i];
    }
}
