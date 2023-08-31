#include <cmath>
#include "Encoder.hpp"
#include "Coeffs.hpp"
#include "stage.h"
#include <numeric>
#include <algorithm>
#include "codebook_fmcq.hpp"
#include "MVSQ.hpp"
Encoder::Encoder()
	: pavg(50),G2p(20),array_hamming({}), envelopes_filter(smooth_num, smooth_den),
	  cheb_filter(dcr_num, dcr_den),
	  butter_filter(butt_1000num, butt_1000den),
	  band_filter({
		BandFilter(butt_bp_num[0],butt_bp_den[0]),
		BandFilter(butt_bp_num[1],butt_bp_den[1]),
		BandFilter(butt_bp_num[2],butt_bp_den[2]),
		BandFilter(butt_bp_num[3],butt_bp_den[3]),
		BandFilter(butt_bp_num[4],butt_bp_den[4])
	  }),
	  hpFilter(butt_1000num, butt_1000den),
	  conv_lpc(),
	  conv_e_resid(),
	  lpc({}), conv({}), qlpc({}), buffer({50,50,50}),
	  koef_lpc({}), temp({}), LSF({}),
	  sig_in({}), sig_1000({}),envelopes({{}}),
	  bands({{}}),
	  tmp1({{}}),
	  e_resid({}),
	  fltd_resid({}),
	  tresid({}),
	  copy_sigin({}),
	  resid2({}),
	  tmp({}){
	hamming.use_hamming(200, array_hamming);
}


void Encoder::encode(std::array<float, 180>& signal){
	int cur_intp = 0, Qpitch = 0;
	std::array<float,10>lsfs;
	std::array<float,10>lpc2;
	mvsq *mv = new mvsq();
	int jitter = 0;
	float peak = 0;
	std::array<float,1> a({});
	std::array<float, 2> G({});
	std::array<float, 5> vp({});
	std::array<float, 4> ls({});
	std::array<float, 11> klpc({});
	std::pair<float,float> tp = {0.0f,0.0f};
	std::pair<float, float> pr3 = {0.0f, 0.0f};

	updating_buffers();

	cheb_filter.step(signal.begin(), sig_in.begin() + 180);

	butter_filter.step(sig_in.begin() + 180, sig_1000.begin() + 180);

	cur_intp = intpitch(sig_1000.begin(),160, 40);

	melp_5b(sig_in.begin() + 180);

	tp = pitch2(bands[0].begin(),cur_intp);

	vp[0] = tp.second;

	melp_bpva(bands, envelopes,tp.first,vp.begin() + 1);

	if(vp[0] < 0.5) jitter = 1;
	else jitter = 0;

	std::copy(sig_in.begin() + 80, sig_in.begin() + 280, temp.begin());

	lpc.lpc(array_hamming, temp, koef_lpc);
    for (int j = 2; j < 12; j++) {
        koef_lpc[j - 2] *= pow(0.994, j);
    }
	std::copy(koef_lpc.begin(), koef_lpc.end(), klpc.begin() + 1);
    klpc[0] = 1.0;
    // lpc_residual(klpc,sig_in,e_resid);
    conv_e_resid.set_coeffs(klpc, a);
    conv_e_resid.step(sig_in.begin(), e_resid.begin());

    peak = calculatePeak(e_resid.begin() + 10,106,265);
    if(peak > 1.34) vp[0] = 1;
    if(peak > 1.6){
    	vp[1] = 1;
    	vp[2] = 1;
    }

    hpFilter.step(e_resid.begin()+10, fltd_resid.begin() + 5);

    pr3 = pitch3(sig_in, fltd_resid,tp.first,pavg);

    melp_gain(sig_in, vp[0], tp.first, G);

    melp_APU(pr3,G[1]);

    conv.lpc2lsf(koef_lpc, LSF);

    lsf_clmp(LSF, LSF);

    // qlpc.msvq(koef_lpc, LSF, ls);
    melp_msvq(mv, koef_lpc, LSF, st1, st2, ls);
    QGain(G,G2p);

    G2p = G[1];

    if(vp[0] > 0.6) Qpitch = melp_Qpitch(pr3.first);

    d_lsf(ls, lsfs);

    conv.lsf2lpc(lsfs, lpc2);
    conv_lpc.set_coeffs(lpc2, a);
    std::copy(sig_in.begin() + 75, sig_in.begin() + 285, copy_sigin.begin());
    conv_lpc.step(copy_sigin.begin() + 75, tresid.begin());

    for(int i = 0; i < 200; i++){
    	resid2[i] = *(tresid.begin() + 10 + i) * array_hamming[i];
    }
    fft<512>(resid2);
    for(int i = 0;i < 512; i++){
    	magf[i] = std::abs(resid2[i]);
    }
    find_harm(magf, pr3.first);
    computeWeights(pr3.first);
    int QFM = melp_FMCQ(FMCQ_CODEBOOK);
    // ParamsFrame(ls,G,vp,QFM,Qpitch,jitter);
    int breakp = 0;
}

void Encoder::updating_buffers(){

	unsigned int i;

	std::move(std::begin(sig_in) + 180, std::end(sig_in), std::begin(sig_in));
	std::move(std::begin(sig_1000) + 180, std::end(sig_1000), std::begin(sig_1000));

	for(i = 0; i < 5; i++){
		std::move(std::begin(bands[i]) + 180, std::end(bands[i]), std::begin(bands[i]));
	}

	for(i = 0; i < 4; i++){
		std::move(std::begin(envelopes[i]) + 180, std::end(envelopes[i]), std::begin(envelopes[i]));
	}
}

int Encoder::intpitch(array_2d::iterator signal, int ipmax, int ipmin){
	float r = 0, r_new = 0, c0_t = 0, c0_0 = 0, ct_t = 0, den = 0;

	int T = 80, k = 0;


	for(int tao = ipmin; tao < ipmax; tao++){
		k = int(tao / 2);

		c0_t = std::inner_product(signal + (100 - k),  signal + (259 - k),
				signal + (100 - k + tao), 0.0);
		c0_0 = std::inner_product(signal + (100 - k), signal + (259 - k),
				signal + (100 - k), 0.0);
		ct_t = std::inner_product(signal + (100 - k + tao), signal + (259 - k + tao),
				signal + (100 - k + tao), 0.0);

		den = c0_0 * ct_t;

		if(den > 0) r_new = c0_t * c0_t / den;

		if(r_new > r) {
			r = r_new;
			T = tao;
		}
	}
	return T;
}


void Encoder::melp_5b(array_1d::iterator sig_in){
	unsigned int i;

	for(i = 0; i < 5; i++){
		band_filter[i].step(sig_in, bands[i].begin() + 180);
	}

	for(i = 1; i < 5; i++){
		std::transform(bands[i].begin() + 180, bands[i].end(), fbands[i - 1].begin(),
				[](float value){return std::fabs(value);});
	}

	for(i = 0; i < 4; i++){
		envelopes_filter.step(fbands[i].begin(), envelopes[i].begin() + 180);
	}

}

std::pair<float, float> Encoder::fpr(const array_2d::iterator signal, int cur_intp){

	int k = int(cur_intp / 2);
	float fp,fr;
	float c0_tm1, c0_t1, c0_t, ct_t, c0_0, ct_t1, ct1_t1, den, delta;

	c0_tm1 = std::inner_product(signal + (100 - k),  signal + (259 - k),
					signal + (100 - k + cur_intp - 1), 0.0f); //c[0,t - 1]

	c0_t1 = std::inner_product(signal + (100 - k), signal + (259 - k),
					signal + (100 - k + cur_intp + 1), 0.0f);

	c0_t = std::inner_product(signal + (100 - k), signal + (259 - k),
					signal + (100 - k + cur_intp), 0.0f);

	if(c0_tm1 > c0_t1){
		c0_t1 = c0_t;
		c0_t = c0_tm1;
		cur_intp--;
	}

	ct_t = std::inner_product(signal + (100 - k + cur_intp), signal + (259 - k + cur_intp),
					signal + (100 - k + cur_intp), 0.0f);

	c0_0 = std::inner_product(signal + (100 - k), signal + (259 - k), signal + (100 - k), 0.0f);
	ct_t1 = std::inner_product(signal + (100 - k + cur_intp), signal + (259 - k + cur_intp),
					signal + (100 - k + cur_intp + 1), 0.0f);
	ct1_t1 = std::inner_product(signal + (100 - k + cur_intp + 1), signal +(259 - k + cur_intp + 1),
					signal + (100 - k + cur_intp + 1), 0.0f);

	den = c0_t1 * (ct_t - ct_t1) + c0_t * (ct1_t1 - ct_t1);

	if(abs(den) > 0.01) delta = (c0_t1 * ct_t - c0_t *ct_t1) / den;
	else delta = 0.5;

	if(delta < -1) delta = -1;
	if(delta > 2) delta = 2;

	fp = cur_intp + delta;
	den = c0_0 * (ct_t * pow((1-delta),2) + 2 * delta *(1 - delta) * ct_t1 + pow(delta, 2) * ct1_t1);
	den = sqrt(den);

	if(den > 0.01) fr = ((1 - delta) * c0_t + delta * c0_t1) / den;
	else fr = 0;

	if (fp < 20) fp = 20;
	if (fp > 160) fp = 160;

	return std::make_pair(fp,fr);
}

std::pair<float,float> Encoder::pitch2(const array_2d::iterator signal, int cur_intp){
	int low = cur_intp - 5;

	if(low < 20) low = 20;

	int up = cur_intp + 5;

	if(up > 160) up = 160;

	int intp = intpitch(signal,up,low);
	std::pair<float,float> temp = fpr(signal,intp);

	return temp;
}

void Encoder::melp_bpva(const std::array<array_2d, 5>& bands, const std::array<array_2d,4>&envelopes,int p2, std::array<float, 4>::iterator vbp){

	unsigned int k = 0;
	float temp = 0;

	std::pair<float,float> pr1 = {0.0f,0.0f}, pr2 = {0.0f,0.0f};

	for(unsigned int j = 1; j < 4; j++){
		k = j + 1;
		pr1 = fpr((std::array<float, 360>::iterator)bands[k].begin(), p2);
		pr2 = fpr((std::array<float, 360>::iterator)envelopes[j].begin(), p2);

		pr2.second -= 0.1;
		if(pr2.second > pr1.second) temp = pr2.second;
		else temp = pr1.second;

		if(temp > 0.6) *(vbp + j) = 1;
		else *(vbp + j) = 0;
	}
	if (std::equal(vbp + 1, vbp + 4, std::vector<int>{0, 0, 0}.begin())) {
	        *(vbp + 4) = 0;
	}
}

float Encoder::calculatePeak(const std::array<float, 350>::iterator e_resid, int start, int end) {
    std::vector<float> squaredResid;
    squaredResid.reserve(end - start + 1);

    // Calculate the squared residuals
    for (int i = start; i <= end; i++) {
        squaredResid.push_back(*(e_resid + i) * *(e_resid + i));
    }

    // Calculate the sum of squared residuals
    float sumSquaredResid = std::accumulate(squaredResid.begin(), squaredResid.end(), 0.0);

    // Calculate the square root of the sum of squared residuals divided by the number of elements
    float rms = std::sqrt(sumSquaredResid / (end - start + 1));

    // Calculate the sum of absolute residuals
    float sumAbsResid = std::accumulate(e_resid + start, e_resid + end + 1, 0.0, [](float sum, float val) {
        return sum + std::abs(val);
    });

    // Calculate the peak value
    float peak = rms / (sumAbsResid / (end - start + 1));

    return peak;
}

std::pair<float,float> Encoder::pitch3(const std::array<float, 360>& sig_in, const std::array<float, 360>& resid,
		float p2, float pavg){
	float Dth = 0;

	p2 = int(p2);
	std::pair<float, float> pr3{0.0f,0.0f};
	pr3 = pitch2((std::array<float, 360>::iterator)resid.begin(),p2);
	if(pr3.second >= 0.6){
		Dth = 0.5;
		if(pr3.first <= 100) Dth = 0.75;
		pr3 = double_ck(resid, pr3.first, Dth);
	}else{
		pr3 = fpr((std::array<float, 360>::iterator)sig_in.begin(), p2);
		if(pr3.second < 0.55) pr3.first = pavg;
		else{
			Dth = 0.7;
			if(pr3.first <= 100) Dth = 0.9;
			pr3 = double_ck(sig_in,pr3.first,Dth);
		}
	}
	if(pr3.second < 0.55)pr3.first =pavg;
	return pr3;
}

std::pair<float, float> Encoder::double_ck(const std::array<float, 360>& sig_in, int p, float Dth){
	int pmin = 20;
	int k, temp_pit;
	std::pair<float, float>npc = fpr((std::array<float, 360>::iterator)sig_in.begin(), round(p));
	std::pair<float, float>temp = {0.0f,0.0f};
	for(int n = 0; n < 7; n++){
		k = 9 - (n + 1);
		temp_pit = round(npc.first / k);
		if(temp_pit >= pmin){
			temp = fpr((std::array<float, 360>::iterator)sig_in.begin(), temp_pit);
			if(temp.first < 30) temp.second = double_ver(sig_in,temp.first, temp.second);
			if(temp.second > Dth * npc.second){
				npc = fpr((std::array<float, 360>::iterator)sig_in.begin(), round(temp.first));
				break;
			}
		}
	}
	if(npc.first < 30) npc.second = double_ver(sig_in, npc.first, npc.second);
	return npc;
}

float Encoder::double_ver(const std::array<float, 360>& sig_in, float pp, float cor_p){
	std::pair<float,float> npc {0.0f,0.0f};
	npc = fpr((std::array<float, 360>::iterator)sig_in.begin(), round(2*pp));
	if(npc.second < cor_p) cor_p = npc.second;
	return cor_p;
}

void Encoder::melp_gain(const std::array<float, 360>& s, float vpb1, float p2, std::array<float,2>& OG){
	int k = 1;
	int Ltmp = int(p2);
	int Lfr = int(p2);

	if(vpb1 > 0.6){
		while(Ltmp < 180){
			k++;
			Lfr = Ltmp;
			Ltmp = int(p2 * k);
		}
	}else{
		Lfr = 120;
	}
	int HL = round(Lfr/2);
	Lfr = HL * 2;
	OG[0] = calculateG(s, HL, Lfr, 91 - HL, 90 + HL);
	OG[1] = calculateG(s, HL, Lfr, 181 - HL, 180 + HL);
	if(OG[0] < 0) OG[0] = 0;
	if(OG[1] < 0) OG[1] = 0;
}

float Encoder::calculateG(const std::array<float, 360>& s, int HL, int Lfr, int start, int end) {
    std::vector<float> squaredS;
    squaredS.reserve(end - start + 1);

    // Calculate the squared values of s
    for (int i = start; i <= end; i++) {
        squaredS.push_back(s[i] * s[i]);
    }

    // Calculate the sum of squared values
    float sumSquaredS = std::accumulate(squaredS.begin(), squaredS.end(), 0.0);

    // Calculate the gain (G)
    float G = 10 * std::log10(0.01 + sumSquaredS / Lfr);

    return G;
}

void Encoder::melp_APU(std::pair<float, float>& pr3,float G2){

	if(pr3.second > 0.8 && G2 > 30){
		std::copy(buffer.begin() + 1, buffer.end(), buffer.begin());
		buffer[2] = pr3.first;
	}
	else std::transform(buffer.begin(), buffer.end(), buffer.begin(),
			[](float b) { return b * 0.95 + 2.5; });
	pavg = std::accumulate(buffer.begin(), buffer.end(), 0.0) / buffer.size();
}

void Encoder::lsf_clmp(const std::array<float, 10>& LSF, std::array<float, 10>& out){

	int dmin = 50;
	float d = 0, s1 = 0, s2 = 0, tmpv = 0;
	for(int i = 0; i < 10; i++){
		tmp[i] = LSF[i] * 4000 / 3.14159265358979323846;
	}

	for(int i = 0; i < 9; i++){
		d = tmp[i + 1] - tmp[i];
		if(d < dmin){
			s1 = (dmin - d) / 2;
			s2 = s1;
			if(i == 0 && (tmp[i] < dmin)) s1 = tmp[i] / 2;
			else if (i > 1){
				tmpv = tmp[i] - tmp[i - 1];
				if(tmpv < dmin) s1 = 0;
				else if(tmpv < 2 * dmin) s1 = (tmpv - dmin) / 2;
			}
			if(i == 8 && (tmp[i + 1] > 4000 - dmin)) s2 = (4000 - tmp[i + 1]) /2;
			else if (i < 8){
				tmpv = tmp[i + 2] - tmp[i + 1];
				if(tmpv < dmin) s2 = 0;
				else if(tmpv < 2 * dmin) s2 = (tmpv - dmin) / 2;
			}
			tmp[i] = tmp[i] - s1;
			tmp[i + 1] = tmp[i + 1] + s2;
		}
	}

	std::copy(tmp.begin(), tmp.end(), out.begin());
}

void Encoder::QGain(std::array<float, 2>& G, int G2p){

	float Q1,Q2, gain_max, gain_min, delta,temp;

	if(G[0] < 10) G[0] = 10;

	if(G[0] > 77) G[0] = 77;

	if(std::abs(G2p - G[1]) < 5 && std::abs(G[0] - 0.5 * (G[1] + G2p)) < 3) Q1 = 0;

	else{
		gain_max = std::max((float)G2p, G[1]) + 6;
		gain_min = std::min((float)G2p, G[1]) - 6;
		if(gain_min < 10) gain_min = 10;
		if(gain_max > 77) gain_max = 77;
		delta = (gain_max - gain_min) / 7;
		temp = G[0] - gain_min;
		Q1 = 1 + fix(temp / delta);
		if(Q1 > 7) Q1 = 7;
	}

	delta = 67 / 32;

	Q2 = fix((G[1] - 10)/delta);

	if(Q2 > 31) Q2 = 31;

	G[0] = Q1;
	G[1] = Q2;
}

int Encoder::melp_Qpitch(int G2p){
    int Q;

    // Perform quantization
    if (G2p >= 0) {
        Q = (G2p + 50) / 100 * 100;
    } else {
        Q = (G2p - 50) / 100 * 100;
    }

    return Q;
}

void Encoder::d_lsf(const std::array<float, 4>& codeword, std::array<float,10>& lsfs){
	for(int i = 0; i < 10; i++){
		tmp1[0][i] = st1[((int)codeword[0] - 1) * 10 + i];
		tmp1[1][i] = st2[0][((int)codeword[1] - 1) * 10 + i];
		tmp1[2][i] = st2[0][((int)codeword[2] - 1) * 10 + i];
		tmp1[3][i] = st2[0][((int)codeword[3] - 1) * 10 + i];
	}
	for(int i = 0; i < 10; i++){
		lsfs[i] = tmp1[0][i] + tmp1[1][i] + tmp1[2][i] + tmp1[3][i];
	}
}
