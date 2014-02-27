/*  FM-Index - Text Index
 *  Copyright (C) 2011  Matthias Petri
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. */
 
#ifndef FM_H
#define	FM_H

#include <stdarg.h>
#include <vector>

#include "util.h"

/* libcds includes */
#include "libcdsBasics.h"
#include "BitString.h"
#include "BitSequence.h"
#include "BitSequenceRG.h"
#include "Sequence.h"
#include "WaveletTreeNoptrs.h"
#include "Mapper.h"
#include "BitSequenceBuilder.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

#define DEFAULT_SAMPLERATE      64
#define RRR_SAMPLERATE			20

struct SA_interval {
  uint32_t sp, ep;

  void setRange(uint32_t sp, uint32_t ep) {
    this->sp = sp;
    this->ep = ep;
  }
};
struct SA_intervals {
  SA_interval ivl, r_ivl;

  SA_intervals(): ivl(), r_ivl() {}
  SA_intervals(SA_interval ivl, SA_interval r_ivl): ivl(ivl), r_ivl(r_ivl) {}
  void setIvls(SA_interval ivl, SA_interval r_ivl) {
    this->ivl = ivl;
    this->r_ivl = r_ivl;
  }
};

class FM {
public:
    FM(uint8_t* T,uint32_t n,uint32_t samplerate);
    void build(uint8_t* T,uint32_t n,uint32_t samplerate);
    static FM* load(char* filename);
    int32_t save(char* filename);
    uint8_t* remap0(uint8_t* T,uint32_t n);
    uint8_t* reverse(uint8_t* T,uint32_t n);
    uint32_t count(uint8_t* pattern,uint32_t m);
    uint32_t* locate(uint8_t* pattern,uint32_t m,uint32_t* matches);
    std::vector<SA_intervals> locate1(uint8_t* pattern, uint32_t m);
    std::vector<SA_intervals> locate2(uint8_t* pattern, uint32_t m);
    uint32_t* getLocations(SA_interval ivl, uint32_t* matches);
    SA_intervals updateForwardBackward(SA_intervals ivls, uint8_t c, WaveletTreeNoptrs* T);
    SA_interval updateBackward(SA_interval ivl, uint8_t c, WaveletTreeNoptrs* T);
    SA_intervals* backwardSearch(uint8_t* pattern, uint32_t s, uint32_t e, SA_intervals ivls);
    SA_intervals* forwardSearch(uint8_t* pattern, uint32_t s, uint32_t e, SA_intervals ivls);
    uint32_t getSize();
    uint8_t* extract(uint32_t start,uint32_t stop);
    uint8_t* reconstructText(uint32_t* n);
    float getSizeN();
    virtual ~FM();
public:
	static void info(const char *format,...)
	{
		if(FM::verbose == 1) {
			va_list vargs;
			va_start (vargs, format);
			vfprintf (stderr, format, vargs);
			fprintf (stderr, "\n");
		}
	}
	static int verbose;
private:
    FM();
private:
  uint32_t sigma;
  uint32_t samplerate;
  int32_t I;
  uint32_t n;
  uint32_t C[size_uchar+1];
  uint8_t remap[size_uchar];
  uint8_t* remap_reverse;
  uint32_t* suffixes;
  uint32_t* positions;
  BitSequence* sampled;
  WaveletTreeNoptrs *T_bwt;
  WaveletTreeNoptrs *T_bwt_reverse;
  Sequence *N;
  Sequence *M;
};


#endif	/* FM_H */

