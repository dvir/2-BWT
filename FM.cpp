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

#include "FM.h"
#include "util.h"
#include "divsufsort.h"

#include <algorithm>
#include <vector>
#include <cmath>

int FM::verbose = 0;

FM::FM(uint8_t* T,uint32_t N,uint32_t samplerate = DEFAULT_SAMPLERATE) {
    this->samplerate = samplerate;
    this->n = N;
    
    /* 0 terminate */
    if(T[N-1] != 0) {
        T = (uint8_t*) safe_realloc(T,(N+1) * sizeof(uint8_t));
        T[N] = 0;
        this->n++;
    } 
    
    build(T,n,samplerate);
}

FM::FM() {
    this->n = 0;
    this->sigma = 0;
    this->I = 0;
    this->remap_reverse = NULL;
    this->T_bwt = NULL;
    this->T_bwt_reverse = NULL;
    this->sampled = NULL;
    this->suffixes = NULL;
}

FM::~FM() {
    free(remap_reverse);
    free(suffixes);
	free(positions);
    delete T_bwt;
    delete T_bwt_reverse;
    delete sampled;
}

uint32_t
FM::getSize() {
    uint32_t bytes = 0;
    
    bytes += sizeof(this->n);
    bytes += sizeof(this->samplerate);
    bytes += sizeof(this->sigma);
    bytes += sizeof(this->I);
    bytes += sizeof(this->remap);
    bytes += sizeof(this->C);
    bytes += this->sigma * sizeof(uint8_t); /* remap_reverse */
    bytes += ((n/samplerate)+1) * sizeof(uint32_t); /* suffixes */
	bytes += ((n/samplerate)+2) * sizeof(uint32_t); /* positions */
    bytes += this->sampled->getSize();
    bytes += this->T_bwt->getSize();
    bytes += this->T_bwt_reverse->getSize();
    
    return bytes;
}

float
FM::getSizeN() {
    uint32_t bytes = getSize();
    return (float)(bytes)/(float)(n);
}

uint8_t*
FM::reverse(uint8_t* X,uint32_t n) {
  uint8_t* X_ = (uint8_t*) malloc(n * sizeof(uint8_t));
  for (uint32_t i = 0; i < n-1; ++i) {
    X_[n-2-i] = X[i];
  }
  X_[n-1] = X[n-1];

  return X_;
}

uint8_t*
FM::remap0(uint8_t* T,uint32_t n) {
    uint8_t* X;
    uint32_t i,j,size=0;
    uint32_t freqs[size_uchar];
    
    for(i=0;i<size_uchar;i++) freqs[i]=0;
    for(i=0;i<n;i++) if(freqs[T[i]]++==0) size++;
    
    this->sigma=size;
    
    // remap alphabet
    if (freqs[0]>1) {i=1;sigma++;} //test if some character of T is zero, we already know that text[n-1]='\0'
    else i=0;

    remap_reverse = (uint8_t*) malloc(size*sizeof(uint8_t));
    for(j=0;j<size_uchar;j++) {
      if(freqs[j]!=0) {
        remap[j]=i;
        remap_reverse[i++]=j;
      }
    }
    // remap text
    X = (uint8_t*) malloc(n * sizeof(uint8_t));
    for(i=0;i<n-1;i++) // the last character must be zero
      X[i]=remap[T[i]];
    
    return X;
}

void
FM::build(uint8_t* T,uint32_t n,uint32_t samplerate) {
    uint8_t* X;
    uint8_t* X_bwt;
    int32_t* SA;
    uint8_t* X_;
    uint8_t* X_bwt_reverse;
    int32_t* SA_;
    uint32_t i,prev,tmp,start,stop;
    float elapsed;
    
    start = gettime();
	
	info("building index.");
    
    /* remap if 0 in text */
    info("- remapping alphabet.");
    X = remap0(T,n);
    free(T);
    
    /* create cumulative counts */
    info("- creating cumulative counts C[].");
    for (i=0;i<size_uchar+1;i++) C[i]=0;
    for (i=0;i<n;++i) C[X[i]]++;
    prev=C[0];C[0]=0;
    for (i=1;i<size_uchar+1;i++) {
      tmp = C[i];
      C[i]=C[i-1]+prev;
      prev = tmp;
    }
    
    /* perform k-BWT */
    info("- performing bwt.");
    SA = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
    if( divsufsort(X,SA,n) != 0 ) {
        fatal("error divsufsort");
    }
    
    /* sample SA for locate() */
    info("- sample SA locations.");
    suffixes = (uint32_t*) safe_malloc( ((n/samplerate)+1) * sizeof(uint32_t));
    BitString B(n);
    tmp = 0;
    for(i=0;i<n;i++) {
        if( SA[i] % samplerate == 0) {
            suffixes[tmp] = SA[i];
            B.setBit(i,true);
            tmp++;
        } else B.setBit(i,false);
    }
    /* enable rank on context vector */
    this->sampled = new BitSequenceRRR(B,RRR_SAMPLERATE);
	
	/* sample SA for display() */
	positions = (uint32_t*) safe_malloc( ((n/samplerate)+2) * sizeof(uint32_t));
    for (i=0;i<this->n;i++) {
        if (SA[i] % samplerate == 0) this->positions[SA[i]/samplerate] = i;
	}
    positions[(this->n-1)/samplerate+1] = positions[0];
	
    info("- creating bwt output.");
    X_bwt = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
    for(i=0;i<n;i++) {
        if(SA[i]==0) { 
            X_bwt[i] = X[n-1];
            this->I = i;
        } else X_bwt[i] = X[SA[i]-1];
    }
    free(SA);
    
    info("- create RRR wavelet tree over bwt.");
    MapperNone * map = new MapperNone();
    BitSequenceBuilder * bsb = new BitSequenceBuilderRRR(RRR_SAMPLERATE);
    T_bwt = new WaveletTreeNoptrs((uint32_t*)X_bwt,n,sizeof(uint8_t)*8,bsb,map,true);
    
    // create T_bwt_reverse
    /* perform k-BWT */
    info("- performing bwt (reverse).");
    X_ = reverse(X, n);

    SA_ = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
    if( divsufsort(X_,SA_,n) != 0 ) {
        fatal("error divsufsort (reverse)");
    }
    
    info("- creating bwt output (reverse).");
    X_bwt_reverse = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
    for(i=0;i<n;i++) {
        if(SA_[i]==0) { 
            X_bwt_reverse[i] = X_[n-1];
        } else X_bwt_reverse[i] = X_[SA_[i]-1];
    }
    free(SA_);
    free(X_);
    
    info("- create RRR wavelet tree over bwt (reverse).");
    MapperNone * map_ = new MapperNone();
    BitSequenceBuilder * bsb_ = new BitSequenceBuilderRRR(RRR_SAMPLERATE);
    T_bwt_reverse = new WaveletTreeNoptrs((uint32_t*)X_bwt_reverse,n,sizeof(uint8_t)*8,bsb_,map_,true);

    stop = gettime();
    elapsed = (float)(stop-start)/1000000;
    
    /* build aux data */
    info("build FM-Index done. (%.3f sec)",elapsed);
    
    uint32_t bytes;
    info("space usage:");
    bytes = sigma * sizeof(uint8_t);
    info("- remap_reverse: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = sizeof(this->C);
    info("- C: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = ((n/samplerate)+1) * sizeof(uint32_t);
    info("- Suffixes: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = ((n/samplerate)+2) * sizeof(uint32_t);
    info("- Positions: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = sampled->getSize();
    info("- Sampled: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = T_bwt->getSize();
    info("- T_bwt: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
    bytes = T_bwt_reverse->getSize();
    info("- T_bwt_reverse: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
	info("input Size n = %lu bytes\n",this->n);
	info("index Size = %lu bytes (%.2f n)",getSize(),getSizeN());
}


int32_t
FM::save(char* filename) {
    std::ofstream f;
    f.open(filename,std::ios::out | std::ios::binary); 
    
	info("writing FM Index to file '%s'",filename);
    if(f.is_open()) {
        f.write(reinterpret_cast<char*>(&samplerate),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&sigma),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&I),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&n),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(C),sizeof(uint32_t)*(size_uchar+1));
        f.write(reinterpret_cast<char*>(remap),sizeof(uint8_t)*size_uchar);
        f.write(reinterpret_cast<char*>(remap_reverse),sizeof(uint8_t)*sigma);
        f.write(reinterpret_cast<char*>(suffixes),sizeof(uint32_t)*((n/samplerate)+1));
		f.write(reinterpret_cast<char*>(positions),sizeof(uint32_t)*((n/samplerate)+2));
        T_bwt->save(f);
        T_bwt_reverse->save(f);
        sampled->save(f);
        f.close();
    } else return 1;
    
    return 0;
}

FM*
FM::load(char* filename) {
    FM* newIdx = new FM();
    std::ifstream f;
    f.open(filename,std::ios::in | std::ios::binary); 
    
    if(f.is_open()) {
        info("loading FM Index from file '%s'",filename);
        f.read(reinterpret_cast<char*>(&newIdx->samplerate),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->sigma),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->I),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->n),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(newIdx->C),sizeof(uint32_t)*(size_uchar+1));
        f.read(reinterpret_cast<char*>(newIdx->remap),sizeof(uint8_t)*size_uchar);
        newIdx->remap_reverse = (uint8_t*) safe_malloc(sizeof(uint8_t)*(newIdx->sigma));
        f.read(reinterpret_cast<char*>(newIdx->remap_reverse),sizeof(uint8_t)*newIdx->sigma);
        newIdx->suffixes = (uint32_t*) safe_malloc(sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+1));
        f.read(reinterpret_cast<char*>(newIdx->suffixes),sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+1));
        newIdx->positions = (uint32_t*) safe_malloc(sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+2));
        f.read(reinterpret_cast<char*>(newIdx->positions),sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+2));
        newIdx->T_bwt = WaveletTreeNoptrs::load(f);
        newIdx->T_bwt_reverse = WaveletTreeNoptrs::load(f);
        newIdx->sampled = BitSequenceRRR::load(f);
        f.close();
        info("samplerate = %d",newIdx->samplerate);
        info("sigma = %d",newIdx->sigma);
        info("I = %d",newIdx->I);
        info("n = %d",newIdx->n);
    } else {
        delete newIdx;
        return NULL;
    }
        
    return newIdx;
}

uint32_t
FM::count(uint8_t* pattern,uint32_t m) {
    uint8_t c = remap[pattern[m-1]]; /* map pattern to our alphabet */
    uint32_t i=m-1;
    uint32_t j = 1;
    
    uint32_t sp = C[c]; /* starting range in M from p[m-1] */
    uint32_t ep = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (sp<=ep && i>=1) { 
      c = remap[pattern[--i]]; /* map pattern to our alphabet */
      sp = C[c] + T_bwt->rank(c, sp-1); /* LF Mapping */
      ep = C[c] + T_bwt->rank(c, ep)-1; /* LF Mapping */
      j++;
    }
    
    if (sp<=ep) {
      return ep-sp+1;
    } else {
      return 0;
    }
}

inline static size_t Occ(uint8_t c, uint32_t range, WaveletTreeNoptrs* T) {
  if (range == 4294967295) return 0;
  if (range == 0) return 0;

  size_t sz = T->rank(c, range);
  return sz;
}
inline static size_t OccLT(uint8_t c, uint32_t range, WaveletTreeNoptrs* T) {
  if (range == 4294967295) return 0;
  if (range == 0) return 0;

  size_t sum = 0;
  for (uint8_t i = 0; i < c; ++i) {
    sum += Occ(i, range, T);
  }

  return sum;
}
inline static SA_intervals reverseIntervals(SA_intervals ivls) {
  return {ivls.r_ivl, ivls.ivl};
}

SA_intervals FM::updateForwardBackward(SA_intervals ivls, uint8_t c, WaveletTreeNoptrs* T) {
    /*
     * l' <- l' + OccLTf(a, u) - OccLTf(a, l - 1)
     * u' <- l' + Occf(a, u) - Occf(a, l - 1) - 1
     */
    SA_interval ivl = ivls.ivl;
    SA_interval r_ivl = ivls.r_ivl;
    r_ivl.sp = r_ivl.sp + OccLT(c, ivl.ep, T) - OccLT(c, ivl.sp - 1, T);
    r_ivl.ep = r_ivl.sp + Occ(c, ivl.ep, T) - Occ(c, ivl.sp - 1, T) - 1;
    ivl = updateBackward(ivl, c, T);
    return {ivl, r_ivl};
}

SA_interval FM::updateBackward(SA_interval ivl, uint8_t c, WaveletTreeNoptrs* T) {
    /*
     * l <- Cx(a) + Occx(a, l-1)
     * u <- Cx(a) + Occx(a, u)-1
     */
    ivl.sp = C[c] + Occ(c, ivl.sp - 1, T);
    ivl.ep = C[c] + Occ(c, ivl.ep, T) - 1;
    return ivl;
}

SA_intervals* FM::backwardSearch(uint8_t* pattern, uint32_t s, uint32_t e, SA_intervals ivls) {
    if (s > e) {
      return new SA_intervals{ivls.ivl, ivls.r_ivl};
    }

    for (uint32_t i = e; 
         ivls.ivl.sp <= ivls.ivl.ep 
         && ivls.r_ivl.sp <= ivls.r_ivl.ep 
         && i >= s
         && i <= e; 
         --i) {
      ivls = updateForwardBackward(ivls, remap[pattern[i]], T_bwt);
    }

    if (ivls.ivl.sp > ivls.ivl.ep || ivls.r_ivl.sp > ivls.r_ivl.ep) {
      return NULL;
    }

   return new SA_intervals{ivls.ivl, ivls.r_ivl};
}

SA_intervals* FM::forwardSearch(uint8_t* pattern, uint32_t s, uint32_t e, SA_intervals ivls) {
    if (s > e) {
      return new SA_intervals{ivls.ivl, ivls.r_ivl};
    }

    ivls = reverseIntervals(ivls);
    for (uint32_t i = s; 
         ivls.ivl.sp <= ivls.ivl.ep 
         && ivls.r_ivl.sp <= ivls.r_ivl.ep 
         && i >= s
         && i <= e; 
         ++i) {
      ivls = updateForwardBackward(ivls, remap[pattern[i]], T_bwt_reverse);
    }

    if (ivls.ivl.sp > ivls.ivl.ep || ivls.r_ivl.sp > ivls.r_ivl.ep) {
      return NULL;
    }

   return new SA_intervals{ivls.r_ivl, ivls.ivl};
}

std::vector<SA_intervals>
FM::locate1(uint8_t* pattern, uint32_t m) {
  uint8_t c;
  std::vector<SA_intervals> matching_intervals;
  SA_interval ivl, r_ivl;
  SA_intervals ivls;

  uint32_t middle = ceil(m/2);

  /*
   * Two possible cases for 1 error:
   * str: abcdef
   * - assume error is in first half.
   *   - backward search def
   *   - try all options for abc errors
   * - assume error is in second half.
   *   - forward search abc
   *   - try all options for def errors
   */
  SA_intervals* res_ptr = NULL; 
  SA_intervals res; 

  c = remap[pattern[m-1]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // backward search [middle + 1, m - 2]
  // m-2 because we already searched for m-1
  res_ptr = backwardSearch(pattern, middle + 1, m - 2, ivls); 
  if (res_ptr) {
    SA_intervals ivls_second_half = *res_ptr;

    // choose a remaining character to be the one with the error
    for (uint32_t i = middle; i >= 0 && i <= m - 1; --i) {
      // backward search up until the character with the error
      // so [i + 1, middle]
      res_ptr = backwardSearch(pattern, i + 1, middle, ivls_second_half); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_error = *res_ptr;

      // choose which character it should actually represent
      for (uint32_t c = 1; c < sigma; ++c) {
        if (c == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(ivls_upto_error, c, T_bwt);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction = res;

        // backward search the remaining of the pattern
        // so [0, i - 1]
        if (i == 0) {
          // nothing to search for. we have a match!
          matching_intervals.push_back(ivls_with_correction);
        } else {
          res_ptr = backwardSearch(pattern, 0, i - 1, ivls_with_correction); 
          if (!res_ptr) continue;
          SA_intervals ivls_full_fixed = *res_ptr;

          // if we got here, we have a match!
          matching_intervals.push_back(ivls_full_fixed);
        }
      }
    }
  }

  /*
   * - assume error is in second half.
   *   - forward search abc
   *   - try all options for def errors
   */

  c = remap[pattern[0]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // forward search [1, middle]
  res_ptr = forwardSearch(pattern, 1, middle, ivls); 
  if (res_ptr) {
    SA_intervals ivls_first_half = *res_ptr;

    // choose a remaining character to be the one with the error
    for (uint32_t i = middle + 1; i < m; ++i) {
      // forward search up until the character with the error
      // so [middle + 1, i - 1]
      res_ptr = forwardSearch(pattern, middle + 1, i - 1, ivls_first_half); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_error = *res_ptr;

      // choose which character it should actually represent
      for (uint32_t c = 1; c < sigma; ++c) {
        if (c == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(reverseIntervals(ivls_upto_error), c, T_bwt_reverse);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction = reverseIntervals(res);
        // forward search the remaining of the pattern
        // so [i + 1, m - 1]
        res_ptr = forwardSearch(pattern, i + 1, m - 1, ivls_with_correction); 
        if (!res_ptr) continue;
        SA_intervals ivls_full_fixed = *res_ptr;

        // if we got here, we have a match!
        matching_intervals.push_back(ivls_full_fixed);
      }
    }
  }

  return matching_intervals;
}

std::vector<SA_intervals>
FM::locate2(uint8_t* pattern, uint32_t m) {
  uint8_t c;
  std::vector<SA_intervals> matching_intervals;
  SA_interval ivl, r_ivl;
  SA_intervals ivls;
  SA_intervals* res_ptr = NULL; 
  SA_intervals res; 

  // split separators
  // note: s1 and s2 are the beginning of their respective parts
  uint32_t s1 = floor(m/3);
  uint32_t s2 = m - s1;

  /*
   * Four possible cases for 2 error:
   * str: abcdefghi
   * 200
   * 020
   * 110
   * A. assume errors are in the first two parts.
   *   - backward search ghi
   *   - try all options for abcdef with 2 errors
   * 002
   * B. assume both errors are in the last part.
   *   - forward search abcdef
   *   - try all options for ghi with 2 errors
   * 011
   * C. assume errors are in the second and third part respectively.
   *   - forward search abc
   *   - try all options for def with 1 error
   *   - try all options for ghi with 1 error
   * 101
   * D. assume errors are in the first and last part respectively.
   *   - forward search def
   *   - try all options for abc with 1 errors (backward search)
   *   - try all options for ghi with 1 errors (forward search)
   */



  /*
   * 200
   * 020
   * 110
   * A. assume errors are in the first two parts.
   *   - backward search 'ghi'
   *   - try all options for abcdef with 2 errors
   */
  c = remap[pattern[m-1]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // backward search [s2, m - 1]
  // m-2 because we already searched for m-1
  res_ptr = backwardSearch(pattern, s2, m - 2, ivls); 
  if (res_ptr) {
    SA_intervals ivls_last_part = *res_ptr;

    // choose a remaining character to be the one with the error
    for (uint32_t i = s2 - 1; i >= 1 && i <= m - 1; --i) {
      // backward search up until the character with the error
      // so [i + 1, s2 - 1]
      res_ptr = backwardSearch(pattern, i + 1, s2 - 1, ivls_last_part); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_i = *res_ptr;

      // choose which character i should actually represent
      for (uint32_t e1 = 1; e1 < sigma; ++e1) {
        if (e1 == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(ivls_upto_i, e1, T_bwt);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction_for_i = res;

        for (uint32_t j = i - 1; j >= 0 && j <= m - 1; --j) {
          // backward search up until the character with the -NEXT- error
          // so [j + 1, i - 1]
          res_ptr = backwardSearch(pattern, j + 1, i - 1, ivls_with_correction_for_i); 
          if (!res_ptr) continue;
          SA_intervals ivls_upto_j = *res_ptr;

          // choose which character j should actually represent
          for (uint32_t e2 = 1; e2 < sigma; ++e2) {
            if (e2 == remap[pattern[j]]) {
              continue;
            }

            // search for the specific character we are trying
            res = updateForwardBackward(ivls_upto_j, e2, T_bwt);
            if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
              continue;
            }
            SA_intervals ivls_with_correction_for_j = res;

            // backward search the remaining of the pattern
            // so [0, j - 1]
            if (j == 0) {
              // nothing to search for. we have a match!
              matching_intervals.push_back(ivls_with_correction_for_j);
            } else {
              res_ptr = backwardSearch(pattern, 0, j - 1, ivls_with_correction_for_j); 
              if (!res_ptr) continue;
              SA_intervals ivls_full_fixed = *res_ptr;

              // if we got here, we have a match!
              matching_intervals.push_back(ivls_full_fixed);
            }
          }
        }
      }
    }
  }

  /*
   * 002
   * B. assume both errors are in the last part.
   *   - forward search abcdef
   *   - try all options for ghi with 2 errors
   */
  c = remap[pattern[0]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // forward search [1, s2 - 1]
  // starting from 1 because we already searched for 0
  res_ptr = forwardSearch(pattern, 1, s2 - 1, ivls); 
  if (res_ptr) {
    SA_intervals ivls_first_and_second = *res_ptr;

    // choose a remaining character to be the one with the error
    for (uint32_t i = s2; i <= m - 2; ++i) {
      // forward search up until the character with the error
      // so [s2, i - 1]
      res_ptr = forwardSearch(pattern, s2, i - 1, ivls_first_and_second); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_i = *res_ptr;

      // choose which character i should actually represent
      for (uint32_t e1 = 1; e1 < sigma; ++e1) {
        if (e1 == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(reverseIntervals(ivls_upto_i), e1, T_bwt_reverse);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction_for_i = reverseIntervals(res);

        for (uint32_t j = i + 1; j <= m - 1; ++j) {
          // forward search up until the character with the -NEXT- error
          // so [i + 1, j - 1]
          res_ptr = forwardSearch(pattern, i + 1, j - 1, ivls_with_correction_for_i); 
          if (!res_ptr) continue;
          SA_intervals ivls_upto_j = *res_ptr;

          // choose which character j should actually represent
          for (uint32_t e2 = 1; e2 < sigma; ++e2) {
            if (e2 == remap[pattern[j]]) {
              continue;
            }

            // search for the specific character we are trying
            res = updateForwardBackward(reverseIntervals(ivls_upto_j), e2, T_bwt_reverse);
            if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
              continue;
            }
            SA_intervals ivls_with_correction_for_j = reverseIntervals(res);

            // forward search the remaining of the pattern
            // so [j + 1, m - 1]
            res_ptr = forwardSearch(pattern, j + 1, m - 1, ivls_with_correction_for_j); 
            if (!res_ptr) continue;
            SA_intervals ivls_full_fixed = *res_ptr;

            // if we got here, we have a match!
            matching_intervals.push_back(ivls_full_fixed);
          }
        }
      }
    }
  }

  /*
   * 011
   * C. assume errors are in the second and third part respectively.
   *   - forward search abc
   *   - try all options for def with 1 error
   *   - try all options for ghi with 1 error
   */
  c = remap[pattern[0]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // forward search [1, s1 - 1]
  // starting from 1 because we already searched for 0
  res_ptr = forwardSearch(pattern, 1, s1 - 1, ivls); 
  if (res_ptr) {
    SA_intervals ivls_first = *res_ptr;

    // choose a character in the second part to be the one with the error
    for (uint32_t i = s1; i <= s2 - 1; ++i) {
      // forward search up until the character with the error
      // so [s1, i - 1]
      res_ptr = forwardSearch(pattern, s1, i - 1, ivls_first); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_i = *res_ptr;

      // choose which character i should actually represent
      for (uint32_t e1 = 1; e1 < sigma; ++e1) {
        if (e1 == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(reverseIntervals(ivls_upto_i), e1, T_bwt_reverse);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction_for_i = reverseIntervals(res);

        // forward search to the end of the second part 
        // so [i + 1, s2 - 1]
        res_ptr = forwardSearch(pattern, i + 1, s2 - 1, ivls_with_correction_for_i); 
        if (!res_ptr) continue;
        SA_intervals ivls_second = *res_ptr;

        for (uint32_t j = s2; j <= m - 1; ++j) {
          // forward search up until the character with the error
          // so [s2, j - 1]
          res_ptr = forwardSearch(pattern, s2, j - 1, ivls_second); 
          if (!res_ptr) continue;
          SA_intervals ivls_upto_j = *res_ptr;

          // choose which character j should actually represent
          for (uint32_t e2 = 1; e2 < sigma; ++e2) {
            if (e2 == remap[pattern[j]]) {
              continue;
            }

            // search for the specific character we are trying
            res = updateForwardBackward(reverseIntervals(ivls_upto_j), e2, T_bwt_reverse);
            if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
              continue;
            }
            SA_intervals ivls_with_correction_for_j = reverseIntervals(res);

            // forward search to the end of the third part 
            // so [j + 1, m - 1]
            res_ptr = forwardSearch(pattern, j + 1, m - 1, ivls_with_correction_for_j); 
            if (!res_ptr) continue;
            SA_intervals ivls_full_fixed = *res_ptr;

            // if we got here, we have a match!
            matching_intervals.push_back(ivls_full_fixed);
          }
        }
      }
    }
  }

  /*
   * 101
   * D. assume errors are in the first and last part respectively.
   *   - forward search def
   *   - try all options for abc with 1 errors (backward search)
   *   - try all options for ghi with 1 errors (forward search)
   */
  c = remap[pattern[s1]];
  ivl = {C[c], C[c+1]-1};
  r_ivl = ivl;
  ivls = {ivl, r_ivl};

  // forward search [s1 + 1, s2 - 1]
  // starting from s1 + 1 because we already searched for s1
  res_ptr = forwardSearch(pattern, s1 + 1, s2 - 1, ivls); 
  if (res_ptr) {
    SA_intervals ivls_second = *res_ptr;

    // choose a character in the third part to be the one with the error
    for (uint32_t i = s2; i <= m - 1; ++i) {
      // forward search up until the character with the error
      // so [s2, i - 1]
      res_ptr = forwardSearch(pattern, s2, i - 1, ivls_second); 
      if (!res_ptr) continue;
      SA_intervals ivls_upto_i = *res_ptr;

      // choose which character i should actually represent
      for (uint32_t e1 = 1; e1 < sigma; ++e1) {
        if (e1 == remap[pattern[i]]) {
          continue;
        }

        // search for the specific character we are trying
        res = updateForwardBackward(reverseIntervals(ivls_upto_i), e1, T_bwt_reverse);
        if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
          continue;
        }
        SA_intervals ivls_with_correction_for_i = reverseIntervals(res);

        // forward search to the end of the third part 
        // so [i + 1, m - 1]
        res_ptr = forwardSearch(pattern, i + 1, m - 1, ivls_with_correction_for_i); 
        if (!res_ptr) continue;
        SA_intervals ivls_second_and_third = *res_ptr;

        // choose a character in the first part to be the one with the error
        for (uint32_t j = s1 - 1; j >= 0 && j <= m - 1; --j) {
          // backward search up until the character with the error
          // so [j + 1, s1 - 1]
          res_ptr = backwardSearch(pattern, j + 1, s1 - 1, ivls_second_and_third); 
          if (!res_ptr) continue;
          SA_intervals ivls_upto_j = *res_ptr;

          // choose which character j should actually represent
          for (uint32_t e2 = 1; e2 < sigma; ++e2) {
            if (e2 == remap[pattern[j]]) {
              continue;
            }

            // search for the specific character we are trying
            // backward search
            res = updateForwardBackward(ivls_upto_j, e2, T_bwt);
            if (res.ivl.sp > res.ivl.ep || res.r_ivl.sp > res.r_ivl.ep) {
              continue;
            }
            SA_intervals ivls_with_correction_for_j = res;

            // backward search to the end of the first part 
            // so [0, j - 1]
            if (j == 0) {
              // nothing to search for. we have a match!
              matching_intervals.push_back(ivls_with_correction_for_j);
            } else {
              res_ptr = backwardSearch(pattern, 0, j - 1, ivls_with_correction_for_j); 
              if (!res_ptr) continue;
              SA_intervals ivls_full_fixed = *res_ptr;

              // if we got here, we have a match!
              matching_intervals.push_back(ivls_full_fixed);
            }
          }
        }
      }
    }
  }

  return matching_intervals;
}

uint32_t*
FM::getLocations(SA_interval ivl, uint32_t* matches) {
  *matches = ivl.ep - ivl.sp + 1;
  uint32_t* locations = (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
  uint32_t locate=0;
  uint32_t i=ivl.sp;
  int32_t j,dist,rank;
  uint8_t c;
  while (i<=ivl.ep) {
      j=i,dist=0;
      while (!sampled->access(j)) {
          c = T_bwt->access(j);
          rank = T_bwt->rank(c,j)-1;
          j = C[c]+rank; // LF-mapping
          ++dist;
      }
      locations[locate]=suffixes[sampled->rank1(j)-1]+dist;
      locate++;
      ++i;
  }
  /* locations are in SA order */
  std::sort(locations,locations+(*matches));
  return locations;
}

uint32_t*
FM::locate(uint8_t* pattern,uint32_t m,uint32_t* matches) {
    uint32_t* locations;
    uint8_t c =  remap[pattern[m-1]];
    uint32_t i=m-1;
    
    /* count occs */
    uint32_t sp = C[c];
    uint32_t ep = C[c+1]-1;
    while (sp<=ep && i>=1) {
      c =  remap[pattern[--i]];
      sp = C[c] + T_bwt->rank(c, sp-1);
      ep = C[c] + T_bwt->rank(c, ep)-1;
    }
    
    if (sp<=ep) {
        /* determine positions */
        *matches = ep-sp+1;
        uint32_t locate=0;
        locations= (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
        i=sp;
        int32_t j,dist,rank;
        while (i<=ep) {
            j=i,dist=0;
            while (!sampled->access(j)) {
                c = T_bwt->access(j);
                rank = T_bwt->rank(c,j)-1;
                j = C[c]+rank; // LF-mapping
                ++dist;
            }
            locations[locate]=suffixes[sampled->rank1(j)-1]+dist;
            locate++;
            ++i;
        }
        /* locations are in SA order */
        std::sort(locations,locations+(*matches));
        return locations;
    } else {
      /* no matches */
      *matches = 0;
      return NULL;
    }
    
    return locations;
}

uint8_t*
FM::extract(uint32_t start,uint32_t stop)
{
    uint8_t* T;
	uint32_t m,j,skip,todo,dist;
	uint8_t c;
	
	/* last text pos is n-2 */
	if(stop > (this->n-1) ) stop = n-2; 
    if(start > stop) {
		return NULL;
	}
	
	m = stop-start+1; /* snippet len */
	T = (uint8_t*) safe_malloc( (m+1) * sizeof(uint8_t)  );
	
	/* determine start pos of backwards search */
	j = positions[(stop/samplerate)+1];
	
	/* determine distance from start pos to the text snippet we want */
	if ((stop/samplerate+1) == ((n-1)/samplerate+1)) 
	   skip = n-2 - stop;
	else 
	   skip = (samplerate-stop)%samplerate-1;
	   
	/* start the backwards search */
	todo = m;
	dist = 0;
	while(todo>0) {
		c = T_bwt->access(j);
		j = C[c] + T_bwt->rank(c,j)-1;
		
		/* check if we are at the snippet */
		if(dist>=skip) {
			c = remap_reverse[c];
			T[todo-1] = c;
			todo--;
		}
		dist++;
	}
	
	/* terminate */
	T[m] = 0;
	
    return T;
}

uint8_t*
FM::reconstructText(uint32_t* size)
{
    uint8_t* T;
    uint8_t c;
    uint32_t j,i;
    
    T = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
    
    j = I; /* I is sa[I] = 0 -> last sym in T */
    for(i=0;i<n;i++) {
        c = T_bwt->access(j); /* L[j] = c */
        T[n-i-1] = remap_reverse[c]; /* undo sym mapping */
        j = C[c]+T_bwt->rank(c,j)-1; /* LF-mapping: j = LF[j] */
    }
    
	if(T[n-1] == 0) *size = n-1;
    else *size = n;
    
    return T;
}

