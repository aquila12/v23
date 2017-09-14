#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#define SAMPLERATE      44100
#define F_MARK_FREQ     1300
#define F_SPACE_FREQ    2100
#define B_MARK_FREQ     390
#define B_SPACE_FREQ    450

#define F_BIT_RATE      1200
#define B_BIT_RATE      75

#define FRAME_MASK      0b11000000001 // Idle, Start, 8, Stop
#define FRAME_PATTERN   0b10000000001 // Frame format (masked)
#define FRAME_OFFSET    1             // Number of bits after data
#define FRAME_SIZE      10            // Number of bits in a frame

struct maf {
  int16_t *buf;
  size_t N;
  size_t p;
  int32_t sum;
};

struct iq {
  int16_t *buf;
  size_t N;
  size_t i;
  size_t q;
};

int16_t* make_buffer(size_t N)
{
  return (int16_t*)calloc(N, sizeof(int16_t));
}

void maf_process(maf& maf, int16_t *samples_in, int16_t *samples_out,
  size_t n_samples, bool nodivide=false)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    maf.sum        -= maf.buf[maf.p];
    maf.buf[maf.p] =  samples_in[i];
    maf.sum        += maf.buf[maf.p];

    if(nodivide)
    {
        if(maf.sum > 32767) samples_out[i]=32767;
        else if(maf.sum < -32767) samples_out[i]=-32767;
        else samples_out[i] = maf.sum;
    }
    else
        samples_out[i] = ( maf.sum + (int32_t)maf.N/2 ) / (int32_t)maf.N;

    ++maf.p;
    if(maf.p >= maf.N) maf.p -= maf.N;
  }
}

void iq_get_samples(iq& iq, int16_t *i_samples_out, int16_t* q_samples_out,
  size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    i_samples_out[i] = iq.buf[iq.i];
    q_samples_out[i] = iq.buf[iq.q];

    ++iq.i;
    ++iq.q;
    if(iq.i >= iq.N) iq.i -= iq.N;
    if(iq.q >= iq.N) iq.q -= iq.N;
  }
}

bool iq_init(iq& iq, size_t N)
{
  iq.buf = make_buffer(N);
  if(!iq.buf) return false;

  iq.N = N;
  iq.i = 0;
  iq.q = N/4;

  for(size_t i=0; i<iq.N; ++i)
  {
    double x = 2.0 * M_PI * (double)i / (double)N;
    iq.buf[i] = (int16_t)(32767.0 * cos(x));
  }

  return true;
}

bool maf_init(maf& maf, size_t N)
{
  maf.buf = make_buffer(N);
  if(!maf.buf) return false;

  maf.N   = N;
  maf.p   = 0;
  maf.sum = 0;

  return true;
}

void mul_samples(int16_t *samples_a, int16_t *samples_b,
  int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    int32_t product = ((int32_t)samples_a[i] * (int32_t)samples_b[i]) / (int32_t)32768;
    if(product >  32767)
    {
        printf("mul: clipped\n");
        product =  32767;
    }
    else if(product < -32767)
    {
        printf("mul: clipped\n");
        product = -32767;
    }
    samples_out[i] = (int16_t)product;
  }
}

//NB: Halves magnitude
void sub_samples(int16_t *samples_a, int16_t *samples_b,
  int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    samples_out[i] = ( samples_a[i]/2 - samples_b[i]/2 );
  }
}

void sgn_samples(int16_t *samples_in, int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    samples_out[i] = ( samples_in[i] > 0 ) - ( samples_in[i] < 0 );
  }
}

void mag_iq_samples(int16_t *samples_i, int16_t *samples_q,
  int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    int32_t sumsqr;
    sumsqr  = ((int32_t)samples_i[i] * (int32_t)samples_i[i]);
    sumsqr += ((int32_t)samples_q[i] * (int32_t)samples_q[i]);
    if(sumsqr >  1073676289)
    {
        printf("mag: clipped\n");
      samples_out[i] = 32767;
    }
    else
      samples_out[i] = (int16_t)sqrt(sumsqr);
  }
}

void dump_buf(int16_t *samples_i, size_t n_samples)
{
  write(2, samples_i, n_samples * sizeof(int16_t));
}

size_t get_input_samples(int16_t *buf, size_t n) {
  return read(0, (char*)buf, n * sizeof(int16_t)) / sizeof(int16_t);
}

int main()
{
  // Set up some working constants
  // By default, input is the backward channel.
  // If reversed, input is the forward channel.
  bool forward = true;
  size_t mark_period_samples  = SAMPLERATE / ( forward ? F_MARK_FREQ  : B_MARK_FREQ );
  size_t space_period_samples = SAMPLERATE / ( forward ? F_SPACE_FREQ : B_SPACE_FREQ );
  size_t bit_period_samples   = SAMPLERATE / ( forward ? F_BIT_RATE   : B_BIT_RATE );
  
  printf("Demodulating the %s channel\n", forward ? "FORWARD" : "BACKWARD");
  printf("Mark period:  %d samples\n", mark_period_samples);
  printf("Space period: %d samples\n", space_period_samples);
  printf("Bit period:   %d samples\n", bit_period_samples);
    
  // Local mark / space oscillators
  iq oscMark, oscSpace;
  
  // Working registers
  int16_t *bufIn, *bufOut;
  int16_t *bufWork;
  int16_t *bufI, *bufQ;
  int16_t *bufMark, *bufSpace;
  int16_t out_shift;  // Raw serial shift-register
  int frame_hold;   // How many bits to hold off for
  
  maf mafMarkI, mafMarkQ, mafSpaceI, mafSpaceQ, mafOut, mafBit;
  size_t N = 1024; // Maximum samples we can take at once
  
  size_t bit_wait = bit_period_samples;  // Samples left until we read a bit
  
  // Set up the IQ and MAF units
  if(! (
    iq_init(oscMark,  mark_period_samples)  &&
    iq_init(oscSpace, space_period_samples) &&
    maf_init(mafMarkI,  bit_period_samples) &&
    maf_init(mafMarkQ,  bit_period_samples) &&
    maf_init(mafSpaceI, bit_period_samples) &&
    maf_init(mafSpaceQ, bit_period_samples) &&
    maf_init(mafOut,    bit_period_samples) &&
    maf_init(mafBit,    bit_period_samples) )) {
    
    printf("Failed to initialize MAFs / IQs\n");
    exit(1);
  }
  
  bufIn    = make_buffer(N);
  bufOut   = make_buffer(N);
  bufWork  = make_buffer(N);
  bufI     = make_buffer(N);
  bufQ     = make_buffer(N);
  bufMark  = make_buffer(N);
  bufSpace = make_buffer(N);
  
  if(!( bufIn && bufOut && bufWork && bufI && bufQ && bufMark && bufSpace))
  {
      printf("Failed to allocate buffers\n");
      exit(1);
  }
  
  printf("Initialized.  Processing samples.\n");

  size_t n;     // Number of samples we have this time
  int state=0;  // What was the last state
  while(n = get_input_samples(bufIn, N))
  {
    //printf("Got %d of %d samples...\n", n, N);
    // Mix and filter the Mark oscillator
    iq_get_samples(oscMark, bufI, bufQ, n);
    mul_samples(bufIn, bufI, bufWork, n);
    maf_process(mafMarkI, bufWork, bufI, n);
    mul_samples(bufIn, bufQ, bufWork, n);
    maf_process(mafMarkQ, bufWork, bufQ, n);
    mag_iq_samples(bufI, bufQ, bufMark, n);

    // Mix and filter the Space oscillator
    iq_get_samples(oscSpace, bufI, bufQ, n);
    mul_samples(bufIn, bufI, bufWork, n);
    maf_process(mafSpaceI, bufWork, bufI, n);
    mul_samples(bufIn, bufQ, bufWork, n);
    maf_process(mafSpaceQ, bufWork, bufQ, n);
    mag_iq_samples(bufI, bufQ, bufSpace, n);

    // Subtract the magnitudes and feed output MAF
    sub_samples(bufMark, bufSpace, bufWork, n);
    maf_process(mafOut, bufWork, bufOut, n);
    
    // Bit sampling
    sgn_samples(bufOut, bufWork, n);
    maf_process(mafBit, bufWork, bufOut, n, true);
    
    // Run through the output samples
    int last;
    for(size_t i=0; i<n; ++i)
    {
        last = state;
        state = (bufOut[i] > 0);
        
        // Edge detected - re-align
        if(last != state) {
            //printf("Realign\n");
            bit_wait = bit_period_samples / 2;
        }
        
        if(--bit_wait <= 0)
        {
            out_shift <<= 1;
            out_shift += state;
            
            if(frame_hold > 0)                                  // Frame Hold-off
            {
                --frame_hold;
            }
            else if((out_shift & FRAME_MASK) == FRAME_PATTERN)    // Frame is valid
            {
                // TODO: Parity checking ?
                
                // Get raw byte
                int32_t b = 0xff & (out_shift >> FRAME_OFFSET);
                
                // Reverse bits in byte (LSB is first transmitted)
                // http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv
                b = (b * 0x0202020202ULL & 0x010884422010ULL) % 1023;
                
                // TODO: Output
                printf("%c", b);
                
                frame_hold = FRAME_SIZE - 1;
            }
            
            bit_wait += bit_period_samples;
        }
    }
  }
  
  free(bufIn);
  free(bufOut);
  free(bufWork);
  free(bufI);
  free(bufQ);
  free(bufMark);
  free(bufSpace);
  free(mafMarkI.buf);
  free(mafMarkQ.buf);
  free(mafSpaceI.buf);
  free(mafSpaceQ.buf);
  free(mafOut.buf);
  free(mafBit.buf);
  free(oscMark.buf);
  free(oscSpace.buf);
  
  return 0;
}
