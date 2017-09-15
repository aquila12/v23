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

#define MAX_SKEW        5

#define DEF_FRAME_FORMAT    "10dddddddp1"
#define FRAME_LSB_FIRST 1
#define FRAME_OVERLAP   1

struct maf {
  int16_t *buf;
  size_t N;
  size_t p;
  int32_t sum;
};

struct osc {
  int freqhz;
  int p;
};

int16_t *sinebuf;
size_t sinelen;

int16_t* make_buffer(size_t N)
{
  return (int16_t*)calloc(N, sizeof(int16_t));
}

bool sine_init(size_t N)
{
  sinebuf = make_buffer(N);
  if(!sinebuf) return false;
  
  sinelen = N;

  for(size_t i=0; i<N; ++i)
  {
    double x = 2.0 * M_PI * (double)i / (double)N;
    sinebuf[i] = (int16_t)(32767.0 * sin(x));
  }

  return true;
}

void sin_get_samples(int& p, int freqhz, int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    samples_out[i] = sinebuf[p];

    p += freqhz;
    while(p >= sinelen) p -= sinelen;
  }
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

void osc_get_samples(osc& o, int16_t *samples_out, size_t n_samples)
{
    // Note: getting samples increments the phase variable
    sin_get_samples(o.p, o.freqhz, samples_out, n_samples);
}

void osc_get_complex_samples(osc& o, int16_t *i_samples_out, int16_t *q_samples_out,
  size_t n_samples)
{
  int& q=o.p;                  // This is the Q phase variable (sine)
  int  i=o.p + sinelen / 4;    // The I phase variable is always a quarter wave ahead of Q (cosine)
  
  // Note: getting samples increments the phase variable
  sin_get_samples(i, o.freqhz, i_samples_out, n_samples);
  sin_get_samples(q, o.freqhz, q_samples_out, n_samples);
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

void mag_complex_samples(int16_t *samples_i, int16_t *samples_q,
  int16_t *samples_out, size_t n_samples)
{
  for(size_t i=0; i<n_samples; ++i)
  {
    // Fast vector magnitude calculation
    // http://www.embedded.com/design/real-time-and-performance/4007218/Digital-Signal-Processing-Tricks--High-speed-vector-magnitude-approximation
    int32_t x, y, max, min, mag;
    x = samples_i[i];
    y = samples_q[i];
    
    x = (x < 0) ? -x : x;
    y = (y < 0) ? -y : y;
    
    max = (x > y) ? x : y;
    min = (x < y) ? x : y;
    
    mag = ( 15 * (max + min / 2) ) / 16;

    if(mag >  32767)
    {
        printf("mag: clipped\n");
        samples_out[i] = 32767;
    }
    else
        samples_out[i] = mag;
  }
}

void output_buf(int16_t *samples_i, size_t n_samples)
{
  write(1, samples_i, n_samples * sizeof(int16_t));
}

size_t get_input_samples(int16_t *buf, size_t n) {
  return read(0, (char*)buf, n * sizeof(int16_t)) / sizeof(int16_t);
}

bool parity(unsigned int v){
    // http://graphics.stanford.edu/~seander/bithacks.html#ParityWith64Bits
    v ^= v >> 1;
    v ^= v >> 2;
    v = (v & 0x11111111U) * 0x11111111U;
    return ((v >> 28) & 1 != 0);
}

int main(int argc, char* argv[])
{
    bool d_forward = false;    // By default, input is the backward channel.
    int debug = 0, quiet = 0;
    char errchar = 0;   // No output for errors
    const char *frame_format = DEF_FRAME_FORMAT;
    bool lsb_first = true;
    
    // Process args
    for(int i=0; i<argc; ++i)
    {
        char* arg = argv[i];
        switch(arg[0])
        {
        case '\0':  // Empty arg
            fprintf(stderr, "Error: argument %d is empty\n", i);
            exit(1);
        case '-':   // Start of option
            switch(arg[1])
            {
                case 'D':   // Set demodulation (input) channel
                    switch(arg[2]) {
                        case 'f': d_forward = true; break;
                        case 'b': d_forward = false; break;
                        default:
                            fprintf(stderr, "Error: use -Df for forward or -Db for backward channel demodulation\n");
                            exit(1);
                    }
                    break;
                case 'M':   // Set modulation (output) channel
                    break;
                case 'd':   // Debugging output
                    ++debug;
                    break;
                case 'q':
                    ++quiet;
                    break;
                case 'e':
                    errchar = arg[2];
                    break;
                case 'f':
                    frame_format = &arg[2];
                    break;
                case 'E':   // Set endianism
                    switch(arg[2]) {
                        case 'l': lsb_first = true; break;
                        case 'b': lsb_first = false; break;
                        default:
                            fprintf(stderr, "Error: use -El for lsb first or -Eb for msb first\n");
                            exit(1);
                    }
                    break;
            }
            break;
        }
    }
    
    // Set up a the sine buffer
    if(!sine_init(SAMPLERATE)) {
        fprintf(stderr, "Failed to initialize sine buffer\n");
        exit(1);
    }
    
    // Set up frame constants
    int32_t frame_pattern = 0;
    int32_t frame_mask    = 0;
    int32_t parity_mask   = 0;
    bool parity_enable  = false;
    bool parity_even    = false;
    int data_offset = 0;
    int data_mask  = 0;
    int data_size  = 0;
    int frame_size = -FRAME_OVERLAP;
    for(size_t i=0; frame_format[i] != '\0'; ++i)
    {
        char c = frame_format[i];
        frame_mask    <<= 1;
        frame_pattern <<= 1;
        parity_mask   <<= 1;
        data_mask     <<= 1;
        ++data_offset;
        ++frame_size;
        
        switch(c)
        {
            case '1':
                frame_mask    |= 1;
                frame_pattern |= 1;
                break;
            case '0':
                frame_mask    |= 1;
                frame_pattern |= 0;
                break;
            case 'd':
                data_mask   |= 1;
                data_offset = 0;
                ++data_size;
                break;
            case 'p':
                parity_mask   |= 1;
                parity_enable = true;
                parity_even = false;
                break;
            case 'P':
                parity_mask   |= 1;
                parity_enable = true;
                parity_even = true;
                break;
            default:
                fprintf(stderr, "Invalid frame format specifier: %c\n", c);
        }
    }

  // Set up some working constants
  size_t bit_period_samples   = SAMPLERATE / ( d_forward ? F_BIT_RATE : B_BIT_RATE );
  
  // Local mark / space demodulation oscillators
  osc oscMarkIn, oscSpaceIn;
  oscMarkIn.p = 0;
  oscSpaceIn.p = 0;
  oscMarkIn.freqhz  = d_forward ? F_MARK_FREQ  : B_MARK_FREQ  ;
  oscSpaceIn.freqhz = d_forward ? F_SPACE_FREQ : B_SPACE_FREQ ;
  
  if(!quiet)
  {
    fprintf(stderr, "Demodulating the %s channel\n", d_forward ? "FORWARD" : "BACKWARD");
    fprintf(stderr, "Mark frequency:  %d Hz\n", oscMarkIn.freqhz);
    fprintf(stderr, "Space frequency: %d Hz\n", oscSpaceIn.freqhz);
    fprintf(stderr, "Bit period:   %d samples\n", bit_period_samples);
    fprintf(stderr, "Frame size: %d, format %s\n", frame_size, frame_format);
    fprintf(stderr, "Data size:  %d, %s first, with %s parity\n",
            data_size, lsb_first ? "lsb":"msb", parity_enable ? (parity_even ? "even" : "odd" ) : "no");
  }
  
  // Working registers
  int16_t *bufIn, *bufOut;
  int16_t *bufWork;
  int16_t *bufI, *bufQ;
  int16_t *bufMark, *bufSpace;
  int32_t out_shift = -1; // Raw serial shift-register
  int frame_hold = 50;    // How many bits to hold off for - preset to stabilize the MAFs
  
  // Quality monitoring
  int num_transitions = 0;
  int total_skew_correction = 0;
  
  maf mafMarkI, mafMarkQ, mafSpaceI, mafSpaceQ, mafOut, mafBit;
  size_t N = 1024; // Maximum samples we can take at once
  
  size_t bit_wait = bit_period_samples;  // Samples left until we read a bit
  
  // Set up the moving average filters
  if(! (
    maf_init(mafMarkI,  bit_period_samples) &&
    maf_init(mafMarkQ,  bit_period_samples) &&
    maf_init(mafSpaceI, bit_period_samples) &&
    maf_init(mafSpaceQ, bit_period_samples) &&
    maf_init(mafOut,    bit_period_samples) &&
    maf_init(mafBit,    bit_period_samples) )) {
    
    fprintf(stderr, "Failed to initialize MAFs\n");
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
      fprintf(stderr, "Failed to allocate buffers\n");
      exit(1);
  }
  
  if(!quiet)
    fprintf(stderr, "Initialized.  Processing samples.\n");

  size_t n;     // Number of samples we have this time
  int state=0;  // What was the last state
  while(n = get_input_samples(bufIn, N))
  {
    //printf("Got %d of %d samples...\n", n, N);
    // Mix and filter the Mark oscillator
    osc_get_complex_samples(oscMarkIn, bufI, bufQ, n);
    mul_samples(bufIn, bufI, bufWork, n);
    maf_process(mafMarkI, bufWork, bufI, n);
    mul_samples(bufIn, bufQ, bufWork, n);
    maf_process(mafMarkQ, bufWork, bufQ, n);
    mag_complex_samples(bufI, bufQ, bufMark, n);

    // Mix and filter the Space oscillator
    osc_get_complex_samples(oscSpaceIn, bufI, bufQ, n);
    mul_samples(bufIn, bufI, bufWork, n);
    maf_process(mafSpaceI, bufWork, bufI, n);
    mul_samples(bufIn, bufQ, bufWork, n);
    maf_process(mafSpaceQ, bufWork, bufQ, n);
    mag_complex_samples(bufI, bufQ, bufSpace, n);

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
            int adj = (bit_period_samples / 2) - bit_wait;
            if(debug > 1)
                fprintf(stderr, "Realigned by %d samples\n", adj);
            
            bit_wait = bit_period_samples / 2;
            
            total_skew_correction += (adj >= 0) ? adj : -adj;
            ++num_transitions;
        }
        
        // If the shift register is all ones or all zeros, dump the skew counter 
        if(out_shift == -1 || out_shift == 0)
        {
            total_skew_correction = 0;
            num_transitions = 0;
            if(debug > 3)
                fprintf(stderr, "Line idle (%04x)\n", out_shift);
        }
        
        if(--bit_wait <= 0)
        {
            if(debug > 2)
                fprintf(stderr, "Read bit '%d'\n", state);
            out_shift <<= 1;
            out_shift += state;
            
            if(frame_hold > 0)                                  // Frame Hold-off
            {
                --frame_hold;
                if(debug > 2)
                    fprintf(stderr, "Frame hold (%d left)\n", frame_hold);
            }
            else if((out_shift & frame_mask) == frame_pattern)    // Frame is valid
            {
                // Check the quality
                if(num_transitions == 0)
                    fprintf(stderr, "Dropping frame with no transitions!\n");
                else
                {
                    int avg_skew = total_skew_correction / num_transitions;
                    
                    // Reset counters
                    total_skew_correction = 0;
                    num_transitions = 0;
                    frame_hold = frame_size - 1;
                    
                    if(avg_skew > MAX_SKEW)
                    {
                        if(debug > 1)
                            fprintf(stderr, "Dropping frame with high skew of %d\n", avg_skew);
                    }
                    else
                    {
                        int32_t frame_data = out_shift & ((1<<frame_size) - 1);
                        if(debug > 1)
                            fprintf(stderr, "Processing frame: 0x%x, skew %d\n",
                                    frame_data, avg_skew);
                        
                        bool parity_bit = (frame_data & parity_mask) != 0;
                        uint32_t data = (out_shift & data_mask) >> data_offset;
                        bool data_parity = parity(data);
                        
                        if(debug > 1)
                            fprintf(stderr, "Data: 0x%02x Parity: %c Data parity: %c\n", (int)data,
                                parity_bit ? '1' : '0', data_parity ? '1' : '0'
                            );
                        
                        // Check parity
                        if(!parity_even) data_parity = !data_parity;
                        
                        // Parity check
                        if(!parity_enable || (data_parity == parity_bit))
                        {
                            // All OK
                            if(lsb_first)
                            {
                                // Assume we're working with no more than 8 data bits!
                                data <<= (8 - data_size);
                                
                                // Reverse bits in byte (LSB is first transmitted)
                                // http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv
                                data = (data * 0x0202020202ULL & 0x010884422010ULL) % 1023;
                            }
                            
                            data &= 0xff;
                            
                            if(debug > 1)
                                fprintf(stderr, "Got byte: 0x%02x\n", data);
                            printf("%c", (char)data);
                            fflush(stdout);
                        }
                        else
                        {
                            if(debug > 1)
                                fprintf(stderr, "Dropping frame with bad parity\n");
                            if(errchar)
                            {
                                printf("%c", errchar);
                                fflush(stdout);
                            }
                        }
                    }
                }
            }
            else if(debug > 2)
                fprintf(stderr, "Waiting for a valid frame\n");
            
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
  free(sinebuf);
  
  return 0;
}
