#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>

#define DEF_SAMPLE_RATE 44100
#define F_MARK_FREQ     1300
#define F_SPACE_FREQ    2100
#define B_MARK_FREQ     390
#define B_SPACE_FREQ    450

#define F_BIT_RATE      1200
#define B_BIT_RATE      75

#define SKEW_LIMIT          0.4
#define SKEW_CORRECT_FACTOR 3

#define DEF_FRAME_FORMAT    "10dddddddp1"

int quiet=0;
int debug=0;
int monit=0;

struct maf {
  int16_t *buf;
  size_t N;
  size_t p;
  int32_t sum;
};

struct differentiator {
  int16_t last;
};

struct osc {
  int freqhz;
  int p;
};

struct framefmt {
    int frame_size;         // Overall size of a frame
    int32_t frame_pattern;  // Pattern to look for, including previous idle / stop bit
    int32_t frame_mask;     // Mask to apply before checking for pattern
    int32_t parity_mask;    // Mask to apply to find parity bit
    bool parity_enable;     // Check parity at all ?
    bool parity_even;       // Set to use even rather than odd parity
    int data_offset;        // Number of bits after data
    int data_mask;          // Mark to apply to get data bits
    int data_size;          // Total number of data bits
    bool lsb_first;         // Does lsb come first or last (endianism)
};

struct modemcfg {
    int sample_rate;
    int mark_freqhz;
    int space_freqhz;
    framefmt ff;
    int samples_per_bit;
    int max_skew;
    char errchar;
};

int16_t *sinebuf;
size_t sinelen;

int16_t* make_buffer(size_t N)
{
  return (int16_t*)calloc(N, sizeof(int16_t));
}

bool sin_init(float amplitude, size_t N)
{
  sinebuf = make_buffer(N);
  if(!sinebuf) return false;
  
  sinelen = N;

  for(size_t i=0; i<N; ++i)
  {
    double x = 2.0 * M_PI * (double)i / (double)N;
    sinebuf[i] = (int16_t)(amplitude * sin(x));
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
        fprintf(stderr, "mul: clipped\n");
        product =  32767;
    }
    else if(product < -32767)
    {
        fprintf(stderr, "mul: clipped\n");
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

void deriv_samples(differentiator& d, int16_t *samples_in, int16_t *samples_out, size_t n_samples)
{
    int16_t last = d.last;
    for(size_t i=0; i<n_samples; ++i)
    {
        samples_out[i] = samples_in[i] - last;
        last = samples_in[i];
    }
    d.last = last;
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
        fprintf(stderr, "mag: clipped\n");
        samples_out[i] = 32767;
    }
    else
        samples_out[i] = mag;
  }
}

void ang_complex_samples(int16_t *samples_i, int16_t *samples_q,
  int16_t *samples_out, size_t n_samples)
{
    for(size_t i=0; i<n_samples; ++i)
    {        
        // Hacky inverse tangent approximation
        // Taken from my Wheeliebot guidance code
        int32_t x, y, abs_x, abs_y, angle;
        x = samples_i[i];
        y = samples_q[i];
        
        if(x==0 && y==0)
        {
            samples_out[i] = 0;
            continue;
        }
        
        // NOTE: Output phase units are 1/65536 revolution,
        // i.e. output variable overflows once per cycle
        abs_x = (x < 0) ? -x : x;
        abs_y = (y < 0) ? -y : y;
        
        if(abs_x > abs_y)
        {
            angle = (8192 * y) / x;
            if(x < 0) angle += 32768;
        }
        else
        {
            angle = 16384 - (8192 * x) / y;
            if(y < 0) angle += 32768;
        }
            
        samples_out[i] = angle;
    }
}

void output_buf(int16_t *samples_i, size_t n_samples)
{
  write(1, samples_i, n_samples * sizeof(int16_t));
}

void output_multi(int16_t *buffers[], size_t n_bufs, size_t n_samples)
{
    int16_t outbuf[n_bufs];
    
    for(size_t i=0; i<n_samples; ++i)
    {
        for(size_t j=0; j<n_bufs; ++j)
        {
            int16_t *b = buffers[j];
            outbuf[j]=b[i];
        }
        write(1, outbuf, sizeof(int16_t) * n_bufs);
    }
}

size_t get_input_samples(int16_t *buf, size_t n) {
  return read(0, (char*)buf, n * sizeof(int16_t)) / sizeof(int16_t);
}

// Parity of some data - returns true if an odd number of bits are set
bool parity(unsigned int v){
    // http://graphics.stanford.edu/~seander/bithacks.html#ParityWith64Bits
    v ^= v >> 1;
    v ^= v >> 2;
    v = (v & 0x11111111U) * 0x11111111U;
    return ((v >> 28) & 1 != 0);
}

// Note: Overlap of 1 allows checking for previous stop / idle bit
bool init_framefmt(framefmt& ff, const char* fmt, int overlap)
{
    // Set up frame constants
    ff.frame_pattern = 0;
    ff.frame_mask    = 0;
    ff.parity_mask   = 0;
    ff.parity_enable = false;
    ff.parity_even   = false;
    ff.data_offset = 0;
    ff.data_mask  = 0;
    ff.data_size  = 0;
    ff.frame_size = -overlap;
    ff.lsb_first  = true;
    for(size_t i=0; fmt[i] != '\0'; ++i)
    {
        char c = fmt[i];
        ff.frame_mask    <<= 1;
        ff.frame_pattern <<= 1;
        ff.parity_mask   <<= 1;
        ff.data_mask     <<= 1;
        ++ff.data_offset;
        ++ff.frame_size;
        
        switch(c)
        {
            case '1':
                ff.frame_mask    |= 1;
                ff.frame_pattern |= 1;
                break;
            case '0':
                ff.frame_mask    |= 1;
                ff.frame_pattern |= 0;
                break;
            case 'd':
            case 'D':
                ff.data_mask   |= 1;
                ff.data_offset = 0;
                ++ff.data_size;
                ff.lsb_first = (c == 'd');
                break;
            case 'p':
            case 'P':
                ff.parity_mask   |= 1;
                ff.parity_enable = true;
                ff.parity_even = (c == 'P');
                break;
            default:
                fprintf(stderr, "Invalid frame format specifier in %s: %c\n", fmt, c);
                return false;
        }
    }
    
    return true;
}

// Note: This is a nasty function
uint64_t bin_as_octal(uint32_t w)
{
    uint64_t d=0;
    
    for(int i=0; i<32; ++i)
    {
        d <<= 3;
        if(w & 0x80000000) ++d;
        w <<= 1;
    }
    
    return d;
}

void init_modemcfg(modemcfg& m, int mark, int space, int samplerate, int baudrate, float skew_limit) {
    m.sample_rate     = samplerate;
    m.mark_freqhz     = mark;
    m.space_freqhz    = space;
    m.samples_per_bit = samplerate / baudrate;
    m.max_skew        = (float)samplerate * skew_limit / (float)baudrate;
    m.errchar         = 0;
}

void v23_demodulate(modemcfg& m) {
    FILE* out = (monit > 0) ? stderr : stdout;    // Output chars to stderr if we're monitoring
    framefmt& f = m.ff;
    // Local oscillator
    osc o;
    o.p = 0;
    differentiator diffAng;
    diffAng.last = 0;
  
    // Working registers
    int16_t *bufIn;
    int16_t *bufI, *bufQ;
    int16_t *bufAng;
    int16_t *bufWork, *bufOut;
    
    int32_t out_shift = -1; // Raw serial shift-register
    int frame_hold = 50;    // How many bits to hold off for - preset to stabilize the MAFs
    
    // Quality monitoring
    int num_transitions = 0;
    int total_skew = 0;
    
    maf mafI, mafQ, mafAng, mafOut, mafBit;
    size_t N = 1024; // Maximum samples we can take at once
    
    int bit_wait = m.samples_per_bit;  // Samples left until we read a bit
    
    // Set up the moving average filters
    // NOTE: delta_freqhz +ve when MARK is higher frequency
    int delta_freqhz  = m.mark_freqhz - m.space_freqhz;
    o.freqhz = (m.mark_freqhz + m.space_freqhz) / 2;
    
    delta_freqhz = (delta_freqhz < 0) ? -delta_freqhz : delta_freqhz;
    int input_maf_samples = m.sample_rate / delta_freqhz;
    if(debug > 0)
    {
        fprintf(stderr, "IQ delta freq:  %d Hz\n", delta_freqhz);
        fprintf(stderr, "LO centre freq: %d Hz\n", o.freqhz);
        fprintf(stderr, "IQ MAF:         %d samples\n", input_maf_samples);
    }
    
    if(! (
        maf_init(mafI,  input_maf_samples)     &&
        maf_init(mafQ,  input_maf_samples)     &&
        maf_init(mafOut,    m.samples_per_bit) &&
        maf_init(mafBit,    m.samples_per_bit) )) {
        
        fprintf(stderr, "Failed to initialize MAFs\n");
        exit(1);
    }
    
    bufIn    = make_buffer(N);
    bufI     = make_buffer(N);
    bufQ     = make_buffer(N);
    bufAng   = make_buffer(N);
    bufWork  = make_buffer(N);
    bufOut   = make_buffer(N);
  
    if(!( bufIn && bufI && bufQ && bufAng && bufWork  && bufOut))
    {
        fprintf(stderr, "Failed to allocate buffers\n");
        exit(1);
    }
    
    if(!quiet)
        fprintf(stderr, "Initialized.  Processing samples.\n");

    size_t n;     // Number of samples we have this time
    int state=0;  // What was the last state
    bool line_idle = true;  // Are we in idle mode?
    while(n = get_input_samples(bufIn, N))
    {
        //printf("Got %d of %d samples...\n", n, N);
        // Mix and filter the local oscillator
        osc_get_complex_samples(o, bufI, bufQ, n);
        mul_samples(bufIn, bufI, bufWork, n);
        maf_process(mafI, bufWork, bufI, n);
        mul_samples(bufIn, bufQ, bufWork, n);
        maf_process(mafQ, bufWork, bufQ, n);
        
        // Determine the phase, phase change, then filter it
        ang_complex_samples(bufI, bufQ, bufAng, n);
        deriv_samples(diffAng, bufAng, bufWork, n);
        maf_process(mafOut, bufWork, bufOut, n);
        
        if(monit > 0)
        {
          int16_t *bufs[] = {bufIn, bufI, bufQ, bufAng, bufWork, bufOut};
          output_multi(bufs, 6, n);
        }
        
        // Bit sampling
        sgn_samples(bufOut, bufWork, n);
        maf_process(mafBit, bufWork, bufOut, n, true);
        
        // Run through the output samples
        int last;
        for(size_t i=0; i<n; ++i)
        {
            last = state;
            state = (bufOut[i] > 0) ? 1 : 0;
            
            // Edge detected - re-align
            if(last != state) {
                int adj = (m.samples_per_bit / 2) - bit_wait;
                if(debug > 2)
                    fprintf(stderr, "Transition, skew: %d samples\n", adj);
                
                // Don't count the first correction, and correct completely
                if(line_idle)
                    line_idle = false;
                else
                {
                    total_skew += (adj >= 0) ? adj : -adj;
                    ++num_transitions;
                    
                    // Figure out the adjustment to make
                    // ALWAYS adjust in the correct direction
                    // ALWAYS correct by at least one, unless the error is zero
                    if(adj > 0) {
                        adj /= SKEW_CORRECT_FACTOR;
                        adj += 1;
                    }
                    else if(adj < 0) {
                        adj /= SKEW_CORRECT_FACTOR;
                        adj -= 1;
                    }
                    bit_wait += adj;
                }
                if(debug > 2)
                    fprintf(stderr, "Adjusting by %d samples\n", adj);
                
                bit_wait += adj;
            }
            
            if(--bit_wait <= 0)
            {
                if(debug > 3)
                    fprintf(stderr, "Read bit '%d'\n", state);
                out_shift <<= 1;
                out_shift += state;
                
                // If the shift register is all ones or all zeros, the line is idle
                if(!line_idle && out_shift == -1 || out_shift == 0)
                {
                    line_idle = true;
                    if(debug > 1)
                        fprintf(stderr, "Line idle (%04x)\n", out_shift);
                }
                
                if(frame_hold > 0)                                  // Frame Hold-off
                {
                    --frame_hold;
                    if(debug > 2)
                        fprintf(stderr, "Frame hold (%d left)\n", frame_hold);
                }
                else if((out_shift & f.frame_mask) == f.frame_pattern)    // Frame is valid
                {                    
                    // Check the quality
                    if(num_transitions == 0)
                        fprintf(stderr, "Dropping frame '%o' with no transitions!\n", out_shift);
                    else
                    {
                        int avg_skew = total_skew / num_transitions;
                        
                        // Set line idle in case we have partial stop bits
                        line_idle = true;
                        
                        if(avg_skew > m.max_skew)
                        {
                            if(debug > 1)
                                fprintf(stderr, "Dropping frame with high skew of %d\n", avg_skew);
                        }
                        else
                        {
                            uint32_t frame_data = out_shift & ((1 << (f.frame_size+1)) - 1);
                            if(debug > 1)
                                fprintf(stderr, "Processing frame: %llo, skew %d\n",
                                        bin_as_octal(frame_data), avg_skew);
                            
                            bool parity_bit = (frame_data & f.parity_mask) != 0;
                            uint32_t data =   (frame_data & f.data_mask  ) >> f.data_offset;
                            bool data_parity = parity(data);
                            
                            if(debug > 1)
                                fprintf(stderr, "Data: 0x%02x Parity: %c Data parity: %c\n", (int)data,
                                    parity_bit ? '1' : '0', data_parity ? '1' : '0'
                                );
                            
                            // Check parity
                            if(!f.parity_even) data_parity = !data_parity;
                            
                            // Parity check
                            if(!f.parity_enable || (data_parity == parity_bit))
                            {
                                // All OK
                                if(f.lsb_first)
                                {
                                    // Assume we're working with no more than 8 data bits!
                                    data <<= (8 - f.data_size);
                                    
                                    // Reverse bits in byte (LSB is first transmitted)
                                    // http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv
                                    data = (data * 0x0202020202ULL & 0x010884422010ULL) % 1023;
                                }
                                
                                data &= 0xff;
                                
                                if(debug > 1)
                                    fprintf(stderr, "Got byte: 0x%02x\n", data);
                                
                                fprintf(out, "%c", (char)data );
                                fflush(out);
                            }
                            else
                            {
                                if(debug > 1)
                                    fprintf(stderr, "Dropping frame with bad parity\n");
                                if(m.errchar)
                                {
                                    fprintf(out, "%c", m.errchar);
                                    fflush(out);
                                }
                            }
                        }
                    }
                }
                else if(!line_idle)
                {
                    if (debug > 2)
                        fprintf(stderr, "Waiting for a valid frame\n");
                }
                
                // If the line is in idle state, reset the skew and transition count
                if(line_idle)
                {
                    total_skew = 0;
                    num_transitions = 0;
                    frame_hold = f.frame_size - 1;
                }
                
                bit_wait += m.samples_per_bit;
            }
        }
    }
    
    free(bufIn);
    free(bufI);
    free(bufQ);
    free(bufAng);
    free(bufWork);
    free(bufOut);
    free(mafI.buf);
    free(mafQ.buf);
    free(mafOut.buf);
    free(mafBit.buf);
}

void v23_modulate(modemcfg& m) {    
    // Set non-blocking input (STDIN)
    int flags = fcntl(0, F_GETFL, 0);
    fcntl(0, F_SETFL, flags | O_NONBLOCK);
    
    framefmt& f = m.ff;
    
    osc o;
    o.p = 0;
    o.freqhz = m.mark_freqhz;
    
    int32_t out_shift = -1;
    int bits_in_buffer = 0;
    int bit_wait = m.sample_rate;   // At least 1s leader tone
    
    bool line_idle = true;
    
    size_t N = m.samples_per_bit;   // Number of samples to handle at once
    
    int16_t *bufOut;
    bufOut = make_buffer(N);
    if(!bufOut)
    {
        fprintf(stderr, "Failed to allocate buffers\n");
        exit(1);
    }
    
    // Not ideal, but it's about what we want
    while(true)
    {
        // Time for another byte?
        if( bits_in_buffer < 1 )
        {
            // Is there a byte available ?
            unsigned char c_in;
            if(read(0, &c_in, 1) > 0)
            {
                // Set up the frame
                out_shift = f.frame_pattern;
                
                // Truncate data if needed
                uint32_t data = (int)c_in & (( 1 << f.data_size ) - 1);
                
                // Work out parity
                if( f.parity_enable && ( parity(data) == f.parity_even ) )
                {
                    // Set all parity bits if needed
                    out_shift |= f.parity_mask;
                }
                
                // Sort out the data order
                if(f.lsb_first)
                {
                    // Assume we're working with no more than 8 data bits!
                    data <<= (8 - f.data_size);
                    
                    // Reverse bits in byte (LSB is first transmitted)
                    // http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv
                    data = (data * 0x0202020202ULL & 0x010884422010ULL) % 1023;
                }
                
                // Manipulate the data bits into the right position
                data <<= f.data_offset;
                data &=  f.data_mask;  // Probably not necessary... just in case
                
                // Set the data bits
                out_shift |= data;
                
                bits_in_buffer = f.frame_size;
                
                if(debug > 1)
                    fprintf(stderr, "Frame for input 0x%02x: %llo\n", (int)c_in, bin_as_octal(out_shift));
                
                // One last manipulation: shift the data to the top of the word
                out_shift <<= (32 - f.frame_size);
            }
        }
        
        if( bits_in_buffer > 0 )
        {
            int i_out = 0;
            // Get next bit
            if(out_shift & 0x80000000) i_out = 1;
                
            if(debug > 2)
                fprintf(stderr, "State '%d'\n", i_out);
            
            if(i_out)
                o.freqhz = m.mark_freqhz;
            else
                o.freqhz = m.space_freqhz;
            
            out_shift <<= 1;
            --bits_in_buffer;
        }
        else
            o.freqhz = m.mark_freqhz;   // Idle
            
        // Get the oscillator samples and dump them to stdout
        osc_get_samples(o, bufOut, N);
        output_buf(bufOut, N);
    }
    
    free(bufOut);
};

int main(int argc, char* argv[])
{
    bool demodulate = true;     // By default, demodulate the backward channel.
    bool forward = false;
    char errchar = 0;           // No output for errors
    const char *frame_format = DEF_FRAME_FORMAT;
    modemcfg modem;
    int sample_rate = DEF_SAMPLE_RATE;
    float amplitude = 32767.0;  // Full-scale
    
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
                case 'A':   // Set amplitude for modulation in dB
                    {
                        float dB=0.0;
                        if(sscanf(&arg[2],"%f",&dB) < 1)
                        {
                            fprintf(stderr, "Error: -A requires a float e.g. -A3 for -3dB FS amplitude\n");
                            exit(1);
                        }
                        amplitude = 32767.0 / pow(10.0, dB / 20.0);
                        fprintf(stderr, "Set amplitude to -%f (amplitude %f)\n", dB, amplitude);
                    }
                    break;
                    
                case 'c':   // Set channel
                    switch(arg[2]) {
                        case 'f': forward = true; break;
                        case 'b': forward = false; break;
                        default:
                            fprintf(stderr, "Error: use -cf for forward or -cb for backward channel\n");
                            exit(1);
                    }
                    break;
                case 'm':   // Set mode
                    switch(arg[2]) {
                        case 'm': demodulate = false; break;
                        case 'd': demodulate = true; break;
                        default:
                            fprintf(stderr, "Error: use -mm to modulate or -md to demodulate\n");
                            exit(1);
                    }
                    break;
                case 'd':   // Debugging output
                    ++debug;
                    break;
                case 'q':   // Quiet mode
                    ++quiet;
                    break;
                case 'r':   // Sample Rate
                    sscanf(&arg[2],"%d",&sample_rate);
                    fprintf(stderr, "Set sample rate to %d\n", sample_rate);
                    break;
                case 'e':   // Character to output for parity errors
                    errchar = arg[2];
                    break;
                case 'f':   // Frame format specifier
                    frame_format = &arg[2];
                    break;
                case 'M':   // Monitor mode
                    ++monit;
                    break;
                default:
                    fprintf(stderr, "Unknown flag: %c\n", arg[1]);
                    exit(1);
            }
            break;
        }
    }
    
    // Demodulation expects the amplitude to be set to this!
    if(demodulate) amplitude = 32767.0;
    
    // Set up a the sine buffer
    if(!sin_init(amplitude, sample_rate)) {
        fprintf(stderr, "Failed to initialize sine buffer\n");
        exit(1);
    }
    
    if( !init_framefmt(modem.ff, frame_format, 1) )
    {
        fprintf(stderr, "Failed to initialize frame format\n");
        exit(1);
    }
    
    if(forward)
        init_modemcfg(modem, F_MARK_FREQ, F_SPACE_FREQ, sample_rate, F_BIT_RATE, SKEW_LIMIT);
    else
        init_modemcfg(modem, B_MARK_FREQ, B_SPACE_FREQ, sample_rate, B_BIT_RATE, SKEW_LIMIT);
    
    modem.errchar = errchar;
  
    if(!quiet)
    {
        framefmt& ff = modem.ff;
        fprintf(stderr, "%s the %s channel\n",
                demodulate ? "Demodulating" : "Modulating", forward ? "FORWARD" : "BACKWARD");
        fprintf(stderr, "Mark frequency:  %d Hz\n", modem.mark_freqhz);
        fprintf(stderr, "Space frequency: %d Hz\n", modem.space_freqhz);
        fprintf(stderr, "Bit period:      %d samples\n", modem.samples_per_bit);
        fprintf(stderr, "Max skew:        %d samples\n", modem.max_skew);
        fprintf(stderr, "Frame size:      %d, format %s\n", ff.frame_size, frame_format);
        fprintf(stderr, "Data size:       %d, %s first, with %s parity\n",
                ff.data_size, ff.lsb_first ? "lsb":"msb",
                ff.parity_enable ? (ff.parity_even ? "even" : "odd" ) : "no");
    }
    
    if(demodulate)
        v23_demodulate(modem);
    else
        v23_modulate(modem);
  
    
    free(sinebuf);
    
    return 0;
}
