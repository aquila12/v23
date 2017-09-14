function b=fskdecode2(signal, baudrate, samplefreq, f_mark, f_space)

samples_per_bit = samplefreq / baudrate;
n_samples = length(signal);
maf = ones(1,ceil(samples_per_bit));
maf = maf / length(maf);

t = (0:(n_samples - 1)) / samplefreq;
p0 = 2*pi*f_space*t;
p1 = 2*pi*f_mark*t;

iq0 = cos(p0) + i * sin(p0);
iq1 = cos(p1) + i * sin(p1);

m0 = abs(conv(maf, iq0 .* signal));
m1 = abs(conv(maf, iq1 .* signal));

output = conv(maf, m1 - m0);

binsignal = (output > 0);

# Need to downsample the signal
last = binsignal(1,1);
output = [];
triggered = 0;
active = 0;
integrator = 0;
i_samples = 0;
bits = [];
toread = 0;
for s = binsignal
  if triggered
    # Integrate
    integrator += sign(s - 0.5);
    if i_samples <= 0
      # Are we a mark or a space?
      vote = integrator > 0;
        
      if !active
        # Either become active or go back to idle
        if(vote == 0)
          active = 1;
          toread = 10;
        else
          triggered = 0;
        endif
      endif

      if active
        # Append to output
        output = [output, vote];
        toread -= 1;
        i_samples += samples_per_bit;
        integrator = 0;

        if(toread==0)
          # Back to idle
          active = 0;
          triggered = 0;
        endif
      endif
    endif
  else
    # Wait for a negative-going edge
    if ((s==0) && (last==1))
      triggered = 1;
      i_samples += samples_per_bit;
      integrator = 0;
    else
      i_samples = 1;
    endif
  endif
  i_samples -= 1;
  last = s;
endfor

b=output;
