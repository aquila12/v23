function b=fskdecode3(signal, baudrate, samplefreq, f_mark, f_space)

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
integrator = conv(maf, binsignal);

# Need to downsample the signal
baud = baudrate;
t=0;
s_prev = binsignal(1,1);

output = [];
offset = 0;
next_sample = samples_per_bit;

for s = binsignal
  next_sample -= 1;
  t += 1;

  if s_prev != s
    t_errbits = (next_sample / samples_per_bit);
    offset = round(t_errbits) - t_errbits;
    offset *= samples_per_bit;
    next_sample += offset;
  end

  if next_sample <= 0
    vote = integrator(t) > 0.5;
    output = [output, vote];
    next_sample += samplefreq/baud;
  end

  s_prev = s;
endfor

b=output;
