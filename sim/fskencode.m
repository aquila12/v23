function y=fskencode(bitstream, baudrate, samplefreq, f_mark, f_space)

samples_per_bit = samplefreq / baudrate;
outlen = length(bitstream) * samples_per_bit;
sample_indices = ceil((1:(outlen - 1)) /samples_per_bit);

samples = bitstream(sample_indices);

base   = ones(1,length(samples)) * f_space;
offset = samples * (f_mark - f_space);

frequencies = base + offset;
phase = cumsum(2 * pi * frequencies * (1/samplefreq));

y = sin(phase);
