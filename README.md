# v23
A v23 softmodem

# How to use
To build the software, make sure you have a suitable `gcc`, then run `make`.  The program will be built in the build directory.

The program can either modulate or demodulate a signal - not both at the same time.
If you want both, you'll need to run the program twice.

The program is set up with good channel isolation, so you can run the forward and backward channels on one audio line without,
for example, needing sidetone cancellation.

`v23` is now able to use ALSA directly, so you don't need to pipe data around yourself.  There are no other audio backends.

## Basics
By default, v23 either:
* Reads audio data from ALSA (for demodulation) and writes characters to STDOUT, and messages to STDERR
* Reads characters from STDIN (for modulation) and writes audio data to ALSA, and messages to STDERR

## Command-line arguments
Note the command-line interpretation is fairly dumb; each argument must be separate,
and any parameters must be specified without spaces in between.

### Default configuration
The default configuration is equivalent to specifying:
```shell
v23 -md -cb -r44100 -f10dddddddp1 -Ddefault -L100
```
This means v23 will demodulate the backward channel, using a 44100Hz sample rate, and no character
will be output when there is a parity error.  The frame format is equivalent to 7o1 - see below for
a description of the frame specifier.  The ALSA "default" device is used, with a 100ms latency.

Basic command-line options:
* `-m` selects whether `v23` should modulate or demodulate a signal.  Use `-mm` to modulate, and `-md` to demodulate.
* `-c` selects the channel `v23` should work on.  Use `-cf` for the forward channel, and `-cb` for the backward channel.
* `-d` increases debugging output.  Use `-d -d -d ...` for more debugging.
* `-q` increases quietness.  This disables some status messages.
* `-r` overrides the default sample rate - e.g. use `-r48000` for 48kHz sampling.
* `-f` overrides the default frame format.  See below for details.
* `-D` overrides the ALSA audio device.
* `-L` overrides the ALSA latency in ms.

The following command-line options are understood by `v23` for _modulation only_:
* `-A` specifies the amplitude of the output, in dB relative to full-scale.  Specify `-A6` for -6dB, for example.

The following command-line options are understood by `v23` for _demodulation only_:
* `-e` specifies the character to be output in the case of a parity error.  Use `-e?` to specify `?` as the error character.
* `-M` puts the program in monitor mode, for debugging the signal operations.  See below for details.

Note that you can't alter the FSK frequencies.  These are set within the code.
If you want to change them, pick the null frequencies for `init_modemcfg` carefully.

### Frame specifiers
The following characters can be used in frame format specifications:
* `1` or `0`: This bit must be in the correct state for the frame to be recognised.  Examples: start / stop bits.
* `d` or `D`: Data bits, LSb or MSb first, respectively.
* `p` or `P`: Parity bit, odd or even, respectively.

The first bit in one frame may overlap the last bit in the previous (assuming they're set the same).  This allows the
program to ensure the line was in the stop/idle state at the beginning and end of the frame, for alignment purposes.

Note that all data bits must be consecutive, or `v23` won't work as expected.

You'll only get the 8 least-significant bits of the data (in case you try to use a format such as 9n1) from stdout,
and you can only control the 8 least-significant bits at stdin.

Some examples of valid frame formats:
* `-f10dddddddp1` (7o1, the default)
* `-f10dddddddd1` (8n1, a common frame format with no parity)
* `-f10dddddddP1` (7e1, as used by viewdata)
* `-f100DDDDD11` (5 data bits, MSb first, with two start and stop bits)

### Monitor mode
If `v23` is put in monitor mode, the following things happen:
* Raw 16-bit signed audio data is written to STDOUT.  It contains one channel per signal monitored (at present, 8).
* Any demodulated data received is sent to STDERR (along with the usual messages).

## Example usage
For demodulation, use a command-line like:
```shell
build/v23 -e'?' -f10dddddddP1
```

For modulation, use a command-line like:
```shell
build/v23 -cf -mm -d -A6 -f10dddddddP1
```

You'll almost certainly want to build a phone line injector or simulator to hook this up to, but you can just
connect audio leads for testing, or to talk to another softmodem at the far end.

If you want to monitor the signals, use something like:
```shell
build/v23 -e'?' -f10dddddddP1 -M > dump.raw
```

The resulting raw audio can be imported into a tool such as Audacity.  The format is signed 16-bit little-endian, at the
same sample rate as the application was running (Default 44100Hz), with 8 channels (at present) containing the monitored
signals.
