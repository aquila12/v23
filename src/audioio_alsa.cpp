#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>

snd_pcm_t *pcm = 0;

bool audioio_alsa_init(const char* device, int rate, int audio_latency, char mode)
{
    snd_pcm_hw_params_t *hw_params;
    snd_pcm_stream_t stream;
    snd_pcm_format_t format;
    int i;
    
    if(pcm)
    {
        fprintf(stderr, "Audio device is already initialized\n");
        return false;
    }
    
    // FIXME: This might not be done at compile time?
    char test_bytes[] = {1,0};
    int16_t *test_i = (int16_t*)test_bytes;
    format = (*test_i == 1) ? SND_PCM_FORMAT_S16_LE : SND_PCM_FORMAT_S16_BE;
    
    switch(mode) {
        case 'r':
            stream = SND_PCM_STREAM_CAPTURE;
            break;
        case 'w':
            stream = SND_PCM_STREAM_PLAYBACK;
            break;
        default:
            fprintf(stderr, "Invalid mode specified (%c)\n", mode);
            return false;
    }
    
    i=snd_pcm_open(&pcm, device, stream, 0);
    if(i<0)
    {
        fprintf(stderr, "Failed to open pcm device for mode '%c': %s (ALSA error: %s)\n",
                mode, device, snd_strerror(i));
        return false;
    }
    
    i=snd_pcm_set_params(pcm, format, SND_PCM_ACCESS_RW_INTERLEAVED, 1, rate, 0, 1000 * audio_latency);
    if(i<0)
    {
        fprintf(stderr, "Failed to set the hardware params (ALSA error: %s)\n",
                snd_strerror (i));
        return false;
    }
    
    i=snd_pcm_prepare(pcm);
    if(i<0)
    {
        fprintf(stderr, "Failed to prepare the pcm device (ALSA error: %s)\n",
                snd_strerror (i));
        return false;
    }
    
    return true;
}

size_t audioio_alsa_getsamples(int16_t *buf, size_t n)
{
    snd_pcm_sframes_t frames_read=0;
    int retry=0;
    
    while(frames_read <= 0 && retry < 2) {
        ++retry;
        frames_read = snd_pcm_readi(pcm, buf, n);
        if (frames_read < 0) {
            if(frames_read == -EPIPE)
            {
                fprintf(stderr, "Buffer overrun\n");
                if(snd_pcm_prepare(pcm) < 0)
                {
                    fprintf(stderr, "  Failed to re-prepare\n");
                    return 0;
                }
                else
                {
                    continue;
                }
            }
            fprintf(stderr, "Failed to read audio frames (%s)\n",
                    snd_strerror (frames_read));
            return 0;
        }
    }
    
    return frames_read;
}

size_t audioio_alsa_putsamples(int16_t *buf, size_t n)
{
    snd_pcm_sframes_t frames_written=0;
    int retry=0;
    
    while(frames_written <= 0 && retry < 2) {
        ++retry;
        frames_written = snd_pcm_writei(pcm, buf, n);
        if (frames_written < 0) {
            if(frames_written == -EPIPE)
            {
                fprintf(stderr, "Buffer underrun\n");
                if(snd_pcm_prepare(pcm) < 0)
                {
                    fprintf(stderr, "  Failed to re-prepare\n");
                    return 0;
                }
                else
                {
                    continue;
                }
            }
            fprintf(stderr, "Failed to write audio frames (%s)\n",
                    snd_strerror (frames_written));
            return 0;
        }
    }
    
    return frames_written;
}

void audioio_alsa_stop()
{
    if(!pcm) return;
    
    snd_pcm_close(pcm);
}

