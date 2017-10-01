#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>

snd_pcm_t *pcm = 0;

bool audioio_alsa_init(const char* device, unsigned int &rate, char mode)
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
    
    i=snd_pcm_hw_params_malloc(&hw_params);
    if(i<0)
    {
        fprintf(stderr, "Failed to malloc hardware params (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
    
    i=snd_pcm_hw_params_any(pcm, hw_params);
    if(i<0)
    {
        fprintf(stderr, "Failed to fill hardware params (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
                        
    i=snd_pcm_hw_params_set_access(pcm, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED);
    if(i<0)
    {
        fprintf(stderr, "Failed to set access (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
    
    // FIXME: This might not be done at compile time?
    char test_bytes[] = {1,0};
    int16_t *test_i = (int16_t*)test_bytes;
    format = (*test_i == 1) ? SND_PCM_FORMAT_S16_LE : SND_PCM_FORMAT_S16_BE;
    i=snd_pcm_hw_params_set_format(pcm, hw_params, format);
    if(i<0)
    {
        fprintf(stderr, "Failed to set format (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
    
    i=snd_pcm_hw_params_set_rate_near (pcm, hw_params, &rate, 0);
    if(i<0)
    {
        fprintf(stderr, "Failed to set sample rate (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
    
    snd_pcm_hw_params_set_channels(pcm, hw_params, 1);
    if(i<0)
    {
        fprintf(stderr, "Failed to set number of channels (ALSA error: %s)\n",
                snd_strerror(i));
        return false;
    }
    
    // Now apply the hardware parameters
    i=snd_pcm_hw_params(pcm, hw_params);
    if(i<0)
    {
        fprintf(stderr, "Failed to apply the hardware params (ALSA error: %s)\n",
                snd_strerror (i));
        return false;
    }
    
    snd_pcm_hw_params_free(hw_params);
    
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
    snd_pcm_sframes_t frames_read = snd_pcm_readi(pcm, buf, n);
    if (frames_read < 0) {
        fprintf(stderr, "Failed to read audio frames (%s)\n",
                snd_strerror (frames_read));
        return 0;
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

