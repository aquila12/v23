#include <alsa/asoundlib.h>

bool audioio_alsa_init(const char* device, unsigned int &rate, char mode);
size_t audioio_alsa_getsamples(int16_t *buf, size_t n);
size_t audioio_alsa_putsamples(int16_t *buf, size_t n);
void audioio_alsa_stop();
