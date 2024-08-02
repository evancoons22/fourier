#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <portaudio.h>
#include <time.h>
#include <string.h>


#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 256
#define BUFFER_SIZE 8192  // Adjust this value as needed
#define N 1024  // FFT size
#define WINDOW_SIZE N
#define OVERLAP (N / 2)  // 50% overlap
#define PI 3.14159265358979323846

float circular_buffer[N * 2];
int buffer_pos = 0;

typedef struct {
    double *buffer;
    int bufferSize;
    int position;
} PlaybackData;


static int paCallbackPlay(const void *inputBuffer, void *outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo* timeInfo,
                      PaStreamCallbackFlags statusFlags,
                      void *userData) {
    PlaybackData *data = (PlaybackData*)userData; // cast back to our Playbackdata struct (because we passed a pointer) 
    float *out = (float*)outputBuffer; // Cast outputBuffer to a float pointer
    (void) inputBuffer; // Prevent unused variable warning

    for (unsigned long i=0; i<framesPerBuffer; i++) {
        if (data->position >= data->bufferSize) {
            data->position = 0; // Loop the sound
        }
        *out++ = (float)data->buffer[data->position++]; // Left channel
        *out++ = (float)data->buffer[data->position-1]; // Right channel (duplicate for stereo)
    }

    return paContinue;
}


void play_sound(double *buffer, int bufferSize, double duration) {
    PaStream *stream;
    PaError err;
    PlaybackData data;

    data.buffer = buffer;
    data.bufferSize = bufferSize;
    data.position = 0;

    err = Pa_Initialize();
    if (err != paNoError) goto error;

    err = Pa_OpenDefaultStream(&stream,
                               0,          // no input channels
                               2,          // stereo output
                               paFloat32,  // 32 bit floating point output
                               SAMPLE_RATE,
                               FRAMES_PER_BUFFER,
                               paCallbackPlay,
                               &data);
    if (err != paNoError) goto error;

    err = Pa_StartStream(stream);
    if (err != paNoError) goto error;

    Pa_Sleep(duration * 1000); // Sleep for the duration of the sound

    err = Pa_StopStream(stream);
    if (err != paNoError) goto error;

    err = Pa_CloseStream(stream);
    if (err != paNoError) goto error;

    Pa_Terminate();
    return;

error:
    Pa_Terminate();
    fprintf(stderr, "An error occurred while using the portaudio stream\n");
    fprintf(stderr, "Error number: %d\n", err);
    fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
}

static int paCallback(const void *inputBuffer, void *outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo* timeInfo,
                      PaStreamCallbackFlags statusFlags,
                      void *userData)
{
    float *in = (float*)inputBuffer;
    (void) outputBuffer; // Prevent unused variable warning

    for (unsigned long i=0; i<framesPerBuffer; i++) {
        circular_buffer[buffer_pos] = in[i];
        circular_buffer[buffer_pos + N] = in[i];  // Duplicate for easy windowing
        buffer_pos = (buffer_pos + 1) % N;
    }

    return paContinue;
}

typedef struct {
    double real;
    double imag;
} complex;

// Utility functions for complex numbers
complex complex_add(complex a, complex b) {
    complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

complex complex_sub(complex a, complex b) {
    complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

complex complex_mul(complex a, complex b) {
    complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

complex complex_from_polar(double r, double theta) {
    complex result;
    result.real = r * cos(theta);
    result.imag = r * sin(theta);
    return result;
}

double test_func(double x) {
    // the period is 2pi/b, the frequency is b/2pi
    // return cos(10 * 2 * PI * x) + 2 * cos(20 * 2 * PI * x) + 3 * cos(30 * 2 * PI * x) + 4 * cos(40 * 2 * PI * x);  // frequency of 10
    return cos(10 * 2 * PI * x) + 2 * cos(20 * 2 * PI * x);  // frequency of 10 and 20
    // return cos(10 * 2 * PI * x);  // frequency of 10 and 20
}

void write_frequencies_csv(complex *X, int K) {
    // write to data.csv
    FILE *f = fopen("data.csv", "w");
    // write columns, index, real, imag
    fprintf(f, "index,real,imag\n");
    for (int k = 0; k < K; k++) {
        // the indices are i * range / N
        fprintf(f, "%d,%f,%f\n", k , sqrt(X[k].real*X[k].real + X[k].imag *  X[k].imag), 0.0);
        // fprintf(f, "%d,%f,%f\n", k, X[k].real, 0.0);
    }
}

unsigned int bit_reversal(unsigned int x,  int bits) {
    // want to go from 00011 to 11000, for example
    unsigned int y = 0;
    for (int i = 0; i < bits; i++) { 
        y = y << 1 | (x & 1);
        x >>= 1;
    } 
    return y;
} 

void bit_reverse_copy(complex *a, complex *b, int l) { 
    // for each number in the a array, reverse the bits
    // copy into array b
    int bits = (int)log2(l);
    for (unsigned int i = 0; i < l; i++) { 
        unsigned int revk = bit_reversal(i, bits);
        b[revk] = a[i];
    } 
} 

void print_binary(int x) { 
    // int is 16 bits? 
    for (int i = 1 << 15; i > 0; i = i >> 1) { 
        if ((i & x) != 0) { 
            printf("1");
        } else { 
            printf("0");
        } 
    } 
} 

void print_array(double *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        printf("%f ", a[i]);
    } 

}

void print_complex_array(complex *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        printf("%d: %f %f \n", i, a[i].real, a[i].imag);
    } 
}


void print_binary_array(int *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        print_binary(a[i]);
        printf(" ");
    } 
} 

int is_power_of_two(int x) {
    return (x > 0) && ((x & (x-1)) == 0);
}  

// Faster random number generator
static unsigned int seed = 1;
inline static unsigned int fast_rand() {
    seed = (214013 * seed + 2531011);
    return (seed >> 16) & 0x7FFF;
}

void generate_random_signal(complex *signal, int n, int num_frequencies) {
    for (int i = 0; i < n; i++) {
        signal[i].real = 0;
        signal[i].imag = 0;
    }

    for (int f = 0; f < num_frequencies; f++) {
        double amplitude = (double)fast_rand() / 0x7FFF * 10;
        int frequency = fast_rand() % (n/2);
        double phase = (double)fast_rand() / 0x7FFF * 2 * PI;

        for (int i = 0; i < n; i++) {
            signal[i].real += amplitude * cos(2 * PI * frequency * i / n + phase);
        }
    }
}

// in this better implementation, the input is not preserved
void cooley_tukey_iterative(complex *x, int n) {
    // bit reversal
    for (int i = 0; i < n; i++) {
        int j = bit_reversal(i, (int)log2(n));
        if (i < j) {
            complex temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    // main loop
    for (int step = 2; step <= n; step *= 2) {
        double angle = -2 * PI / step;
        complex w = complex_from_polar(1.0, angle);
        
        for (int k = 0; k < n; k += step) {
            complex w_k = {1, 0}; // Start with w^0 = 1
            
            for (int j = 0; j < step / 2; j++) {
                complex t = complex_mul(w_k, x[k + j + step/2]);
                complex u = x[k + j];
                // even and odd parts
                x[k + j] = complex_add(u, t);
                x[k + j + step/2] = complex_sub(u, t);
                w_k = complex_mul(w_k, w);
            }
        }
    }
}

void i_cooley_tukey_iterative(complex *X, int n) {
    // Bit-reversal permutation
    for (int i = 0; i < n; i++) {
        int j = 0;
        for (int k = 0; k < log2(n); k++) {
            j = (j << 1) | ((i >> k) & 1);
        }
        if (i < j) {
            complex temp = X[i];
            X[i] = X[j];
            X[j] = temp;
        }
    }

    // iterative inverse cooley-tukey
    for (int step = 2; step <= n; step *= 2) {
        double angle = 2 * PI / step;  // Note the positive angle for inverse FFT
        complex w = {cos(angle), sin(angle)};
        
        for (int k = 0; k < n; k += step) {
            complex w_k = {1, 0}; // Start with w^0 = 1
            
            for (int j = 0; j < step / 2; j++) {
                complex t = complex_mul(w_k, X[k + j + step/2]);
                complex u = X[k + j];
                X[k + j] = complex_add(u, t);
                X[k + j + step/2] = complex_sub(u, t);
                w_k = complex_mul(w_k, w);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        X[i].real /= n;
        X[i].imag /= n;
    }
}


void print_to_terminal(complex *x, int n) {
    static double magnitudes[N/2];
    static int d[N/2];
    static char buffer[N/2][N+3];  // Pre-allocate buffer for each line

    // Calculate magnitudes and scale to 0-50 range
    double max_magnitude = 0;
    for (int i = 0; i < n/2; i++) {
        magnitudes[i] = sqrt(x[i].real * x[i].real + x[i].imag * x[i].imag);
        if (magnitudes[i] > max_magnitude) max_magnitude = magnitudes[i];
    }

    for (int i = 0; i < n/2; i++) {
        d[i] = (int)(magnitudes[i] / max_magnitude * 50);
    }

    // Prepare output buffer
    for (int i = 0; i < n/2; i++) {
        int pos = sprintf(buffer[i], "%2d ", i);
        for (int j = 0; j < d[i]; j++) {
            buffer[i][pos++] = '*';
        }
        buffer[i][pos] = '\n';
        buffer[i][pos+1] = '\0';
    }

    // Clear terminal and print
    printf("\033[2J\033[1;1H");
    for (int i = 0; i < n/2; i++) {
        printf("%s", buffer[i]);
    }
}

// Generate a single random sample
double generate_random_sample() {
    double sample = 0;
    int num_frequencies = 3;
    for (int f = 0; f < num_frequencies; f++) {
        double amplitude = (double)fast_rand() / 0x7FFF * 10;
        int frequency = fast_rand() % (N/2);
        double phase = (double)fast_rand() / 0x7FFF * 2 * PI;
        sample += amplitude * sin(2 * PI * frequency * buffer_pos / N + phase);
    }
    return sample;
}

// Add a new sample to the circular buffer
void add_sample_to_buffer(double sample) {
    circular_buffer[buffer_pos] = sample;
    circular_buffer[buffer_pos + N] = sample;  // Duplicate for easy windowing
    buffer_pos = (buffer_pos + 1) % N;
}

// Get a window of samples from the circular buffer
void get_window(complex *window) {
    int start = (buffer_pos - WINDOW_SIZE + N) % N;
    for (int i = 0; i < WINDOW_SIZE; i++) {
        window[i].real = circular_buffer[start + i] * (0.54 - 0.46 * cos(2 * PI * i / (WINDOW_SIZE - 1)));  // Hamming window
        window[i].imag = 0;
    }
}

// Generate sine wave
void generate_sine(double *buffer, int n) {
    for (int i = 0; i < n; i++) {
        buffer[i] = sin(2 * PI * i * 2 / n);
    }
}

// Generate square wave
void generate_square(double *buffer, int n) {
    for (int i = 0; i < n; i++) {
        buffer[i] = (i < n / 2) ? 1.0 : -1.0;
    }
}

// Generate sawtooth wave
void generate_sawtooth(double *buffer, int n) {
    for (int i = 0; i < n; i++) {
        buffer[i] = 2.0 * (i / (double)n) - 1.0;
    }
}

// Generate triangle wave
void generate_triangle(double *buffer, int n) {
    for (int i = 0; i < n; i++) {
        if (i < n / 2) {
            buffer[i] = 4.0 * i / n - 1.0;
        } else {
            buffer[i] = 3.0 - 4.0 * i / n;
        }
    }
}

void low_pass_filter(complex *X, int n, double cutoff_frequency, double sample_rate) {
    double cutoff = cutoff_frequency / sample_rate;
    
    for (int i = 0; i < n; i++) {
        double frequency = (double)i / n;
        
        if (i > n/2) {
            frequency = (double)(n-i) / n;  // Mirror for negative frequencies
        }
        
        if (frequency > cutoff) {
            // Attenuate frequencies above the cutoff
            double attenuation = 1.0 - ((frequency - cutoff) / (0.5 - cutoff));
            if (attenuation < 0) attenuation = 0;
            
            X[i].real *= attenuation;
            X[i].imag *= attenuation;
        }
    }
}



typedef struct {
    double frequency;
    double amplitude;
    int waveform_type;
    double filter_cutoff;
    
    // LFO parameters
    double lfo_frequency;
    double lfo_depth;
    
    // LFO modulation parameters
    double lfo_freq_mod_rate;
    double lfo_freq_mod_depth;
    
    // LFO envelope parameters
    double lfo_envelope_attack;
    double lfo_envelope_decay;
    double lfo_envelope_sustain;
    double lfo_envelope_release;
} SynthParams;


typedef struct {
    SynthParams params;
    double phase;
    
    // New fields for LFO and envelope
    double lfo_phase;
    double lfo_envelope_value;
    // double note_on_time;
    double start_time;
    float buffer[BUFFER_SIZE];
    int buffer_index;
    int samples_generated;
} SynthData;


double apply_envelope(SynthData *data, double current_time) {
    double elapsed = current_time - data->start_time;
    double attack = data->params.lfo_envelope_attack;
    double decay = data->params.lfo_envelope_decay;
    double sustain = data->params.lfo_envelope_sustain;
    // double release = data->params.lfo_envelope_release;
    //
    if (elapsed < 0) {
        elapsed = 0; // Ensure elapsed time is not negative
    }

    // printf("envelope data. current_time =  %f: elapsed = %f, attack = %f, decay = %f, sustain = %f\n", current_time, elapsed, attack, decay, sustain);
    // printf("note on time = %f\n", data->note_on_time);
    if (elapsed < attack) {
        return elapsed / attack;
    } else if (elapsed < attack + decay) {
        return 1.0 - ((elapsed - attack) / decay) * (1.0 - sustain);
    } else {
        return sustain;
    }
    // Note: Release is not implemented here as it requires note-off information
}

// clock_t start_time;

double generate_sample(SynthData *data, double current_time) {
    double sample = 0.0;
    // double t = data->phase / SAMPLE_RATE;
    double t = current_time - data->start_time;
    // double t = Pa_GetStreamTime(NULL) - data->note_on_time;  // Calculate relative time, returns unix time? This doesn't work
    // double t = ((double)(clock() - start_time)) / CLOCKS_PER_SEC;  // Calculate relative time, THIS WORKS!

    // Apply LFO with separate phase, modulation, and envelope
    double lfo_freq_mod = sin(2 * M_PI * data->params.lfo_freq_mod_rate * t) * data->params.lfo_freq_mod_depth;
    double lfo_freq = data->params.lfo_frequency + lfo_freq_mod;
    double lfo = sin(data->lfo_phase);
    // Apply envelope to LFO
    double envelope = apply_envelope(data, current_time);
    double lfo_with_envelope = lfo * data->params.lfo_depth * envelope;
    
    // Calculate modulated frequency
   double modulated_freq = data->params.frequency + lfo_with_envelope;
    // printf("\tmodulated freq: %f = freq (%f) + lfo (%f)\n", modulated_freq, data->params.frequency, lfo_with_envelope);
    
    // Generate waveform using modulated frequency
   switch(data->params.waveform_type) {
       case 0:  // Sine
           sample = sin(2 * M_PI * modulated_freq * t);
           break;
       case 1:  // Square
           sample = (sin(2 * M_PI * modulated_freq * t) > 0) ? 1.0 : -1.0;
           break;
       case 2:  // Sawtooth
           sample = 2.0 * (modulated_freq * t - floor(0.5 + modulated_freq * t));
           break;
       case 3:  // Triangle
           sample = fabs(4.0 * (modulated_freq * t - floor(0.5 + modulated_freq * t))) - 1.0;
           break;
   }
    
    // Update phases
    data->phase += 2 * M_PI * modulated_freq / SAMPLE_RATE;
    data->lfo_phase += 2 * M_PI * lfo_freq / SAMPLE_RATE;
    
    // Wrap phases
    if (data->phase >= 2 * M_PI) data->phase -= 2 * M_PI;
    if (data->lfo_phase >= 2 * M_PI) data->lfo_phase -= 2 * M_PI;


    // Apply filter
    // return sin(2 * M_PI * data->params.frequency * t) * data->params.amplitude;

    
      return sample * data->params.amplitude;
}

void fill_buffer(SynthData *data, double current_time) {
    int samples_to_generate = BUFFER_SIZE - data->samples_generated;
    
    for (int i = 0; i < samples_to_generate; i++) {
        double t = current_time + (double)(data->samples_generated + i) / SAMPLE_RATE;
        double sample = generate_sample(data, t);
        data->buffer[(data->buffer_index + data->samples_generated + i) % BUFFER_SIZE] = (float)sample;
    }
    
    data->samples_generated += samples_to_generate;
}

static int paCallbackSynth(const void *inputBuffer, void *outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo* timeInfo,
                      PaStreamCallbackFlags statusFlags,
                      void *userData) {
    SynthData *data = (SynthData*)userData;
    float *out = (float*)outputBuffer;
    (void) inputBuffer; // Prevent unused variable warning

    //memset(out, 0, framesPerBuffer * 2 * sizeof(float));

    double current_time = timeInfo->currentTime;

    if (data->samples_generated < framesPerBuffer) {
        fill_buffer(data, current_time);
    }

    // Copy samples from the buffer to the output
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        float sample = data->buffer[data->buffer_index];
        *out++ = sample;  // Left channel
        *out++ = sample;  // Right channel

        data->buffer_index = (data->buffer_index + 1) % BUFFER_SIZE;
        data->samples_generated--;
    }

    // Refill the buffer if it's running low
    if (data->samples_generated < BUFFER_SIZE / 2) {
        fill_buffer(data, current_time + (double)framesPerBuffer / SAMPLE_RATE);
    }

    return paContinue;

   // static int call_count = 0;
   // if (++call_count % 1000 == 0) {
   //     printf("Time: %f, Frequency: %f, Amplitude: %f\n", 
   //             current_time, data->params.frequency, data->params.amplitude);
   // }

   // for (unsigned long i=0; i<framesPerBuffer; i++) {
   //     double sample = generate_sample(data, current_time + (double)i / SAMPLE_RATE);  
   //     *out++ = (float)sample;  // Left channel
   //     *out++ = (float)sample;  // Right channel
 
   // }

   // return paContinue;
}


void play_sound_stream() {
    PaStream *stream;
    PaError err;


    SynthData data = {
        {
            140.0,  // frequency
            2.0,    // amplitude
            3,      // waveform_type
            1000.0, // filter_cutoff
            2.0,    // lfo_frequency
            0.0,   // lfo_depth
            0.1,    // lfo_freq_mod_rate
            0.5,    // lfo_freq_mod_depth
            0.1,    // lfo_envelope_attack
            0.2,    // lfo_envelope_decay
            0.7,    // lfo_envelope_sustain
            0.3     // lfo_envelope_release
        },
        0.0,  // phase
        0.0,  // lfo_phase
        0.0,  // lfo_envelope_value
        0.0   // start_time
    };


    err = Pa_Initialize();
    if (err != paNoError) goto error;

    err = Pa_OpenDefaultStream(&stream,
                               0,          // no input channels
                               2,          // stereo output
                               paFloat32,  // 32 bit floating point output
                               SAMPLE_RATE,
                               FRAMES_PER_BUFFER,
                               paCallbackSynth,
                               &data);

    if (err != paNoError) goto error;

    // data.note_on_time = Pa_GetStreamTime(stream);
    //start_time = clock();
    data.start_time = Pa_GetStreamTime(stream);


    err = Pa_StartStream(stream);
    if (err != paNoError) goto error;

    printf("Playing. Ctrl + C to stop.\n");
    getchar();

    err = Pa_StopStream(stream);
    if (err != paNoError) goto error;

    err = Pa_CloseStream(stream);
    if (err != paNoError) goto error;

    Pa_Terminate();

error:
    Pa_Terminate();
    fprintf(stderr, "An error occurred while using the portaudio stream\n");
    fprintf(stderr, "Error number: %d\n", err);
    fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
}




int main() { 
    // --------------------------- play stream, no fft filters yet ------------------------------
    play_sound_stream();
    return 0;
    
    // --------------------------- Signal generation and play ------------------------------
    double time_domain[N];
    complex freq_domain[N];

    // Generate a sawtooth wave
    generate_sine(time_domain, N);

    // Convert to complex numbers for FFT
    for (int i = 0; i < N; i++) {
        freq_domain[i].real = time_domain[i];
        freq_domain[i].imag = 0;
    }

    // Apply FFT
    cooley_tukey_iterative(freq_domain, N);

    // Here you could modify freq_domain for various effects
    // print_to_terminal(freq_domain, N);
    low_pass_filter(freq_domain, N, 0.5, SAMPLE_RATE);

    // Apply inverse FFT
    i_cooley_tukey_iterative(freq_domain, N);

    // Convert back to real numbers
    for (int i = 0; i < N; i++) {
        time_domain[i] = freq_domain[i].real;
    }

    // play_sound(time_domain, N, (double)N / SAMPLE_RATE);
    play_sound(time_domain, N, 3);
    
    
    return 1;

    //--------------------------- Port audio io ------------------------------------
    // listens to audio input
     PaStream *stream;
    PaError err;

    // portaudio boiler plate
    err = Pa_Initialize();
    if (err != paNoError) goto error;
    err = Pa_OpenDefaultStream(&stream,
                               1,           // input channel count
                               0,           // output channel count
                               paFloat32,   // sample format
                               SAMPLE_RATE,
                               FRAMES_PER_BUFFER,
                               paCallback,
                               NULL);
    if (err != paNoError) goto error;
    err = Pa_StartStream(stream);
    if (err != paNoError) goto error;
    complex window[WINDOW_SIZE];

    // translate signals from audio input using port audio
    while(1) {
        // Get window of samples
        get_window(window);

        // Perform FFT
        cooley_tukey_iterative(window, WINDOW_SIZE);

        // Print results
        // print_to_terminal(window, WINDOW_SIZE);

        usleep(100000);  // 0.1 second delay
    }

    err = Pa_StopStream(stream);
    if (err != paNoError) goto error;

    err = Pa_CloseStream(stream);
    if (err != paNoError) goto error;

    Pa_Terminate();
    return 0;

    error:
        Pa_Terminate();
        fprintf(stderr, "An error occurred while using the portaudio stream\n");
        fprintf(stderr, "Error number: %d\n", err);
        fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
        return -1;


    //--------------------------- Random signals ------------------------------------
    // translate random signals, version 1
    complex signal[N];
    while(0) {
        generate_random_signal(signal, N, 3);
        cooley_tukey_iterative(signal, N);
        print_to_terminal(signal, N);
        usleep(25000);
    }
    return 0;

}
