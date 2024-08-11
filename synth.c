#include <stdio.h> 
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <portaudio.h>
#include <time.h>
#include <string.h>
#include <ncurses.h>
#include <unistd.h>

FILE *log_file; 

#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 256
#define BUFFER_SIZE 8192  // Adjust this value as needed
#define PI 3.14159265358979323846
#define MAX_SOUNDS 4


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
    float buffer1[BUFFER_SIZE];
    float buffer2[BUFFER_SIZE];
    float *current_buffer;
    float *next_buffer;

    int buffer_index;
    int samples_generated;

    pthread_mutex_t buffer_mutex;
    pthread_cond_t buffer_cond;
    int buffer_ready;
} SynthData;



typedef struct { 
    SynthData sounds[MAX_SOUNDS];
    int num_sounds;
    pthread_mutex_t mutex;
} MultiSynthData;


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


void fill_buffer(SynthData *sound, double current_time) {
    for (unsigned long i = 0; i < BUFFER_SIZE; i++) {
        sound->next_buffer[i] = generate_sample(sound, current_time + i / (double)SAMPLE_RATE);
    }
    pthread_mutex_lock(&sound->buffer_mutex);
    sound->buffer_ready = 1;
    pthread_cond_signal(&sound->buffer_cond);
    pthread_mutex_unlock(&sound->buffer_mutex);
}


double current_time;

void* generate_thread(void *arg) {
    MultiSynthData *data = (MultiSynthData*)arg;
    current_time = 0.0;

    while (1) {
        // sleep for a bit
        usleep(100);
        for (int i = 0; i < data->num_sounds; i++) {
            SynthData *sound = &data->sounds[i];
            pthread_mutex_lock(&sound->buffer_mutex);
            while (sound->buffer_ready) {
                pthread_cond_wait(&sound->buffer_cond, &sound->buffer_mutex);
            }
            // fprintf(log_file, "pthread condition hit, unlocking mutex for next sound to generate\n");
            // fflush(log_file);
            pthread_mutex_unlock(&sound->buffer_mutex);

            fill_buffer(sound, current_time);
            current_time += (double)BUFFER_SIZE / SAMPLE_RATE;
        }
    }

    return NULL;
}


static int paCallbackSynth(const void *inputBuffer, void *outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo* timeInfo,
                      PaStreamCallbackFlags statusFlags,
                      void *userData) {
    MultiSynthData *data = (MultiSynthData*)userData;
    float *out = (float*)outputBuffer;
    (void) inputBuffer; // Prevent unused variable warning
                        //
                        
    memset(out, 0, framesPerBuffer * 2 * sizeof(float));


    // Copy samples from the buffer to the output
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        float sample = 0.0f;
        // cycle through each sound to generate a sample and send to the output
        for (int j = 0; j < data->num_sounds; j++) {
            SynthData *sound = &data->sounds[j];
            pthread_mutex_lock(&sound->buffer_mutex);
            if (sound->samples_generated == 0) {
                // Switch buffers if current buffer is exhausted
                float *temp = sound->current_buffer;
                sound->current_buffer = sound->next_buffer;
                sound->next_buffer = temp;
                sound->buffer_index = 0;
                sound->samples_generated = BUFFER_SIZE;
                sound->buffer_ready = 0;
                pthread_cond_signal(&sound->buffer_cond);
            }

            sample += sound->current_buffer[sound->buffer_index];
            sound->buffer_index++;
            sound->samples_generated--;
            pthread_mutex_unlock(&sound->buffer_mutex);

        }

        // Normalize the mixed sample to avoid clipping
        sample = sample / data->num_sounds;

        // Output the mixed sample to both left and right channels
        *out++ = sample;  // Left channel
        *out++ = sample;  // Right channel
    }


    return paContinue;
}



void initialize_buffer(SynthData *data) {
    mvprintw(10, 0, "Initializing buffer with frequency: %.2f", data->params.frequency); // Debug print
    refresh();
    data->buffer_index = 0;
    data->samples_generated = BUFFER_SIZE;
    data->current_buffer = data->buffer1;
    data->next_buffer = data->buffer2;
    pthread_mutex_init(&data->buffer_mutex, NULL);
    pthread_cond_init(&data->buffer_cond, NULL);
    data->buffer_ready = 0;
    fill_buffer(data, 0.0);
}

void initialize_multi_buffer(MultiSynthData *data) {
    // data->num_sounds = 0;
    for (int i = 0; i < data->num_sounds; i++) {
        initialize_buffer(&data->sounds[i]);
    }
}


void add_synth(MultiSynthData *multi_data, double frequency) {
    if (multi_data->num_sounds >= MAX_SOUNDS) {
        mvprintw(12, 0, "Max sounds reached. Cannot add more.");
        refresh();
        return;
    }

    SynthData *new_sound = (SynthData *)malloc(sizeof(SynthData));
    if (new_sound == NULL) {
        mvprintw(12, 0, "Memory allocation failed!");
        refresh();
        return;
    }

    *new_sound = (SynthData){
        // params
        {
            frequency,  // frequency
            1.0,    // amplitude
            0,      // waveform_type
            1000.0, // filter_cutoff
            0.5,    // lfo_frequency
            0.5,    // lfo_depth
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
        current_time  // start_time
    };

    initialize_buffer(new_sound);

    pthread_mutex_lock(&multi_data->mutex);
    multi_data->sounds[multi_data->num_sounds] = *new_sound;
    multi_data->num_sounds++;
    pthread_mutex_unlock(&multi_data->mutex);

    mvprintw(11, 0, "Added sound with frequency: %.2f", frequency);
    refresh();
    fprintf(log_file, "Added sound with frequency: %.2f\n", frequency);
    fflush(log_file);

    free(new_sound);  // Free the dynamically allocated memory
}

// this is a bad way to do this right now
// really should be using a linked list or give each synth an id
void remove_synth(MultiSynthData *multi_data, SynthData *synth_data) {
    pthread_mutex_lock(&multi_data->mutex);
    for (int i = 0; i < multi_data->num_sounds; i++) {
        if (multi_data->sounds[i].params.frequency == synth_data->params.frequency) {
            multi_data->num_sounds--;
            for (int j = i; j < multi_data->num_sounds; j++) {
                multi_data->sounds[j] = multi_data->sounds[j+1];
            }
            break;
        }
    }
    pthread_mutex_unlock(&multi_data->mutex);
}

void initialize_multi_syth_data(MultiSynthData *data) {
    data->num_sounds = 0;
    pthread_mutex_init(&data->mutex, NULL);
}

void fill_multi_buffer(MultiSynthData *data, double current_time) {
    for (int i = 0; i < data->num_sounds; i++) {
        fill_buffer(&data->sounds[i], current_time);
    }
}

void print_synth_info(SynthData *data) {
    fprintf(log_file, "Here is how bad your new synth is:\n");
    fprintf(log_file, "\tfrequency: %.2f\n", data->params.frequency);
    fprintf(log_file, "\tamplitude: %.2f\n", data->params.amplitude);
    fprintf(log_file, "\twaveform_type: %d\n", data->params.waveform_type);
    fprintf(log_file, "\tfilter_cutoff: %.2f\n", data->params.filter_cutoff);
    fprintf(log_file, "\tlfo_frequency: %.2f\n", data->params.lfo_frequency);
    fprintf(log_file, "\tlfo_depth: %.2f\n", data->params.lfo_depth);
    fprintf(log_file, "\tlfo_freq_mod_rate: %.2f\n", data->params.lfo_freq_mod_rate);
    fprintf(log_file, "\tlfo_freq_mod_depth: %.2f\n", data->params.lfo_freq_mod_depth);
    fprintf(log_file, "\tlfo_envelope_attack: %.2f\n", data->params.lfo_envelope_attack);
    fprintf(log_file, "\tlfo_envelope_decay: %.2f\n", data->params.lfo_envelope_decay);
    fprintf(log_file, "\tlfo_envelope_sustain: %.2f\n", data->params.lfo_envelope_sustain);
}





void* ncurses_thread(void *arg) {
    MultiSynthData *multi_data = (MultiSynthData*)arg;
    initscr();
    noecho();
    cbreak();
    timeout(100); // Non-blocking input with a timeout

    int line = 5;
    while (1) {
        int ch = getch();
        if (ch == 'q') {
            break;
        }
        // SynthData new_synth;
        switch (ch) {
            case 'a':
                // add_default_sound(multi_data, 440.0);
                add_synth(multi_data, 164.0);
                mvprintw(0, 0, "A pressed: 440Hz");
                break;
            case 's':
                // add_default_sound(multi_data, 494.0);
                add_synth(multi_data, 196.0);
                mvprintw(1, 0, "S pressed: 494Hz");
                break;
            case 'd':
                add_synth(multi_data, 200.0);
                mvprintw(2, 0, "D pressed: 523Hz");
                break;
            case 'f':
                add_synth(multi_data, 587.0);
                mvprintw(3, 0, "F pressed: 587Hz");
                break;
        }
        // mvprintw(line++, 0, "Key %c pressed, num_sounds: %d", ch, multi_data->num_sounds);
        refresh();
    }

    endwin();
    return NULL;
}

void play_sound_stream() {
    PaStream *stream;
    PaError err;


    SynthData data = {
        {
            // C
            440.81,  // frequency
            1.0,    // amplitude
            0,      // waveform_type
            1000.0, // filter_cutoff
            0.2,    // lfo_frequency
            1.0,   // lfo_depth
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


    SynthData data2 = {
        {
            // E
            164.81,  // frequency
            1.0,    // amplitude
            0,      // waveform_type
            1000.0, // filter_cutoff
            0.5,    // lfo_frequency
            0.5,   // lfo_depth
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
    

    SynthData data3 = {
        {
            // G
            196.00,  // frequency
            1.0,    // amplitude
            0,      // waveform_type
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

    MultiSynthData multi_data;
    initialize_multi_syth_data(&multi_data);
    // print number of sounds
    fprintf(log_file, "Initializing program... num_sounds: %d\n", multi_data.num_sounds);
    fflush(log_file);
    // add_synth(&multi_data, &data);
    // add_synth(&multi_data, &data2);
    // add_synth(&multi_data, 196.00);

    //multi_data.sounds[0] = data;
    //multi_data.sounds[1] = data2;
    //multi_data.sounds[2] = data3;

    // print some multidata parameters


    err = Pa_Initialize();
    initialize_multi_buffer(&multi_data);
    fprintf(log_file, "Initializing buffer ... num_sounds: %d\n", multi_data.num_sounds);
    fflush(log_file);
    if (err != paNoError) goto error;

    pthread_t thread;
    pthread_create(&thread, NULL, generate_thread, &multi_data);

    pthread_t ui_thread;
    pthread_create(&ui_thread, NULL, ncurses_thread, &multi_data);


    err = Pa_OpenDefaultStream(&stream,
                               0,          // no input channels
                               2,          // stereo output
                               paFloat32,  // 32 bit floating point output
                               SAMPLE_RATE,
                               FRAMES_PER_BUFFER,
                               paCallbackSynth,
                               &multi_data);

    if (err != paNoError) goto error;

    // data.note_on_time = Pa_GetStreamTime(stream);
    //start_time = clock();
    // data.start_time = Pa_GetStreamTime(stream);


    err = Pa_StartStream(stream);
    //Pa_Sleep(500);
    if (err != paNoError) goto error;
    
    pthread_join(ui_thread, NULL);
    fclose(log_file);

    //printf("Playing. Ctrl + C to stop.\n");
    //getchar();

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


void print_params(MultiSynthData *data) {
    // using ncurser this will be my simple terminal interface for now
    printf("Frequency: %f\n", data->sounds[0].params.lfo_frequency);
}





int main() { 
    // --------------------------- play stream, no fft filters yet ------------------------------
    
    log_file = fopen("debug.log", "w");
    if (!log_file) {
        fprintf(stderr, "Failed to open log file\n");
        return 1;
    }

    play_sound_stream();

    // figure out how to get ncurses working with portaudio callback still going
    
    return 0;
} 
