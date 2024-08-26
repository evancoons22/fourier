#include <math.h>
#include <ncurses.h>
#include <portaudio.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdatomic.h>

FILE *log_file;

#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 256
#define BUFFER_SIZE 8192 // Adjust this value as needed
#define PI 3.14159265358979323846
#define MAX_SOUNDS 4
#define NUM_BUFFERS 3


#if defined(__GNUC__) || defined(__clang__)
    #define atomic_store_double(ptr, val) __atomic_store_n((ptr), (val), __ATOMIC_SEQ_CST)
    #define atomic_load_double(ptr) __atomic_load_n((ptr), __ATOMIC_SEQ_CST)
    #define atomic_fetch_add_double(ptr, val) __atomic_fetch_add((ptr), (val), __ATOMIC_SEQ_CST)
    #define atomic_double_t double _Atomic
#else
    #error "Compiler does not support atomic operations on doubles"
#endif

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
    float *buffers[NUM_BUFFERS];
    atomic_int current_buffer;
    atomic_int next_buffer;
    atomic_int buffer_ready[NUM_BUFFERS];
    pthread_mutex_t mutexes[NUM_BUFFERS];
    pthread_cond_t conds[NUM_BUFFERS];
    atomic_int buffer_index;
} BufferManager;

typedef struct {
    SynthParams params;
    BufferManager buffer_manager;
    double phase;
    double lfo_phase;
    double lfo_envelope_value;
    // atomic_double_t start_time;
    double start_time;
} SynthData;

typedef struct {
    SynthData sounds[MAX_SOUNDS];
    int num_sounds;
    pthread_mutex_t mutex;
    // atomic_double_t current_time;
    double current_time;
} MultiSynthData;

// Function prototypes
void initialize_buffer_manager(BufferManager *bm);
void fill_buffer(SynthData *sound, int buffer_index, double current_time);
void *generate_thread(void *arg);
static int paCallbackSynth(const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo *timeInfo,
                           PaStreamCallbackFlags statusFlags, void *userData);
void add_synth(MultiSynthData *multi_data, double frequency);
double generate_sample(SynthData *data, double current_time);
double apply_envelope(SynthData *data, double current_time);



//intialize
void initialize_buffer_manager(BufferManager *bm) {
    for (int i = 0; i < NUM_BUFFERS; i++) {
        bm->buffers[i] = (float *)calloc(BUFFER_SIZE, sizeof(float));
        atomic_init(&bm->buffer_ready[i], 0);
        pthread_mutex_init(&bm->mutexes[i], NULL);
        pthread_cond_init(&bm->conds[i], NULL);
    }
    atomic_init(&bm->current_buffer, 0);
    atomic_init(&bm->next_buffer, 1);
    atomic_init(&bm->buffer_index, 0);
}

// Fill buffer function
void fill_buffer(SynthData *sound, int buffer_index, double current_time) {
    BufferManager *bm = &sound->buffer_manager;
    double sample_duration = 1.0 / SAMPLE_RATE;
    
    // this was previously the generate_sample function
    for (int i = 0; i < BUFFER_SIZE; i++) {
        double t = current_time + i * sample_duration;
        
        // Calculate LFO
        double lfo_freq_mod = sin(2 * M_PI * sound->params.lfo_freq_mod_rate * t) * 
                              sound->params.lfo_freq_mod_depth;
        double lfo_freq = sound->params.lfo_frequency + lfo_freq_mod;
        double lfo = sin(sound->lfo_phase);
        
        // Apply envelope to LFO
        double envelope = apply_envelope(sound, t);
        double lfo_with_envelope = lfo * sound->params.lfo_depth * envelope;
        
        // Calculate modulated frequency
        double modulated_freq = sound->params.frequency + lfo_with_envelope;
        
        // Generate waveform
        double sample = 0.0;
        switch (sound->params.waveform_type) {
            case 0: // Sine
                sample = sin(2 * M_PI * modulated_freq * t);
                break;
            case 1: // Square
                sample = (sin(2 * M_PI * modulated_freq * t) > 0) ? 1.0 : -1.0;
                break;
            case 2: // Sawtooth
                sample = 2.0 * (modulated_freq * t - floor(0.5 + modulated_freq * t));
                break;
            case 3: // Triangle
                sample = fabs(4.0 * (modulated_freq * t - floor(0.5 + modulated_freq * t))) - 1.0;
                break;
        }
        
        // Apply amplitude
        sample *= sound->params.amplitude;
        
        // Store in buffer
        bm->buffers[buffer_index][i] = (float)sample;
        
        // Update phases
        sound->phase += 2 * M_PI * modulated_freq * sample_duration;
        sound->lfo_phase += 2 * M_PI * lfo_freq * sample_duration;
        
        // Wrap phases
        if (sound->phase >= 2 * M_PI) sound->phase -= 2 * M_PI;
        if (sound->lfo_phase >= 2 * M_PI) sound->lfo_phase -= 2 * M_PI;
    }
}

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

double current_time;

void *generate_thread(void *arg) {
    MultiSynthData *data = (MultiSynthData *)arg;
    while (1) {
        for (int i = 0; i < data->num_sounds; i++) {
            SynthData *sound = &data->sounds[i];
            BufferManager *bm = &sound->buffer_manager;
            int buffer_to_fill = atomic_load(&bm->next_buffer);
            
            pthread_mutex_lock(&bm->mutexes[buffer_to_fill]);
            while (atomic_load(&bm->buffer_ready[buffer_to_fill])) {
                pthread_cond_wait(&bm->conds[buffer_to_fill], 
                                  &bm->mutexes[buffer_to_fill]);
            }
            
            // double current_time = atomic_load_double(&data->current_time);
            double current_time = data->current_time;
            fill_buffer(sound, buffer_to_fill, current_time);
            
            atomic_store(&bm->buffer_ready[buffer_to_fill], 1);
            pthread_cond_signal(&bm->conds[buffer_to_fill]);
            pthread_mutex_unlock(&bm->mutexes[buffer_to_fill]);
            
            atomic_store(&bm->next_buffer, (buffer_to_fill + 1) % NUM_BUFFERS);
        }
        // atomic_fetch_add_double(&data->current_time, (double)BUFFER_SIZE / SAMPLE_RATE);
        data->current_time += (double)BUFFER_SIZE / SAMPLE_RATE;
    }
    return NULL;
}

static int paCallbackSynth(const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo *timeInfo,
                           PaStreamCallbackFlags statusFlags, void *userData) {
    MultiSynthData *data = (MultiSynthData *)userData;
    float *out = (float *)outputBuffer;
    
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        float sample = 0.0f;
        for (int j = 0; j < data->num_sounds; j++) {
            SynthData *sound = &data->sounds[j];
            BufferManager *bm = &sound->buffer_manager;
            
            if (bm->buffer_index >= BUFFER_SIZE) {
                pthread_mutex_lock(&bm->mutexes[bm->current_buffer]);
                atomic_store(&bm->buffer_ready[bm->current_buffer], 0);
                pthread_cond_signal(&bm->conds[bm->current_buffer]);
                pthread_mutex_unlock(&bm->mutexes[bm->current_buffer]);
                
                bm->current_buffer = (bm->current_buffer + 1) % NUM_BUFFERS;
                bm->buffer_index = 0;
            }
            
            sample += bm->buffers[bm->current_buffer][bm->buffer_index++];
        }
        
        sample /= data->num_sounds;  // Normalize
        *out++ = sample;  // Left channel
        *out++ = sample;  // Right channel
    }
    
    return paContinue;
}



// Add synth function

void add_synth(MultiSynthData *multi_data, double frequency) {
    if (multi_data->num_sounds >= MAX_SOUNDS) {
        mvprintw(12, 0, "Max sounds reached. Cannot add more.");
        fprintf(log_file, "Max sounds reached. Cannot add more.\n");
        fflush(log_file);
        refresh();
        return;
    }

    pthread_mutex_lock(&multi_data->mutex);
    // new sound is a pointer to the next available slot
    SynthData *new_sound = &multi_data->sounds[multi_data->num_sounds];
    
    // Initialize SynthData fields
    new_sound->params = (SynthParams){
        frequency, // frequency
        1.0,       // amplitude
        0,         // waveform_type
        1000.0,    // filter_cutoff
        0.5,       // lfo_frequency
        0.5,       // lfo_depth
        0.1,       // lfo_freq_mod_rate
        0.5,       // lfo_freq_mod_depth
        0.1,       // lfo_envelope_attack
        0.2,       // lfo_envelope_decay
        0.7,       // lfo_envelope_sustain
        0.3        // lfo_envelope_release
    };
    new_sound->phase = 0.0;
    new_sound->lfo_phase = 0.0;
    new_sound->lfo_envelope_value = 0.0;
    
    // double current_time = atomic_load_double(&multi_data->current_time);
    // atomic_store_double(&new_sound->start_time, current_time);
    double current_time = multi_data->current_time;
    new_sound->start_time = current_time;
    // memcpy(new_sound->start_time, multi_data->current_time, sizeof(current_time));

    initialize_buffer_manager(&new_sound->buffer_manager);

    // need to fill here, or for some reason, we get more clicking noises
    fill_buffer(new_sound, 0, current_time);
    fill_buffer(new_sound, 1, current_time);
    fill_buffer(new_sound, 2, current_time);
    atomic_store(&new_sound->buffer_manager.buffer_ready[0], 1);

    multi_data->num_sounds++;
    pthread_mutex_unlock(&multi_data->mutex);

    mvprintw(11, 0, "Added sound with frequency: %.2f", frequency);
    // refresh();

    fprintf(log_file, "Added sound with frequency: %.2f\n", frequency);
    fflush(log_file);
}

//
//
// this is a bad way to do this right now
// really should be using a linked list or give each synth an id
void remove_synth(MultiSynthData *multi_data, SynthData *synth_data) {
  pthread_mutex_lock(&multi_data->mutex);
  for (int i = 0; i < multi_data->num_sounds; i++) {
    if (multi_data->sounds[i].params.frequency ==
        synth_data->params.frequency) {
      multi_data->num_sounds--;
      for (int j = i; j < multi_data->num_sounds; j++) {
        multi_data->sounds[j] = multi_data->sounds[j + 1];
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

//void fill_multi_buffer(MultiSynthData *data, double current_time) {
//  for (int i = 0; i < data->num_sounds; i++) {
//    fill_buffer(&data->sounds[i], current_time);
//  }
//}

void draw_piano2(WINDOW *win) {
    int startx = 2;
    int starty = 2;
    int key_width = 4;
    int key_height = 3;
    int white_keys = 7;

    // Clear the window
    wclear(win);

    // Draw the outline of the keyboard
    box(win, 0, 0);
    mvwprintw(win, 0, 1, " Piano Keyboard ");

    // Draw white keys
    for (int i = 0; i < white_keys; i++) {
        int x = startx + i * key_width;
        mvwprintw(win, starty + key_height, x + 1, " %c ", 'A' + i);
    }

    mvwhline(win, starty + key_height + 1, startx, ACS_HLINE, white_keys * key_width);

    // Draw black keys
    mvwaddch(win, starty, startx + key_width - 2, 'S');
    mvwaddch(win, starty, startx + 2 * key_width - 1, 'D');
    mvwaddch(win, starty, startx + 4 * key_width - 3, 'F');
    mvwaddch(win, starty, startx + 5 * key_width - 2, 'G');
    mvwaddch(win, starty, startx + 6 * key_width - 1, 'A');

    // Draw black key symbols (use custom or available characters)
    mvwaddch(win, starty + 1, startx + key_width - 2, ACS_BLOCK);
    mvwaddch(win, starty + 1, startx + 2 * key_width - 1, ACS_BLOCK);
    mvwaddch(win, starty + 1, startx + 4 * key_width - 3, ACS_BLOCK);
    mvwaddch(win, starty + 1, startx + 5 * key_width - 2, ACS_BLOCK);
    mvwaddch(win, starty + 1, startx + 6 * key_width - 1, ACS_BLOCK);

    // Refresh the window
    wrefresh(win);
}



void *ncurses_thread(void *arg) {
    MultiSynthData *multi_data = (MultiSynthData *)arg;
    initscr();
    start_color();
    cbreak();
    noecho();
    keypad(stdscr, TRUE);
    timeout(100); // Non-blocking input with a timeout

    int max_x, max_y;
    getmaxyx(stdscr, max_y, max_x);
    fprintf(log_file, "max_x: %d, max_y: %d\n", max_x, max_y);
    fflush(log_file);

    // Create windows for the keyboard and the keypress info
    WINDOW *keyboard_win = newwin(10, 33, 1, 1);
    WINDOW *keypress_win = newwin(3, max_x - 2, 12, 1);

    // Draw the piano
    draw_piano2(keyboard_win);
    
    // Draw the keypress info box
    box(keypress_win, 0, 0);
    mvwprintw(keypress_win, 0, 1, " Key Press Info ");

    // Refresh both windows initially
    wrefresh(keyboard_win);
    wrefresh(keypress_win);

    while (1) {
        int ch = getch();
        if (ch == 'q') {
            break;
        }

        // Clear only the content area of the keypress window, not the whole window
        werase(keypress_win);
        box(keypress_win, 0, 0);

        switch (ch) {
            case 'a':
                add_synth(multi_data, 164.0);
                mvwprintw(keypress_win, 1, 1, "A pressed: 164Hz");
                break;
            case 's':
                add_synth(multi_data, 196.0);
                mvwprintw(keypress_win, 1, 1, "S pressed: 196Hz");
                break;
            case 'd':
                add_synth(multi_data, 200.0);
                mvwprintw(keypress_win, 1, 1, "D pressed: 200Hz");
                break;
            case 'f':
                add_synth(multi_data, 587.0);
                mvwprintw(keypress_win, 1, 1, "F pressed: 587Hz");
                break;
            default:
                mvwprintw(keypress_win, 1, 1, "Invalid key pressed");
                break;
        }

        // Refresh both windows to ensure they are displayed correctly
        wrefresh(keypress_win);
    }

    delwin(keyboard_win);
    delwin(keypress_win);
    endwin();
    return NULL;
}

void play_sound_stream() {
  PaStream *stream;
  PaError err;

  MultiSynthData multi_data = {0};
  // atomic_init(&multi_data.current_time, 0.0);
  // set the current time to zero but without the atomic part
  multi_data.current_time = 0.0;
  pthread_mutex_init(&multi_data.mutex, NULL);

  // add_synth(&multi_data, 200.0);
  // add_synth(&multi_data, 193.0);
  add_synth(&multi_data, 150.0);

  fprintf(log_file, "Initializing program... num_sounds: %d\n",
          multi_data.num_sounds);
  fflush(log_file);

  err = Pa_Initialize();
  // initialize_multi_buffer(&multi_data);
  if (err != paNoError)
    goto error;

  fprintf(log_file, "Initializing buffer ... num_sounds: %d\n", multi_data.num_sounds);
  fflush(log_file);

  pthread_t thread;
  pthread_create(&thread, NULL, generate_thread, &multi_data);

  pthread_t ui_thread;
  pthread_create(&ui_thread, NULL, ncurses_thread, &multi_data);

  err = Pa_OpenDefaultStream(&stream,
                             0,         // no input channels
                             2,         // stereo output
                             paFloat32, // 32 bit floating point output
                             SAMPLE_RATE, FRAMES_PER_BUFFER, paCallbackSynth,
                             &multi_data);

  if (err != paNoError)
    goto error;

  err = Pa_StartStream(stream);
  // Pa_Sleep(500);
  if (err != paNoError)
    goto error;

  pthread_join(ui_thread, NULL);
  fclose(log_file);

  // printf("Playing. Ctrl + C to stop.\n");
  // getchar();

  err = Pa_StopStream(stream);
  if (err != paNoError)
    goto error;

  err = Pa_CloseStream(stream);
  if (err != paNoError)
    goto error;

  Pa_Terminate();

error:
  Pa_Terminate();
  fprintf(stderr, "An error occurred while using the portaudio stream\n");
  fprintf(stderr, "Error number: %d\n", err);
  fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
}


int main() {
  log_file = fopen("debug.log", "w");
  if (!log_file) {
    fprintf(stderr, "Failed to open log file\n");
    return 1;
  }

  play_sound_stream();


  return 0;
}
