// define termbox2 first for some macros
#define TB_IMPL
#include "termbox2.h"

#include <math.h>
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
#define MAX_SOUNDS 12
#define NUM_BUFFERS 3
#define RECENT_ACTIONS_MAX 10

int global_waveform_type = 0;  // Default to sine wave

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
    double start_time;
    double release_start_time;
    int is_active;
} SynthData;

typedef struct {
    SynthData sounds[MAX_SOUNDS];
    int num_sounds;
    pthread_mutex_t mutex;
    // atomic_double_t current_time;
    double current_time;
} MultiSynthData;

typedef struct {
    char actions[RECENT_ACTIONS_MAX][50];
    int action_count;
} RecentActions;

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
            case 4: // Pulse
                sample = (fmod(t * modulated_freq, 1.0) < 0.5) ? 1.0 : -1.0;
                break;
            case 5: // Noise
                sample = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
                break;
            case 6: // Sine Squared
                sample = pow(sin(M_PI * modulated_freq * t), 2);
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
    double release = data->params.lfo_envelope_release;

    if (elapsed < 0) {
        elapsed = 0; // Ensure elapsed time is not negative
    }

    if (data->is_active) {
        if (elapsed < attack) {
            return elapsed / attack;
        } else if (elapsed < attack + decay) {
            return 1.0 - ((elapsed - attack) / decay) * (1.0 - sustain);
        } else {
            return sustain;
        }
    } else {
        fprintf(log_file, "Note is not active\n");
        fflush(log_file);
        double release_elapsed = current_time - data->release_start_time;
        if (release_elapsed < release) {
            return sustain * (1.0 - release_elapsed / release);
        } else {
            fprintf(log_file, "Note completely off\n");
            fflush(log_file);
            return 0.0;
        }
    }
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
                pthread_cond_wait(&bm->conds[buffer_to_fill], &bm->mutexes[buffer_to_fill]);
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
        fprintf(log_file, "Max sounds reached. Cannot add more.\n");
        fflush(log_file);
        return;
    }

    pthread_mutex_lock(&multi_data->mutex);
    // new sound is a pointer to the next available slot
    SynthData *new_sound = &multi_data->sounds[multi_data->num_sounds];

    // Initialize SynthData fields
    new_sound->params = (SynthParams){
        frequency, // frequency
        1.0,       // amplitude
        global_waveform_type, // waveform_type (now using the global variable)
        1000.0,    // filter_cutoff
        1.2,       // lfo_frequency
        0.1,       // lfo_depth
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
    new_sound->is_active = 1;
    double current_time = multi_data->current_time;
    new_sound->start_time = current_time;
    new_sound->release_start_time = 0.0;

    initialize_buffer_manager(&new_sound->buffer_manager);

    // need to fill here, or for some reason, we get more clicking noises
    fill_buffer(new_sound, 0, current_time);
    fill_buffer(new_sound, 1, current_time);
    fill_buffer(new_sound, 2, current_time);
    atomic_store(&new_sound->buffer_manager.buffer_ready[0], 1);

    multi_data->num_sounds++;
    pthread_mutex_unlock(&multi_data->mutex);


    fprintf(log_file, "Added sound with frequency: %.2f\n", frequency);
    fflush(log_file);
}

void remove_synth(MultiSynthData *multi_data, double frequency) {
    pthread_mutex_lock(&multi_data->mutex);
    for (int i = 0; i < multi_data->num_sounds; i++) {
        if (multi_data->sounds[i].params.frequency == frequency) {
            // Free any resources associated with this sound
            // (e.g., buffer memory)
            for (int j = 0; j < NUM_BUFFERS; j++) {
                free(multi_data->sounds[i].buffer_manager.buffers[j]);
            }
            
            // Shift remaining sounds
            for (int j = i; j < multi_data->num_sounds - 1; j++) {
                multi_data->sounds[j] = multi_data->sounds[j + 1];
            }
            
            multi_data->num_sounds--;
            break;
        }
    }
    pthread_mutex_unlock(&multi_data->mutex);
}

void initialize_multi_syth_data(MultiSynthData *data) {
    data->num_sounds = 0;
    pthread_mutex_init(&data->mutex, NULL);
}


//void tb_print(int x, int y, uint16_t fg, uint16_t bg, const char *str) {
//    while (*str) {
//        tb_set_cell(x++, y, *str++, fg, bg);
//    }
//}

// Function to add recent action
void add_recent_action(RecentActions *recent_actions, const char *action) {
    if (recent_actions->action_count < RECENT_ACTIONS_MAX) {
        strcpy(recent_actions->actions[recent_actions->action_count], action);
        recent_actions->action_count++;
    } else {
        // Shift actions up and add new action at the end
        for (int i = 1; i < RECENT_ACTIONS_MAX; i++) {
            strcpy(recent_actions->actions[i - 1], recent_actions->actions[i]);
        }
        strcpy(recent_actions->actions[RECENT_ACTIONS_MAX - 1], action);
    }
}

void draw_box(int startx, int starty, int width, int height, const char *title) {
    // Draw top and bottom borders
    for (int x = startx; x <= startx + width; x++) {
        tb_set_cell(x, starty, '-', TB_WHITE, TB_DEFAULT);
        tb_set_cell(x, starty + height, '-', TB_WHITE, TB_DEFAULT);
    }
    // Draw left and right borders
    for (int y = starty; y <= starty + height; y++) {
        tb_set_cell(startx, y, '|', TB_WHITE, TB_DEFAULT);
        tb_set_cell(startx + width, y, '|', TB_WHITE, TB_DEFAULT);
    }
    // Draw corners
    tb_set_cell(startx, starty, '+', TB_WHITE, TB_DEFAULT);
    tb_set_cell(startx + width, starty, '+', TB_WHITE, TB_DEFAULT);
    tb_set_cell(startx, starty + height, '+', TB_WHITE, TB_DEFAULT);
    tb_set_cell(startx + width, starty + height, '+', TB_WHITE, TB_DEFAULT);
    // Draw title
    tb_print(startx + 2, starty, TB_WHITE, TB_DEFAULT, title);
}

void draw_recent_actions(int startx, int starty, RecentActions *recent_actions) {
    int width = 40;
    int height = RECENT_ACTIONS_MAX + 2;
    draw_box(startx, starty, width, height, " Recent Actions ");
    for (int i = 0; i < recent_actions->action_count; i++) {
        tb_print(startx + 2, starty + 2 + i, TB_WHITE, TB_DEFAULT, recent_actions->actions[i]);
    }
}

void draw_constants(int startx, int starty) {
    int width = 30;
    int height = 6;
    draw_box(startx, starty, width, height, " Constants ");
    char buffer[50];
    snprintf(buffer, sizeof(buffer), "SAMPLE_RATE: %d", SAMPLE_RATE);
    tb_print(startx + 2, starty + 2, TB_WHITE, TB_DEFAULT, buffer);
    snprintf(buffer, sizeof(buffer), "FRAMES_PER_BUFFER: %d", FRAMES_PER_BUFFER);
    tb_print(startx + 2, starty + 3, TB_WHITE, TB_DEFAULT, buffer);
    snprintf(buffer, sizeof(buffer), "BUFFER_SIZE: %d", BUFFER_SIZE);
    tb_print(startx + 2, starty + 4, TB_WHITE, TB_DEFAULT, buffer);

    switch (global_waveform_type) {
        case 0: // Sine
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: SINE ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 1: // Square
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: SQUARE ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 2: // Sawtooth
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: SAWTOOTH ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 3: // Triangle
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: TRIANGLE ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 4: // Pulse
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: PULSE ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 5: // Noise
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: NOISE ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
        case 6: // Sine Squared
            snprintf(buffer, sizeof(buffer), "WAVEFORM_TYPE: SINE SQUARED ");
            tb_print(startx + 2, starty + 5, TB_WHITE, TB_DEFAULT, buffer);
            
    }
}


void draw_synth_params(int startx, int starty, MultiSynthData *multi_data) {
    int width = 150;  // Adjust width based on the number of parameters
    int height = multi_data->num_sounds * 2 + 4;  // Height based on the number of sounds + header row + extra space
    draw_box(startx, starty, width, height, " Parameters ");

    // Header row
    tb_print(startx + 2, starty + 2, TB_WHITE, TB_DEFAULT, "Key ");
    tb_print(startx + 6, starty + 2, TB_WHITE, TB_DEFAULT, "Freq ");
    tb_print(startx + 16, starty + 2, TB_WHITE, TB_DEFAULT, "Amp  ");
    tb_print(startx + 26, starty + 2, TB_WHITE, TB_DEFAULT, "Waveform ");
    tb_print(startx + 36, starty + 2, TB_WHITE, TB_DEFAULT, "Cutoff ");
    tb_print(startx + 46, starty + 2, TB_WHITE, TB_DEFAULT, "LFO Freq ");
    tb_print(startx + 58, starty + 2, TB_WHITE, TB_DEFAULT, "LFO Depth ");
    tb_print(startx + 70, starty + 2, TB_WHITE, TB_DEFAULT, "LFO Mod Rate ");
    tb_print(startx + 84, starty + 2, TB_WHITE, TB_DEFAULT, "LFO Mod Depth ");
    tb_print(startx + 98, starty + 2, TB_WHITE, TB_DEFAULT, "Env Attack ");
    tb_print(startx + 110, starty + 2, TB_WHITE, TB_DEFAULT, "Env Decay ");
    tb_print(startx + 122, starty + 2, TB_WHITE, TB_DEFAULT, "Env Sustain ");
    tb_print(startx + 134, starty + 2, TB_WHITE, TB_DEFAULT, "Env Release ");

    // Print each sound's parameters
    for (int i = 0; i < multi_data->num_sounds; i++) {
        SynthParams *params = &multi_data->sounds[i].params;
        int y_offset = starty + 4 + i * 2;

        char buffer[20];
        snprintf(buffer, sizeof(buffer), "%2d", i);
        tb_print(startx + 2, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->frequency);
        tb_print(startx + 6, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->amplitude);
        tb_print(startx + 16, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%d", params->waveform_type);
        tb_print(startx + 26, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->filter_cutoff);
        tb_print(startx + 36, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_frequency);
        tb_print(startx + 46, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_depth);
        tb_print(startx + 58, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_freq_mod_rate);
        tb_print(startx + 70, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_freq_mod_depth);
        tb_print(startx + 84, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_envelope_attack);
        tb_print(startx + 98, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_envelope_decay);
        tb_print(startx + 110, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_envelope_sustain);
        tb_print(startx + 122, y_offset, TB_WHITE, TB_DEFAULT, buffer);

        snprintf(buffer, sizeof(buffer), "%.1f", params->lfo_envelope_release);
        tb_print(startx + 134, y_offset, TB_WHITE, TB_DEFAULT, buffer);
    }
}


void draw_piano_keyboard(int startx, int starty, int *active_keys) {
    int white_key_width = 4;
    int white_key_height = 7;
    int black_key_width = 3;
    int black_key_height = 4;

    const char *white_keys = "SDFGHJK";
    const char *black_keys = "ERYUI";
    int black_key_positions[] = {1, 2, 4, 5, 6};

    // Draw white keys
    for (int i = 0; i < 7; i++) {
        int x = startx + i * white_key_width;
        int y = starty;
        uint16_t bg = (active_keys[i]) ? TB_GREEN : TB_WHITE;
        // Draw key body
        for (int dx = 0; dx < white_key_width; dx++) {
            for (int dy = 0; dy < white_key_height; dy++) {
                tb_set_cell(x + dx, y + dy, ' ', TB_BLACK, bg);
            }
        }
        // Draw key label
        tb_set_cell(x + white_key_width / 2, y + white_key_height - 1, white_keys[i], TB_BLACK, bg);
    }

    // Draw black keys
    for (int i = 0; i < 5; i++) {
        int x = startx + black_key_positions[i] * white_key_width - black_key_width / 2;
        int y = starty;
        uint16_t bg = (active_keys[i + 7]) ? TB_RED : TB_BLACK;
        // Draw key body
        for (int dx = 0; dx < black_key_width; dx++) {
            for (int dy = 0; dy < black_key_height; dy++) {
                tb_set_cell(x + dx, y + dy, ' ', TB_WHITE, bg);
            }
        }
        // Draw key label
        tb_set_cell(x + black_key_width / 2, y + black_key_height - 1, black_keys[i], TB_WHITE, bg);
    }
}

void draw_interface(MultiSynthData *multi_data, RecentActions *recent_actions, int highlighted_key, int *active_keys) {
    tb_clear();

    // Define positions
    int recent_actions_x = 2;
    int recent_actions_y = 2;

    int piano_x = 45;
    int piano_y = 5;

    int constants_x = 75;
    int constants_y = 2;

    int synth_params_x = 2;
    int synth_params_y = 18; 

    // Draw components
    draw_recent_actions(recent_actions_x, recent_actions_y, recent_actions);
    draw_piano_keyboard(piano_x, piano_y, active_keys);
    draw_constants(constants_x, constants_y);
    draw_synth_params(synth_params_x, synth_params_y, multi_data);

    tb_present();
}


void *termbox_thread(void *arg) {
    MultiSynthData *multi_data = (MultiSynthData *)arg;
    RecentActions recent_actions = {.action_count = 0};
    struct tb_event ev;
    int highlighted_key = -1;
    int key_to_freq[] = {
        130, // S
        146, // D
        164, // F
        174, // G
        196, // H
        220, // J
        246, // K
        138, // E
        155, // R
        185, // Y
        207, // U
        233  // I
    };
    char key_chars[] = {'s','d','f','g','h','j','k','e','r','y','u','i'};
    int active_keys[12] = {0};

    if (tb_init() != 0) {
        fprintf(stderr, "tb_init() failed\n");
        return NULL;
    }

    tb_set_input_mode(TB_INPUT_ESC | TB_INPUT_MOUSE);
    tb_set_output_mode(TB_OUTPUT_NORMAL);

    draw_interface(multi_data, &recent_actions, highlighted_key, active_keys);

    while (1) {
        int tb_result = tb_poll_event(&ev);
        if (tb_result == -1) break;
        if (ev.type == TB_EVENT_KEY) {
            if (ev.key == TB_KEY_ESC) {
                break;
            }
            char action_text[50];
            highlighted_key = -1;
            for (int i = 0; i < 12; i++) {
                if (ev.ch == key_chars[i]) {
                    highlighted_key = i;
                    if (!active_keys[i]) {
                        add_synth(multi_data, key_to_freq[i]); // waveform_type set to 0
                        active_keys[i] = 1;
                        snprintf(action_text, sizeof(action_text), "Key %c pressed: %dHz", ev.ch, key_to_freq[i]);
                        add_recent_action(&recent_actions, action_text);
                    } else {
                        remove_synth(multi_data, key_to_freq[i]);
                        active_keys[i] = 0;
                        snprintf(action_text, sizeof(action_text), "Key %c released", ev.ch);
                        add_recent_action(&recent_actions, action_text);
                    }
                    break;
                }
            }
            draw_interface(multi_data, &recent_actions, highlighted_key, active_keys);
        }
    }

    tb_shutdown();
    return NULL;
}


int main(int argc, char *argv[]) {
    if (argc != 3 || strcmp(argv[1], "-type") != 0) {
        fprintf(stderr, "Usage: %s -type <waveform_type>\n", argv[0]);
        return 1;
    }

    int waveform_type = atoi(argv[2]);
    if (waveform_type < 0 || waveform_type > 6) {
        fprintf(stderr, "Invalid waveform type. Must be between 0 and 6.\n");
        return 1;
    }

    global_waveform_type = waveform_type;


    log_file = fopen("debug.log", "w");
    if (!log_file) {
        fprintf(stderr, "Failed to open log file\n");
        return 1;
    }

    PaStream *stream;
    PaError err;

    MultiSynthData multi_data = {0};
    multi_data.current_time = 0.0;
    pthread_mutex_init(&multi_data.mutex, NULL);

    add_synth(&multi_data, 150.0);
    multi_data.sounds[0].params.amplitude = 0.0;
    multi_data.sounds[0].params.waveform_type = waveform_type;

    fprintf(log_file, "Initializing program... num_sounds: %d\n", multi_data.num_sounds);
    fflush(log_file);

    err = Pa_Initialize();
    if (err != paNoError)
        goto error;

    fprintf(log_file, "Initializing buffer ... num_sounds: %d\n", multi_data.num_sounds);
    fflush(log_file);

    pthread_t thread;
    pthread_create(&thread, NULL, generate_thread, &multi_data);

    pthread_t ui_thread;
    pthread_create(&ui_thread, NULL, termbox_thread, &multi_data);

    err = Pa_OpenDefaultStream(&stream,
            0,         // no input channels
            2,         // stereo output
            paFloat32, // 32 bit floating point output
            SAMPLE_RATE, FRAMES_PER_BUFFER, paCallbackSynth,
            &multi_data);

    if (err != paNoError)
        goto error;

    err = Pa_StartStream(stream);
    if (err != paNoError)
        goto error;

    pthread_join(ui_thread, NULL);
    fclose(log_file);

    err = Pa_StopStream(stream);
    if (err != paNoError)
        goto error;

    err = Pa_CloseStream(stream);
    if (err != paNoError)
        goto error;

    Pa_Terminate();

    return 0;

error:
    Pa_Terminate();
    fprintf(stderr, "An error occurred while using the portaudio stream\n");
    fprintf(stderr, "Error number: %d\n", err);
    fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));

    return 1;
}
