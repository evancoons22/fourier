## fast fourier transform
- initial project:  audio filters with fast fourier transforms (in fft.c) 
- this has turned into a digital synth project to generate waveforms first

Project path: 
1. finish digital synth with midi input and adjustable parameters
2. build concurrent i/o program to transform sound with dynamic filters (fft, apply filter, ifft)
3. embedded system? 


Usage: `$ make && ./synth -type <waveform_type> -octave <octave_number>`

## example

<img src="https://github.com/evancoons22/FT/blob/main/output.gif">



To dos: 
- [X] fix black key highlight.
- [X] new waveform types (sine squared, noise, pulse)
- [X] key params in UI
- [X] Larger keyboard in UI
- [X] display waveform type in constants
- [X] Key highlights removed when pressed again
- [X] Remove sound functionality. (on second press)
- [X] add args to the program for waveform type and octave
    - [X] waveform type
    - [X ] octave
- [X] remove dummy sound (initialized but silent)
- [ ] implement decay for keys and then naturally stop them
- [ ] properly remove sounds after release time is up
- [ ] instead of toggle, press and release for sustain and release



