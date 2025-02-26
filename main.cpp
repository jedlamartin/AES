#include <getopt.h>
#include <math.h>
#include <papi.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "K.h"
#include "euler.h"

#define _USE_MATH_DEFINES

std::vector<double> generateInput(const double t,
                                  const double f,
                                  const double ampl,
                                  const double fs);

void print_usage(const char* prog_name);

void write_to_file(std::vector<double>& input, std::vector<double>& output);

int main(int argc, char* argv[]) {
    const option longopts[] = {{"method", required_argument, 0, 'm'},
                               {"N", required_argument, 0, 'n'},
                               {"input-length", required_argument, 0, 't'},
                               {"input-amplitude", required_argument, 0, 'a'},
                               {"input-frequency", required_argument, 0, 'f'},
                               {"samplerate", required_argument, 0, 's'},
                               {"help", no_argument, 0, 'h'},
                               {0, 0, 0, 0}};
    int opt;
    int N;
    std::string method;
    bool hasN = false, hasMethod = false;

    double t = 1.f;
    double f = 440.f;
    double fs = 48000.f;
    double ampl = 5;

    while((opt = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(opt) {
            case 'n':
                N = std::stoi(optarg);
                if(N < 1) {
                    std::cerr << "Error: N must be greater than 0!\n";
                    exit(1);
                }
                hasN = true;
                break;
            case 'm':
                method = std::string(optarg);
                if(method != "Kmethod" && method != "Euler") {
                    std::cerr
                        << "Error: Invalid method. Use 'Kmethod' or 'Euler'.\n";
                    exit(1);
                }
                hasMethod = true;
                break;
            case 't':
                t = std::stod(optarg);
                if(t <= 0) {
                    std::cerr << "Error: The length of the input signal must "
                                 "be greater than 0!\n";
                    exit(1);
                }
                break;
            case 'a':
                ampl = std::stod(optarg);
                if(ampl <= 0) {
                    std::cerr << "Error: The amplitude of the input signal "
                                 "must be greater than 0!\n";
                    exit(1);
                }
                break;
            case 'f':
                f = std::stod(optarg);
                if(f <= 0) {
                    std::cerr << "Error: The frequency if the input signal "
                                 "must be greater than 0!\n";
                    exit(1);
                }
                break;
            case 's':
                fs = std::stod(optarg);
                if(fs <= 0) {
                    std::cerr
                        << "Error: Sampling rate must be greater than 0!\n";
                    exit(1);
                }
                break;
            case 'h':
                print_usage(argv[0]);
                break;
            default:
                std::cerr << "Error: Invalid argument." << std::endl;
                print_usage(argv[0]);
                break;
        }
    }

    if(!hasN || !hasMethod) {
        std::cerr << "Error: Missing mandatory parameters." << std::endl;
        print_usage(argv[0]);
    }

    std::cout << "N=" << N << "\n";
    std::cout << "Method=" << method << "\n";
    std::cout << "t=" << t << "\n";
    std::cout << "Amplitude=" << ampl << "\n";
    std::cout << "Frequency=" << f << "\n";
    std::cout << "Sampling rate=" << fs << "\n";

    std::vector<double> input(generateInput(t, f, ampl, fs));
    std::vector<double> output(input.size(), 0.f);

    // Initialize PAPI
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if(retval != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI Library initialization error!\n");
        exit(1);
    }

    // Create an event set
    int event_set = PAPI_NULL;
    retval = PAPI_create_eventset(&event_set);
    if(retval != PAPI_OK) {
        fprintf(stderr, "PAPI EventSet creation error!\n");
        exit(1);
    }

    // Add the floating-point operation event
    retval = PAPI_add_event(event_set, PAPI_TOT_CYC);    // PAPI_FP_OPS);
    if(retval != PAPI_OK) {
        fprintf(stderr, "PAPI Add Event error!\n");
        exit(1);
    }

    if(method == "Kmethod") {
        Kmethod kmethod(N, fs, 10e-6);

        // Start counting the events
        retval = PAPI_start(event_set);
        if(retval != PAPI_OK) {
            fprintf(stderr, "PAPI Start error!\n");
            exit(1);
        }

        kmethod.process(input, output);
    } else if(method == "Euler") {
        Euler euler(N, fs);

        // Start counting the events
        retval = PAPI_start(event_set);
        if(retval != PAPI_OK) {
            fprintf(stderr, "PAPI Start error!\n");
            exit(1);
        }

        euler.process(input, output);
    }

    // Stop counting and get the result
    long long int fp_ops;
    retval = PAPI_stop(event_set, &fp_ops);
    if(retval != PAPI_OK) {
        fprintf(stderr, "PAPI Stop error!\n");
        exit(1);
    }

    // Output the number of floating-point operations
    printf("Floating-point operations: %lld\n", fp_ops);

    // Cleanup PAPI resources
    retval = PAPI_cleanup_eventset(event_set);
    if(retval != PAPI_OK) {
        fprintf(stderr, "PAPI Cleanup error!\n");
    }
    retval = PAPI_destroy_eventset(&event_set);
    if(retval != PAPI_OK) {
        fprintf(stderr, "PAPI Destroy error!\n");
    }

    write_to_file(input, output);

    return 0;
}

std::vector<double> generateInput(const double t,
                                  const double f,
                                  const double ampl,
                                  const double fs) {
    int i = 0;
    int length = static_cast<int>(t * fs);
    std::vector<double> input(length);
    double theta = 2 * M_PI * f / fs;
    while(i < length) {
        input[i] = ampl * std::sin(theta * i);
        i++;
    }
    return input;
}

void print_usage(const char* prog_name) {
    std::cerr << "Usage: " << prog_name
              << " --N <value> --method <Kmethod|Euler> [options]\n";
    std::cerr
        << "Options:\n"
        << "  --input-length <value>   Length of input signal (default: 1 s)\n"
        << "  --amplitude <value>      Amplitude (default: 5)\n"
        << "  --frequency <value>      Frequency (default: 440 Hz)\n"
        << "  --fs <value>             Sampling rate (default: 48000 Hz)\n"
        << "  --help                   Show this message and exit\n";
    exit(1);
}

void write_to_file(std::vector<double>& input, std::vector<double>& output) {
    std::ofstream iFile("input.txt");
    std::ofstream oFile("output.txt");

    for(int i = 0; i < input.size(); i++) {
        iFile << input[i] << std::endl;
        oFile << output[i] << std::endl;
    }

    iFile.close();
    oFile.close();
}
