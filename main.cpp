#include <getopt.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Kmethod.h"
#include "PAPI.h"
#include "euler.h"
#include "eulerSimplified.h"

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
                               {"increase", no_argument, 0, 'i'},
                               {"constant", no_argument, 0, 'c'},
                               {0, 0, 0, 0}};

    int opt;
    int N;
    std::string method;
    bool hasN = false, hasMethod = false;
    bool increase = true;

    double t = 0.01f;
    double f = 2000.f;
    double fs = 48000.f;
    double ampl = 0.5;

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
                if(method != "Kmethod" && method != "Euler" &&
                   method != "EulerSimplified") {
                    std::cerr << "Error: Invalid method. Use 'Kmethod', "
                                 "'Euler' or 'EulerSimplified'.\n";
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
                    std::cerr << "Error: The frequency of the input signal "
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
            case 'i':
                increase = true;
                break;
            case 'c':
                increase = false;
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
    } else if(method == "EulerSimplified" && increase) {
        std::cerr << "Error: The simplified Euler method can only be used with "
                     "a single nonlinear component."
                  << std::endl;
        print_usage(argv[0]);
    }

    std::cout << "N=" << N << "\n";
    std::cout << "increase=" << increase << "\n";
    std::cout << "Method=" << method << "\n";
    std::cout << "t=" << t << "\n";
    std::cout << "Amplitude=" << ampl << "\n";
    std::cout << "Frequency=" << f << "\n";
    std::cout << "Sampling rate=" << fs << "\n";

    std::vector<double> input(generateInput(t, f, ampl, fs));
    std::vector<double> output(input.size(), 0.f);

    int event_set = PAPI::setup(PAPI_TOT_CYC);

    if(method == "Kmethod") {
        Kmethod kmethod(N, fs, 0.5f, increase);

        // Start counting the events
        PAPI::start(event_set);

        kmethod.process(input, output);
    } else if(method == "Euler") {
        Euler euler(N, fs, increase);

        // Start counting the events
        PAPI::start(event_set);

        euler.process(input, output);
    } else if(method == "EulerSimplified") {
        EulerSimplified eulerSimplified(N, fs, increase);

        // Start counting the events
        PAPI::start(event_set);
        eulerSimplified.process(input, output);
    }

    // Output the number of floating-point operations
    printf("Floating-point operations: %lld\n", PAPI::close(event_set));

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
    std::cerr
        << "Usage: " << prog_name
        << " --N <value> --method <Kmethod|Euler|EulerSimplified> [options]\n"
        << "Options:\n"
        << "  --input-length <value>   Length of input signal (default: 1 s)\n"
        << "  --amplitude <value>      Amplitude (default: 5)\n"
        << "  --frequency <value>      Frequency (default: 440 Hz)\n"
        << "  --fs <value>             Sampling rate (default: 48000 Hz)\n"
        << "  --increase               Simultaneously increase nonlinear "
           "elements with N (default)\n"
        << "  --constant               Use a single nonlinear element\n"
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