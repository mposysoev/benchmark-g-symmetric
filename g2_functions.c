#include <math.h>

// C function: distance_cutoff
double distance_cutoff(double distance, double rcutoff) {
    if (distance > rcutoff) {
        return 0.0;
    }
    else {
        return 0.5 * (cos(M_PI * distance / rcutoff) + 1.0);
    }
}

// C function: compute_g2_element
double compute_g2_element(double distance, double eta, double rcutoff, double rshift) {
    if (distance > 0.0) {
        return exp(-eta * pow(distance - rshift, 2)) * distance_cutoff(distance, rcutoff);
    }
    else {
        return 0.0;
    }
}

