#pragma once

int interleave_even(int n, int i) {
    assert(n % 2 == 0);
    if (i % 2 == 0) {
        return i;
    } else {
        return n - i;
    }
}

// Apply interleaving order
int interleave(int n, int i) {
    if (n % 2 == 0) {
        return interleave_even(n, i);
    } else {
        if (i == 0) {
            return 0;
        } else {
            return 1 + interleave_even(n - 1, i - 1);
        }
    }
}
