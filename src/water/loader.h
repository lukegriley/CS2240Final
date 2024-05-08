#pragma once

#include "plant.h"

namespace water {

class Loader {
public:
    // Read a line formatted as <header name> <header value> and return the value.
    template <class T>
    static T read_header(std::istream &stream, const std::string &name) {
        if (!stream.good()) {
            throw std::runtime_error("bad stream");
        }
        std::string header_name;
        T value;
        stream >> header_name >> value;
        // Check that the correct header and value were read
        if (!stream.good()) {
            throw std::runtime_error("not a header line");
        }
        if (header_name != name) {
            throw std::runtime_error(
                "invalid header name " + header_name + ". Expected " + name);
        }
        // Return the header value
        return value;
    }

    // Expect a certain character in a stream
    struct expect {
        explicit expect(int c) : c(c) {}

        friend std::istream &operator>>(std::istream &stream, expect e) {
            using namespace std::literals;
            int x = stream.get();
            if (x != e.c) {
                throw std::runtime_error((std::ostringstream()
                                          << "expected char '" << e.c << "', received '" << x << "'"
                                          ).str());
            }
            return stream;
        }
    private:
        int c;
    };
    static void load_S_decomp(const std::string &plant_path, Plant &plant);
};

}

