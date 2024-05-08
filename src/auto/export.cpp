#include "export.h"

#include <iomanip>

std::ostream &operator<<(std::ostream &out, const rod::Tree &tree) {
    out << "rods " << tree.rods.size() << std::endl;
    for (const rod::Rod &rod : tree.rods) {
        out << std::setprecision(0) << rod.index << ","
            << std::setprecision(3) << rod.direction(tree).norm() << ","
            << std::setprecision(3) << rod.orientation.w() << ","
            << std::setprecision(3) << rod.orientation.x() << ","
            << std::setprecision(3) << rod.orientation.y() << ","
            << std::setprecision(3) << rod.orientation.z() << std::endl;
    }
    return out;
}
