#ifndef EECS_FAST_USINGS_H
#define EECS_FAST_USINGS_H

namespace fastEEC{
    using namespace boost;
    using namespace std;

    using uvec = std::vector<unsigned>;

    using axis_t = histogram::axis::variable<double>;
    using axisptr = std::shared_ptr<axis_t>;

    using umat = multi_array<unsigned, 2>;
    using umatptr = std::shared_ptr<umat>;
}

#endif
