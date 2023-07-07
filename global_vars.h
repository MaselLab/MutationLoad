#ifndef GLOBAL_VARS_H_INCLUDED
#define GLOBAL_VARS_H_INCLUDED 1

#define VERBOSE 1
#define VERYVERBOSE 1
#define MISCELLANEOUS 1
#define INDIVIDUALWIDATA 1
#define OVERLAPPINGGENERATIONS 1 //probably should be an input argument at some point.
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }//error checking for tree sequence recording

#endif // GLOBAL_VARS_H_INCLUDED

