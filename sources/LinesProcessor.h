#pragma once

#include "LinesTaskQueue11.h"

class CLinesTaskQueue;
class CagmVectorField;
class CagmRKF45;

class CLinesProcessor
{
public:
    uint32_t queueID;

protected:
    CagmRKF45 *rkf45;
    CagmVectorField *v;
    int dir;
    double step;
    double boundAchieve, boundAchieveBottom;
    int maxLength;
    int *passed;
    double *coord;
    double *linesteps;

    double point[3];

    LQPSupervisor *supervisor;

public:
    CLinesProcessor(LQPSupervisor *, CagmVectorField *, int, double, double, double, double, double, int, int *);
    virtual ~CLinesProcessor();

    virtual uint32_t setTaskParams(void * params);
    virtual uint32_t ActionCore();
};
