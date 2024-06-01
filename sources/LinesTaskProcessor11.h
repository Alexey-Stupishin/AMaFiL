#pragma once

#include "LinesTaskQueue11.h"
#include "LinesProcessor.h"

class LQPProcessor : public ATQPProcessor
{
protected:
    LQPTask *this_task;
    CLinesProcessor *w;

public:
    LQPProcessor(LQPSupervisor *supervisor, int _id, CagmVectorField *_v, int _dir, double _step, double _relErr, double _absErr
        , double _boundAchieve, double _boundAchieveBottom, int _maxLength, int *_passed, ATQPSynchonizer *_sync) : ATQPProcessor(_id, _sync)
    {
        w = new CLinesProcessor(supervisor, _v, _dir, _step, _relErr, _absErr
            , _boundAchieve, _boundAchieveBottom, _maxLength, _passed);
    }
    virtual ~LQPProcessor()
    {
        delete w;
    }

    virtual bool setTask(ATQPTask * _task)
    {
        bool finish = ATQPProcessor::setTask(_task);
        if (!finish)
        {   // cast *_task to the real object:
            this_task = (LQPTask *)_task;
        }
        return finish;
    }

    virtual bool proceed()
    {
        double *params = this_task->getData();
        w->setTaskParams(params);
        w->queueID = this_task->getID();
        w->ActionCore();

        delete this_task;

        return true;
    }
};
