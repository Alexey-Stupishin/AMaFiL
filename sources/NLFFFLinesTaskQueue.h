#pragma once

#include "LinesTaskQueue.h"

class CNLFFFLinesTaskQueue : public CLinesTaskQueue
{
protected:
    double *seeds;

    double *params;

private:
    void parseGlobalIdx(int idx, int *kx, int *ky, int *kz)
    {
        *kx = idx % NF[0];
        int kyz = (idx-*kx) / NF[0];
        *ky = kyz % NF[1];
        *kz = (kyz-*ky) / NF[1];
    }

    void InitOutput(int id)
    {
        if (physLength)
            physLength[id] = 0;
        if (avField)
            avField[id] = 0;
        if (startIdx)
            startIdx[id] = 0;
        if (endIdx)
            endIdx[id] = 0;
        if (voxelStatus)
            voxelStatus[id] = Status::None;
        if (passed)
            passed[id] = false;
    }

public:
    CNLFFFLinesTaskQueue(CagmVectorFieldOps *_field, double *_seeds, int _Nseeds,
        uint32_t _cond = 0x3, double _chromoLevel = 0,
        double *_physLength = nullptr, double *_avField = nullptr,
        int *_voxelStatus = nullptr, int *_codes = nullptr,
        int *_startIdx = nullptr, int *_endIdx = nullptr, int *_apexIdx = nullptr,
        int maxResult = 50000,
        uint64_t _maxCoordLength = 0, int *_linesLength = nullptr, double *_coords = nullptr, uint64_t *_linesStart = nullptr, int *_linesIndex = nullptr, int *seedIdx = nullptr)
        : CLinesTaskQueue(_field, 
            _cond, _chromoLevel,
            _physLength, _avField,
            _voxelStatus, _codes, 
            _startIdx, _endIdx, _apexIdx, maxResult,
            _maxCoordLength, _linesLength, _coords, _linesStart, _linesIndex, seedIdx)
      , seeds(_seeds)
    {
        Nseeds = _Nseeds;
        globalIdx = nullptr;

        autoParams = (Nseeds <= 0 || !seeds);
        field->dimensions(NF);
        //memcpy(NF, field->GetDimensions(), 3*sizeof(int));
        int nVox = (autoParams ? NF[0]*NF[1]*NF[2] : Nseeds);

        InitQueue(nVox);

        if (cond != Conditions::NoCond)
            passed = new int[nVox];

        if (!autoParams)
        {
            globalIdx = new int[nVox];
            params = new double[3 * sizeof(double) * nVox];
            for (int id = 0; id < Nseeds; id++)
            {
                InitOutput(id);
                params[3*id + 0] = seeds[3*id + 0];
                params[3*id + 1] = seeds[3*id + 1];
                params[3*id + 2] = seeds[3*id + 2];
                globalIdx[id] = getGlobalID(seeds + 3*id);
            }
        }
        else
        {
            params = new double[3 * sizeof(double)];

            for (int kz = 0; kz < NF[2]; kz++)
                for (int ky = 0; ky < NF[1]; ky++)
                    for (int kx = 0; kx < NF[0]; kx++)
                        InitOutput(getGlobalID(kx, ky, kz));
        }
    }

    virtual ~CNLFFFLinesTaskQueue()
    {
        delete[] params;
        delete[] globalIdx;
        delete[] passed;
    }

    virtual void * GetParams(uint32_t _queueID)
    {
        if (!autoParams)
            return params + _queueID*3;
        else
        {
            int kx, ky, kz;
            parseGlobalIdx(_queueID, &kx, &ky, &kz);
            params[0] = kx;
            params[1] = ky;
            params[2] = kz;
            return params;
        }
    }
};
