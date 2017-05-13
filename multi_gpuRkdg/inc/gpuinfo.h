/**
 * gpuinfo.h
 *
 *  Created on: 2017Äê3ÔÂ28ÈÕ
 *      Author: xyzhou
 */

#include "cudarkdgsolver.h"
#include <string.h>

class GpuInfo {
private:
	int deviceId;
	CCUDARkdgSolver solver;
public:

	GpuInfo(string config_file);

	void setDeviceId(int id);

	void runGpu();

	~GpuInfo();
};