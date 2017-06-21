#include "../inc/gpuinfo.h"

GpuInfo::GpuInfo(string config_file):
deviceId(-1)
{
	solver.config_file = config_file;
}

void GpuInfo::setDeviceId(int id) {
	deviceId = id;
}

void GpuInfo::runGpu() {
	solver.run();
}

GpuInfo::~GpuInfo() {

}
