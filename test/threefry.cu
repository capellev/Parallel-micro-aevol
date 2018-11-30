#include "../Threefry.h"
#include "../ThreefryGPU.h"

#include <cuda.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

struct TestStateDevice
{
	std::vector<unsigned long long> rndDevice_host_;

	unsigned long long *threefryDevice_;
	unsigned long long *rndDevice_;

	int nNum_, nIndiv_;

		TestStateDevice(int nNum, int nIndiv, Threefry& host);
		~TestStateDevice();

	bool compare(const std::vector<unsigned long long>& host, const char* prefix = "");
	template<class Kernel>
	void getNum(Kernel kernel, int blocks, int threads);
};

__global__ void genThreefryDevice(unsigned long long* state, unsigned long long* rnd, int n, int nIndiv)
{
	const int indiv = blockIdx.x*blockDim.x+threadIdx.x;
	Threefry::Device gen(state, indiv, Threefry::MUTATION, nIndiv);

	const auto rndLocal = rnd+n*indiv;
	for(int a = 0; a < n; ++a)
		rndLocal[a] = gen.random(); // not coalesced

	__syncthreads(); // make sure block is synced, shared state ptr is pvisible to all
}

__global__ void genThreefryDeviceBlock(unsigned long long* state, unsigned long long* rnd, int n, int nIndiv)
{
	const int indiv = blockIdx.x;
	Threefry::DeviceCollectiveBlock gen(state, indiv, Threefry::MUTATION, nIndiv);

	const auto rndLocal = rnd+n*indiv;
	for(int a = threadIdx.x; a < n; a += blockDim.x)
		rndLocal[a] = gen.random(); // coalesced

	__syncthreads(); // make sure block is synced, shared state ptr is pvisible to all
}

void genCPU(std::vector<unsigned long long>& rndHost, Threefry& master, int nNum, int nIndiv)
{
	for(int a = 0; a < nIndiv; ++a)
	{
		Threefry::Gen gen = std::move(master.gen(a, Threefry::MUTATION));
		for(int s = 0; s < nNum; ++s)
		{
			rndHost[a*nNum+s] = gen.random_raw()[0];
		}
	}
}

int main(int argc, char* argv[])
{
	const int seed = [argc, argv](int argi=1){return argc < argi+1 ? 23 : atoi(argv[argi]);}();
	const int nIndiv = [argc, argv](int argi=2){return argc < argi+1 ? 128 : atoi(argv[argi]);}();
	const int nNum = [argc, argv](int argi=3){return argc < argi+1 ? 128 : atoi(argv[argi]);}();

	std::vector<unsigned long long> rndHost(nNum*nIndiv);

	Threefry threefry(nIndiv, 1, seed);
	threefry.initDevice();
	TestStateDevice stateSingle(nNum, nIndiv, threefry), stateBlock(nNum, nIndiv, threefry);

	bool fail = false;
	for(int a = 0; a < 2; ++a)
	{
		genCPU(rndHost, threefry, nNum, nIndiv);

		stateSingle.getNum(genThreefryDevice, nIndiv/64, 64);
		fail = !stateSingle.compare(rndHost);

		stateBlock.getNum(genThreefryDeviceBlock, nIndiv, 128);
		fail &= !stateBlock.compare(rndHost, "Block: ");

		if(fail)
		{
			std::cout << "Failed in round " << a << std::endl;
			break;
		}
	}

	return fail;
}

TestStateDevice::TestStateDevice(int nNum, int nIndiv, Threefry& host)
	: rndDevice_host_(nNum*nIndiv), nNum_(nNum), nIndiv_(nIndiv)
{
	cudaMalloc(&threefryDevice_, host.counters().size()*sizeof(unsigned long long));
	cudaMemset(threefryDevice_, 0, host.counters().size()*sizeof(unsigned long long));
	cudaMalloc(&rndDevice_, rndDevice_host_.size()*sizeof(unsigned long long));
}

bool TestStateDevice::compare(const std::vector<unsigned long long>& host, const char* prefix)
{
	assert(host.size() == rndDevice_host_.size());
	bool ret = true;
	for(int a = 0; a < nIndiv_ && ret; ++a)
	{
		for(int s = 0; s < nNum_; ++s)
		{
			const size_t idx = a*nNum_+s;
			if(host[idx] != rndDevice_host_[idx])
			{
				std::cout << prefix << "Different number at (" << a << "," << s << "): "
					<< rndDevice_host_[idx] << " vs. " << host[idx] << std::endl;
				ret = false;
			}
		}
	}
	return ret;
}

template<class Kernel>
void TestStateDevice::getNum(Kernel kernel, int blocks, int threads)
{
	kernel<<<blocks, threads>>>(threefryDevice_, rndDevice_, nNum_, nIndiv_);
	cudaMemcpy(rndDevice_host_.data(), rndDevice_, rndDevice_host_.size()*sizeof(unsigned long long), cudaMemcpyDeviceToHost);
}

TestStateDevice::~TestStateDevice()
{
	cudaFree(rndDevice_);
	cudaFree(threefryDevice_);
}
