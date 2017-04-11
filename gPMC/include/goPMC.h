/*
 *
 * file  goPMC.h
 * goPMC interfaces declaration.
 *
 * author Nan Qin
 *
 * last update on 9/21/2016
 *
 */

#pragma once

#ifdef GOPMC_EXPORTS
#define GOPMC_API __declspec(dllexport) 
#else
#define GOPMC_API __declspec(dllimport) 
#endif

#include "cl.hpp"
#include <string>
#include <ctime>
#include <vector>
#include <map>

class OpenCLStuff;
class Phantom;
class ParticleStatus;
class Proton;
class Secondary;
class DensCorrection;
class MacroCrossSection;
class MSPR;
class RSPW;

namespace goPMC{
	class GOPMC_API MCEngine
	{
	public:
		MCEngine();
		virtual ~MCEngine();

		void initializeComputation(cl::Platform &, cl::Device &);
		void initializePhysics(const std::string & dir = "input");
		void initializePhantom(const std::string & dicomDir);
		void simulate(cl_float *, cl_float3 *, cl_float3 *, cl_float *, cl_uint, std::string &);
		void getResult(std::vector<cl_float> &, std::vector<cl_float> &);
		void clearCounter();
	private:
		ParticleStatus * primary, *secondary;
		Phantom * phantom;
		MacroCrossSection * macroSigma;
		MSPR * mspr;
		RSPW * rspw;
		DensCorrection * densCorrection;
		OpenCLStuff * stuff;
		std::clock_t start;

		static const std::map<std::string, cl_int> types;
		std::string quantity;
		cl_uint nPaths;
		cl_double totalWeight;
	};

}

