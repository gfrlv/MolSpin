/////////////////////////////////////////////////////////////////////////
// TaskMultiStaticSS(RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
// -- When applicable uses a Tridagonal Block Solver instead of the Armadillo solver --
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiStaticSS
#define MOD_RunSection_TaskMultiStaticSS

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskMultiStaticSS : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		bool productYieldsOnly;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiStaticSS(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiStaticSS();												   // Destructor
	};
}

#endif
