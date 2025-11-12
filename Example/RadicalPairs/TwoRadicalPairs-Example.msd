SpinSystem RPC
{
	Spin electron1
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin electron2
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin nucleus1
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}
	
	Spin nucleus2
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}

	Spin nucleus3
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}
	
	Interaction Hyperfine1
 	{
		type = Hyperfine;
 		group1 = electron1;
 		group2 = nucleus1;
		tensor= matrix(" 2.13814981 3.19255832 -2.48895215;
			3.19255832 15.45032887 -12.44778343;
			-2.48895215 -12.44778343 12.49532827");
		prefactor =3.568245455e-5;
 	}

 	Interaction Hyperfine2
 	{
 		type = Hyperfine;
 		group1 = electron2;
 		group2 = nucleus2;
		tensor=matrix(" -0.36082693 -0.0702137 -1.41518116;
						-0.0702137 -0.60153649 0.32312139;
						-1.41518116 0.32312139 50.80213093");
		prefactor =3.568245455e-5;
 	}

	Interaction Zeeman
	{
		type = Zeeman;
		field = "5e-5 0 0";
		spins = electron1, electron2;
	}


	State Singlet
	{
		spins(electron1, electron2) = |1/2,-1/2> - |-1/2,1/2>;
	}

	State T0
	{
		spins(electron1, electron2) = |1/2,-1/2> + |-1/2,1/2>;
	}

	State TPlus
	{
		spins(electron1) = |1/2>;
		spins(electron2) = |1/2>;
	}

	State TMinus
	{
		spins(electron1) = |-1/2>;
		spins(electron2) = |-1/2>;
	}
	
	State Identity
	{
	}

	Transition SingletDecay
	{
		rate = 10;
		source = Singlet;
	}
	
	Transition spinindependent_decay
	{
		rate = 1;
		source = Identity;
	}

	Transition RPtransitionSinglet 	{rate = 1; source = Singlet; targetsystem = RPD; targetstate = Singlet;}
	Transition RPtransitionTriplet 	{rate = 1; source = T0; targetsystem = RPD; targetstate = T0;}
	Transition RPtransitionTPlus 	{rate = 1; source = TPlus; targetsystem = RPD; targetstate = TPlus;}
	Transition RPtransitionTMinus 	{rate = 1; source = TMinus; targetsystem = RPD; targetstate = TMinus;}

	Properties prop
	{
		initialstate = Singlet;
	}	
		
}

SpinSystem RPD
{
	Spin electron1
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin electron2
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin nucleus1
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}
	
	Spin nucleus2
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}

	Spin nucleus3
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	

	Interaction Zeeman
	{
		type = Zeeman;
		field = "0 0 5e-5";
		spins = electron1, electron2;
	}

	Interaction Hyperfine3
	{
		type = Hyperfine;
		group1 = electron1;
		group2 = nucleus1;
		tensor=matrix(" 0.98491908 3.28010265 -0.53784491;
						3.28010265 25.88547678 -1.6335986;
						-0.53784491 -1.6335986 1.41368001 ");
		prefactor =3.568245455e-5;
	}

	Interaction Hyperfine4
	{
		type = Hyperfine;
		group1 = electron2;
		group2 = nucleus3;
		tensor=matrix(" -0.294412424 -0.0568059200 -1.02860888;
 						-0.0568059200 -0.540578469  -0.0267686240;
 						-1.02860888 -0.0267686240 50.5815320 ");
		prefactor =3.568245455e-5;
	}

	State Singlet
	{
		spins(electron1, electron2) = |1/2,-1/2> - |-1/2,1/2>;
	}

	State T0
	{
		spins(electron1, electron2) = |1/2,-1/2> + |-1/2,1/2>;
	}

	State TPlus
	{
		spins(electron1) = |1/2>;
		spins(electron2) = |1/2>;
	}

	State TMinus
	{
		spins(electron1) = |-1/2>;
		spins(electron2) = |-1/2>;
	}
	
	State Identity
	{
	}

	Transition RPtransitionSinglet 	{rate = 1; source = Singlet; 	targetsystem = RPC; targetstate = Singlet;}
	Transition RPtransitionTriplet 	{rate = 1; source = T0; 	targetsystem = RPC; targetstate = T0;}
	Transition RPtransitionTPlus 	{rate = 1; source = TPlus; 		targetsystem = RPC; targetstate = TPlus;}
	Transition RPtransitionTMinus 	{rate = 1; source = TMinus; 	targetsystem = RPC; targetstate = TMinus;}


	Transition spinindependent_decay2
	{
		rate = 1;
		source = Identity;
	}

	Properties prop
	{
		initialstate = zero;
	}

}

Run
{
	Task CalculateQuantumYeild
	{
		type = MultiStaticSS-timeevolution; 
		//MultiStaticSS
		logfile = "logfile.log";
		datafile = "result.dat";
		timestep = 1e-4;
		totaltime = 12;
		propagator = exp; 
		//RK45
	}
}

Settings
{
	Settings general
	{
		steps = 1;
	}
}

