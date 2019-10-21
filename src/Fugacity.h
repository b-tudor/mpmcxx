#pragma once


class Fugacity
{
	Fugacity() {};
	~Fugacity() {};
	//constants for peng_robinson equation of state
	typedef struct _peng_robinson_constants {
		double Tc;
		double Pc;
		double w;
	} peng_robinson_constants;

public:

	static double h2_fugacity(double temperature, double pressure);

	// use the semi-empirical BACK equation of state 
	// Tomas Boublik, "The BACK equation of state for hydrogen and related compounds", Fluid Phase Equilibria, 240, 96-100 (2005) 
	static double h2_fugacity_back(double temperature, double pressure);
	
	static double h2_comp_back(double temperature, double pressure);

	// calculate the fugacity correction for H2 for 0 C and higher this empirical relation follows from: 
	// H.R. Shaw, D.F. Wones, American Journal of Science, 262, 918-929 (1964) 
	static double h2_fugacity_shaw(double temperature, double pressure);
	
	// fugacity for low temperature and up to 200 atm: Zhou, Zhou, Int. J. Hydrogen Energy, 26, 597-601 (2001) 
	static double h2_fugacity_zhou( double pressure );
	

	


	static double ch4_fugacity(double temperature, double pressure);

	// Incorporate BACK EOS 
	static double ch4_fugacity_back(double temperature, double pressure);
	
	static double ch4_comp_back(double temperature, double pressure);

	// Apply the Peng-Robinson EoS to methane 
	static double ch4_fugacity_PR(double temperature, double pressure);
	
	



	static double n2_fugacity(double temperature, double pressure);

	// Incorporate BACK EOS 
	static double n2_fugacity_back(double temperature, double pressure);
	
	static double n2_comp_back(double temperature, double pressure);

	// Apply the Peng-Robinson EoS to N2 
	static double n2_fugacity_PR(double temperature, double pressure);
	
	// Apply the Zhou function to N2 
	static double n2_fugacity_zhou( double pressure );
	


	
	static void get_peng_robinson_constants(_peng_robinson_constants peng_robinson_constants, std::string species);
	// reads in temperature in K, and pressure (of the ideal gas in the resevoir) in atm 
	// return the CO2 fugacity via the Peng-Robinson equation of state
	// else return 0.0 on error - I don't have an error statement
	// units are K, atmstatic
	static double get_peng_robinson_fugacity(double temperature, double pressure, std::string species);
	
};