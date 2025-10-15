#include <math.h>

inline double Q(int l, int m)
{
	return pow(-1, (m+fabs(m))/2)* sqrt((2*l+1)*tgamma(l+fabs(m)+1)/(2*tgamma(l-fabs(m)+1)));
}

inline double w_adk(const double B_m, const double kappa, const double Z_c, const double F, const int m)
{
	const double first_part= pow(B_m, 2)/(pow(2, fabs(m))*tgamma(fabs(m)+1)*pow(kappa,2*Z_c/kappa-1));
	const double second_part= pow(2*pow(kappa,3)/F, 2*Z_c/kappa-fabs(m)-1);
	const double third_part= exp(-2*pow(kappa,3)/(3*F));
	return first_part*second_part*third_part;
}

inline double Probab_ioniz_per_time_unit_ADK_N2 (const double E_t) 
{
	const double Ip=15.58; // first ionization potential for N2 molecule
	const double I_H=13.59; // ionization potential for hydrogen atom
    const double kappa = sqrt(Ip/I_H);
    const int Z_c = 1;
	const int m = 0; // magnetic quantum number: for sigma orbitals, m=0; for pi orbitals, m=1
	// From Phys. Rev. A 81, 033423 (2010)
	const double C_0_m = 2.68;
	const double C_2_m = 1.1;
	const double C_4_m = 0.06;
	// From Tong et al., Phys. Rev. A 66, 033402 (2002)
	//const double C_0_m = 2.02;
	//const double C_2_m = 0.78;
	//const double C_4_m = 0.04;
	const double B_m = C_0_m*Q(0,m)+C_2_m*Q(2,m)+C_4_m*Q(4,m);
	return 1/3*w_adk(B_m, kappa, Z_c, E_t, m); // 1/3 is due to averaging over the angles, Phys. Rev. A 66, 033402 (2002)
}


inline double Probab_ioniz_per_time_unit_ADK_O2 (const double E_t) 
{
	const double Ip=12.06; // first ionization potential for O2 molecule
    const double I_H=13.59; // ionization potential for hydrogen atom
    const double kappa = sqrt(Ip/I_H);
    const int Z_c = 1;
	const int m = 1; // magnetic quantum number: for sigma orbitals, m=0; for pi orbitals, m=1
	// From Phys. Rev. A 81, 033423 (2010)
	const double C_2_m = 0.52;
	const double C_4_m = 0.03;
	// From Tong et al., Phys. Rev. A 66, 033402 (2002)
	//const double C_2_m = 0.62;
	//const double C_4_m = 0.03;
	const double B_m = C_2_m*Q(2,m)+C_4_m*Q(4,m);
	return 2*w_adk(B_m, kappa, Z_c, E_t, m); // 2.0 is due to averaging over the angles, Phys. Rev. A 66, 033402 (2002)	    
}
