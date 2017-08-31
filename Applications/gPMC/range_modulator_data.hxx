#include <vector>

/*
####################################################################
# Range Modulator for UFPTI&KNCC Gantry1 (modified from MGH nozzle by dongho shin)
####################################################################
# RM_number    RM_ID         Wheel     track        
# 2            1             SW1       Track1
# 2            2             SW1       Track2
# 3            3             SW1       Track3
# 8            4             SW2       Track1
# 5            5             SW2       Track2
# 4            6             SW2       Track3
# 6            7             SW3       Track1
# 7            8             SW3       Track2
# 1            9             SW3       Track3 (for uniform scanning)
*/


// Angles in deg
const std::vector<const double> W1T1_angles{
	5.00,   114.13, 146.54, 175.66, 197.17,
	216.15, 231.94, 246.00, 258.64, 270.60,
	282.24, 294.56, 306.17, 324.00
};

// Thickness in cm
const std::vector<const double> W1T1_heigths{
	7.7717, 8.2501, 8.7268, 9.2017, 9.6747,
	10.1455, 10.6142, 11.0811, 11.5326, 11.9494,
	12.3663, 12.7832, 13.2003, NULL // NULL >= 6.5 cm brass
};

/* d: Ge / RMW1_TR1_Offset = 0.0 deg + Ge / RM_SW / InnerTrackOffset
d: Ge / RMW1_TR1_Lexan_Brass / Offset = Ge / RMW1_TR1_Offset deg */
const char* W1T1_material = "lexan";


/* # RMW1_TR1_th_lead
d: Ge / RMW1_TR1_Lead / Offset = Ge / RMW1_TR1_Offset deg */
// Lead Angles in deg
const std::vector<const double> W1T1_lead_angles{
	3.00, 114.13, 146.54, 175.66, 197.17,
	216.15, 231.94, 246.00, 258.64 
};

// Lead Thickness in cm
const std::vector<const double> W1T1_lead_heigths{
	0.0893, 0.0766, 0.0642, 0.0522, 0.0405,
	0.0291, 0.0182, 0.0075, 0.0000
};

/* ###### Wheel 1 Track 2: 15 steps ######
# RMW1_TR2_weights
# RMW1_TR2_th_lexan[]
# Stop block(Brass)
# 10.5*degree - (-2.67)
d:Ge / RMW1_TR2_Offset = 0.0 deg #13.17 deg + Ge / RM_SW / MiddleTrackOffset
d : Ge / RMW1_TR2_Lexan_Brass / Offset = Ge / RMW1_TR2_Offset deg */
// Angles in deg
const std::vector<const double> W1T2_angles{
	33.09, 152.94, 188.22, 217.83, 240.90,
	260.15, 276.74, 291.03, 303.82, 315.03,
	325.73, 334.44, 344.26, 349.99, 357.06
};

// Thickness in cm
const std::vector<const double> W1T2_heigths{
	4.4355, 5.1476, 5.8601, 6.5729, 7.2858,
	7.9988, 8.7116, 9.4242, 10.1364, 10.8481,
	11.5590, 12.2690, 12.9778, 13.6853, NULL // NULL >= 6.5 cm brass
};

/* d: Ge / RMW1_TR1_Offset = 0.0 deg + Ge / RM_SW / InnerTrackOffset
d: Ge / RMW1_TR1_Lexan_Brass / Offset = Ge / RMW1_TR1_Offset deg */
const char* W1T2_material = "lexan";


/*
# RMW1_TR2_th_lead
# RM1_start_angle = -2.67
# start0 = 10.5 - (start_angle) = > offset
# In order to calculate appropriate values.
d:Ge / RMW1_TR2_Lead / Offset = Ge / RMW1_TR2_Offset deg */
// Lead Angles in deg
const std::vector<const double> W1T2_lead_angles{
	31.09, 152.94, 188.22, 217.83, 240.90,
	260.15, 276.74, 291.03, 303.82, 315.03,
	325.73, 334.44, 344.26, 349.99, 357.06
};

// Lead Thickness in cm
const std::vector<const double> W1T2_lead_heigths{
	0.8395, 0.8046, 0.7696, 0.7345, 0.6993,
	0.6640, 0.6288, 0.5935, 0.5582, 0.5230,
	0.4879, 0.4529, 0.4181, 0.3834, 0.0000
};


/*###### Wheel 1 Track 3: 22 steps ######
# RMW1_TR3_weights
# RMW1_TR3_th_lexan[]  Stop block(Brass)
d:Ge / RMW1_TR3_Offset = 0.0 deg #13.22 deg + Ge / RM_SW / OuterTrackOffset
d : Ge / RMW1_TR3_Lexan_Brass / Offset = Ge / RMW1_TR3_Offset deg */
// Angles in deg
const std::vector<const double> W1T3_angles{
	4.12, 37.12, 136.76, 168.87, 194.80, 214.48,
	231.46, 246.30, 259.15, 270.43, 280.70,
	290.22, 298.98, 307.06, 314.51, 321.45,
	327.99, 334.07, 339.92, 345.28, 351.66,
	356.83
};

// Thickness in cm
const std::vector<const double> W1T3_heigths{
	NULL, 4.4437, 4.9485, 5.4528, 5.9561, 6.4587,
	6.9604, 7.4609, 7.9602, 8.4580, 8.9544,
	9.4488, 9.9413, 10.4317, 10.9196, 11.4047,
	11.8870, 12.3671, 12.8286, 13.2539, 13.6793,
	14.1050 // NULL >= 6.5 cm brass
};

/* d: Ge / RMW1_TR1_Offset = 0.0 deg + Ge / RM_SW / InnerTrackOffset
d: Ge / RMW1_TR1_Lexan_Brass / Offset = Ge / RMW1_TR1_Offset deg */
const char* W1T3_material = "lexan";


/*
# RMW1_TR3_th_lead
d : Ge / RMW1_TR3_Lead / Offset = Ge / RMW1_TR3_Offset deg */
// Lead Angles in deg
const std::vector<const double> W1T3_lead_angles{
	35.12, 136.76, 168.87, 194.80, 214.48,
	231.46, 246.30, 259.15, 270.43, 280.70,
	290.22, 298.98, 307.06, 314.51, 321.45,
	327.99, 334.07, 339.92
};

// Lead Thickness in cm
const std::vector<const double> W1T3_lead_heigths{
	0.2379, 0.2218, 0.2059, 0.1901, 0.1745,
	0.1589, 0.1436, 0.1285, 0.1137, 0.0991,
	0.0848, 0.0709, 0.0574, 0.0443, 0.0317,
	0.0196, 0.0079, 0.0000
};

/*###########	Wheel 2 Track 1: 31 steps ####################
### RM track data final check by dongho Shin at May 13 2014 ###
# RMW2_TR1_weights
# RMW2_TR1_th_lexan[]
# Stop block(Brass)
d:Ge / RMW2_TR1_Offset = 0 deg + Ge / RM_SW / InnerTrackOffset
d: Ge / RMW2_TR1_Lexan_Brass / Offset = Ge / RMW2_TR1_Offset deg */
// Angles in deg
const std::vector<const double> W2T1_angles{
	5.00,  75.51, 101.42, 120.91, 135.92, 149.12, 160.53, 170.61, 179.71, 188.42,
	196.69, 204.63, 212.15, 219.50, 226.64, 233.34, 239.54, 245.30, 250.99, 256.77,
	262.92, 268.90, 275.47, 282.31, 288.99, 295.57, 302.44, 308.99, 315.15, 320.93,
	326.33
};

// Thickness in cm
const std::vector<const double> W2T1_heigths{
	0.071,  0.613,  1.153,  1.719,  2.231,  2.791,  3.303,  3.857,  4.368,  4.914,
	5.423,  5.962,  6.467,  6.998,  7.512,  8.018,  8.522,  9.014,  9.494,  9.938,
	10.467, 10.929, 11.394, 11.868, 12.330, 12.790, 13.250, 13.710, 14.170, 14.630,
	NULL // NULL >= 6.5 cm brass
};

/* d: Ge / RMW1_TR1_Offset = 0.0 deg + Ge / RM_SW / InnerTrackOffset
d: Ge / RMW1_TR1_Lexan_Brass / Offset = Ge / RMW1_TR1_Offset deg */
const char* W2T1_material = "lexan";


/*
# RMW2_TR1_th_lead
d : Ge / RMW2_TR1_Lead / Offset = Ge / RMW2_TR1_Offset deg */
// Lead Angles in deg
const std::vector<const double> W2T1_lead_angles{
	3.00,  75.51, 101.42, 120.91, 135.92, 149.12, 160.53, 170.61, 179.71, 188.42,
	196.69, 204.63, 212.50, 219.50, 226.64, 233.34, 239.54, 245.30, 250.99, 256.77
};

// Lead Thickness in cm
const std::vector<const double> W2T1_lead_heigths{
	0.244, 0.227, 0.211, 0.189, 0.179, 0.159, 0.148, 0.129, 0.119, 0.101,
	0.091, 0.075, 0.065, 0.051, 0.040, 0.030, 0.021, 0.015, 0.010, 0.014
};

				 

/*###########	Wheel 2 Track 2: 23 steps ####################
# RMW2_TR2_weights
# RMW2_TR2_th_lexan[]
# Stop block(Brass)
d:Ge / RMW2_TR2_Offset = 0.0 deg # 13.17 deg + Ge / RM_SW / MiddleTrackOffset
d:Ge / RMW2_TR2_Lexan_Brass / Offset = Ge / RMW2_TR2_Offset deg*/
// Angles in deg
const std::vector<const double> W2T2_angles{
	33.09, 127.42, 158.91, 182.78, 202.07, 217.94, 231.98, 244.46, 255.57, 265.87,
	275.28, 284.26, 292.74, 300.96, 308.47, 315.87, 322.99, 329.23, 335.22, 341.06,
	346.47, 351.66, 357.09
};

/*#NOTICE : non zero height appeared first!
#see : start, RMW2_TR2_weights[], RMW2_Tr2_start_angle ..
#may requires divide angles!*/
// Thickness in cm
const std::vector<const double> W2T2_heigths{
	0.106, 0.791, 1.520, 2.164, 2.887, 3.531, 4.209, 4.918, 5.562, 6.262,
	6.903, 7.594, 8.254, 8.913, 9.567, 10.235, 10.883, 11.524, 12.165, 12.741,
	13.317, 13.893, NULL // NULL >= 5 cm brass
};

const char* W2T2_material = "lexan";

/*
# RMW2_TR2_th_lead
d : Ge / RMW2_TR2_Lead / Offset = Ge / RMW2_TR2_Offset deg */
// Lead Angles in deg
const std::vector<const double> W2T2_lead_angles{
	31.09, 127.42, 158.91, 182.78, 202.07, 217.94, 231.98, 244.46, 255.57, 265.87,
	275.28, 284.26, 292.74, 300.96, 308.47, 315.87, 322.99, 329.23, 335.22
};

// Lead Thickness in cm
const std::vector<const double> W2T2_lead_heigths{
	0.317, 0.297, 0.269, 0.256, 0.229, 0.216, 0.197, 0.172, 0.160, 0.136,
	0.124, 0.102, 0.087, 0.071, 0.056, 0.038, 0.025, 0.012, 0.000
};


/* ###########	Wheel 2 Track 3: 21 steps ####################
# RMW2_TR3_weights
# RMW2_TR3_th_lexan[]
# Stop block(Brass)
d:Ge / RMW2_TR3_Offset = 0.0 deg # 13.22 deg + Ge / RM_SW / OuterTrackOffset
d : Ge / RMW2_TR3_Lexan_Brass / Offset = Ge / RMW2_TR3_Offset deg*/
// Angles in deg
const std::vector<const double> W2T3_angles{
	4.12, 37.12, 146.00, 183.23, 210.46, 231.72, 249.27, 264.15, 276.80, 288.18, 298.33,
	307.29, 315.50, 322.81, 329.52, 335.56, 341.18, 346.35, 351.18, 355.71, 360.02
};

// Thickness in cm
const std::vector<const double> W2T3_heigths{
	NULL, 0.108, 0.856, 1.603, 2.353, 3.096, 3.835, 4.619, 5.313, 6.044, 6.812,
	7.504, 8.260, 8.978, 9.719, 10.405, 11.135, 11.837, 12.538, 13.233, 13.943 // NULL >= 5 cm brass
};

const char* W2T3_material = "lexan";

/*
# RMW2_TR3_th_lead
d : Ge / RMW2_TR3_Lead / Offset = Ge / RMW2_TR3_Offset deg */
// Lead Angles in deg
const std::vector<const double> W2T3_lead_angles{
	35.12, 146.00, 183.23, 210.46, 231.72, 249.27, 264.15, 276.80, 288.18, 298.33,
	307.29, 315.50, 322.81, 329.52, 335.56, 341.18, 346.35, 351.18, 355.71, 360.02
};

// Lead Thickness in cm
const std::vector<const double> W2T3_lead_heigths{
	0.425, 0.399, 0.374, 0.347, 0.322, 0.297, 0.264, 0.249, 0.226, 0.195,
	0.180, 0.152, 0.132, 0.107, 0.093, 0.070, 0.052, 0.035, 0.019, 0.000
};

/*###########	Wheel 3 Track 1: 14 steps ####################
d : Ge / RMW3_TR1_Offset = 0.0 deg #13.0 deg + Ge / RM_SW / InnerTrackOffset
d : Ge / RMW3_TR1_Carbon_Brass / Offset = Ge / RMW3_TR1_Offset deg*/
// Angles in deg
const std::vector<const double> W3T1_angles{
	5.00,  90.86, 123.91, 150.89, 172.83, 191.84, 209.99, 227.36, 244.30, 261.37,
	283.82, 296.92, 310.93, 324.00
};

// Thickness in cm
const std::vector<const double> W3T1_heigths{
	2.318, 2.746, 3.174, 3.602, 4.030, 4.458, 4.886, 5.314, 5.743, 6.174,
	6.621, 7.020, 7.432, NULL // NULL >= 5 cm brass
};

const char* W3T1_material = "carbon";

/*
# RMW3_TR1_th_lead
d : Ge / RMW3_TR1_Lead / Offset = Ge / RMW3_TR1_Offset deg */
// Lead Angles in deg
const std::vector<const double> W3T1_lead_angles{
	3.00, 90.86, 123.91, 150.89, 172.83, 191.84, 209.99, 227.35, 244.30, 261.37,
	283.82
};

// Lead Thickness in cm
const std::vector<const double> W3T1_lead_heigths{
	0.123, 0.111, 0.099, 0.086, 0.074, 0.062, 0.050, 0.038, 0.027, 0.015,
	0.000
};

/*###########	Wheel 3 Track 2: 24 steps ####################
d : Ge / RMW3_TR2_Offset = 0.0 deg #13.17 deg + Ge / RM_SW / MiddleTrackOffset
d : Ge / RMW3_TR2_Carbon_Brass / Offset = Ge / RMW3_TR2_Offset deg*/
// Angles in deg
const std::vector<const double> W3T2_angles{
	140.64, 171.76, 200.39, 220.19, 236.84, 250.46, 262.24, 272.21, 281.06, 288.82,
	295.63, 301.79, 307.18, 311.98, 316.59, 321.10, 325.60, 330.08, 334.49, 338.96,
	343.40, 347.85, 352.39, 357.09
};

// Thickness in cm
const std::vector<const double> W3T2_heigths{
	0.517,  1.041,  1.562, 2.085, 2.601, 3.118, 3.636, 4.151, 4.665, 5.175,
	5.687,  6.197,  6.708, 7.184, 7.628, 8.073, 8.517, 8.961, 9.405, 9.850,
	10.294, 10.739, 11.184, NULL // NULL >= 5 cm brass
};

const char* W3T2_material = "carbon";

/*
# RMW3_TR2_th_lead
d : Ge / RMW3_TR2_Lead / Offset = Ge / RMW3_TR2_Offset deg */
// Lead Angles in deg
const std::vector<const double> W3T2_lead_angles{
	26.09, 140.64, 171.76, 200.39, 220.19, 236.84, 250.46, 262.24, 272.21, 281.06,
	288.82, 295.63, 301.79, 307.18, 311.98
};

// Lead Thickness in cm
const std::vector<const double> W3T2_lead_heigths{
	0.274, 0.253, 0.231, 0.209, 0.187, 0.167, 0.146, 0.125, 0.105, 0.085,
	0.066, 0.047, 0.029, 0.009, 0.000
};

/*###########	Wheel 3 Track 3: 25 steps ####################
#this is track for wobbling mode!
d:Ge / RMW3_TR3_Offset = 0.0 deg # 19.16 deg + Ge / RM_SW / OuterTrackOffset
d : Ge / RMW3_TR3_Carbon / Offset = Ge / RMW3_TR3_Offset deg*/
// Angles in deg
const std::vector<const double> W3T3_angles{
	9.38,  26.18,  40.48,  54.78,  69.08,  83.38,
	97.68, 111.98, 126.28, 140.58, 154.88, 169.18, 183.48, 197.78, 212.08, 226.38,
	240.68, 254.98, 269.28, 283.58, 297.88, 312.18, 326.48, 340.78, 355.08
};

// Thickness in cm
const std::vector<const double> W3T3_heigths{
	12.680, 11.374, 8.758, 6.144, 3.531,  0.000,
	0.662, 1.455, 2.247, 3.039, 3.833, 4.625, 5.417, 6.210, 7.002, 7.794,
	8.587, 9.379, 10.171, 10.964, 11.757, 12.549, 13.341, 14.134, 13.987 // NULL >= 5 cm brass
};

const char* W3T3_material = "carbon2";

/*
d : Ge / RMW3_TR3_Aluminum / Offset = Ge / RMW3_TR3_Offset deg */
// Lead Angles in deg
const std::vector<const double> W3T3_aluminium_angles{
	9.38, 26.18, 40.48, 54.78, 69.08, 83.38, 355.08
};

// Lead Thickness in cm
const std::vector<const double> W3T3_aluminium_heigths{
	0.000, 0.001, 0.001, 0.001, 0.001, 0.000, 0.001
};


// FUNCTIONS
std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
get_range_modulator(size_t ID){
	switch (ID){
	case 1:
		return std::make_tuple(W1T1_angles, W1T1_heigths, W1T1_material, W1T1_lead_angles, W1T1_lead_heigths);
	case 2:
		return std::make_tuple(W1T2_angles, W1T2_heigths, W1T1_material, W1T2_lead_angles, W1T2_lead_heigths);
	case 3:
		return std::make_tuple(W1T3_angles, W1T3_heigths, W1T1_material, W1T3_lead_angles, W1T3_lead_heigths);
	case 4:
		return std::make_tuple(W2T1_angles, W2T1_heigths, W2T1_material, W2T1_lead_angles, W2T2_lead_heigths);
	case 5:
		return std::make_tuple(W2T2_angles, W2T2_heigths, W2T2_material, W2T2_lead_angles, W2T3_lead_heigths);
	case 6:
		return std::make_tuple(W2T3_angles, W2T3_heigths, W2T3_material, W2T3_lead_angles, W2T1_lead_heigths);
	case 7:
		return std::make_tuple(W3T1_angles, W3T1_heigths, W3T1_material, W3T1_lead_angles, W3T1_lead_heigths);
	case 8:
		return std::make_tuple(W3T2_angles, W3T2_heigths, W3T2_material, W3T2_lead_angles, W3T2_lead_heigths);
	case 9:
		return std::make_tuple(W3T3_angles, W3T3_heigths, W3T3_material, W3T3_aluminium_angles, W3T3_aluminium_heigths);
	}
	std::cout << "\a" << "Range modulator ID undefined!! ID: " << ID << std::endl;
	std::cout << "returning RM_8 as default" << std::endl;
	return std::make_tuple(W3T2_angles, W3T2_heigths, W3T2_material, W3T2_lead_angles, W3T2_lead_heigths);
}

// bethe for water and water equivalent
double bethe(const double T) {
	// const double pi = 3.1415926535897932384626433832795028841971693993751;
	// const double N_A = 6.022140857e23; // mol^-1
	// const double mass_to_energy = 931.4940954; //MeV
	// const double Z = 0.555086707; // Z/A of water
	// const double I = 67.2e-6; // Exitation potential of water [MeV]

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                  pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     Z/A of water
	const double second = 6.022140857e23 * 0.555086707 / b_2; // / M_m=1;

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 67.2e-6) - b_2;

	return first * second * third * 0.1; // MeV/mm
}

double bethe_carbon(const double T) {
	// "Carbon, 6": Material(A=12.0, Z=6.0, mass=12.011, I=81.0, density=2.0,

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                         pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     rho[g/cm3], Z   &  M_m of Carbon
	const double second = 6.022140857e23 * 2.0 *       6.0 / (12.011 * b_2); // / M_u = 1 g/cm3

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 81.0e-6) - b_2;

	return first * second * third; // MeV/cm
}

double bethe_carbon2(const double T) {
	// "Carbon, 6": Material(A=12.0, Z=6.0, mass=12.011, I=81.0, density=2.0,

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                         pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     rho[g/cm3], Z   &  M_m of Carbon
	const double second = 6.022140857e23 * 1.816 *     6.0 / (12.011 * b_2); // / M_u = 1 g/cm3

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 78.0e-6) - b_2;

	return first * second * third; // MeV/cm
}

double bethe_aluminium(const double T) {
	// "Aluminum, 13": Material(A=27.0, Z=13.0, mass=26.9815385, I=166.0, density=2.6989

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                         pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     rho[g/cm3], Z   &   M_m of Aluminium
	const double second = 6.022140857e23 * 2.6989 *    13.0 / (26.9815385 * b_2); // / M_u = 1 g/cm3

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 166.0e-6) - b_2;

	return first * second * third; // MeV/cm
}

double bethe_lead(const double T) {
	// "Lead, 82": Material(A=207.0, Z=82.0, mass=207.2, I=823.0, density=11.35

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                         pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     rho[g/cm3], Z   &   M_m of Lead
	const double second = 6.022140857e23 * 11.35 *     81.0 / (207.2 * b_2); // / M_u = 1 g/cm3

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 823.0e-6) - b_2;

	return first * second * third; // MeV/cm
}

double bethe_lexan(const double T) {
	// Polycarbonate

	// classical electron radius = e^2/(4*pi*eps * m_e)
	const double r = 2.8179403227e-15 * 100.0; // convert to cm
	const double m_e_x2 = 0.5109989461 * 2.0; // MeV times two to shave off some instructions

	// Convert beam mass from AMU to MeV
	const double M_b = 1.0072766 * 931.4940954;

	// Beta^2
	const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b)); //<- reduce this if possible if bethe is too time consuming
	// Beta^2 * Gamma^2
	const double b_2_x_g_2 = b_2 / (1.0 - b_2);

	//                         pi
	const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

	//                    N_A [mol^-1]     rho[g/cm3],            Z                       &   M_m of Polycarbonate
	const double second = 6.022140857e23 * 1.2 * (0.055491 + 6.0 * 0.755751 + 8.0 * 0.188758) / ((1.008*0.055491 + 12.011 * 0.755751 + 15.999 * 0.188758) * b_2); // / M_u = 1 g/cm3

	//                                            I [MeV]
	const double third = log(m_e_x2 * b_2_x_g_2 / 73.1e-6) - b_2;

	return first * second * third; // MeV/cm
}

static const double cdf_xs[] = { 0.0, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 2.00E-05, 3.00E-05, 3.00E-05, 5.00E-05, 5.00E-05, 6.00E-05, 7.00E-05,
1.00E-04, 0.00011, 0.00012, 0.00012, 0.00017, 0.00017, 0.00019, 2.00E-04, 0.00021, 0.00021, 0.00021, 0.00023, 0.00025, 0.00028, 0.00031, 0.00036, 4.00E-04, 0.00043,
0.00049, 0.00056, 6.00E-04, 0.00064, 0.00071, 0.00079, 0.00087, 0.00093, 0.001, 0.00108, 0.00121, 0.00132, 0.00147, 0.00161, 0.00174, 0.0019, 0.00208, 0.00222, 0.00234,
0.00256, 0.00277, 0.00294, 0.00318, 0.00339, 0.0036, 0.00386, 0.00419, 0.00447, 0.00469, 0.0051, 0.00538, 0.00572, 0.00614, 0.00655, 0.00701, 0.00734, 0.00776, 0.0083,
0.00886, 0.00953, 0.01017, 0.01082, 0.01149, 0.01226, 0.013, 0.01371, 0.01469, 0.0157, 0.01688, 0.01778, 0.01908, 0.02032, 0.02124, 0.02242, 0.02354, 0.02493, 0.02634,
0.02781, 0.02927, 0.03087, 0.03264, 0.03437, 0.03611, 0.03799, 0.04022, 0.04233, 0.04457, 0.04674, 0.04922, 0.05162, 0.05416, 0.05665, 0.05927, 0.06213, 0.06501, 0.0683,
0.07143, 0.07445, 0.07792, 0.0819, 0.08562, 0.08946, 0.09345, 0.09728, 0.101, 0.10525, 0.1095, 0.11401, 0.11851, 0.12344, 0.12861, 0.13349, 0.13879, 0.14444, 0.14966,
0.15552, 0.16136, 0.16793, 0.17379, 0.17996, 0.18626, 0.19268, 0.19956, 0.20597, 0.21333, 0.22052, 0.22712, 0.23422, 0.2417, 0.24874, 0.25595, 0.26383, 0.27163, 0.27921,
0.28714, 0.29521, 0.30349, 0.31165, 0.31974, 0.32844, 0.33683, 0.3454, 0.35449, 0.36317, 0.37233, 0.38124, 0.38985, 0.39872, 0.40802, 0.41726, 0.4266, 0.43564, 0.44476,
0.45438, 0.46376, 0.4729, 0.48204, 0.4915, 0.50134, 0.51052, 0.52039, 0.53014, 0.53899, 0.54843, 0.5579, 0.56746, 0.5767, 0.58548, 0.59521, 0.60492, 0.61417, 0.62274,
0.63167, 0.64041, 0.65011, 0.65854, 0.66755, 0.67579, 0.68426, 0.69258, 0.70105, 0.70923, 0.71715, 0.7244, 0.73231, 0.73993, 0.74732, 0.75476, 0.76239, 0.76944, 0.77586,
0.78261, 0.78904, 0.79591, 0.80282, 0.80943, 0.81584, 0.82177, 0.82799, 0.83405, 0.84022, 0.84617, 0.85182, 0.85714, 0.86243, 0.8675, 0.87258, 0.87747, 0.88225, 0.88698,
0.89136, 0.8957, 0.89957, 0.90388, 0.90797, 0.91174, 0.91536, 0.9189, 0.92206, 0.92523, 0.92878, 0.9319, 0.93495, 0.93811, 0.9409, 0.94344, 0.9459, 0.94824, 0.95064,
0.95319, 0.95543, 0.95762, 0.95988, 0.96184, 0.96368, 0.96543, 0.96726, 0.969, 0.97052, 0.97191, 0.97341, 0.9748, 0.97617, 0.97757, 0.9789, 0.97996, 0.98111, 0.98227,
0.98333, 0.98428, 0.98531, 0.98616, 0.98694, 0.98768, 0.98837, 0.98914, 0.9898, 0.99045, 0.99095, 0.99161, 0.9921, 0.99256, 0.99296, 0.99336, 0.99382, 0.99424, 0.99452,
0.99492, 0.99532, 0.99559, 0.99594, 0.99624, 0.99645, 0.99671, 0.99697, 0.99715, 0.99741, 0.99763, 0.99779, 0.99796, 0.99807, 0.99825, 0.99836, 0.99849, 0.9986, 0.99873,
0.99877, 0.99887, 0.99897, 0.99902, 0.99909, 0.9992, 0.99925, 0.99929, 0.99935, 0.9994, 0.99945, 0.9995, 0.99957, 0.9996, 0.99963, 0.99966, 0.99968, 0.99971, 0.99975,
0.99977, 0.99978, 0.99981, 0.99983, 0.99984, 0.99986, 0.99987, 0.99987, 0.99987, 0.99989, 0.9999, 0.9999, 0.99992, 0.99993, 0.99993, 0.99993, 0.99993, 0.99993, 0.99993,
0.99993, 0.99994, 0.99995, 0.99995, 0.99996, 0.99997, 0.99997, 0.99997, 0.99997, 0.99998, 0.99999, 0.99999, 1.0 };

static const double cdf_val[] = { -1.7, -1.69, -1.68, -1.67, -1.66, -1.65, -1.64, -1.63, -1.62, -1.61, -1.6, -1.59, -1.58, -1.57, -1.56, -1.55, -1.54, -1.53, -1.52, -1.51, -1.5, -1.49, -1.48,
-1.47, -1.46, -1.45, -1.44, -1.43, -1.42, -1.41, -1.4, -1.39, -1.38, -1.37, -1.36, -1.35, -1.34, -1.33, -1.32, -1.31, -1.3, -1.29, -1.28, -1.27, -1.26, -1.25, -1.24, -1.23, -1.22, -1.21, -1.2, -1.19,
-1.18, -1.17, -1.16, -1.15, -1.14, -1.13, -1.12, -1.11, -1.1, -1.09, -1.08, -1.07, -1.06, -1.05, -1.04, -1.03, -1.02, -1.01, -1, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9,
-0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62,
-0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4, -0.39, -0.38, -0.37, -0.36, -0.35, -0.34,
-0.33, -0.32, -0.31, -0.3, -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23, -0.22, -0.21, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.1, -0.09, -0.08, -0.07, -0.06,
-0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27,
0.28, 0.29, 0.3, .31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61,
0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95,
0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29,
1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63,
1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78 };