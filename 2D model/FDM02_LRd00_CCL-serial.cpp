/* The Luo-Rudy Dynamic (LRd) Model of the Mammalian Ventricular Myocyte */
/* Gregory Faber */
/* This code requires a C++ compiler */
/* Detailed list of equations and model description are provided in */

/* Circ Res 1991;68:1501-1526 */
/* Circ Res 1994;74:1071-1096 */
/* Circ Res 1994;74:1097-1113 */
/* Circ Res 1995;77:140-152 */
/* Biophys J 1995;68:949-964 */
/* Cardiovasc Res 1997;35:256-272 */
/* Circulation 1999;99:2466-2474 */
/* Cardiovas Res 1999;42:530-542 */
/* Nature 1999;400:566-569 */
/* Circulation 2000;101:1192-1198 */
/* Biophy J 2000;78:2392-2404 */


/* IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USER SHOULD PACE THE MODEL UNTIL STEADY-STATE IS REACHED.
THIS REQUIRES APPROXIMATELY 5-20 MINUTES OF PACING DEPENDING ON THE RATE.*/

#include <iostream> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 
#include <stdlib.h> 
#include <stdio.h> 
#include <sys/timeb.h>

#define bcl 1000 // Basic Cycle Length (ms) 
#define beats 30 // Number of Beats 

//*******FDM parameters for LRd00 *******
int const nx = 5, ny = 5;//grid numbers
double dx = 0.015, dy = 0.015;//space step, 3cm*3cm
double D = 0.001;//D: diffusion coefficient cm^2/ms
double dt_max = 0.01; //the time step for the PDEs in opreator splitting method (ms)

/* Numerical variables */
int sti_flag = -1;// turn to 1 if stimulus is applied
const double Voffset = 0.1;//mv, a given membrane potential offset
double D1V[nx + 2][nx + 2] = { 0 }, D1V_old[nx + 2][nx + 2] = { 0 }, D2V[nx + 2][nx + 2] = { 0 };
double D2V_nai, D2V_ki, D2V_cai, D2V_nsr, D2V_djsr;
double v_old[nx + 2][nx + 2] = { 0 };
double dt_rem; //the left time of the current cycle
int b_n = 0;// beat count
const double dt_univ = dt_max;// universal maximum time-step size (ms) 
const double dt_min = 0.001;//ms
double epsilon = 1e-10;

double dV2[nx + 2][nx + 2]; // second order derivatives of Voltage (mv)

//performance compared
double Vmax=0, V_left = 0, V_right = 0, left_peak, right_peak, conduction_t = 0;
double APD90; // Time of 90% Repolarization 
double Vold, v_onset;
double dvdtmax, APD90_start, APD90_end;

//*******FDM parameters for LRd00 *******
int cutcount = 40 / dt_max;
int stim_step = (int)(0.5 / dt_max + 0.6); //Time period during which stimulus is applied 

/*****************************************/

/* List of variables and parameters (this code uses all global variables) */

/* Creation of Data File */
FILE *ap;
FILE *fmaxs;
FILE *fpara;
FILE *CPU_time_CCL;

/* Cell Geometry */
const double l = 0.01; // Length of the cell (cm) 
const double a = 0.0011; // Radius of the cell (cm) 
const double pi = 3.141592; // Pi 
double vcell; // Cell volume (uL) 
double ageo; // Geometric membrane area (cm^2) 
double acap; // Capacitive membrane area (cm^2) 
double vmyo; // Myoplasm volume (uL) 
double vmito; // Mitochondria volume (uL) 
double vsr; // SR volume (uL) 
double vnsr; // NSR volume (uL) 
double vjsr; // JSR volume (uL) 
double vcleft; // Cleft volume (uL) 

/* Voltage */
double v[nx + 2][nx + 2]; // Membrane voltage (mV) 
double vnew[nx + 2][nx + 2]; // New Voltage (mV) 
double dvdt[nx + 2][nx + 2]; // Change in Voltage / Change in Time (mV/ms) 
double dvdtnew[nx + 2][nx + 2]; // New dv/dt (mV/ms) 
double flag[nx + 1][nx + 1]; // Flag condition to test for dvdtmax

/* Time Step */
double dt[nx + 1][nx + 1]; // Time step (ms) 
double t=0; // Time (ms) 
int steps; // Number of Steps 
int ncount=0; // Loop Control Variable

/* Action Potential Duration and Max. Info */
//double vmax[nx + 1][nx + 1][beats]; // Max. Voltage (mV) 
//double dvdtmax[nx + 1][nx + 1][beats]; // Max. dv/dt (mV/ms) 
//double apd[nx + 1][nx + 1][beats]; // Action Potential Duration 
//double toneapd[nx + 1][nx + 1][beats]; // Time of dv/dt Max. 
//double ttwoapd[nx + 1][nx + 1][beats]; // Time of 90% Repolarization 
double rmbp[nx + 1][nx + 1]; // Resting Membrane Potential 
double nair[nx + 1][nx + 1]; // Intracellular Na At Rest 
double cair[nx + 1][nx + 1]; // Intracellular Ca At Rest 
double kir[nx + 1][nx + 1]; // Intracellular K At Rest 
double caimax; // Peak Intracellular Ca 
int b_i; // Stimulation Counter

/* Total Current and Stimulus */
double stim = 0; // Stimulus (uA/cm^2) 
const double amp=-80.0; // const stimulus amplitude (uA/cm^2) 
double tstim; // Time Stimulus is Applied (ms) 
double stimtime; // Time period during which stimulus is applied (ms) 
double it[nx + 1][nx + 1]; // Total current (uA/cm^2) 
double duration; // duration of the stimulus (ms) 

/* Terms for Solution of Conductance and Reversal Potential */
const double Rgas = 8314; // Universal Gas Constant (J/kmol*K) 
const double frdy = 96485; // Faraday's Constant (C/mol) 
const double temp = 310; // Temperature (K) 

/* Ion Valences */
const double zna = 1; // Na valence 
const double zk = 1; // K valence 
const double zca = 2; // Ca valence 

/* Ion Concentrations */
double nai[nx + 1][nx + 1]; // Intracellular Na Concentration (mM) 
double nao; // Extracellular Na Concentration (mM) 
double nabm; // Bulk Medium Na Concentration (mM) 
double dnao; // Change in Cleft Na Concentration (mM) 
double ki[nx + 1][nx + 1]; // Intracellular K Concentration (mM) 
double ko; // Extracellular K Concentration (mM) 
double kbm; // Bulk Medium K Concentration (mM) 
double dko; // Change in Cleft K Concentration (mM) 
double cai[nx + 1][nx + 1]; // Intracellular Ca Concentration (mM) 
double cao; // Extracellular Ca Concentration (mM) 
double cabm; // Bulk Medium Ca Concentration (mM) 
double dcao; // Change in Cleft Ca Concentration (mM) 
double cmdn[nx + 1][nx + 1]; // Calmodulin Buffered Ca Concentration (mM) 
double trpn[nx + 1][nx + 1]; // Troponin Buffered Ca Concentration (mM) 
double nsr[nx + 1][nx + 1]; // NSR Ca Concentration (mM) 
double jsr[nx + 1][nx + 1]; // JSR Ca Concentration (mM) 
double csqn[nx + 1][nx + 1]; // Calsequestrin Buffered Ca Concentration (mM) 
const double taudiff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

/* Myoplasmic Na Ion Concentration Changes */
double naiont[nx + 1][nx + 1]; // Total Na Ion Flow (uA/uF) 
double dnai; // Change in Intracellular Na Concentration (mM) 

/* Myoplasmic K Ion Concentration Changes */
double kiont[nx + 1][nx + 1]; // Total K Ion Flow (uA/uF) 
double dki; // Change in Intracellular K Concentration (mM) 

/* NSR Ca Ion Concentration Changes */
double dnsr; // Change in [Ca] in the NSR (mM) 
double iup[nx + 1][nx + 1]; // Ca uptake from myo. to NSR (mM/ms) 
double ileak[nx + 1][nx + 1]; // Ca leakage from NSR to myo. (mM/ms) 
double kleak; // Rate constant of Ca leakage from NSR (ms^-1) 
const double kmup = 0.00092; // Half-saturation concentration of iup (mM) 
const double iupbar = 0.00875; // Max. current through iup channel (mM/ms) 
const double nsrbar = 15; // Max. [Ca] in NSR (mM)

/* JSR Ca Ion Concentration Changes */
double djsr; // Change in [Ca] in the JSR (mM) 
const double tauon = 0.5; // Time constant of activation of Ca release from JSR (ms) 
const double tauoff = 0.5; // Time constant of deactivation of Ca release from JSR (ms) 
double tcicr[nx + 1][nx + 1]; // t=0 at time of CICR (ms) 
double irelcicr[nx + 1][nx + 1]; // Ca release from JSR to myo. due to CICR (mM/ms) 
const double csqnth = 8.75; // Threshold for release of Ca from CSQN due to JSR overload (mM) 
const double gmaxrel = 150; // Max. rate constant of Ca release from JSR due to overload (ms^-1) 
double grelbarjsrol[nx + 1][nx + 1]; // Rate constant of Ca release from JSR due to overload (ms^-1) 
double greljsrol; // Rate constant of Ca release from JSR due to CICR (ms^-1) 
double tjsrol[nx + 1][nx + 1]; // t=0 at time of JSR overload (ms) 
double ireljsrol[nx + 1][nx + 1]; // Ca release from JSR to myo. due to JSR overload (mM/ms) 
const double csqnbar = 10; // Max. [Ca] buffered in CSQN (mM) 
const double kmcsqn = 0.8; // Equilibrium constant of buffering for CSQN (mM) 
double bjsr; // b Variable for analytical computation of [Ca] in JSR (mM) 
double cjsr; // c Variable for analytical computation of [Ca] in JSR (mM) 
double on; // Time constant of activation of Ca release from JSR (ms) 
double off; // Time constant of deactivation of Ca release from JSR (ms) 
double magrel; // Magnitude of Ca release 
double dcaiont[nx + 1][nx + 1]; // Rate of change of Ca entry 
double dcaiontnew[nx + 1][nx + 1]; // New rate of change of Ca entry 
double caiontold[nx + 1][nx + 1]; // Old rate of change of Ca entry 

/* Translocation of Ca Ions from NSR to JSR */
double itr[nx + 1][nx + 1]; // Translocation current of Ca ions from NSR to JSR (mM/ms) 
const double tautr = 180; // Time constant of Ca transfer from NSR to JSR (ms) 

/* Myoplasmic Ca Ion Concentration Changes */
double caiont[nx + 1][nx + 1]; // Total Ca Ion Flow (uA/uF) 
double dcai; // Change in myoplasmic Ca concentration (mM) 
double catotal[nx + 1][nx + 1]; // Total myoplasmic Ca concentration (mM) 
double bmyo; // b Variable for analytical computation of [Ca] in myoplasm (mM) 
double cmyo; // c Variable for analytical computation of [Ca] in myoplasm (mM) 
double dmyo; // d Variable for analytical computation of [Ca] in myoplasm (mM) 
double gpig; // Tribute to all the guinea pigs killed for the advancement of knowledge 
const double cmdnbar = 0.050; // Max. [Ca] buffered in CMDN (mM) 
const double trpnbar = 0.070; // Max. [Ca] buffered in TRPN (mM) 
const double kmcmdn = 0.00238; // Equilibrium constant of buffering for CMDN (mM) 
const double kmtrpn = 0.0005; // Equilibrium constant of buffering for TRPN (mM) 

/* Fast Sodium Current (time dependant) */
double ina[nx + 1][nx + 1]; // Fast Na Current (uA/uF) 
double gna; // Max. Conductance of the Na Channel (mS/uF) 
double ena[nx + 1][nx + 1]; // Reversal Potential of Na (mV) 
double am; // Na alpha-m rate constant (ms^-1) 
double bm; // Na beta-m rate constant (ms^-1) 
double ah; // Na alpha-h rate constant (ms^-1) 
double bh; // Na beta-h rate constant (ms^-1) 
double aj; // Na alpha-j rate constant (ms^-1) 
double bj; // Na beta-j rate constant (ms^-1) 
double mtau; // Na activation 
double htau; // Na inactivation 
double jtau; // Na inactivation 
double mss; // Na activation 
double hss; // Na inactivation 
double jss; // Na inactivation 
double m[nx + 1][nx + 1]; // Na activation 
double h[nx + 1][nx + 1]; // Na inactivation 
double jj[nx + 1][nx + 1]; // Na inactivation 

/* Current through L-type Ca Channel */
double ilca[nx + 1][nx + 1]; // Ca current through L-type Ca channel (uA/uF) 
double ilcana[nx + 1][nx + 1]; // Na current through L-type Ca channel (uA/uF) 
double ilcak[nx + 1][nx + 1]; // K current through L-type Ca channel (uA/uF) 
double ilcatot[nx + 1][nx + 1]; // Total current through the L-type Ca channel (uA/uF) 
double ibarca; // Max. Ca current through Ca channel (uA/uF) 
double ibarna; // Max. Na current through Ca channel (uA/uF) 
double ibark; // Max. K current through Ca channel (uA/uF) 
double d[nx + 1][nx + 1]; // Voltage dependant activation gate 
double dss; // Steady-state value of activation gate d 
double taud; // Time constant of gate d (ms^-1) 
double f[nx + 1][nx + 1]; // Voltage dependant inactivation gate 
double fss; // Steady-state value of inactivation gate f 
double tauf; // Time constant of gate f (ms^-1) 
double fca; // Ca dependant inactivation gate 
const double kmca = 0.0006; // Half-saturation concentration of Ca channel (mM) 
const double pca = 0.00054; // Permeability of membrane to Ca (cm/s) 
const double gacai = 1; // Activity coefficient of Ca 
const double gacao = 0.341; // Activity coefficient of Ca 
const double pna = 0.000000675; // Permeability of membrane to Na (cm/s) 
const double ganai = 0.75; // Activity coefficient of Na 
const double ganao = 0.75; // Activity coefficient of Na 
const double pk = 0.000000193; // Permeability of membrane to K (cm/s) 
const double gaki = 0.75; // Activity coefficient of K 
const double gako = 0.75; // Activity coefficient of K 

/* Current through T-type Ca Channel */
double icat[nx + 1][nx + 1]; // Ca current through T-type Ca channel (uA/uF) 
double gcat; // Max. Conductance of the T-type Ca channel (mS/uF) 
double eca[nx + 1][nx + 1]; // Reversal Potential of the T-type Ca channel (mV) 
double b[nx + 1][nx + 1]; // Voltage dependant activation gate 
double bss; // Steady-state value of activation gate b 
double taub; // Time constant of gate b (ms^-1) 
double g[nx + 1][nx + 1]; // Voltage dependant inactivation gate 
double gss; // Steady-state value of inactivation gate g 
double taug; // Time constant of gate g (ms^-1) 

/* Potassium Current (time-dependent) in LRd94*/
//double x;
//double gk;
//double ek;
//double xi;
//double ax;
//double bx;
//double taux;
//double xss;
//double ik;
//const double prnak = 0.01833;

/* Rapidly Activating Potassium Current */
double ikr[nx + 1][nx + 1]; // Rapidly Activating K Current (uA/uF) 
double gkr; // Channel Conductance of Rapidly Activating K Current (mS/uF) 
double ekr[nx + 1][nx + 1]; // Reversal Potential of Rapidly Activating K Current (mV) 
double xr[nx + 1][nx + 1]; // Rapidly Activating K time-dependant activation 
double xrss; // Steady-state value of inactivation gate xr 
double tauxr; // Time constant of gate xr (ms^-1) 
double r; // K time-independent inactivation

/* Slowly Activating Potassium Current */
double iks[nx + 1][nx + 1]; // Slowly Activating K Current (uA/uF) 
double gks; // Channel Conductance of Slowly Activating K Current (mS/uF) 
double eks[nx + 1][nx + 1]; // Reversal Potential of Slowly Activating K Current (mV) 
double xs1[nx + 1][nx + 1]; // Slowly Activating K time-dependant activation 
double xs1ss; // Steady-state value of inactivation gate xs1 
double tauxs1; // Time constant of gate xs1 (ms^-1) 
double xs2[nx + 1][nx + 1]; // Slowly Activating K time-dependant activation 
double xs2ss; // Steady-state value of inactivation gate xs2 
double tauxs2; // Time constant of gate xs2 (ms^-1) 
const double prnak = 0.01833; // Na/K Permeability Ratio

/* Potassium Current (time-independent) */
double iki[nx + 1][nx + 1]; // Time-independent K current (uA/uF) 
double gki; // Channel Conductance of Time Independant K Current (mS/uF) 
double eki[nx + 1][nx + 1]; // Reversal Potential of Time Independant K Current (mV) 
double aki; // K alpha-ki rate constant (ms^-1) 
double bki; // K beta-ki rate constant (ms^-1) 
double kin; // K inactivation 

/* Plateau Potassium Current */
double ikp[nx + 1][nx + 1]; // Plateau K current (uA/uF) 
double gkp; // Channel Conductance of Plateau K Current (mS/uF) 
double ekp; // Reversal Potential of Plateau K Current (mV) 
double kp; // K plateau factor 

/* Na-Activated K Channel */
double ikna[nx + 1][nx + 1]; // Na activated K channel 
double pona; // Open probability dependant on Nai 
double pov; // Open probability dependant on Voltage 
double ekna[nx + 1][nx + 1]; // Reversal potential 
const double gkna = 0.12848; // Maximum conductance (mS/uF) 
const double nkna = 2.8; // Hill coefficient for Na dependance 
const double kdkna = 66; // Dissociation constant for Na dependance(mM)

/* ATP-Sensitive K Channel */
double ikatp[nx + 1][nx + 1]; // ATP-sensitive K current (uA/uF) 
double ekatp[nx + 1][nx + 1]; // K reversal potential (mV) 
double gkbaratp; // Conductance of the ATP-sensitive K channel (mS/uF) 
double gkatp; // Maximum conductance of the ATP-sensitive K channel (mS/uF) 
double patp; // Percentage availability of open channels 
const double natp = 0.24; // K dependence of ATP-sensitive K current 
const double nicholsarea = 0.00005; // Nichol's area (cm^2) 
const double atpi = 3; // Intracellular ATP concentraion (mM) 
const double hatp = 2; // Hill coefficient 
const double katp = 0.250; // Half-maximal saturation point of ATP-sensitive K current (mM)

/* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */
double ito[nx + 1][nx + 1]; // Transient outward current 
double gitodv; // Maximum conductance of Ito 
double ekdv[nx + 1][nx + 1]; // Reversal Potential of Ito 
double rvdv; // Time independent voltage dependence of Ito 
double zdv[nx + 1][nx + 1]; // Ito activation 
double azdv; // Ito alpha-z rate constant 
double bzdv; // Ito beta-z rate constant 
double tauzdv; // Time constant of z gate 
double zssdv; // Steady-state value of z gate 
double ydv[nx + 1][nx + 1]; // Ito inactivation 
double aydv; // Ito alpha-y rate constant 
double bydv; // Ito beta-y rate constant 
double tauydv; // Time constant of y gate 
double yssdv; // Steady-state value of y gate

/* Sodium-Calcium Exchanger V-S */
double inaca[nx + 1][nx + 1]; // NaCa exchanger current (uA/uF) 
const double c1 = 0.00025; // Scaling factor for inaca (uA/uF) 
const double c2 = 0.0001; // Half-saturation concentration of NaCa exhanger (mM) 
const double gammas = 0.15; // Position of energy barrier controlling voltage dependence of inaca 

///* Sodium-Calcium Exchanger V-S in LRd94 be replaced by LRd00*/  
//double inaca; // NaCa exchanger current (uA/uF) 
//const double knaca = 2000; // Scaling factor for inaca (uA/uF) 
//const double ksat = 0.1; // Half-saturation concentration of NaCa exhanger (mM) 
//const double gamma = 0.35; //eta=0.35 in the paper// Position of energy barrier controlling voltage dependence of inaca 
//const double kmna = 87.5; // Half-saturation concentration of Na channel (mM)
//const double kmca2 = 1.38; // Half-saturation concentration of Ca channel (mM)

/* Sodium-Potassium Pump */
double inak[nx + 1][nx + 1]; // NaK pump current (uA/uF) 
double fnak; // Voltage-dependence parameter of inak 
double sigma; // [Na] o dependence factor of fnak 
/* In unpublished changes on rudy's website(ref: History of LRd Model Development), ibarnak is changed from 1.5(LRd94) to 2.*/
const double ibarnak = 2.25; // Max. current through Na-K pump (uA/uF) 
const double kmnai = 10; // Half-saturation concentration of NaK pump (mM) 
const double kmko = 1.5; // Half-saturation concentration of NaK pump (mM)

/* Nonspecific Ca-activated Current */
double insna[nx + 1][nx + 1]; // Non-specific Na current (uA/uF) 
double insk[nx + 1][nx + 1]; // Non-specific K current (uA/uF) 
double ibarnsna; // Max. Na current through NSCa channel (uA/uF) 
double ibarnsk; // Max. K current through NSCa channel (uA/uF) 
const double pnsca = 0.000000175; // Permeability of channel to Na and K (cm/s) 
const double kmnsca = 0.0012; // Half-saturation concentration of NSCa channel (mM) 

/* Sarcolemmal Ca Pump */
double ipca[nx + 1][nx + 1]; // Sarcolemmal Ca pump current (uA/uF) 
const double ibarpca = 1.15; // Max. Ca current through sarcolemmal Ca pump (uA/uF) 
const double kmpca = 0.0005; // Half-saturation concentration of sarcolemmal Ca pump (mM)

/* Ca Background Current */
double icab[nx + 1][nx + 1]; // Ca background current (uA/uF) 
double gcab; // Max. conductance of Ca background (mS/uF) 
double ecan[nx + 1][nx + 1]; // Nernst potential for Ca (mV) 

/* Na Background Current */
double inab[nx + 1][nx + 1]; // Na background current (uA/uF) 
double gnab; // Max. conductance of Na background (mS/uF) 
double enan; // Nernst potential for Na (mV) 

/* Ion Current Functions */
void comp_ina(int i, int j); // Calculates Fast Na Current 
void comp_ical(int i, int j); // Calculates Currents through L-Type Ca Channel 
void comp_icat(int i, int j); // Calculates Currents through T-Type Ca Channel 
//comp_ik();
void comp_ikr(int i, int j); // Calculates Rapidly Activating K Current 
void comp_iks(int i, int j); // Calculates Slowly Activating K Current 
void comp_iki(int i, int j); // Calculates Time-Independent K Current 
void comp_ikp(int i, int j); // Calculates Plateau K Current 
void comp_ikna(int i, int j); // Calculates Na-activated K Current 
void comp_ikatp(int i, int j); // Calculates ATP-Sensitive K Current 
void comp_ito(int i, int j); // Calculates Transient Outward Current 
void comp_inaca(int i, int j); // Calculates Na-Ca Exchanger Current 
void comp_inak(int i, int j); // Calculates Na-K Pump Current 
void comp_insca(int i, int j); // Calculates Non-Specific ca-Activated Current 
void comp_ipca(int i, int j); // Calculates Sarcolemmal Ca Pump Current 
void comp_icab(int i, int j); // Calculates Ca Background Current 
void comp_inab(int i, int j); // Calculates Na Background Current 
void comp_it(int i, int j); // Calculates Total Current 

/* Ion Concentration Functions */
void conc_nai(int i, int j, double dt); // Calculates new myoplasmic Na ion concentration 
void conc_ki(int i, int j, double dt); // Calculates new myoplasmic K ion concentration 
void conc_nsr(int i, int j, double dt); // Calculates new NSR Ca ion concentration 
void conc_jsr(int i, int j, double dt); // Calculates new JSR Ca ion concentration 
void calc_itr(int i, int j, double dt); // Calculates Translocation of Ca from NSR to JSR 
void conc_cai(int i, int j, double dt); // Calculates new myoplasmic Ca ion concentration 
void conc_cleft(int i, int j, double dt); // Calculates new cleft ion concentrations 

/* Additive funtions */
void comp_revpots(int i, int j);//Calculating new Reversal Potential
double comp_currents(int i, int j, double dt);//a package of current functions for updating new currents
void comp_concentrations(int i, int j, double dt);//a package for calaculating ion concentrations
double comp_CCL(int i, int j, double dt); // Calculating new time-step size
void Rush_Larsen(int i, int j, double dt);//Calculating new Gatings
void comp_voltage(int i, int j, double dt);//Calculating new voltage
void performance();

int main()
{
	/* Opening of Datafiles */
	ap = fopen("ap", "w");
	fpara = fopen("fpara", "w");
	fmaxs = fopen("fmaxs", "w");
	CPU_time_CCL = fopen("CPU_time_CCL.dat", "w");

	/* Cell Geometry */
	vcell = 1000 * pi*a*a*l; // 3.801e-5 uL 
	ageo = 2 * pi*a*a + 2 * pi*a*l; // 7.671e-5 cm^2 
	acap = ageo * 2; // 1.534e-4 cm^2 
	vmyo = vcell*0.68;
	vmito = vcell*0.26;
	vsr = vcell*0.06;
	vnsr = vcell*0.0552;
	vjsr = vcell*0.0048;
	vcleft = vcell*0.12 / 0.88;

	/* Time Loop Conditions */
	//t = 0.0; // Time (ms) 
	//dt = 0.01; // Time step (ms) 
	//steps = (bcl*beats) / dt; // Number of ms 
	//st = -80.0; // Stimulus 
	tstim = 10.0; // Time to begin stimulus 
	stimtime = 10.0; // Initial Condition for Stimulus 
	duration = 0.5;// duration of the stimulus (ms) 
	int i, j;
	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			v[i][j] = -88.654973; // Initial Voltage (mv)
		}
	}
	for (i = 1; i < nx + 1; i++){
		for (j = 1; j < ny + 1; j++){
			nai[i][j] = 12.236437; // Initial Intracellular Na (mM) 
			ki[i][j] = 136.89149; // Initial Intracellular K (mM) 
			cai[i][j] = 0.000079; // Initial Intracellular Ca (mM) 
			/* Initial Gate Conditions */
			m[i][j] = 0.000838;
			h[i][j] = 0.993336;
			jj[i][j] = 0.995484;
			d[i][j] = 0.000003;
			f[i][j] = 0.999745;
			xs1[i][j] = 0.004503;
			xs2[i][j] = 0.004503;
			xr[i][j] = 0.000129;
			b[i][j] = 0.000994;
			g[i][j] = 0.994041;
			zdv[i][j] = 0.0120892;
			ydv[i][j] = 0.999978;
			/* Initial Conditions */
			grelbarjsrol[i][j] = 0;
			tjsrol[i][j] = 1000;
			tcicr[i][j] = 1000;
			jsr[i][j] = 1.179991;
			nsr[i][j] = 1.179991;
			trpn[i][j] = 0.0143923;
			cmdn[i][j] = 0.00257849;
			csqn[i][j] = 6.97978;
			flag[i][j] = 0;
			dcaiont[i][j] = 0;
			dt[i][j] = dt_max;
		}
	}

	/* Beginning Ion Concentrations */
	//nai = 12.236437; // Initial Intracellular Na (mM) 
	nao = 140; // Initial Extracellular Na (mM) 
	nabm = 140; // Initial Bulk Medium Na (mM) 
	//ki = 136.89149; // Initial Intracellular K (mM) 
	ko = 4.5; // Initial Extracellular K (mM) 
	kbm = 4.5; // Initial Bulk Medium K (mM) 
	//cai = 0.000079; // Initial Intracellular Ca (mM) 
	cao = 1.8; // Initial Extracellular Ca (mM) 
	cabm = 1.8; // Initial Bulk Medium Ca (mM) 
	
	b_i = -1;

	int nstep = bcl / dt_max; // snapshot interval 10 ms to save data files
	int index = 0;// filename index
	char filename[100];

	struct timeb start, end;
	int diff;
	ftime(&start);

	/* Beginning of Time Loop */
	int b_n_temp = 1;
	for (ncount = 1; t <= bcl*beats; ncount++)
	{
		b_n = t / bcl+1;
		if (b_n>b_n_temp){
			APD90 = APD90_end - APD90_start;
			fprintf(fmaxs, "%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", b_n_temp, Vmax, dvdtmax, APD90,
				APD90_start, APD90_end, nair[nx / 2][1], kir[nx / 2][1], cair[nx / 2][1], caimax, rmbp[nx / 2][1]);
			Vmax = 0;//reset
			dvdtmax = 0;
			caimax = 0;
			b_n_temp = b_n;
		}

		for (i = 1; i < nx + 1; i++){
			//****no flux boundary conditions*****
			v[i][0] = v[i][1];
			v[i][ny + 1] = v[i][ny];
		}
		for (j = 1; j < ny + 1; j++){
			v[0][j] = v[1][j];
			v[nx + 1][j] = v[nx][j];
		}
		performance();
		//**** save data in file "ap"
		int fileflag = 0;

		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				if (ncount%nstep == 0){//
					if (fileflag == 0){
						sprintf(filename, "ap%d", index);
						ap = fopen(filename, "w");
						fileflag = 1;
						index++;
					}
					fprintf(ap, "%g\t", v[i][j]);
					//printf("%g\t", v[i][j]);
					if (j == ny){
						fprintf(ap, "\n");
						//printf("\n");
					}
				}
			}
		}
		if (fileflag == 1){
			fclose(ap);
		}
		//********** save data in file "ap"

		//*********** step 1 *******
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				dV2[i][j] = D*((v[i + 1][j] + v[i - 1][j] - 2 * v[i][j]) / (dx*dx) + (v[i][j + 1] + v[i][j - 1] - 2 * v[i][j]) / (dy*dy));
			}
		}
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				//Forward Euler
				vnew[i][j] = v[i][j] + dt_max / 2 * dV2[i][j];
				v[i][j] = vnew[i][j];
			}
		}
		//*********** step 1 *******

		//*********** step 2 *******
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){				
				it[i][j] = comp_currents(i, j, dt_max);
				D1V[i][j] = -it[i][j];
				//stim = 0;//reset the stimulus to 0
				//sti_flag = -1;//for CCL, reset the stimulus flag to -1
			}
		}
		//*****stimulation with a plane waves****
		//if (ncount >= 1000 && ncount <= 1000+stim_step) { //stimulus is hold with 0.5 ms, 0.02*15 = 0.3 ms			
		//if (t >= tstim + (b_n - 1)*bcl && t<=(tstim + (b_n - 1)*bcl + duration)){
		if (stimtime >= 0 && stimtime <= duration){
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j <= 5; j++){
					//stim = amp;
					//sti_flag = 1;//for CCL
					D1V[i][j] = D1V[i][j] + (-stim);
				}
			}
		}

		if (ncount == 1){// in order to get D1V[i][j], for computing D2V[i][j] in CCL(i, j, dt_max);
			/* The first time step*/
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					D1V_old[i][j] = D1V[i][j];
					comp_concentrations(i, j, dt_max);
					Rush_Larsen(i, j, dt_max);
					//v[i][j] = v[i][j] + dt_max * D1V[i][j];
					comp_voltage(i, j, dt_max);
				}
			}
		}
		else{
			double dt_sum;
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					/******************* adjust or correct time step---CCL method ***********************/
					dt_sum = 0;
					dt[i][j] = comp_CCL(i, j, dt_max);
					do
					{
						dt_sum = dt_sum + dt[i][j];
						if (dt_sum<dt_max){
							D1V_old[i][j] = D1V[i][j];
							comp_concentrations(i, j, dt[i][j]);
							Rush_Larsen(i, j, dt[i][j]);
							comp_voltage(i, j, dt[i][j]);
						}
						else{
							dt[i][j] = dt_max - (dt_sum - dt[i][j]);// here is a new dt  !!!
							D1V_old[i][j] = D1V[i][j];
							comp_concentrations(i, j, dt[i][j]);
							Rush_Larsen(i, j, dt[i][j]);
							comp_voltage(i, j, dt[i][j]);
							break;
						}
						it[i][j] = comp_currents(i, j, dt[i][j]);
						D1V[i][j] = -it[i][j] + (-stim);//here stim has obtained values
						dt[i][j] = comp_CCL(i, j, dt[i][j]);
					} while (true);
				}
			}
		}
		//*********** step 2 *******

		//*********** step 3 *******
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				dV2[i][j] = D*((v[i + 1][j] + v[i - 1][j] - 2 * v[i][j]) / (dx*dx) + (v[i][j + 1] + v[i][j - 1] - 2 * v[i][j]) / (dy*dy));
			}
		}

		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				//Forward Euler
				vnew[i][j] = v[i][j] + dt_max / 2 * dV2[i][j];
				v[i][j] = vnew[i][j];
			}
		}
		//*********** step 3 *******
		stimtime = stimtime + dt_max;
		t = t + dt_max;
		//if (b_i >= 0){
		//	if (vnew>vmax[b_i])
		//		vmax[b_i] = vnew;
		//	if (cai[i][j]>caimax[b_i])
		//		caimax[b_i] = cai[i][j];
		//	if (dvdtnew>dvdtmax[b_i]){
		//		dvdtmax[b_i] = dvdtnew;
		//		toneapd[b_i] = t;
		//	}
		//	if (vnew >= (vmax[b_i] - 0.9*(vmax[b_i] - rmbp[b_i])))
		//		ttwoapd[b_i] = t;
		//}
	}
	ftime(&end);
	diff = (int)(1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	conduction_t = (right_peak - left_peak)*0.001; //condution time from left side to right side
	//APD90 = APD90_end - APD90_start;
	//fprintf(fmaxs, "%d\t%g\t%g\t%g\t%g\t", diff / 1000, Vmax, dvdtmax, APD90, dx*nx / conduction_t);
	APD90 = APD90_end - APD90_start;
	fprintf(fmaxs, "%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", b_n_temp, Vmax, dvdtmax, APD90,
		APD90_start, APD90_end, nair[nx / 2][1], kir[nx / 2][1], cair[nx / 2][1], caimax, rmbp[nx / 2][1]);
	fclose(ap);
	fclose(fmaxs);

	//fprintf(fpara, "%.3f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t"
	//	"%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
	//	t, v, nai, ki, cai, jsr, nsr, nao, ko, cao, m, h, jj, d, f, xs1, xs2, xr, b, g, tcicr, flag);
	//for (i = 0; i<beats; i++)
	//{
	//	apd[i] = ttwoapd[i] - toneapd[i];
	//	fprintf(fmaxs, "%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i, vmax[i], dvdtmax[i], apd[i],
	//		toneapd[i], ttwoapd[i], nair[i], kir[i], cair[i], caimax[i], rmbp[i]);
	//}

	//fprintf(CPU_time_CCL, "%g\n", diff);
	//cout << "\a\a";

	return (1);
}

/********************************************************/

//performance 
void performance(){
	if (v[nx / 2][1] - V_left > 0){
		left_peak = t; // peak time at j=1
		V_left = v[nx / 2][1];
	}
	if (v[nx / 2][ny] - V_right > 0){
		right_peak = t; // peak time at j=ny
		V_right = v[nx / 2][ny];
	}
	if (v[nx / 2][1]>Vmax)
		Vmax = v[nx / 2][1];
	if (dvdtnew[nx / 2][1] > dvdtmax){
		dvdtmax = dvdtnew[nx / 2][1];
		APD90_start = t;
	}
	if (v[nx / 2][1] >= (Vmax - 0.9*(Vmax - (-88.654973))))
		APD90_end = t; //  Time of 90% Repolarization 
	if (cai[nx / 2][1]>caimax)
		caimax = cai[nx / 2][1];
}

////calaculate Reversal Potential
//void comp_revpots(int i, int j){
//	ena[i][j] = ((Rgas*temp) / frdy)*log(nao / nai[i][j]);
//	eca[i][j] = (Rgas*temp / (2 * frdy))*log(cao / cai[i][j]);
//	ekr[i][j] = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
//	eks[i][j] = ((Rgas*temp) / frdy)*log((ko + prnak*nao) / (ki[i][j] + prnak*nai[i][j]));
//	eki[i][j] = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
//	ekna[i][j] = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
//	ekatp[i][j] = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
//	//ekdv = ((Rgas*temp) / frdy)*log((ko + prnak*nao) / (ki[i][j] + prnak*nai[i][j])); //LRd94
//	ekdv[i][j] = ((Rgas*temp) / frdy)*log(ko / ki[i][j]); //LRd99
//	ecan[i][j] = ((Rgas*temp) / (2 * frdy))*log(cao / cai[i][j]);
//}

/*  calaculate new Gatings */
void Rush_Larsen(int i, int j, double dt){
	//Fast sodium current
	double am = 0.32*(v[i][j] + 47.13) / (1 - exp(-0.1*(v[i][j] + 47.13)));
	double bm = 0.08*exp(-v[i][j] / 11);
	double ah, bh, aj, bj;
	if (v[i][j] < -40){
		ah = 0.135*exp((80 + v[i][j]) / -6.8);
		bh = 3.56*exp(0.079*v[i][j]) + 310000 * exp(0.35*v[i][j]);
		aj = (-127140 * exp(0.2444*v[i][j]) - 0.00003474*exp(-0.04391*v[i][j]))*((v[i][j] + 37.78) / (1 + exp(0.311*(v[i][j] + 79.23))));
		bj = (0.1212*exp(-0.01052*v[i][j])) / (1 + exp(-0.1378*(v[i][j] + 40.14)));
	}
	else{
		ah = 0;
		bh = 1 / (0.13*(1 + exp((v[i][j] + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*v[i][j])) / (1 + exp(-0.1*(v[i][j] + 32)));
	}
	double mtau = 1 / (am + bm);
	double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

	double mss = am*mtau;
	double hss = ah*htau;
	double jss = aj*jtau;
	m[i][j] = mss - (mss - m[i][j])*exp(-dt / mtau);
	h[i][j] = hss - (hss - h[i][j])*exp(-dt / htau);
	jj[i][j] = jss - (jss - jj[i][j])*exp(-dt / jtau);

	//L-Type Ca Channel
	double dss = 1 / (1 + exp(-(v[i][j] + 10) / 6.24));
	double taud = dss*(1 - exp(-(v[i][j] + 10) / 6.24)) / (0.035*(v[i][j] + 10));
	double fss = (1 / (1 + exp((v[i][j] + 32) / 8))) + (0.6 / (1 + exp((50 - v[i][j]) / 20)));
	double tauf = 1 / (0.0197*exp(-pow(0.0337*(v[i][j] + 10), 2)) + 0.02);
	d[i][j] = dss - (dss - d[i][j])*exp(-dt / taud);
	f[i][j] = fss - (fss - f[i][j])*exp(-dt / tauf);

	//T-Type Ca Channel in LRd95
	double bss = 1 / (1 + exp(-(v[i][j] + 14.0) / 10.8));
	double taub = 3.7 + 6.1 / (1 + exp((v[i][j] + 25.0) / 4.5));
	double gss = 1 / (1 + exp((v[i][j] + 60.0) / 5.6));
	double taug;
	if (v[i][j] <= 0)
		taug = -0.875*v[i][j] + 12.0;
	else
		taug = 12.0;
	b[i][j] = bss - (bss - b[i][j])*exp(-dt / taub);
	g[i][j] = gss - (gss - g[i][j])*exp(-dt / taug);

	//Rapidly Activating K Current in LRd95
	double xrss = 1 / (1 + exp(-(v[i][j] + 21.5) / 7.5));
	double tauxr = 1 / (0.00138*(v[i][j] + 14.2) / (1 - exp(-0.123*(v[i][j] + 14.2))) + 0.00061*(v[i][j] + 38.9) / (exp(0.145*(v[i][j] + 38.9)) - 1));
	xr[i][j] = xrss - (xrss - xr[i][j])*exp(-dt / tauxr);

	//Slowly Activating K Current in LRd99 Viswanathan
	double xs1ss = 1 / (1 + exp(-(v[i][j] - 1.5) / 16.7));
	double xs2ss = xs1ss;
	double tauxs1 = 1 / (0.0000719*(v[i][j] + 30) / (1 - exp(-0.148*(v[i][j] + 30))) + 0.000131*(v[i][j] + 30) / (exp(0.0687*(v[i][j] + 30)) - 1));
	double tauxs2 = 4 * tauxs1;
	xs1[i][j] = xs1ss - (xs1ss - xs1[i][j])*exp(-dt / tauxs1);
	xs2[i][j] = xs2ss - (xs2ss - xs2[i][j])*exp(-dt / tauxs2);

	//Transient Outward Current in LRd99(Dumaine99)
	tauzdv = 1 / (azdv + bzdv);
	zssdv = azdv / (azdv + bzdv);
	zdv[i][j] = zssdv - (zssdv - zdv[i][j])*exp(-dt / tauzdv);
	aydv = 0.015 / (1 + exp((v[i][j] + 60) / 5));
	bydv = (0.1*exp((v[i][j] + 25) / 5)) / (1 + exp((v[i][j] + 25) / 5));
	tauydv = 1 / (aydv + bydv);
	yssdv = aydv / (aydv + bydv);
	ydv[i][j] = yssdv - (yssdv - ydv[i][j])*exp(-dt / tauydv);
}

//Call This Subroutine After Stimulus
//the new dt should not be larger than the left time of the current cycle
void Comp_Time_Left(){
	if (sti_flag < 0)
		dt_rem = tstim + b_n*bcl - t;
	else
		dt_rem = tstim + (1 + b_n)*bcl - t;
}

void comp_voltage(int i, int j, double dt){
	dvdtnew[i][j] = D1V[i][j];
	if (csqn[i][j] >= csqnth && tjsrol[i][j]>50){
		grelbarjsrol[i][j] = 4;
		tjsrol[i][j] = 0;
		printf("Spontaneous Release occurred at time: %g\n", t);
		//cout << "Spontaneous Release occurred at time " << t << endl;
	}
	caiontold[i][j] = caiont[i][j];
	dcaiont[i][j] = dcaiontnew[i][j];
	v[i][j] = v[i][j] + dt*D1V[i][j] + dt*dt * D2V[i][j] / 2;
	//!WRITE(6, *) "[D1V]", D1V
	//	v = v + dt*D1V + dt**2 * D2V / 2d0
	//	!WRITE(6, *)"{voltage_WU}", "[voltage]", v, "[D1V]", D1V, "[D2V]", D2V

	//	!!nai = nai + dt*nai_dot + dt**2 * D2V_nai / 2d0
	//	!WRITE(6, *)"{nai}", nai

	//	!!ki = ki + dt*ki_dot + dt**2 * D2V_ki / 2d0
	//	!WRITE(6, *)"{ki}", ki

	//	!!!djsr = djsr + dt*djsr_dot + dt**2 * D2V_djsr / 2d0  !!***!!
	//	!WRITE(6, *)"{jsr}", jsr

	//	!!nsr = nsr + dt*nsr_dot + dt**2 * D2V_nsr / 2d0
	//	!WRITE(6, *)"{nsr}", nsr

	//	!!!cai = cai + dt*cai_dot + dt**2 * D2V_cai / 2d0
	//	!WRITE(6, *)"{cai}", cai
}

// Calculating new time-step size
double comp_CCL(int i, int j, double dt){
	double dt_range;
	if (dt_univ > dt * 2){//Min(2*dt,dt_univ)
		dt_range = dt * 2;
	}
	else{
		dt_range = dt_univ;
	}

	//the new dt should not be larger than the left time of the current cycle
	//Comp_Time_Left();
	//if (dt_range > dt_rem){
	//	dt_range = dt_rem;
	//}

	D2V[i][j] = (D1V[i][j] - D1V_old[i][j]) / dt;
	double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
	if ((sti_flag >0) && (t + dt_min > bcl*(b_n)-epsilon)){
		//WRITE(6, *) n, t, dt_loc, start + (n)*CL, start + (n)*CL - t
		dt = tstim + bcl*b_n - t;/// not larger than the left time of the cycle
		printf("error !!!!");
		exit(1);
	}
	else{
		if (D1V[i][j] >= 0){
			DiscriminantP = D1V[i][j] * D1V[i][j] + 2 * D2V[i][j] * Voffset;
			if (D2V[i][j]>0){
				dt = (-D1V[i][j] + sqrt(DiscriminantP)) / D2V[i][j];
			}
			else if (D2V[i][j]<0){
				dtz = -D1V[i][j] / D2V[i][j];
				if (DiscriminantP >= 0){
					dt = (-D1V[i][j] + sqrt(DiscriminantP)) / D2V[i][j];
				}
				else{
					dt = dtz;
				}
			}
		}
		else{
			DiscriminantN = D1V[i][j] * D1V[i][j] - 2 * D2V[i][j] * Voffset;
			if (D2V[i][j]>0){
				dtz = -D1V[i][j] / D2V[i][j];
				if (DiscriminantN >= 0){
					dt = (-D1V[i][j] - sqrt(DiscriminantN)) / D2V[i][j];
				}
				else{
					dt = dtz;
				}
			}
			else if (D2V[i][j]<0){
				dt = (-D1V[i][j] - sqrt(DiscriminantN)) / D2V[i][j];
			}
		}
		if (dt>dt_range){
			dt = dt_range;
		}
		if (dt<dt_min){
			dt = dt_min;
		}
	}
	return dt;
}

//a package of current functions for updating new currents
double comp_currents(int i, int j, double dt){
	/* List of functions called for each timestep, currents commented
	out are only used when modeling pathological conditions */
	comp_ina(i, j);
	comp_ical(i, j);
	comp_icat(i, j);
	//comp_ik();
	comp_ikr(i, j);
	comp_iks(i, j);
	comp_iki(i, j);
	comp_ikp(i, j);
	//comp_ikna (); 
	//comp_ikatp (); 
	//comp_ito (); 
	comp_inaca(i, j);
	comp_inak(i, j);
	//comp_insca (); 
	comp_ipca(i, j);
	comp_icab(i, j);
	comp_inab(i, j);
	comp_it(i, j);
	return it[i][j];
}

//a package for calaculating ion concentrations
void comp_concentrations(int i, int j, double dt){
	conc_nai(i, j, dt);
	conc_ki(i, j, dt);
	calc_itr(i, j, dt);
	conc_jsr(i, j, dt);
	conc_nsr(i, j, dt);
	conc_cai(i, j, dt);
	//conc_cleft (); 
	/* Cleft Space disabled, if you want to use cleft space, make sure the initial conditions
	of ion concentrations in the bulk medium are the same as the extracellular concentrations */
}

/* Functions that describe the currents begin here */

/* Calculates fast sodium current */
void comp_ina(int i, int j)
{
	double gna = 16;
	double ena = ((Rgas*temp) / frdy)*log(nao / nai[i][j]);

	//am = 0.32*(v[i][j] + 47.13) / (1 - exp(-0.1*(v[i][j] + 47.13)));
	//bm = 0.08*exp(-v[i][j] / 11);
	//if (v[i][j] < -40){
	//	ah = 0.135*exp((80 + v[i][j]) / -6.8);
	//	bh = 3.56*exp(0.079*v[i][j]) + 310000 * exp(0.35*v[i][j]);
	//	aj = (-127140 * exp(0.2444*v[i][j]) - 0.00003474*exp(-0.04391*v[i][j]))*((v[i][j] + 37.78) / (1 + exp(0.311*(v[i][j] + 79.23))));
	//	bj = (0.1212*exp(-0.01052*v[i][j])) / (1 + exp(-0.1378*(v[i][j] + 40.14)));
	//}else{
	//	ah = 0;
	//	bh = 1 / (0.13*(1 + exp((v[i][j] + 10.66) / -11.1)));
	//	aj = 0;
	//	bj = (0.3*exp(-0.0000002535*v[i][j])) / (1 + exp(-0.1*(v[i][j] + 32)));
	//}
	//mtau = 1 / (am + bm);
	//htau = 1 / (ah + bh);
	//jtau = 1 / (aj + bj);

	//mss = am*mtau;
	//hss = ah*htau;
	//jss = aj*jtau;

	////m = mss - (mss - m)*exp(-dt / mtau);
	////h = hss - (hss - h)*exp(-dt / htau);
	////jj = jss - (jss - jj)*exp(-dt / jtau);

	ina[i][j] = gna*m[i][j] * m[i][j] * m[i][j] * h[i][j] * jj[i][j] * (v[i][j] - ena);
}

/* Calculates Currents through L-Type Ca Channel */
void comp_ical(int i, int j)
{
	//dss = 1 / (1 + exp(-(v[i][j] + 10) / 6.24));
	//taud = dss*(1 - exp(-(v[i][j] + 10) / 6.24)) / (0.035*(v[i][j] + 10));

	//fss = (1 / (1 + exp((v[i][j] + 32) / 8))) + (0.6 / (1 + exp((50 - v[i][j]) / 20)));
	//tauf = 1 / (0.0197*exp(-pow(0.0337*(v[i][j] + 10), 2)) + 0.02);

	//d = dss - (dss - d)*exp(-dt / taud);
	//f = fss - (fss - f)*exp(-dt / tauf);
	double ibarca = pca*zca*zca*((v[i][j]*frdy*frdy) / (Rgas*temp))
		*((gacai*cai[i][j]*exp((zca*v[i][j]*frdy) / (Rgas*temp)) - gacao*cao) / (exp((zca*v[i][j]*frdy) / (Rgas*temp)) - 1));
	double ibarna = pna*zna*zna*((v[i][j]*frdy*frdy) / (Rgas*temp))
		*((ganai*nai[i][j]*exp((zna*v[i][j]*frdy) / (Rgas*temp)) - ganao*nao) / (exp((zna*v[i][j]*frdy) / (Rgas*temp)) - 1));
	double ibark = pk*zk*zk*((v[i][j]*frdy*frdy) / (Rgas*temp))
		*((gaki*ki[i][j]*exp((zk*v[i][j]*frdy) / (Rgas*temp)) - gako*ko) / (exp((zk*v[i][j]*frdy) / (Rgas*temp)) - 1));
	double fca = 1 / (1 + cai[i][j] / kmca);//Hill coefficient=1 in LRd95, Hill coefficient=2 in LRd94
	ilca[i][j] = d[i][j] * f[i][j] * fca*ibarca;
	ilcana[i][j] = d[i][j] * f[i][j] * fca*ibarna;
	ilcak[i][j] = d[i][j] * f[i][j] * fca*ibark;
	ilcatot[i][j] = ilca[i][j] + ilcana[i][j] + ilcak[i][j];
}

/* Calculates Currents through T-Type Ca Channel in LRd95 */
void comp_icat(int i, int j)
{
	//bss = 1 / (1 + exp(-(v[i][j] + 14.0) / 10.8));
	//taub = 3.7 + 6.1 / (1 + exp((v[i][j] + 25.0) / 4.5));

	//gss = 1 / (1 + exp((v[i][j] + 60.0) / 5.6));
	//if (v[i][j] <= 0)
	//	taug = -0.875*v[i][j] + 12.0;
	//else
	//	taug = 12.0;

	//b = bss - (bss - b)*exp(-dt / taub);
	//g = gss - (gss - g)*exp(-dt / taug);
	gcat = 0.05;
	double eca = (Rgas*temp / (2 * frdy))*log(cao / cai[i][j]);
	icat[i][j] = gcat*b[i][j] * b[i][j] * g[i][j] * (v[i][j] - eca);
}

/* Calculates Potassium Current (time-dependent) in LRd94*/
//void comp_ik()
//{
//	gk = 0.282*sqrt(ko / 5.4);
//	xi = 1 / (1 + exp((v[i][j] - 56.26) / 32.1));
//	ax = 0.0000719*(v[i][j] + 30) / (1 - exp(-0.148*(v[i][j] + 30)));
//	bx = 0.000131*(v[i][j] + 30) / (-1 + exp(0.0687*(v[i][j] + 30)));
//	taux = 1 / (ax + bx);
//	xss = ax*taux;
//	x0 = xss - (xss - x)*exp(-dt / taux);
//	ik = gk*xi*x0*x0*(v[i][j] - ek);
////}

/* Calculates Rapidly Activating K Current in LRd95 */
void comp_ikr(int i, int j)
{
	double gkr = 0.02614*sqrt(ko / 5.4);
	double ekr = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);

	//xrss = 1 / (1 + exp(-(v[i][j] + 21.5) / 7.5));
	//tauxr = 1 / (0.00138*(v[i][j] + 14.2) / (1 - exp(-0.123*(v[i][j] + 14.2))) + 0.00061*(v[i][j] + 38.9) / (exp(0.145*(v[i][j] + 38.9)) - 1));
	//xr = xrss - (xrss - xr)*exp(-dt / tauxr);
	double r = 1 / (1 + exp((v[i][j] + 9) / 22.4));
	ikr[i][j] = gkr*xr[i][j] * r*(v[i][j] - ekr);
}

/* Calculates Slowly Activating K Current in LRd99 Viswanathan */
void comp_iks(int i, int j)
{
	double gks = 0.433*(1 + 0.6 / (1 + pow((0.000038 / cai[i][j]), 1.4)));
	double eks = ((Rgas*temp) / frdy)*log((ko + prnak*nao) / (ki[i][j] + prnak*nai[i][j]));

	//xs1ss = 1 / (1 + exp(-(v[i][j] - 1.5) / 16.7));
	//xs2ss = xs1ss;
	//tauxs1 = 1 / (0.0000719*(v[i][j] + 30) / (1 - exp(-0.148*(v[i][j] + 30))) + 0.000131*(v[i][j] + 30) / (exp(0.0687*(v[i][j] + 30)) - 1));
	//tauxs2 = 4 * tauxs1;
	//xs1 = xs1ss - (xs1ss - xs1)*exp(-dt / tauxs1);
	//xs2 = xs2ss - (xs2ss - xs2)*exp(-dt / tauxs2);
	iks[i][j] = gks*xs1[i][j] * xs2[i][j] * (v[i][j] - eks);
}

/* Calculates Potassium Current (time-independent) in LRd94*/
void comp_iki(int i, int j)
{
	double gki = 0.75*(sqrt(ko / 5.4));
	double eki = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);

	double aki = 1.02 / (1 + exp(0.2385*(v[i][j] - eki - 59.215)));
	double bki = (0.49124*exp(0.08032*(v[i][j] - eki + 5.476)) + exp(0.06175*(v[i][j] - eki - 594.31))) / (1 + exp(-0.5143*(v[i][j] - eki + 4.753)));
	double kin = aki / (aki + bki);
	iki[i][j] = gki*kin*(v[i][j] - eki);
}

/* Calculates Plateau K Current */
void comp_ikp(int i, int j)
{
	double gkp = 0.00552; //the data in LRd94 paper is 0.0183 ??????????????
	//gkp = 0.0183;
	//ekp = eki;
	double ekp = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);

	double kp = 1 / (1 + exp((7.488 - v[i][j]) / 5.98));

	ikp[i][j] = gkp*kp*(v[i][j] - ekp);
}

/* Calculates Na-activated K Current in LRd00*/
void comp_ikna(int i, int j)
{
	double ekna = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
	double pona = 0.85 / (1 + pow((kdkna / nai[i][j]), 2.8));
	double pov = 0.8 - 0.65 / (1 + exp((v[i][j] + 125) / 15));
	ikna[i][j] = gkna*pona*pov*(v[i][j] - ekna);
}

/* Calculates ATP-Sensitive K Current in LRd97 */
void comp_ikatp(int i, int j)
{
	/* Note: If you wish to use this current in your simulations, there are additional */
	/* changes which must be made to the code as detailed in Cardiovasc Res 1997;35:256-272 */

	double ekatp = ((Rgas*temp) / frdy)*log(ko / ki[i][j]);
	double gkatp = 0.000195 / nicholsarea;
	double patp = 1 / (1 + (pow((atpi / katp), hatp)));
	double gkbaratp = gkatp*patp*(pow((ko / 4), natp));

	ikatp[i][j] = gkbaratp*(v[i][j] - ekatp);
}

/* Calculates Transient Outward Current in LRd99(Dumaine99) */
void comp_ito(int i, int j)
{
	double gitodv = 0.5;
	double ekdv = ((Rgas*temp) / frdy)*log((ko) / (ki[i][j]));
	//ekdv = ((Rgas*temp) / frdy)*log((ko + prnak*nao) / (ki[i][j] + prnak*nai[i][j])); //LRd94
	double rvdv = exp(v[i][j] / 100);
	double azdv = (10 * exp((v[i][j] - 40) / 25)) / (1 + exp((v[i][j] - 40) / 25));
	double bzdv = (10 * exp(-(v[i][j] + 90) / 25)) / (1 + exp(-(v[i][j] + 90) / 25));
	//tauzdv = 1 / (azdv + bzdv);
	//zssdv = azdv / (azdv + bzdv);
	//zdv = zssdv - (zssdv - zdv)*exp(-dt / tauzdv);

	//aydv = 0.015 / (1 + exp((v[i][j] + 60) / 5));
	//bydv = (0.1*exp((v[i][j] + 25) / 5)) / (1 + exp((v[i][j] + 25) / 5));
	//tauydv = 1 / (aydv + bydv);
	//yssdv = aydv / (aydv + bydv);
	//ydv = yssdv - (yssdv - ydv)*exp(-dt / tauydv);
	ito[i][j] = gitodv*zdv[i][j] * zdv[i][j] * zdv[i][j] * ydv[i][j] * rvdv*(v[i][j] - ekdv);
}

/* Calculates Na-Ca Exchanger Current in LRd00*/
void comp_inaca(int i, int j)
{
	inaca[i][j] = c1*exp((gammas - 1)*v[i][j] * frdy / (Rgas*temp))
		*((exp(v[i][j]*frdy / (Rgas*temp))*nai[i][j]*nai[i][j]*nai[i][j]*cao - nao*nao*nao*cai[i][j])
		/ (1 + c2*exp((gammas - 1)*v[i][j]*frdy / (Rgas*temp))*(exp(v[i][j]*frdy / (Rgas*temp))*nai[i][j]*nai[i][j]*nai[i][j]*cao + nao*nao*nao*cai[i][j])));

	//inaca = knaca / (kmna*kmna*kmna + nao*nao*nao) / (kmca2 + cao) / (1 
	//+ ksat*exp((gamma - 1)*v[i][j]*Frdy / ( Rgas*Temp)))*(exp(gamma*v[i][j]*Frdy / ( Rgas*Temp))*nai[i][j]*nai[i][j]*nai[i][j]*cao 
	//- exp((gamma - 1)*v[i][j]*Frdy / ( Rgas*Temp))*nao*nao*nao*cai[i][j]);
}

/* Calculates Sodium-Potassium Pump Current */
void comp_inak(int i, int j)
{
	double sigma = (exp(nao / 67.3) - 1) / 7;
	double fnak = 1 / (1 + 0.1245*exp((-0.1*v[i][j]*frdy) / (Rgas*temp)) + 0.0365*sigma*exp((-v[i][j]*frdy) / (Rgas*temp)));

	inak[i][j] = ibarnak*fnak*(1 / (1 + pow(kmnai / nai[i][j], 2)))*(ko / (ko + kmko));//in LRd94 paper, index is 1.5, not 2 !!!
}

//LRd94, here insk and insna are not included
/* Calculates Non-Specific ca-Activated Current in LRd94*/
void comp_insca(int i, int j)
{
	//pnsca Permeability of channel to Na and K (cm/s) (equally permeable)
	double ibarnsna = pnsca*zna*zna*((v[i][j]*frdy*frdy) / (Rgas*temp))
		*((ganai*nai[i][j]*exp((zna*v[i][j]*frdy) / (Rgas*temp)) - ganao*nao) / (exp((zna*v[i][j]*frdy) / (Rgas*temp)) - 1));
	double ibarnsk = pnsca*zk*zk*((v[i][j]*frdy*frdy) / (Rgas*temp))
		*((gaki*ki[i][j]*exp((zk*v[i][j]*frdy) / (Rgas*temp)) - gako*ko) / (exp((zk*v[i][j]*frdy) / (Rgas*temp)) - 1));

	insna[i][j] = ibarnsna / (1 + pow(kmnsca / cai[i][j], 3));
	insk[i][j] = ibarnsk / (1 + pow(kmnsca / cai[i][j], 3));
}

/* Calculates Sarcolemmal Ca Pump Current */
void comp_ipca(int i, int j)
{
	ipca[i][j] = (ibarpca*cai[i][j]) / (kmpca + cai[i][j]);
}

/* Calculates Ca Background Current */
void comp_icab(int i, int j)
{
	double gcab = 0.003016;
	double ecan = ((Rgas*temp) / (2 * frdy))*log(cao / cai[i][j]);
	icab[i][j] = gcab*(v[i][j] - ecan);
}

/* Calculates Na Background Current */
void comp_inab(int i, int j)
{
	double gnab = 0.004;//Increase of gNa,b max from 0.00141 in LRd94 to 0.004 in LRd99.
	double ena = ((Rgas*temp) / frdy)*log(nao / nai[i][j]);
	double enan = ena;
	inab[i][j] = gnab*(v[i][j] - enan);
}

/* Total sum of currents is calculated here, if the time is between
stimtime = 0 and stimtime = duration, a stimulus is applied */
void comp_it(int i, int j)
{
	naiont[i][j] = ina[i][j] + inab[i][j] + ilcana[i][j] + 3 * inak[i][j] + 3 * inaca[i][j];//insna is removed(LRd94)
	//kiont below: ik is replaced by ikr + iks(LRd95), ikna(LRd00)+ikatp(LRd97)+insk(LRd94) are removed
	kiont[i][j] = ikr[i][j] + iks[i][j] + iki[i][j] + ikp[i][j] + ilcak[i][j] - 2 * inak[i][j] + ito[i][j];
	caiont[i][j] = ilca[i][j] + icab[i][j] + ipca[i][j] - 2 * inaca[i][j] + icat[i][j];
	if (dvdtnew[i][j] > 10 && tcicr[i][j] > 10 && flag[i][j] == 1){
		flag[i][j] = 0;
	}

	if (t >= tstim && t<(tstim + dt[i][j])){
		stimtime = 0;
		//b_i = b_i + 1;
		tstim = tstim + bcl;
		printf("t:%g\n", t);
		//cout << t << endl;
	}
	if (stimtime == 0){
		rmbp[i][j] = v[i][j];
		nair[i][j] = nai[i][j];
		kir[i][j] = ki[i][j];
		cair[i][j] = cai[i][j];
	}

	it[i][j] = naiont[i][j] + kiont[i][j] + caiont[i][j];
	//if (ncount == 1035){
	//	double aa, bb, cc;
	//	aa = it[i][j];
	//	printf("%.18f\n", aa);
	//	bb = (stim[i][j]);
	//	printf("%.18f\n", bb);
	//	cc = aa + bb;
	//	printf("%.18f\n", cc);
	//	printf("%.18f\n", bb + it[i][j]);
	//}
	//it[i][j] = ina[i][j] + inab[i][j] + ilcana[i][j] + inak[i][j] + inaca[i][j]
	//	+ ikr[i][j] + iks[i][j] + iki[i][j] + ikp[i][j] + ilcak[i][j] + ito[i][j]
	//	+ ilca[i][j] + icab[i][j] + ipca[i][j] + icat[i][j];
	if (stimtime >= 0 && stimtime <= duration){
		stim = amp;
		sti_flag = 1;//for CCL
	}
	else{
		stim = 0;
		sti_flag = -1;//for CCL
	}
}

/* Functions that calculate intracellular ion concentrations begins here */

/* Calculates new myoplasmic Na ion concentration */
void conc_nai(int i, int j, double dt)
{
	// The units of dnai is in mM. Note that naiont should be multiplied by the 
	// cell capacitance to get the correct units. Since cell capacitance = 1 uF/cm^2, 
	// it doesn't explicitly appear in the equation below. 
	// This holds true for the calculation of dki and dcai. */ 

	double dnai = -dt*(naiont[i][j] * acap) / (vmyo*zna*frdy);
	nai[i][j] = dnai + nai[i][j];
}

/* Calculates new myoplasmic K ion concentration */
void conc_ki(int i, int j, double dt)
{
	double dki;
	//if (stimtime >= 0 && stimtime <= duration){
	dki = -dt*((kiont[i][j] + stim)*acap) / (vmyo*zk*frdy);
	//}
	//else{
	//	dki = -dt*(kiont[i][j] * acap) / (vmyo*zk*frdy);
	//}

	ki[i][j] = dki + ki[i][j];
}

/* Translocation of Ca2+ ions from NSR to JSR: */
void calc_itr(int i, int j, double dt)
{
	//Itr, Translocation current of Ca ions from NSR to JSR (mM/ms) 
	itr[i][j] = (nsr[i][j] - jsr[i][j]) / tautr;
}

/* Calculates new JSR Ca ion concentration */
void conc_jsr(int i, int j, double dt)
{
	// there is always calcium release triggered by calcium entry, 
	// and the release becomes significant as calcium entry increases.
	//Ca2+ uptake and leakage of NSR: Iup and Ileak
	double kleak = iupbar / nsrbar;
	ileak[i][j] = kleak*nsr[i][j];

	iup[i][j] = iupbar*cai[i][j] / (cai[i][j] + kmup);

	dcaiontnew[i][j] = (caiont[i][j] - caiontold[i][j]) / dt;	// Total Ca Ion Flow (uA/uF) 

	//dcaiont; // Rate of change of Ca entry 
	// half-saturation potential(Vh),-35mV.
	if (v[i][j]>-35 && dcaiontnew[i][j]>dcaiont[i][j] && flag[i][j] == 0){
		flag[i][j] = 1;// Flag condition to test for dvdtmax
		tcicr[i][j] = 0;// t=0 at time of CICR (ms).initial value is 1000, make exponential term approached 0
	}

	// Time constant of activation of Ca release from JSR (ms)
	// tcicr, t=0 at time of CICR (ms)
	// ref to Biophys J 78:2392-404, 2000
	double on = 1 / (1 + exp((-tcicr[i][j] + 4) / tauon));
	double off = (1 - 1 / (1 + exp((-tcicr[i][j] + 4) / tauoff)));
	double magrel = 1 / (1 + exp(((ilca[i][j] + icab[i][j] + ipca[i][j] - 2 * inaca[i][j] + icat[i][j]) + 5) / 0.9)); // Magnitude of Ca release
	irelcicr[i][j] = gmaxrel*on*off*magrel*(jsr[i][j] - cai[i][j]);// Ca release from JSR to myo. due to CICR (mM/ms)

	tcicr[i][j] = tcicr[i][j] + dt;

	/********* irelcicr *************/
	////dcai2 not consider Ca released from SR
	//dcai2 = -dt*(((caiont*acap) / (vmyo*zca*Frdy)) + ((iup - ileak)*vnsr / vmyo));
	////vnew = v[i][j] - it*dt;
	////dvdtnew = (vnew - v[i][j]) / dt;

	///*Once it is triggered(at 2 milliseconds), the JSR releases almost all of its free Ca2+ contents.
	//In the model, we compute the increase in [Ca2+]	2 milliseconds (rather than 10 milliseconds)
	//from the time of Vmax and denote it as d[Ca2+]i,2 ......or from the time of stimulus..*/
	//if (vnew > v[i][j] && vnew > vmax_jsr){//vmax_jsr's initial value is -90
	//	vmax_jsr = vnew;
	//	tvmax = 0;
	//	dcaisum = 0;
	//	tcicr = 0;
	//}
	//else if (firstvmax == 1){// find Vmax!!
	//	flagvmax = 1;
	//	firstvmax = 0;
	//}

	//if (flagvmax == 1){
	//	tvmax = tvmax + dt;  //cumulative amount of time, not larger than 2 msec
	//	dcaisum = dcaisum + dcai2; //cumulative amount of Ca2+
	//	if ((tvmax <= 2 && dcaisum > dcaith && firstcicr == 1)){//firstcicr ensure only one CICR in one action potential duration
	//		tcicr = 0;
	//		flagcicr = 1;
	//		firstcicr = 0;
	//	}
	//}

	//if (flagcicr == 1){
	//	grelcicr = grelcicrbar*(dcaisum - dcaith) / (kmrel + dcaisum - dcaith)*(1 - exp(-tcicr / tauon))*exp(-tcicr / tauoff);
	//}
	//else{
	//	grelcicr = 0;
	//}

	//irelcicr = grelcicr*(jsr - cai[i][j]);
	//tcicr = tcicr + dt;
	/*********** irelcicr ***********/

	/*********** ireljsrol ***********/
	// Rate constant of Ca release from JSR due to overload (ms^-1) 
	// The threshold of overload is set, can be founded by searching "grelbarjsrol"
	double greljsrol = grelbarjsrol[i][j] * (1 - exp(-tjsrol[i][j] / tauon))*exp(-tjsrol[i][j] / tauoff);
	ireljsrol[i][j] = greljsrol*(jsr[i][j] - cai[i][j]);
	tjsrol[i][j] = tjsrol[i][j] + dt;
	/*********** ireljsrol ***********/

	/*********** jsr ***********/
	//Method 1. Analytical Computation  ref to Circ Res 1995;77:140-152
	csqn[i][j] = csqnbar*(jsr[i][j] / (jsr[i][j] + kmcsqn));//[CSQN],Calsequestrin Buffered Ca Concentration
	double djsr = dt*(itr[i][j] - irelcicr[i][j] - ireljsrol[i][j]);//dCa2+_JSR
	double bjsr = csqnbar - csqn[i][j] - djsr - jsr[i][j] + kmcsqn;
	double cjsr = kmcsqn*(csqn[i][j] + djsr + jsr[i][j]);

	jsr[i][j] = (sqrt(bjsr*bjsr + 4 * cjsr) - bjsr) / 2;// Ca buffering in JSR, new
	/*********** jsr ***********/

	/*********** jsr ***********/
	//Method 2. Steffensen's iterative method
	/*djsr = dt*(itr-irelcicr-ireljsrol);
	jsr1 = jsr+djsr*0.05;

	ii=1;
	while (ii <= 20) // Number of iteration
	{
	csqn1 = csqnbar/(1+kmcsqn/jsr1);
	jsr2 = jsr+djsr-(csqn1-csqn);
	csqn2 = csqnbar/(1+kmcsqn/jsr2);
	jsr3 = jsr+djsr-(csqn2-csqn);
	jsrk1 = jsr1-pow((jsr1-jsr2), 2)/(jsr1- 2*jsr2+jsr3);

	if (abs((jsrk1-jsr1)/jsr1)<0.01){
	jsr = jsrk1;
	csqn = csqnbar/(1+kmcsqn/jsr);
	break;
	} else if (((jsr1-2*jsr2+jsr3)/jsr1)<=0.0001){
	jsr = jsr3;
	csqn = csqnbar/(1+kmcsqn/jsr);
	break;
	} else{
	jsr1 = jsrk1;
	}
	ii++;
	}*/
	/*********** jsr ***********/
}

/* Calculates new NSR Ca ion concentration */
void conc_nsr(int i, int j, double dt)
{
	//itr;  Translocation current of Ca ions from NSR to JSR (mM/ms) 
	//vjsr; JSR volume, vnsr; NSR volume 
	double dnsr = dt*(iup[i][j] - ileak[i][j] - itr[i][j] * vjsr / vnsr);
	nsr[i][j] = nsr[i][j] + dnsr;
}

//take care the total Ca2+ concentration,catotal
//Steffensen Iterative Method can also be used here
/* Calculates total myoplasmic Ca concentration (mM) */
void conc_cai(int i, int j, double dt)
{
	//Method 1. Steffensen Iterative Method .....

	/*
	//Method 2.  Forward Euler Method
	B1 = -(((caiont*acap) / (vmyo*zca*Frdy)) + ((iup - ileak)*vnsr / vmyo) - (irelcicr*vjsr / vmyo) - (ireljsrol*vjsr / vmyo));
	B2 = 1 + (trpnbar*kmtrpn / pow((cai[i][j] + kmtrpn), 2)) + (cmdnbar*kmcmdn / pow((cai[i][j] + kmcmdn), 2));
	dcai = (B1 / B2)*dt; // Change in myoplasmic Ca concentration (mM)
	cai[i][j] = cai[i][j] + dcai;
	*/

	//Myoplasmic Ca Ion Concentration Changes
	double dcai = -dt*(((caiont[i][j] * acap) / (vmyo*zca*frdy)) + ((iup[i][j] - ileak[i][j])*vnsr / vmyo)
		- (irelcicr[i][j] * vjsr / vmyo) - (ireljsrol[i][j] * vjsr / vmyo));
	trpn[i][j] = trpnbar*(cai[i][j] / (cai[i][j] + kmtrpn));
	cmdn[i][j] = cmdnbar*(cai[i][j] / (cai[i][j] + kmcmdn));

	//Method 3. Analytical Computation  ref to Circ Res 1995;77:140-152
	//Total myoplasmic Ca concentration (mM)
	catotal[i][j] = trpn[i][j] + cmdn[i][j] + dcai + cai[i][j];
	double bmyo = cmdnbar + trpnbar - catotal[i][j] + kmtrpn + kmcmdn;
	double cmyo = (kmcmdn*kmtrpn) - (catotal[i][j] * (kmtrpn + kmcmdn)) + (trpnbar*kmcmdn) + (cmdnbar*kmtrpn);
	double dmyo = -kmtrpn*kmcmdn*catotal[i][j];
	double gpig = sqrt(bmyo*bmyo - 3 * cmyo);

	cai[i][j] = (2 * gpig / 3)*cos(acos((9 * bmyo*cmyo - 2 * bmyo*bmyo*bmyo - 27 * dmyo)
		/ (2 * pow((bmyo*bmyo - 3 * cmyo), 1.5))) / 3) - (bmyo / 3);
}

/* Cleft Space disabled, if you want to use cleft space, make sure the initial conditions
of ion concentrations in the bulk medium are the same as the extracellular concentrations */
/* Calculates new cleft ion concentrations  */
//void conc_cleft()
//{
//	dnao = dt*((nabm - nao) / taudiff + naiont*acap / (vcleft*frdy));
//	nao = dnao + nao;
//
//	if (stimtime >= 0 && stimtime <= duration){
//		dko = dt*((kbm - ko) / taudiff + (kiont + stim)*acap / (vcleft*frdy));
//	}
//	else{
//		dko = dt*((kbm - ko) / taudiff + kiont*acap / (vcleft*frdy));
//	}
//	ko = dko + ko;
//
//	dcao = dt*((cabm - cao) / taudiff + caiont*acap / (vcleft*frdy * 2));
//	cao = dcao + cao;
//}

///* Values are printed to a file called ap. The voltage and
//currents can be plotted versus time using graphing software. */
//void prttofile()
//{
//	if (t>(0) && t<(bcl*beats))
//		fprintf(ap, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
//		t, v, nai, ki, cai, jsr, nsr, ina, ikr, iks, iki, ilca, icat, inab, icab, ipca, inaca, inak);
//
//}