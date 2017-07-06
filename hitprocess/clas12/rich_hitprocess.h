#ifndef RICH_HITPROCESS_H
#define RICH_HITPROCESS_H 1

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>

// gemc headers
#include "HitProcess.h"

// constants to be used in the digitization routine
class richConstants
{
public:

	// database
	int    runNo;
	string variation;
	string date;
	string connection;

	// translation table
	TranslationTable TT;

	// H8500 quantum efficiency table
	static const int nH8500=85;
	double energy_H8500[nH8500];
	double QE_H8500[nH8500];

	// H8500C quantum efficiency table
	static const int nH8500C=91;
	double energy_H8500C[nH8500C];
	double QE_H8500C[nH8500C];


	/* Pixel geometry constants */
	static const int nrows = 8;
	static const int ncols = 8;
	double PixelBig = 6.26;
	double PixelSmall = 6.08;
	double PixelDeadSpace_x = 0.28;
	double PixelDeadSpace_y = 0.28;

	double PhCsize_x;
	double PhCsize_y;
	double PXDX[ncols];
	double PXDY[nrows];


	/* total number of MAPMT (it should be taken from the geometry) */
	double nMAPMT = 391;
	/* Average Gain at 1000 V (default value) */
	double Gain = 3e6;
	/* Average Dark Current at 1000 V (A) (default value) */
	double_t DarkCurrent = 0.5e-9;

	/* Time window (ns), for the noise generation (it sould be taken from rich__hit.txt) */
	double TimeWindow = 250;
};



// Class definition
class rich_HitProcess : public HitProcess
{
public:

	~rich_HitProcess(){;}

	// - integrateDgt: returns digitized information integrated over the hit
	map<string, double> integrateDgt(MHit*, int);

	// - multiDgt: returns multiple digitized information / hit
	map< string, vector <int> > multiDgt(MHit*, int);

	// - charge: returns charge/time digitized information / step
	virtual map< int, vector <double> > chargeTime(MHit*, int);

	// - voltage: returns a voltage value for a given time. The input are charge value, time
	virtual double voltage(double, double, double);
	
	// The pure virtual method processID returns a (new) identifier
	// containing hit sharing information
	vector<identifier> processID(vector<identifier>, G4Step*, detector);

	// creates the HitProcess
	static HitProcess *createHitClass() {return new rich_HitProcess;}

private:

	// constants initialized with initWithRunNumber
	static richConstants richc;

	void initWithRunNumber(int runno);

	// - electronicNoise: returns a vector of hits generated / by electronics.
	vector<MHit*> electronicNoise();

	// - GetCharge: returns the charge on the MAPMT. Inputs are number of p.e., gain
	virtual double GetCharge(int, double);

	// - Charge2ADC: translate the MAPMT charge to MAROC DAC units 
	virtual int Charge2ADC(double charge);

	vector<int> Time2TDC(double time, double charge);
};

#endif
