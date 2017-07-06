// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// gemc headers
#include "rich_hitprocess.h"

/* Electron charge in fC */
double ElectronCharge = 1.6e-4;


static richConstants initializeRICHConstants(int runno)
{
	cout << "Entering initializeRICHConstants" << endl;

	// all these constants should be read from CCDB
	richConstants richc;


	// database
	richc.runNo = runno;
	richc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		richc.connection = (string) getenv("CCDB_CONNECTION");
	else
		richc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	richc.variation  = "main";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(richc.connection));


	/*****************************************************************/
	/* Location of the temporary txt files for the calibration constants */
	char *JLAB_ROOT = getenv("JLAB_ROOT");
	char TxtFileName[200];
	
	/* Loading H8500 and H8500C quantum efficiency tables (must add H12700)
	   NOTE that the table has decreasing energy ordering */
	FILE *fIn = NULL;

	sprintf(TxtFileName, "%s/database/QE_H8500.dat", JLAB_ROOT);
	fIn = fopen(&TxtFileName[0], "r");
	if (fIn) {
	  cout << "  Loading H8500 QE table" << endl;
	  for (unsigned i = 0; i < richc.nH8500; i++) {
	    fscanf(fIn, "%lf %lf \n", &richc.energy_H8500[i], &richc.QE_H8500[i]);
	  }

	  fclose(fIn);
	}
	else { 
	  cout << " ==>> WARNING: Cannot read RICH H8500 QE table, from file: " << TxtFileName << endl;
	  cout << "               QE set to 0" << endl;
	}

	sprintf(TxtFileName, "%s/database/QE_H8500C.dat", JLAB_ROOT);
	fIn = fopen(&TxtFileName[0], "r");
	if (fIn) {
	  cout << "  Loading H8500C QE table" << endl;
	  for (unsigned i = 0; i < richc.nH8500C; i++) {
	    fscanf(fIn, "%lf %lf \n", &richc.energy_H8500C[i], &richc.QE_H8500C[i]);
	  }

	  fclose(fIn);
	}
	else { 
	  cout << " ==>> WARNING: Cannot read RICH H8500C QE table, from file: " << TxtFileName << endl;
	  cout << "               QE set to 0" << endl;
	}

	/* Loading the PMT gain and dark current */
	sprintf(TxtFileName, "%s/database/PMT.dat", JLAB_ROOT);
	fIn = fopen(&TxtFileName[0], "r");
	if (fIn) {
	  cout << "  Loading PMT data" << endl;

	  fscanf(fIn, "%lf %lf \n", &richc.Gain, &richc.DarkCurrent);


	  fclose(fIn);
	}
	else { 
	  cout << " ==>> WARNING: Cannot read RICH H8500 QE table, from file: " << TxtFileName << endl;
	  cout << "               QE set to 0" << endl;
	}
	


	/* Definition of the photocathode pixelization (H8500 and H8500C and H12700) */
	richc.PhCsize_x = 2*richc.PixelBig + 6*richc.PixelSmall;
	richc.PhCsize_y = 2*richc.PixelBig + 6*richc.PixelSmall;


	richc.PXDX[0] = richc.PixelBig;
	richc.PXDX[richc.ncols-1] = richc.PixelBig;
	for (unsigned i = 1; i < richc.ncols-1; i++) {
	  richc.PXDX[i] = richc.PixelSmall;
	}
	richc.PXDY[0] = richc.PixelBig;
	richc.PXDY[richc.nrows-1] = richc.PixelBig;
	for (unsigned i = 1; i < richc.nrows-1; i++) {
	  richc.PXDY[i] = richc.PixelSmall;
	}



	return richc;
}


/* ============================================================================================== */

map<string, double> rich_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
  /* This function is called at the end of the RICH processing 
     Here compute all the final digitized information:
     - check of the MAPMT dead area
     - check of the QE 
     - translation of number of pe into ADC
     - translation of time into TDC for leading and trailing edge

*/

  //cout << "=====>>> Starting integrateDgt "<< endl;

  map<string, double> dgtz;
  vector<identifier> identity = aHit->GetId();

  /* Current hit infos */
  int sector = identity[0].id;
  int pmt    = identity[1].id;
  int pixel  = identity[2].id;
  int thisPid   = aHit->GetPID();


  /* With this here there is segmentation fault for noise hits */
  //trueInfos tInfos(aHit);


  /******************************************************/
  /* Electronic noise hits */
  /******************************************************/
  if(aHit->isElectronicNoise) {
    //cout << "--->>> DGT of noise hit  pmt=" << pmt << "  pixel=" << pixel << endl;
    dgtz["sector"] = sector;
    dgtz["pmt"]    = pmt;
    dgtz["pixel"]  = pixel;

    dgtz["nphotons"]   = 0;
    dgtz["npe"]   = 1;
    dgtz["nphotonsU"]   = 0;
    dgtz["nphotonsO"]   = 0;
    dgtz["nphotonsD"]   = 0;

    dgtz["hitn"]   = -hitn;

    dgtz["TDC1"]   = 0.;
    dgtz["TDC2"]   = 0.;
    dgtz["ADC"]   = 0.;

    return dgtz;
  }
  else {
    //cout << "  No noise hits " << endl;
  }


  /*****************************************************/
  /* Digitization following htcc_hitprocess.cc */
  /*****************************************************/

  /* True information of the current hit 
     If the hit is due to more than 1 tracks, it is the average of all the tracks */
  trueInfos tInfos(aHit);


  if (thisPid != 0) {

    // if anything else than a photon hits the PMT
    // the nphe is the particle id
    // and identifiers are negative
    // this should be changed, what if we still have a photon later?

    dgtz["sector"] = -sector;
    dgtz["pmt"]    = -pmt;
    dgtz["pixel"]  = -pixel;
    dgtz["nphotons"]   = thisPid;
    dgtz["npe"]   = 0;
    dgtz["nphotonsU"]   = 0;
    dgtz["nphotonsO"]   = 0;
    dgtz["nphotonsD"]   = 0;

    dgtz["hitn"]    = hitn;
    dgtz["TDC1"]   = 0.;
    dgtz["TDC2"]   = 0.;
    dgtz["ADC"]   = 0.;
    
    return dgtz;
  }

  
  /************************************************/
  /* LOOKING AT PHOTONS NOW */

  vector<int> tids = aHit->GetTIds();      // track ID at EACH STEP
  vector<int> pids = aHit->GetPIDs();      // particle ID at EACH STEP
  vector<double> Energies = aHit->GetEs(); // energy of the photon as it reach the pmt
  vector<double> Times = aHit->GetTime(); // time of the photon as it reach the pmt
  vector<G4ThreeVector> Lpos = aHit->GetLPos(); // local coordinates of the hit

  int narrived  = 0;
  int noutofrange = 0;
  int ndetected = 0;
  int nunconverted = 0;
  int ndeadarea = 0;
  double firstphotontime = 10000;
	

  /* Loading information from the first step of the hit */
  map<int, double> penergy;  // key is track id
  map<int, G4ThreeVector> lposition;  // key is track id
  map<int, double> ptime;  // key is track id

  for(unsigned int s=0; s<tids.size(); s++) {
    // only insert the first step of each track
    // (i.e. if the map is empty
    if(penergy.find(tids[s]) == penergy.end()) {
      penergy[tids[s]] = Energies[s];

      ptime[tids[s]] = Times[s];

      lposition[tids[s]] = Lpos[s];
    }

  }



  /*****************************************************/
  /* loop over the photons of the current hit */
  for( unsigned int iphoton = 0; iphoton<penergy.size(); iphoton++ )	{
    narrived++;

    double PhotonEnergy = penergy[tids[iphoton]]/eV;
    G4ThreeVector LocalPosition = lposition[tids[iphoton]];
    double PhotonTime = ptime[tids[iphoton]];
    if (PhotonTime < firstphotontime) firstphotontime = PhotonTime;
    //cout << "  t1=" << firstphotontime << endl;

    /* Checking the dead area of the pixel */
    bool _iDeadSpace = false ;
    const double new_xpos = -LocalPosition.y() + richc.PhCsize_x/2. ;
    const double new_ypos =  LocalPosition.y() + richc.PhCsize_y/2. ;

    /* Verify dead area along x */
    double_t _x_progress = 0.0 ;
    for( int _ir = 0 ; _ir <richc.nrows ; _ir++ ){
      if( new_xpos > _x_progress && new_xpos < _x_progress + richc.PXDX[ _ir ] ) {
	/* position of the hit in the pixel */
	double _x_pos_in_the_pixel = new_xpos - _x_progress ;
	
	/* verify it is not in the deadspace */
	if( ( _x_pos_in_the_pixel > ( richc.PXDX[ _ir ] - richc.PixelDeadSpace_x/2.) ) || ( _x_pos_in_the_pixel <  richc.PixelDeadSpace_x/2. ) ) _iDeadSpace = true ;
	
	break ;
	
      }
      _x_progress += richc.PXDX[ _ir ] ;
    }


    /* Verify dead area along y if not already dead */
    if (!_iDeadSpace) {
      
      double_t _y_progress = 0.0 ;
      for( int _ic = 0 ; _ic < richc.ncols; _ic++ ){
	
	if( new_ypos > _y_progress && new_ypos < _y_progress + richc.PXDY[ _ic ] ){
	  /* position of the hit in the pixel */
	  double _y_pos_in_the_pixel = new_ypos - _y_progress ;
	  
	  /* verify it is not in the deadspace */
	  if( ( _y_pos_in_the_pixel > ( richc.PXDY[ _ic ] - richc.PixelDeadSpace_y/2.) ) || ( _y_pos_in_the_pixel <  richc.PixelDeadSpace_y/2. ) ) _iDeadSpace = true ;
	  
	  break ;
	  
	}
	
	_y_progress += richc.PXDY[ _ic ] ;
	
      }
    }

  





    if (_iDeadSpace) {
      ndeadarea++;
    }
    else {/* photon not in dead area, trying to detect it */
    
    
      /* Checking QE for normal glass window MAPMT */
      if ( (PhotonEnergy < richc.energy_H8500[richc.nH8500-1]) || (PhotonEnergy > richc.energy_H8500[0]) ) {
	/* Photon energy outside the QE range */
	noutofrange++;
      }
      else {
	
	/* QE test */
	double QE = -1;
	for(int ie=0; ie<richc.nH8500-1; ie++) {
	  if(PhotonEnergy < richc.energy_H8500[ie] && PhotonEnergy >= richc.energy_H8500[ie+1]) {
	    QE = richc.QE_H8500[ie] + (PhotonEnergy - richc.energy_H8500[ie]) * 
	      (richc.QE_H8500[ie+1] - richc.QE_H8500[ie]) / (richc.energy_H8500[ie+1] - richc.energy_H8500[ie]);
	    break;
	  }
	}
	
	double QEtest = G4UniformRand(); 
	if(QEtest > QE) {
	  nunconverted++;
	}
	else {
	  ndetected++;
	}
	
      }


    }


  }


  /* ============================================================ */
  /* flagging the hits with nodetected photons */
  if (!ndetected) {

    if (ndeadarea == narrived) {
      /* all in dead area */
      pixel = -pixel;
      //cout << "--->> changing sign" << endl;
    }
    else if (nunconverted == narrived) {
      /* all not converted */
      pixel = pixel + 100;
    }
    else if (noutofrange == narrived) {
      /* all not out of range */
      pixel = pixel + 200;
    }
    else pixel = pixel + 300;

  }

  if (pixel<0 && (ndetected==1) ) {
    //cout << "pixel=" << pixel << "  nA=" << narrived << "  nD=" << ndetected << "  nO=" << noutofrange << "  nU=" << nunconverted << "  nK=" << ndeadarea << endl;
    }


  /******************************************************/
  /* Generating the MAPMT charge response
     Using one single mean Gain value for all the pixels
     Individual gain readout from CCDB must be implemented
  */
  double charge = GetCharge(ndetected, richc.Gain);

  /* MAROC charge response, for the moment is just a conversion from fC to DAC units */
  int ADC = Charge2ADC(charge);


  /*****************************************************/
  /* Generating the MAPMT time response */
  vector<int> TDC = Time2TDC(firstphotontime, charge);



  /*****************************************************/
  /* Output of digitized info */

  dgtz["sector"] = sector;
  dgtz["pmt"]    = pmt;
  dgtz["pixel"]  = pixel;
  dgtz["hitn"]   = hitn;
  dgtz["nphotons"]   = narrived;
  dgtz["npe"]   = ndetected;
  dgtz["nphotonsU"]   = nunconverted;
  dgtz["nphotonsO"]   = noutofrange;
  dgtz["nphotonsD"]   = ndeadarea;
  dgtz["TDC1"]   = TDC[0];
  dgtz["TDC2"]   = TDC[1];
  dgtz["ADC"]   = ADC;

  return dgtz;
}

/* ============================================================================================== */

vector<identifier> rich_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
  /* This function is called at each step of the particles entering the RICH sensitive volume
     (-> MAPMT photocathode without pixelization) 
     The pixelization MUST be here in order to add all the steps belonging to the same pixel 

     The check of the dead area and of the QE is in integrateDgt

  */

  /**********************************************/
  /* This part is copied from Ahmed */
  /**********************************************/


  /* Getting the photon information */
  G4StepPoint *prestep = aStep->GetPreStepPoint();
  G4StepPoint *poststep = aStep->GetPostStepPoint();
  G4ThreeVector  xyz   = poststep->GetPosition();
  G4ThreeVector Lxyz = prestep->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz);


  /* MC x axis is mirrored with respect the Lab reference system
     shift the local reference center to up left corner (following configuration files) */
  const double new_xpos = -Lxyz.x() + richc.PhCsize_x/2. ;
  const double new_ypos =  Lxyz.y() + richc.PhCsize_y/2. ;

  /*********************************************/
  /* PIXELIZATION ROUTINE FROM AHMED */
  /* It seems that rows and columns are switched, 
     but they are switched back in the final calculation */
  /*********************************************/

  
  /* ---------------------------------- */
  /* find the pixel row */
  int _irow = -1 ;

  if (new_xpos == richc.PhCsize_x) {
    _irow = richc.nrows;
  }
  else {
    double_t _x_progress = 0.0 ;
    for( int _ir = 0 ; _ir < richc.nrows; _ir++ ){

      if( new_xpos >= _x_progress && new_xpos < _x_progress + richc.PXDX[ _ir ] ){
	
	// cout << " x new pos is " << new_xpos << " and x_progress is " << _x_progress << " for ir = " << _ir << endl ;
	_irow = _ir + 1 ;
	
	break ;
      }
      
      _x_progress += richc.PXDX[ _ir ];
    }

  }



  /* ---------------------------------- */
  /* find the pixel column */
  int _icol = -1 ;

  if (new_ypos == richc.PhCsize_y) {
    _icol = 1;
  }
  else {
    double_t _y_progress = 0.0 ;
    for( int _ic = 0 ; _ic < richc.ncols; _ic++ ){

      if( new_ypos >= _y_progress && new_ypos < _y_progress + richc.PXDY[ _ic ] ){
	
	// cout << " y new pos is " << new_ypos << " and y_progress is " << _y_progress << " for ir = " << _ic << endl ;
	_icol = richc.ncols - _ic ; 
	
	break ;
	
      }

      _y_progress += richc.PXDY[ _ic ] ;
      
    }
  }

  
  
  /* ---------------------------------------------- */
  /* final pixel number */
  int pixel = richc.nrows*( _icol - 1 ) + _irow ;

  /* correcting for photon hits out of the photocathode */
  if (pixel < 0) {
    //fprintf(stdout, "processID:  x=%12.9f   y=%12.9f   c=%d  r=%d  pixel=%d\n", new_xpos, new_ypos, _icol, _irow, pixel);
    pixel = 0;
  }
  id[2].id = pixel;
 
  



  /* ---------------------------------------------- */
  id[id.size()-1].id_sharing = 1;
  
  return id;
}


/* ============================================================================================== */

map< string, vector <int> >  rich_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;


	return MH;
}


void rich_HitProcess::initWithRunNumber(int runno)
{
	if(richc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		richc = initializeRICHConstants(runno);
		richc.runNo = runno;
	}
}


/* ============================================================================================== */
// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> rich_HitProcess :: electronicNoise()
{
  /* The noise is generated uniformly on all the MAPMTS using average values of gain and dark currents.
     To have individual pmt dark counts, a loop over all the MAPMT must be done,
     using the MAPMT gain and dark current values 
  */

  /* Average dark noise rate per MAPMT (e is in fC) */
  double DarkNoiseRate1 = 1e15 * richc.DarkCurrent / (richc.Gain * ElectronCharge);

  /* Average dark noise hits per MAPMT (time window is in ns) */
  double AverageDarkNoiseCounts1 = 1e-9 * DarkNoiseRate1 * richc.TimeWindow;

  /* Expected total dark counts */
  double AverageDarkNoiseCounts = AverageDarkNoiseCounts1 * richc.nMAPMT;


  /* Generating the noise counts */
  unsigned int NoiseNHits = G4Poisson(AverageDarkNoiseCounts);
  //NoiseNHits = 1;
  //cout << "----->>> ELECTRONIC NOISE: <N>=" << AverageDarkNoiseCounts << "  N=" << NoiseNHits << endl;
  

  vector<MHit*> noiseHits;

  // first, identify the cells that would have electronic noise
  // then instantiate hit with energy E, time T, identifier IDF:
  //
  // MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

  // push to noiseHits collection:
  // noiseHits.push_back(thisNoiseHit)
  
  
  /***************************************/
  /* TEST */
  for(unsigned int p=0; p<NoiseNHits; p++) {
    //cout << "   p=" << p << endl;
	  
    vector<identifier> thisID;
	  
    identifier thisIdentifier;

    /* Sector */
    thisIdentifier.id = 4;
    thisIdentifier.name = "richNoise";
    thisID.push_back(thisIdentifier);
	  
    /* PMT */
    int pmtID = richc.nMAPMT * G4UniformRand();
    thisIdentifier.id = pmtID;
    thisIdentifier.name = "richNoise";
    thisID.push_back(thisIdentifier);
	  
    /* Pixel*/
    int pixelID = 1 + richc.nrows*richc.ncols * G4UniformRand();
    thisIdentifier.id = pixelID;
    thisIdentifier.name = "richNoise";
    thisID.push_back(thisIdentifier);
	  

    //cout << " thisID filled nID=" << thisID.size() << endl;



    double energy = 1;
    double time = richc.TimeWindow * G4UniformRand();
    MHit* thisNoiseHit = new MHit(energy, time, thisID, p);
    noiseHits.push_back(thisNoiseHit);
    
    //cout << " noise hit added" << endl;
  }
  /**************************************/


  //cout << "-----> END electronicNoise " << endl;


  return noiseHits;
}

/* ============================================================================================== */


// - charge: returns charge/time digitized information / step
map< int, vector <double> > rich_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The input are charge value, time
double rich_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}


/* ============================================================================================== */

// - GetCharge: generate the charge collected by the MAPMT. Inputs are number of p.e., Gain
double rich_HitProcess :: GetCharge(int npe, double Gain)
{

  /* Returns charge collected on the MAPMT anode in fC

     Assuming the sigma of the charge is proportional to the fluctuations 
     on the first dinode. The H12700 has 10 stages, thus:
     G = g^10
     where g is tha gain of the first stage. 
     g = G^(1/10)
     sigma_g = sqrt(g)
     sigma_G/G = sigma_g/g = 1/sqrt(g) = 1/sqrt(G^(1/10)) 
  */

  double charge = 0;
  if (npe) {
    double Gain = 3.0e6;
    int nstages = 10;
    double gain = pow(Gain, (1./nstages));

    double ne1mean = gain * npe;
    double nemean = Gain * npe;

    double nesigma = nemean / sqrt(ne1mean);
    
    double necollected = G4RandGauss::shoot(nemean, nesigma);
    charge = ElectronCharge * necollected;

  }

  return charge;
}
/* ---------------------------------------------------------------- */
// - Charge2ADC: translate the MAPMT charge to MAROC DAC units 
int rich_HitProcess :: Charge2ADC(double charge)
{
  /* Generating the MAROC response, for the moment is just a conversion */
  double dac2charge = 0.82;//from DAC units to fC
  /* arbitrary ADC pedestal to avoid large events with ADC<0 */
  double AdcPedestal = 200;

  double ADC = AdcPedestal + charge * dac2charge;
  if (ADC > 4095) ADC = 4095;
  if (ADC < 0) ADC = 0;

  return ADC;
}
/* ---------------------------------------------------------------- */
vector<int> rich_HitProcess :: Time2TDC(double time, double charge)
{
  /* 
     For the moment there is no conversion to TDC counts, the output is time in ns

     Conversion from charge to time:
     From Matteo, rich meeting 6/nov/2015
     With 5 pe signal (Q~2100 fC) there is a 70 ns difference between
     fall and raise time
     
     Fall Time is taken from the true time
     Rise time is computed as 
     Tr = Tf + cost*charge
     cost=70/2100 ~0.03 ns/fC 
  */
  double charge2time = 0.03;

  /* RMS of the time resolution */
  double_t timesigma = 1;

  vector<int> TDC;
  int FallTime = 0;
  int RiseTime = 0;

  if (charge) {
    double time0;

    /* Fall time -> TDC1 */
    time0 = G4RandGauss::shoot(time, timesigma);
    FallTime = (int)time0;

    /* Rise time -> TDC2 */
    time0 = G4RandGauss::shoot(time + charge2time * charge, timesigma);
    RiseTime = (int)time0;

  }
  //cout << "Q=" << charge << "  t1=" << FallTime << "  t2=" << RiseTime << endl;

  TDC.push_back(FallTime);
  TDC.push_back(RiseTime);

  return TDC;
}


/* ============================================================================================== */


// this static function will be loaded first thing by the executable
richConstants rich_HitProcess::richc = initializeRICHConstants(-1);





