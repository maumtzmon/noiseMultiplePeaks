//-----This script uses the root file (*_hepeaks.root) generated from hepeaks_*.C
//-----and performs a fit to the electron peaks
#include <string.h>

using namespace std;
#define ROWNUM 1 // 4
#define COLNUM 16 // 4
#include <math.h>



double singleGauss(double *x, double *par){
  double q0 	 =	x[0];
	double norm0peak  =  par[0];
	double mu0    =  par[1];
	double sigma0 =  par[2];
	double val=0;

	double cPi = TMath::Pi();

	for(int w=0; w<1; w++){
		val+=(norm0peak/sqrt(2*cPi*pow(sigma0,2)))*(exp(-0.5*pow(q0-mu0,2)/pow(sigma0,2)));
	}
	return val;
}

double twoGaussians(double *x, double *par){
	double q 	 =	x[0];
	double norm  =  par[0];
	double mu    =  par[1];
	double sigma =  par[2];
	double gain  =  par[3];
	double con   =  par[4];
	double val=0;

	double cPi = TMath::Pi();

	for(int w=0; w<2; w++){
		val+=(norm/sqrt(2*cPi*pow(sigma,2)))*(exp(-0.5*pow(q-mu-w*gain,2)/pow(sigma,2)));
	}
	val+=con;
	return val;
}

double multGaussians(double *x, double *par){
	double q 	 =	x[0];
	double norm  =  par[0];
	double mu    =  par[1];
	double sigma =  par[2];
	double gain  =  par[3];
	double con   =  par[4];
	double val=0;

	double cPi = TMath::Pi();

	for(int w=0; w<3; w++){
		val+=(norm/sqrt(2*cPi*pow(sigma,2)))*(exp(-0.5*pow(q-mu-w*gain,2)/pow(sigma,2)));
	}
	val+=con;
	return val;
}




//--------VARIABLES TO FIT HISTOGRAMS--------
                    //41 para probar solo con el pico 40y 41
//const int numpeaks = 740;		// Number of peaks to fit starting from the 0e peak (numpeaks=2 fits the 0e and the 1e peak)
const int numpeaks = 800;
double xPeak[numpeaks-1], exPeak[numpeaks-1]; 
const int numext = 16;		// Number of working extensions
double gainPeak[numext][numpeaks-1]; // ganancia medida entre picos  gainPeak[0] = peak2mean - peak1mean
double egainPeak[numext][numpeaks-1];
double gainDiv[numext][numpeaks-1];
double egainDiv[numext][numpeaks-1];
double meanPeak[numext][numpeaks-1];
double emeanPeak[numext][numpeaks-1];
double peakGain[numext][numpeaks-1];
double epeakGain[numext][numpeaks-1];
double peakGainDiv[numext][numpeaks-1];

int goodext[numext] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
// float expgain[numext] = {200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};			// Expected gain in ADU/e-

const int fitopt = 0;			// 1 for fitting with 2 gaussians, 2 for fitting with convolution  
float nsmp = 300;
float noise1smp = 5;		// Expected noise for 1 sample

int peakLimit_1=30;
int peakLimit_2=200;

int doFit=0; // 1 do fit, 0 doesn't



void noisedcmm_3peaks_combined (char const* file, char const* file2){
  // File to save output
  std::ofstream outputFile("output_data.txt");
  // --------Retrieve expgain from root file
  TFile *input = new TFile(file2);
  TTree *tree = (TTree*)input->Get("tree");
  float Gain;
  float expgain[16];
  float expgain2;
  tree->SetBranchAddress("Gain",&expgain2);
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i); 
    expgain[i]=expgain2;
  }
  input->Close();

//--------Style--------

  gROOT->Reset();

  TGaxis::SetMaxDigits(3);

//--------Retrieve filename without extension--------

  int length = strlen(file);
  char fileroot[length+1];
  strcpy(fileroot, file);		// Copy the input string to fileroot
  fileroot[length-5] = '\0';		// Throw ".root" from filename (last 5 characters)



//--------Retrieve histograms from the root file--------

  TFile filehist(Form("%s", file));

  TH1F *hpix[numext];
  for (int next=5; next<6; next++){

	hpix[next] = (TH1F*)filehist.Get(Form("ext%i", goodext[next]));
	hpix[next]->SetDirectory(0);
  hpix[next]->SetAxisRange(0.0,5.0);

	}
	
  filehist.Close();
	
//--------Find peaks to fit and fit them--------

  float pixmin, pixmax, xmin, xmax, norm0, offset0, sigma0;
  int nbin, binmin, binmax;
  float fitrange, xminpeak, xmaxpeak, offsetpeak, offsetaux, offsetfit, gainfit;
  int binminpeak, binmaxpeak;
  double chisquare;

  TF1 *fitfun[numpeaks];
  int fitfunstatus[numext][numpeaks], iniPeak=500;
  double pars[numpeaks*5], epars[numpeaks*5], sigma[numext], esigma[numext];
  double pars0[numext][numpeaks], epars0[numext][numpeaks], pars1[numext][numpeaks], epars1[numext][numpeaks], pars2[numext][numpeaks], epars2[numext][numpeaks], xgraph[numpeaks], exgraph[numpeaks]={0};
  double cociente[numpeaks-1];

  TF1 *f2gaus[numext];
  double parsf2gaus[6], pval[numext][numpeaks];

  int fitstatus[numext], reps=2, ndf, tempFit; // reps, si va a hacer una o varias iteraciones para mejorar la calidad del ajuste p/gaussiana
  double gain[numext],offset[numext], lambda[numext], noise[numext], norm[numext];
  double egain[numext],eoffset[numext], elambda[numext], enoise[numext], enorm[numext];
  double gain_fit[numext], noise_fit[numext], lambda_fit[numext], offset_fit[numext], norm_fit[numext];
  double egain_fit[numext], enoise_fit[numext], elambda_fit[numext], eoffset_fit[numext], enorm_fit[numext];

  float mean1, mean_1A, mean2, mean_2A, gain_1, sigma1,const1;

  offsetpeak=15000;
  iniPeak=0;
  int loq=iniPeak*250-100,hiq=iniPeak*250+100;//15900; 
  
  int maxq = iniPeak*250+100; //loq+300; 
  fitrange=1.5*250;

  // ---------------------------------------------
  // Data Extraction, iteracion por cada extension
  for (int next=5; next<6; next++){			// Loop of exts

    //histogram labels
  hpix[next]->SetTitle(Form("ext%i", goodext[next]+1));
  hpix[next]->GetXaxis()->SetTitle("Pixel value [ADU]");
  hpix[next]->GetYaxis()->SetTitle("Number of counts");

  //valores iniciales para el primer ajuste
  binmin = hpix[next]->GetXaxis()->FindBin(loq); 
  binmax = hpix[next]->GetXaxis()->FindBin(maxq);
  offsetaux=loq+fitrange;
  maxq = hpix[next]->GetXaxis()->GetXmax(); 

  printf("\n Peak 0 \n loq = %i, hiq = %i, maxq = %i\n", loq, hiq, maxq);


  //loq=offsetaux-0.5*expgain[next];
  binmin = hpix[next]->GetXaxis()->FindBin(loq); 
  binmax = hpix[next]->GetXaxis()->FindBin(maxq);
  hpix[next]->GetXaxis()->SetRange(binmin, binmax);

	// nbin = (hpix[next]->GetSize())-2;
	


	//--------Gaussian fit of e- peaks--------
  float coef=1;
  int Npar = 5;
  
 
  mean1=loq+80;
  norm0=1;
  gainfit=250;
  const1=4;
  sigma1=60;
  int numGaussians=1;
  
  //double xPeak[numpeaks-1], exPeak[numpeaks-1]; 
  
  // Iteracion para cada pico dentro de la extension

  for (int npeak=iniPeak, increment=numGaussians; npeak<numpeaks; npeak+=numGaussians){ //numpeaks; npeak++){
    exPeak[npeak]=0;
    xPeak[npeak]=npeak;
    
    if (npeak==0){ // fit on 0e- peak
      printf("\n\n Peak %i \n loq = %i, hiq = %i, Expected mean = %f\n\n", npeak, loq-100, hiq+70, mean1);
      fitfun[npeak] = new TF1("fitfun", singleGauss, loq-80, hiq+30 ,3); // Npar = numero de parametros de la funcion multiGaussians necesita 5 parametros  
      
      // set fit parameters

      fitfun[npeak]->SetParameter(0,norm0);        //norm
      fitfun[npeak]->SetParameter(1,0);    //mu
      fitfun[npeak]->SetParameter(2,sigma1);       //Sigma


      // do the fit
    
      hpix[next]->Fit("fitfun","R+");
      
  
      // get fit parameters
      norm0= fitfun[npeak]->GetParameter(0);
      mean1 = fitfun[npeak]->GetParameter(1);
      sigma1= fitfun[npeak]->GetParameter(2);

      fitfun[npeak]->GetParameters(&pars[(npeak+1)*3]);
      ndf=fitfun[npeak]->GetNDF();
      chisquare = fitfun[npeak]->GetChisquare();
      cociente[npeak]=chisquare/ndf;

      if (abs(sigma1)>100){
        loq=mean1+250-50;
        mean1=mean1+250;
        hiq=mean1+50;
      }
      else{
        loq=mean1+250-1.5*abs(sigma1);
        mean1=mean1+250;
        hiq=mean1+1.5*abs(sigma1);
      }
      printf("chisquare=%f; ndf=%i, chisquare/ndf=%f \n\n", chisquare,ndf,cociente[npeak]);
     
      // guardar el valor de la media
      meanPeak[next][npeak]=fitfun[npeak]->GetParameter(1);
      emeanPeak[next][npeak]=fitfun[npeak]->GetParError(1);

    }

    else{ // fit from 1e- to 20e- Peak
      if (0 < npeak && npeak <= peakLimit_1){
        printf("\n\n Peak %i \n loq = %i, hiq = %i, Expected mean = %f\n\n", npeak, loq, hiq, mean1);
        fitfun[npeak] = new TF1("fitfun", singleGauss, loq, hiq ,3); // Npar = numero de parametros de la funcion multiGaussians necesita 5 parametros  
        
        // set fit parameters

        fitfun[npeak]->SetParameter(0,norm0);        //norm
        fitfun[npeak]->SetParameter(1,mean1);    //mu
        fitfun[npeak]->SetParameter(2,sigma1);       //Sigma


        // do the fit
      
        hpix[next]->Fit("fitfun","R+");
        
    
        // get fit parameters
        norm0= fitfun[npeak]->GetParameter(0);
        mean1 = fitfun[npeak]->GetParameter(1);
        sigma1= fitfun[npeak]->GetParameter(2);

        fitfun[npeak]->GetParameters(&pars[(npeak+1)*3]);
        if (abs(sigma1)>100){
          loq=mean1+250-50;
          mean1=mean1+250;
          hiq=mean1+50;
        }
        else{
          loq=mean1+250-1.5*abs(sigma1);
          mean1=mean1+250;
          hiq=mean1+1.5*abs(sigma1);
        }

        chisquare = fitfun[npeak]->GetChisquare();
        ndf=fitfun[npeak]->GetNDF();
        pval[next][npeak]=TMath::Prob(chisquare,ndf);

        cociente[npeak]=chisquare/ndf;

        printf("\nchisquare=%f; ndf=%i, chisquare/ndf=%f ", chisquare,ndf,cociente[npeak]);

  

        // outputFile << "\n Peak " << npeak << "loq = "<< loq << "; " << "hiq=  " << hiq << "; " << "mean =" << mean1 << "; " << "gain ="<< gainfit << "; " << "const = " << const1 << "; "<< "xPeak= "<< xPeak[npeak]<< std::endl;

        // guardar el valor de la media
        meanPeak[next][npeak]=fitfun[npeak]->GetParameter(1);
        emeanPeak[next][npeak]=fitfun[npeak]->GetParError(1);

        //guardar ganancia y errores cuando se calcula en cada pico
        // gain=media actal - media anterior ; error Gain = raiz(cuadrado(error mediaa ctual)-cuadrado(error media previa))
        gainPeak[next][npeak]=meanPeak[next][npeak]-meanPeak[next][npeak-1];
        egainPeak[next][npeak]=sqrt(pow(emeanPeak[next][npeak],2)+pow(emeanPeak[next][npeak-1],2));
        
        gainDiv[next][npeak]=(gainPeak[next][npeak]-245)/245;
        egainDiv[next][npeak]=egainPeak[next][npeak]/245;


        peakGain[next][npeak]=meanPeak[next][npeak]/npeak;
        peakGainDiv[next][npeak]=peakGain[next][npeak]/240;


      }

      if (peakLimit_1 < npeak && npeak <= peakLimit_2){
        printf("\n\n Peak %i \n loq = %i, hiq = %i, mean = %f, offsetaux = %f, xPeak= %f\n\n", npeak, loq, hiq, mean1, offsetaux, xPeak[npeak]);
        if (npeak==peakLimit_1+1){
          numGaussians+=1;
          gainfit=240;
          const1=40;
          hiq=loq+gainfit+3.2*abs(sigma1);
          loq=loq-0.7*abs(sigma1);
          printf("\n\n initial parameters for 2 gaussians in Peak %i \n loq = %i, hiq = %i, mean = %f, sigma= %f, \n initial gain= %f, initial const= %f \n increments= %i \n", npeak, loq, hiq, mean1, sigma1, gainfit,const1, numGaussians);

        }

        fitfun[npeak] = new TF1("fitfun", twoGaussians, loq, hiq ,Npar); // Npar = numero de parametros de la funcion multiGaussians necesita 5 parametros  
        
        // set fit parameters

        fitfun[npeak]->SetParameter(0,norm0);        //norm
        fitfun[npeak]->SetParameter(1,mean1);    //mu
        fitfun[npeak]->SetParameter(2,sigma1);       //Sigma
        fitfun[npeak]->SetParameter(3,gainfit);      //gain
        fitfun[npeak]->SetParameter(4,const1);        //con, fondo en el que estan montados las gaussianas

        // do the fit
      
        hpix[next]->Fit("fitfun","R+");
        
    
        // get fit parameters
        norm0= fitfun[npeak]->GetParameter(0);
        mean1 = fitfun[npeak]->GetParameter(1);
        sigma1= fitfun[npeak]->GetParameter(2);
        gainfit=fitfun[npeak]->GetParameter(3);
        const1=fitfun[npeak]->GetParameter(4);

        fitfun[npeak]->GetParameters(&pars[npeak*5]);
          
        if (abs(sigma1)>100){
          sigma1=60;
          mean1=mean1+2*gainfit;
          loq=mean1-0.4*gainfit;
          
          hiq=mean1+gainfit+0.4*gainfit;
        }
        else{
          mean1=mean1+2*gainfit;
          loq=mean1-0.4*gainfit;
          
          hiq=mean1+gainfit+0.4*gainfit;
          // loq=mean1+250-1.5*abs(sigma1);
          // mean1=mean1+250;
          // hiq=loq+gainfit+3.2*abs(sigma1);
        }

        //guardar ganancia y errores cuando se calcula en cada pico
        gainPeak[next][npeak]=fitfun[npeak]->GetParameter(3);
        egainPeak[next][npeak]=fitfun[npeak]->GetParError(3);

        gainDiv[next][npeak]=(gainPeak[next][npeak]-245)/245;
        egainDiv[next][npeak]=egainPeak[next][npeak]/245;

        // guardar el valor de la media
        meanPeak[next][npeak]=fitfun[npeak]->GetParameter(1);
        emeanPeak[next][npeak]=fitfun[npeak]->GetParError(1);

        peakGain[next][npeak]=meanPeak[next][npeak]/npeak;
        peakGainDiv[next][npeak]=peakGain[next][npeak]/240;

        chisquare = fitfun[npeak]->GetChisquare();
        ndf=fitfun[npeak]->GetNDF();
        pval[next][npeak]=TMath::Prob(chisquare,ndf);

        cociente[npeak]=chisquare/ndf;
        
        printf("\nchisquare=%f; ndf=%i, chisquare/ndf=%f \n\n", chisquare,ndf,cociente[npeak]);
        


        outputFile << "\n Peak " << npeak << "loq = "<< loq << "; " << "hiq=  " << hiq << "; " << "mean =" << mean1 << "; " << "gain ="<< gainfit << "; " << "const = " << const1 << "; "<< "xPeak= "<< xPeak[npeak]<< std::endl;
      
        //guardar ganancia y errores cuando se calcula en cada pico
        gainPeak[next][npeak]=fitfun[npeak]->GetParameter(3);
        egainPeak[next][npeak]=fitfun[npeak]->GetParError(3);

      }
      
      if (peakLimit_2 < npeak && npeak <= 780){
        printf("\n\n Peak %i \n loq = %i, hiq = %i, mean = %f, offsetaux = %f, xPeak= %f\n Actual increment=%i\n", npeak, loq, hiq, mean1, offsetaux, xPeak[npeak], numGaussians);
        if (npeak==peakLimit_2+1){
            numGaussians+=1;
            gainfit=240;
            const1=5;
            hiq=loq+3*gainfit;//+3.2*abs(sigma1);
            loq=loq-0.7*abs(sigma1);
            //printf("\n\n initial parameters for %i gaussians in Peak %i \n loq = %i, hiq = %f, mean = %f, sigma= %f, \n initial gain= %f, initial const= %f \n\n", numGaussians, loq, hiq, mean1, sigma1, gainfit,const1);

        }
        fitfun[npeak] = new TF1("fitfun", multGaussians, loq, hiq ,Npar); // Npar = numero de parametros de la funcion multiGaussians necesita 5 parametros  
          
        // set fit parameters

        fitfun[npeak]->SetParameter(0,norm0);        //norm
        fitfun[npeak]->SetParameter(1,mean1);    //mu
        fitfun[npeak]->SetParameter(2,sigma1);       //Sigma
        fitfun[npeak]->SetParameter(3,gainfit);      //gain
        fitfun[npeak]->SetParameter(4,const1);        //con, fondo en el que estan montados las gaussianas

        // do the fit
      
        hpix[next]->Fit("fitfun","R+");
        
    
        // get fit parameters
        norm0= fitfun[npeak]->GetParameter(0);
        mean1 = fitfun[npeak]->GetParameter(1);
        sigma1= fitfun[npeak]->GetParameter(2);
        if (sigma1<0){sigma1=abs(sigma1);}
        gainfit=fitfun[npeak]->GetParameter(3);
        const1=fitfun[npeak]->GetParameter(4);
        if (const1<0){const1=abs(const1);}


        fitfun[npeak]->GetParameters(&pars[npeak*5]);
          
        if (gainfit<200){
          gainfit=230;
        }

        if (abs(sigma1)>100){
          sigma1=60;
          mean1=mean1+3*gainfit;
          loq=mean1-0.6*gainfit;
          
          hiq=mean1+2*gainfit+0.6*gainfit;
        }
        
        else{
          mean1=mean1+3*gainfit;
          loq=mean1-0.6*gainfit;
          
          hiq=mean1+2*gainfit+0.6*gainfit;
          // loq=mean1+250-1.5*abs(sigma1);
          // mean1=mean1+250;
          // hiq=loq+gainfit+3.2*abs(sigma1);
        }

        

        chisquare = fitfun[npeak]->GetChisquare();
        ndf=fitfun[npeak]->GetNDF();
        pval[next][npeak]=TMath::Prob(chisquare,ndf);

        cociente[npeak]=chisquare/ndf;
        
        printf("\nchisquare=%f; ndf=%i, chisquare/ndf=%f \n\n", chisquare,ndf,cociente[npeak]);
        

        
        //outputFile << "\n Peak " << npeak << "loq = "<< loq << "; " << "hiq=  " << hiq << "; " << "mean =" << mean1 << "; " << "gain ="<< gainfit << "; " << "const = " << const1 << "; "<< "xPeak= "<< xPeak[npeak]<< std::endl;
      
        //guardar ganancia y errores cuando se calcula en cada pico
        gainPeak[next][npeak]=fitfun[npeak]->GetParameter(3);
        egainPeak[next][npeak]=fitfun[npeak]->GetParError(3);

        gainDiv[next][npeak]=(gainPeak[next][npeak]-245)/245;
        egainDiv[next][npeak]=egainPeak[next][npeak]/245;

                // guardar el valor de la media
        meanPeak[next][npeak]=fitfun[npeak]->GetParameter(1);
        emeanPeak[next][npeak]=fitfun[npeak]->GetParError(1);

        peakGain[next][npeak]=meanPeak[next][npeak]/npeak;
        peakGainDiv[next][npeak]=peakGain[next][npeak]/240;        

      }

      else{
        printf("peak %i\n", npeak);
      }
  }

  }


  // plot Histogram and fitting

  TCanvas *c0 = new TCanvas("c0","Draw fitted histograms", 2000, 500);//*ceil(numext/4.));
	gPad->SetLogy();
  hpix[next]->Draw("same");
  }    

  // ---------------------------------------------
  // Plot Section
  double fitPar[numext*2], efitPar[numext*2];

  // ---------------------------------------------
  //Plot gain vs Peak
  printf("\nGain behavior across the peaks\n");
  TGraphErrors *egraph[numext];
  TCanvas *c1 = new TCanvas("c1","Draw gain Behavior", 2000, 500);//*ceil(numext/4.)); //propiedades del canvas

  TF1 *ajusteLin[numext];

  float sizeXpeak=sizeof(xPeak);
  
  int n=numpeaks-iniPeak;
  printf("\n\n numpeaks= %i; iniPeak= %i;xpeak len= %f, n puntos=%i \n\n", numpeaks, iniPeak, sizeXpeak, n);

  for (int next=5; next<6; next++){

    egraph[next] = new TGraphErrors(n, xPeak, gainPeak[next],exPeak, egainPeak[next] ); // Grafica de cada extension
    egraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    egraph[next]->SetMinimum(210);
    egraph[next]->SetMaximum(280);

    ajusteLin[next]=new TF1("f", "[0]+[1]*x",0,numpeaks);  //Fit linear model
    ajusteLin[next]->SetParameters(200,0);
    

    if(doFit == 1){
      egraph[next]->Fit(ajusteLin[next], "R+");

      ajusteLin[next]->GetParameters(&fitPar[next*2]);
      efitPar[0+next*2]=ajusteLin[next]->GetParError(1);
      printf("\nmean: %f +/- %f", fitPar[1+2*next], efitPar[1+2*next]);
    }

    egraph[next]->Draw("AP"); 
   

  }

  // ---------------------------------------------
  //plot Chisquare/ndf
  TCanvas *c2 = new TCanvas("c2","chisquare/ndf", 2000, 500);
  TGraph *pValGraph[numext];  

  for (int next=5; next<6; next++){

    pValGraph[next] = new TGraph(n, xPeak, cociente); // Draw  cociente[npeak]=chisquare/ndf;
    pValGraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    pValGraph[next]->SetMinimum(0);
    pValGraph[next]->SetMaximum(25);
    pValGraph[next]->SetLineColor(4);
    pValGraph[next]->SetMarkerStyle(20); // Filled circle
    pValGraph[next]->SetMarkerSize(1.0);
    pValGraph[next]->Draw("AP");

  }

  
  // ---------------------------------------------
  //plot Mean vs Peak
  printf("\nmean behavior across the peaks\n");
  //TGraphErrors *egraph[numext];
  TCanvas *c3 = new TCanvas("c3","Draw mean Behavior", 2000, 500);//*ceil(numext/4.)); //propiedades del canvas
  for (int next=5; next<6; next++){
    float sizeXpeak=sizeof(xPeak);
  
 
    egraph[next] = new TGraphErrors(numpeaks, xPeak, meanPeak[next],exPeak, emeanPeak[next] ); // Grafica de cada extension
    egraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    egraph[next]->Draw("AP"); 
  }

    // ---------------------------------------------
  //Plot gain vs Peak 2
  printf("\nGain behavior across the peaks\n");
  TCanvas *c4 = new TCanvas("c4","Draw gain Behavior", 2000, 500);//*ceil(numext/4.)); //propiedades del canvas

  
  printf("\n\n numpeaks= %i; iniPeak= %i;xpeak len= %f, n puntos=%i \n\n", numpeaks, iniPeak, sizeXpeak, n);

  for (int next=5; next<6; next++){

    egraph[next] = new TGraphErrors(n, xPeak, gainDiv[next],exPeak, egainDiv[next] ); // Grafica de cada extension
    egraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    egraph[next]->SetMarkerStyle(8);
    egraph[next]->SetMarkerSize(0.5);
    // egraph[next]->SetMinimum(210);
    // egraph[next]->SetMaximum(280);

    ajusteLin[next]=new TF1("f", "[0]+[1]*x",0,numpeaks);  //Fit linear model
    ajusteLin[next]->SetParameters(200,0);
    

    
    if(doFit == 1){
      egraph[next]->Fit(ajusteLin[next], "R+");
      ajusteLin[next]->GetParameters(&fitPar[next*2]);
      efitPar[0+next*2]=ajusteLin[next]->GetParError(1);
      printf("\nmean: %f +/- %f", fitPar[1+2*next], efitPar[1+2*next]);
    }
    egraph[next]->Draw("AP"); 

  }


  // ---------------------------------------------
  //Plot gain vs Peak 3
  printf("\nmean/peakNumber \n");
  TCanvas *c5 = new TCanvas("c5","Gain=mean/peak", 2000, 500);//*ceil(numext/4.)); //propiedades del canvas

  
  printf("\n\n numpeaks= %i; iniPeak= %i;xpeak len= %f, n puntos=%i \n\n", numpeaks, iniPeak, sizeXpeak, n);

  for (int next=5; next<6; next++){

    egraph[next] = new TGraphErrors(n, xPeak, peakGain[next],exPeak, 0 ); // Grafica de cada extension
    egraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    egraph[next]->SetMarkerStyle(8);
    egraph[next]->SetMarkerSize(0.5);
    // egraph[next]->SetMinimum(210);
    // egraph[next]->SetMaximum(280);

    ajusteLin[next]=new TF1("f", "[0]+[1]*x",0,numpeaks);  //Fit linear model
    ajusteLin[next]->SetParameters(200,0);
    

    
    if(doFit == 1){
      egraph[next]->Fit(ajusteLin[next], "R+");
      ajusteLin[next]->GetParameters(&fitPar[next*2]);
      efitPar[0+next*2]=ajusteLin[next]->GetParError(1);
      printf("\nmean: %f +/- %f", fitPar[1+2*next], efitPar[1+2*next]);
    }
    egraph[next]->Draw("AP"); 

  }

  // ---------------------------------------------
  //Plot gain vs Peak 4
  printf("\nDesviacion Fraccionaria mean/peakNumber \n");
  TCanvas *c6 = new TCanvas("c6","Gain=mean/peak", 2000, 500);//*ceil(numext/4.)); //propiedades del canvas

  
  printf("\n\n numpeaks= %i; iniPeak= %i;xpeak len= %f, n puntos=%i \n\n", numpeaks, iniPeak, sizeXpeak, n);

  for (int next=5; next<6; next++){

    egraph[next] = new TGraphErrors(n, xPeak, peakGainDiv[next],exPeak, 0 ); // Grafica de cada extension
    egraph[next]->SetTitle(Form("ext%i",goodext[next]+1));
    egraph[next]->SetMarkerStyle(8);
    egraph[next]->SetMarkerSize(0.5);
    // egraph[next]->SetMinimum(210);
    // egraph[next]->SetMaximum(280);

    ajusteLin[next]=new TF1("f", "[0]+[1]*x",0,numpeaks);  //Fit linear model
    ajusteLin[next]->SetParameters(200,0);
    

    
    if(doFit == 1){
      egraph[next]->Fit(ajusteLin[next], "R+");
      ajusteLin[next]->GetParameters(&fitPar[next*2]);
      efitPar[0+next*2]=ajusteLin[next]->GetParError(1);
      printf("\nmean: %f +/- %f", fitPar[1+2*next], efitPar[1+2*next]);
    }
    egraph[next]->Draw("AP"); 

  }

}
    
  //}
