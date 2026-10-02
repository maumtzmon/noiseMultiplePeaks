//-----This script uses the root file (*_hepeaks.root) generated from hepeaks_*.C
//-----and performs a fit to the electron peaks using TSpectrum to find the peaks 
//-----and then a single Gaussian fit to each peak to extract the mean and sigma. 

#include <string.h>
#include <fstream>
#include <sstream>
#include <regex>
#include <map>
#include <vector>
#include <iostream>

using namespace std;
#include <math.h>

TH1F *hist_arriba;
TGraph *g_abajo;
TVirtualPad *pad_abajo;

double singleGauss(double *x, double *par){
  double q0 	 =	x[0];
	double norm0peak  =  par[0];
	double mu0    =  par[1];
	double sigma0 =  par[2];
  // double con =  par[3];
	double val=0;

	double cPi = TMath::Pi();

		// val+=(norm0peak/sqrt(2.*cPi*pow(sigma0,2.0)))*(exp(-0.5*pow(q0-mu0,2.0)/pow(sigma0,2.0)));
    val+=(norm0peak)*(exp(-0.5*pow(q0-mu0,2.0)/pow(sigma0,2.0)));
	
  // val+=con;
	return val;
}

double twoGaussians(double *x, double *par){
	double q 	 =	x[0];
	double norm  =  par[0];
	double mu1    =  par[1];
	double sigma =  par[2];
	double mu2  =  par[3];
	//double con   =  par[4];
  double norm2   =  par[4];

	double val=0;

	double cPi = TMath::Pi();

	//for(int w=0; w<2; w++){
		// val=(norm/sqrt(2*cPi*pow(sigma,2)))*(exp(-0.5*pow(q-mu1,2)/pow(sigma,2)))+(norm2/sqrt(2*cPi*pow(sigma,2)))*(exp(-0.5*pow(q-mu2,2)/pow(sigma,2)));
    val=(norm)*(exp(-0.5*pow(q-mu1,2)/pow(sigma,2)))+(norm2)*(exp(-0.5*pow(q-mu2,2)/pow(sigma,2)));
    // Agregar una segunda Sigma para hacerlas completamente independientes.
    //}
	//val+=con*exp(-0.01*q);
	return val;
}

double multGaussians(double *x, double *par){
	double q 	 =	x[0];
	double norm  =  par[0];
	double mu    =  par[1];
	double sigma =  par[2];
	double gain  =  par[3];
	double con   =  par[4];
  double expoE =  par[5];
  double sigma2 = par[6];
 	double val=0;

	double cPi = TMath::Pi();

	for(int w=0; w<3; w++){
		val+=((norm/sqrt(2*cPi*pow(sigma,2)))*(exp(-0.5*pow(q-mu-w*gain,2)/pow(sigma2,2))))*exp(-expoE*sqrt(q));
	}
	val+=con;
	return val;
}

Double_t gauss_poisson_fit(Double_t *x, Double_t *par) {
    Int_t k =1000;
    Double_t xval = x[0];
    Double_t a     = par[0];
    Double_t mu    = par[1];
    Double_t sigma = par[2];    
    Double_t lambda_poisson = par[3];
    Double_t gain  = par[4];
        
    Double_t fitval = 0.0;
    for (Int_t p = 0; p <= k; p++){
          // Asi esta
          //fitval += a * TMath::Gaus(xval*gain,p+mu,sigma,1) * TMath::PoissonI(p,lambda_poisson); 
          // Asi la debo cambiar
          fitval += a * TMath::Gaus(xval,(p+mu)*gain,sigma,1) * TMath::PoissonI(p,lambda_poisson);

    }//for Int_t p

    return fitval;

}

void SincronizarZoom() {
    if (!hist_arriba || !g_abajo || !pad_abajo) return;

    // Obtener rango visible del histograma en coordenadas del eje
    TAxis *axis_src = hist_arriba->GetXaxis();
    Double_t x_min_hist = axis_src->GetBinLowEdge(axis_src->GetFirst());
    Double_t x_max_hist = axis_src->GetBinLowEdge(axis_src->GetLast()) + 
                          axis_src->GetBinWidth(axis_src->GetLast());
    
    // --- IMPORTANTE: Transformar coordenadas si los rangos son diferentes ---
    // Obtener rangos originales completos
    Double_t hist_min = hist_arriba->GetXaxis()->GetXmin();
    Double_t hist_max = hist_arriba->GetXaxis()->GetXmax();
    Double_t graph_min = g_abajo->GetXaxis()->GetXmin();
    Double_t graph_max = g_abajo->GetXaxis()->GetXmax();
    
    // Mapeo lineal de [hist_min, hist_max] a [graph_min, graph_max]
    // Fórmula: valor_en_graph = graph_min + (valor_en_hist - hist_min) * (graph_max - graph_min) / (hist_max - hist_min)
    Double_t x_min_graph = graph_min + (x_min_hist - hist_min) * (graph_max - graph_min) / (hist_max - hist_min);
    Double_t x_max_graph = graph_min + (x_max_hist - hist_min) * (graph_max - graph_min) / (hist_max - hist_min);
    
    // Asegurar que no nos salimos del rango del gráfico
    x_min_graph = TMath::Max(x_min_graph, graph_min);
    x_max_graph = TMath::Min(x_max_graph, graph_max);
    
    // Aplicar zoom al gráfico
    g_abajo->GetXaxis()->SetRangeUser(x_min_graph, x_max_graph);
    
    // Actualizar pad inferior
    pad_abajo->Modified();
    pad_abajo->Update();
    
    // Información de depuración (opcional)
    printf("Zoom: Hist=[%.2f, %.2f] -> Graph=[%.2f, %.2f]\n", 
    x_min_hist, x_max_hist, x_min_graph, x_max_graph);
}

TH1F* CloneHistogramRange(TH1F* original, Double_t x_min, Double_t x_max, const char* new_name = nullptr) {
    if (!original) return nullptr;
    
    // Encontrar los bins de límite
    int bin_min = original->GetXaxis()->FindBin(x_min);
    int bin_max = original->GetXaxis()->FindBin(x_max);
    
    // Asegurar que bin_min y bin_max están dentro del rango
    bin_min = TMath::Max(1, bin_min);
    bin_max = TMath::Min(original->GetNbinsX(), bin_max);
    
    if (bin_min >= bin_max) {
        std::cerr << "Error: Rango inválido para clonar" << std::endl;
        return nullptr;
    }
    
    // Crear el nuevo histograma
    int n_bins = bin_max - bin_min + 1;
    Double_t x_min_new = original->GetXaxis()->GetBinLowEdge(bin_min);
    Double_t x_max_new = original->GetXaxis()->GetBinUpEdge(bin_max);
    
    const char* name = new_name ? new_name : Form("%s_range", original->GetName());
    TH1F *cloned = new TH1F(name, original->GetTitle(), n_bins, x_min_new, x_max_new);
    
    // Copiar contenidos
    for (int i = bin_min; i <= bin_max; i++) {
        int new_bin = i - bin_min + 1;
        cloned->SetBinContent(new_bin, original->GetBinContent(i));
        cloned->SetBinError(new_bin, original->GetBinError(i));
    }
    
    // Copiar estadísticas
    cloned->SetEntries(original->GetEntries());
    
    // Copiar estilo y opciones (opcional)
    cloned->SetLineColor(original->GetLineColor());
    cloned->SetLineWidth(original->GetLineWidth());
    cloned->SetFillColor(original->GetFillColor());
    cloned->SetFillStyle(original->GetFillStyle());
    
    return cloned;
}


//--------VARIABLES TO FIT HISTOGRAMS--------
const int numpeaks = 1000;
const int tenFirst = 1;   // modificar el numero de picos que usara TSpactrum. Se sugiere 10 o menos
const int numext = 16;        // Number of working extensions

// char outputFilename[100];
// only mean variables are actually used in this file
double meanPeak[numext][numpeaks-1], meanPeak_2[numext][numpeaks-1], fitGain[numext][numpeaks-1], constFit[numext][numpeaks-1], expoFit[numext][numpeaks-1];
std::array<std::array<double, numext>, numpeaks-1> meanPeak_sort;
double meanPeakErr[numext][numpeaks-1], sigmaFit[numext][numpeaks-1], Gain_Method1[numext], Gain_Method2[numext], egainPeak[numext][numpeaks-1];
double mediana=0;

double pendiente[numext];
double error_pendiente[numext];

int goodext[numext] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
// float expgain[numext] = {200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};			// Expected gain in ADU/e-

double sigma[numext][numpeaks-1], lambda[numext];

int next_ini = 0, pointIndex = 0; // Esta es la extension que va a mostrar 0-15, pointIndex es para la C12
int next_end = next_ini+1; // Esta es la extension que va a mostrar 0-15
Int_t nfound;
Double_t *xpeaks_sorted;
Double_t *ypeaks_sorted;
//--------Expected values----------

// Expected noise         // Valores de prueba, estos se usaron para probar el algoritmo. 
                          // para cada MCM se utilizaran los datos del archivo .root 
double exp_noise[16];      // = {55.21, 53.36, 52.25, 54.81, 56.42, 53.39, 56.04, 52.95, 53.26, 56.1, 55.01, 54.5, 54.06, 54.67, 56.12, 56.57};
// Expected gain
double exp_gain[16];       // = {249.06, 183.66, 243.44, 240.54, 247.95, 239.42, 178.88, 238.12, 233.77, 154.59, 248.49, 250.97, 246.36, 241.86, 247.06, 255.65};
// Expected lamda poisson
double exp_lambda[16]  = {0.152,0.116,0.121,0.106,0.131,0.106,0.01,0.01,0.0144,0.083,0.094,0.0103,0.0157,0.0127,0.107,0.119};

double fitrange, offsetaux, chisquare[numext][numpeaks-1], cociente[numext][numpeaks-1];
char path[1000]="/home/oem/datosFits/DarkBeats/Brenda/outputs_noisemultipeaks/ANSAMP_300_01APR25/";
int Gain_OK[16];

//--------FUNCTIONS--------

void noisedcmm_TSpectrum(char const* file){
  // File to save output
  std::ofstream outputFile("output_data.txt");
  // write header for CSV (extension,peak,index,mean)
  outputFile << "ext,peak,xpeak,meanPeak" << std::endl;
  // --------Retrieve expgain from root file
  TFile *input = new TFile(file);
  TTree *metadata = (TTree*)input->Get("metadata");
  // float Gain;
  // float expgain[16];
  float expgain2;
  float noise;
  int gain_ok;
  metadata->SetBranchAddress("Gain",&expgain2);
  //Long64_t nEntries = tree->GetEntries();
  const char* dirPath = "";

  for (int i = 0; i < metadata->GetEntries(); ++i) {
    metadata->GetEntry(i); 
    exp_gain[i]=expgain2;
  }

  metadata->SetBranchAddress("Noise", &noise);
  for (int i = 0; i < metadata->GetEntries(); ++i) {
    metadata->GetEntry(i); 
    exp_noise[i]=noise;
  }

  metadata->SetBranchAddress("Gain_OK", &gain_ok);
  for (int i = 0; i < metadata->GetEntries(); ++i) {
    metadata->GetEntry(i); 
    Gain_OK[i]=gain_ok;
  }

  //input->Close();

//--------Style--------

  gROOT->Reset();

  TGaxis::SetMaxDigits(3);

//--------Retrieve filename without extension--------

  int length = strlen(file);
  char fileroot[length+1]; //
  
  strcpy(fileroot, file);		// Copy the input string to fileroot
  fileroot[length-5] = '\0';		// Throw ".root" from filename (last 5 characters)

  char* lastSlash = strrchr(fileroot, '/');
  if (lastSlash != nullptr) {
    // Desplazar el nombre al inicio del buffer
    char* nombre = lastSlash + 1;
    // Copiar el nombre al inicio (se puede hacer con memmove por solapamiento)
    memmove(fileroot, nombre, strlen(nombre) + 1);
  }

  

//--------Retrieve histograms from the root file--------

  TFile filehist(Form("%s", file));

  TH1F *hpix[numext];
  for (int next=next_ini; next<next_end; next++){

	hpix[next] = (TH1F*)filehist.Get(Form("ext%i", goodext[next]));
	hpix[next]->SetDirectory(0);
  hpix[next]->SetAxisRange(0.0,5.0);

	}
	
  filehist.Close();
	input->Close();
//--------Find peaks to fit and fit them--------

  // keep only the variables that are actually referenced later
  double norm0, sigma0;
  int binmin, binmax, ndf[numext][numpeaks-1];
  

  TF1 *fitfun[10];
  std::vector<std::vector<TF1*>> fits_by_ext(numext);
  //int iniPeak=500;

  float mean1;  // used only for printing

  int iniPeak=0;
  int loq=-200, hiq=400;//15900; 
  
  int maxq = iniPeak*250+100; //loq+300; 
  fitrange=1.5*250;

  // ---------------------------------------------
  // Data Extraction, iteracion por cada extension
  TCanvas *c1= new TCanvas("c1","Plot #1 histograma",                   2000, 500);//*ceil(numext/4.))
  TCanvas *c2= new TCanvas("c2","Plot #2 Histograma y Fondo asociado",  2000, 500);//*ceil(numext/4.))
  TCanvas *c3= new TCanvas("c3","Plot #3 Histograma sin fondo",         2000, 500);//*ceil(numext/4.))
  TCanvas *c4= new TCanvas("c4","Plot #4 fit two Gaussians",            2000, 500);//*ceil(numext/4.))
  TCanvas *c5= new TCanvas("c5","Plot #5 fit poisson",                  2000, 500);//*ceil(numext/4.))
  TCanvas *c6= new TCanvas("c6","Plot #6 xpeak vs #Pico from TSpectrum",2000, 500);//*ceil(numext/4.))
  TCanvas *c7= new TCanvas("c7","Plot #7 xpeak Fit vs #Pico",           2000, 500);//*ceil(numext/4.))
  TCanvas *c8= new TCanvas("c8","Plot #8 Sigma vs #Pico",               2000, 500);//*ceil(numext/4.))
  TCanvas *c9= new TCanvas("c9","Plot #9 Chi-square vs #Pico",               2000, 500);//*ceil(numext/4.))
  TCanvas *c10= new TCanvas("c10","Plot #10 ajuste lineal de G=media/#pico",               2000, 500);//*ceil(numext/4.))
  TCanvas *c11= new TCanvas("c11","Plot #11 ajuste lineal de G=mu-mu2",               2000, 500);//*ceil(numext/4.))
  TCanvas *c12= new TCanvas("c12","Plot #12 gain vs histogram",               2000, 500);//*ceil(numext/4.))

  std::vector<std::pair<Double_t, Double_t>> peaks;
  int peak2measure = 0;
  for (int next=next_ini; next<next_end; next++){			// Loop of exts
    if(Gain_OK[next]==1){				// Si el gain de la extension es correcto, hacer el análisis, sino no hacer nada
          //histogram labels
      hpix[next]->SetTitle(Form("ext%i", goodext[next]));
      hpix[next]->GetXaxis()->SetTitle("Pixel value [ADU]");
      hpix[next]->GetYaxis()->SetTitle("Number of counts");

      //valores iniciales para el primer ajuste
      binmin = hpix[next]->GetXaxis()->FindBin(loq); 
      binmax = hpix[next]->GetXaxis()->FindBin(maxq);
      offsetaux=loq+fitrange;
      maxq = hpix[next]->GetXaxis()->GetXmax();
      printf("maxq: %i\n", maxq);
      // peak2measure = static_cast<int>(0.80*maxq/exp_gain[next]); // number of peaks to measure, based on the expected gain
      peak2measure = 50;
      printf("peak2measure: %i\n", peak2measure);
      //loq=offsetaux-0.5*expgain[next];
      binmin = hpix[next]->GetXaxis()->FindBin(loq); 
      binmax = hpix[next]->GetXaxis()->FindBin(maxq);
      hpix[next]->GetXaxis()->SetRange(binmin, binmax);

        // int Npar = 6; // 3 for singleGauss, 6 for twoGaussians
        // mean1=0;
        // norm0=1;

            
      //}


      // ---------------------------------------------
      // Backgraund subtraction using TSpectrum
      // Peak finding using TSpectrum
      int min1, max1;

      min1 = hpix[next_ini]->GetXaxis()->FindBin(loq); 
      max1 = hpix[next_ini]->GetXaxis()->FindBin(maxq);
      // hpix[next_ini]->GetXaxis()->SetRange(min1, max1);


      // for (int next=next_ini; next<next_end; next++){
      // Plot Histopgrama original
      c1->Clear();
      c1->cd();
      gPad->SetLogy();
      hpix[next]->GetXaxis()->SetRangeUser(-exp_gain[next], maxq);//numpeaks*250); // Cambia estos valores al rango que quieras mostrar
      hpix[next]->Draw("hist");
      // c1->SaveAs(Form("%s/MCM15_ANSAMP500/C1_espectro/%s_ext%i_hist.png", fileroot, next));
      dirPath=Form("%s/%s/C1_espectro", path,fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      c1->SaveAs(Form("%s/%s/C1_espectro/%s_ext%i_hist.png",path, fileroot, fileroot, next));
      printf("%s/%s/C1_espectro/%s_ext%i_hist.png\n", path, fileroot, fileroot, next);
      c1->Update();

      //histograma en un Double_t para usar el metodo Background de TSpectrum
      // Extraer contenidos del histograma a un array Double_t

      // Plot del histograma original y el fondo estimado por TSpectrum
      TLegend *legend_c2 = new TLegend(0.6,0.7,0.9,0.9);
      c2->Clear();
      c2->cd();
      gPad->SetLogy();

      Int_t nbins = hpix[next]->GetNbinsX();
      Double_t *spectrum = new Double_t[nbins];
      TSpectrum *s = new TSpectrum(peak2measure);
      nfound = s->Search(hpix[next], 3, "", 1.e-09);
      Double_t *xpeaks;
      Double_t *ypeaks;

      xpeaks = s->GetPositionX();
      ypeaks = s->GetPositionY();

      std::vector<std::pair<Double_t, Double_t>> peaks;
      for (Int_t i = 0; i < nfound; i++) {
        peaks.push_back(std::make_pair(xpeaks[i], ypeaks[i]));
      }

      // Ordenar por el primer elemento (x) de menor a mayor
      std::sort(peaks.begin(), peaks.end());

      // Si necesitas los arreglos ordenados
      xpeaks_sorted = new Double_t[nfound];
      ypeaks_sorted = new Double_t[nfound];
      for (Int_t i = 0; i < nfound; i++) {
        xpeaks_sorted[i] = peaks[i].first;
        ypeaks_sorted[i] = peaks[i].second;
      }

      TH1F *d4 = (TH1F*)hpix[next]->Clone("d4");
    
      d4->Reset(); // Limpiar contenidos

      for(Int_t i = 1; i <= nbins; i++){
        spectrum[i-1] = hpix[next]->GetBinContent(i);
      }
      s->Background(spectrum, nbins, 40, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2 , kFALSE, TSpectrum::kBackSmoothing3, kFALSE);
      
      for(Int_t i = 1; i <= nbins; i++){
        d4->SetBinContent(i, spectrum[i-1]);
      }
      
      d4->SetLineColor(kMagenta);
      // hpix[next]->Draw("hist");
      d4->Draw("same");
      
      dirPath=Form("%s/%s/C2_fondo_Tspectrum", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      // c2->SaveAs(Form("%s/MCM15_ANSAMP300/C2_fondo_Tspectrum/%s_ext%i_hist_fondo.png", fileroot, next));
      c2->SaveAs(Form("%s/%s/C2_fondo_Tspectrum/%s_ext%i_hist_fondo.png", path, fileroot, fileroot, next));
      printf("%s/%s/C2_fondo_Tspectrum/%s_ext%i_hist_fondo.png\n", path, fileroot, fileroot, next);
      
      c2->Update();
      // ---------------------------------------------
      // Plot del histograma original con el fondo substraido
      c3->Clear();
      c3->cd();
      TH1F *hist_noBkgd = (TH1F*)hpix[next]->Clone("hist_noBkgd");
      
      //hist_noBkgd->Reset();
      
      hist_noBkgd-> Draw("hist");

      hist_noBkgd->Add(d4,-1);
      gPad->SetLogy();
      hist_noBkgd->Draw();
      
      dirPath=Form("%s/%s/C3_Bkgd_substracted", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      // c3->SaveAs(Form("%s/MCM15_ANSAMP300/C3_Bkgd_substracted/%s_ext%i_hist_nobkgd.png", fileroot, next));
      c3->SaveAs(Form("%s/%s/C3_Bkgd_substracted/%s_ext%i_hist_nobkgd.png", path, fileroot, fileroot, next));
      printf("%s/%s/C3_Bkgd_substracted/%s_ext%i_hist_nobkgd.png\n", path, fileroot, fileroot, next);
      c3->Update();



      // ---------------------------------------------
      // Fit de los picos encontrados por TSpectrum
      // Hace el fit de dos gausianas a cada pareja de picos encontrados por TSpectrum.

      c4->Clear();
      c4->cd();

      gPad->SetLogy();
      TH1F *hist_c4 = (TH1F*)hpix[next]->Clone("hist_c4");
      hist_c4->GetXaxis()->SetRangeUser(-1*exp_gain[next], peak2measure*exp_gain[next]);//numpeaks*250);
      hist_c4->Draw("hist");
      //d4->Draw("same");
      double norm_mu=0,norm_mu2=0;
      // TLegend *legend_c4 = new TLegend(0.6,0.7,0.9,0.9);
      for(int ii=0; ii<10; ii=ii+2){ // Loop de picos encontrados por TSpectrum, se sugiere usar 10 o menos
        // iterar por cada pico encontrado por TSpectrum, hacer un fit con twoGaussians
        // extraer el mean del fit y guardarlo en un CSV junto con la extension y el numero de pico
        printf("\nContador de picos: %d\nposicion en X de los picos: %1.3f, %1.3f\nExpected gain: %1.3f\n", ii, xpeaks_sorted[ii], xpeaks_sorted[ii+1], exp_gain[next]);
        if(ii==0){ // Si es el primer o quinto pico, hacer un fit con twoGaussians, sino con singleGauss
          TF1 *fitfun = new TF1(Form("fitfun_%d", ii), twoGaussians, xpeaks_sorted[ii]-1.05*exp_noise[next], xpeaks_sorted[ii+1]+1.05*exp_noise[next], 6);
          fitfun->SetParameter(0, ypeaks_sorted[ii]);                            // norm
          fitfun->SetParLimits(0, 0.5*ypeaks_sorted[ii], 1.5*ypeaks_sorted[ii]);
          //hacer Setparameter(de todos los parametros) tambien
          fitfun->SetParameter(1, xpeaks_sorted[ii]);                            // mu
          fitfun->SetParLimits(1, xpeaks_sorted[ii]-exp_noise[next], xpeaks_sorted[ii]+exp_noise[next]);
          //fitfun->SetParLimits(1, xpeaks_sorted[ii]-10, xpeaks_sorted[ii]+10);          // mu 

          fitfun->SetParameter(2, exp_noise[next]);                       // sigma  
          //fitfun->SetParLimits(2, 0.5*exp_noise[next], 1.5*exp_noise[next]);                                 // sigma  

          //fitfun->SetParLimits(3, xpeaks_sorted[ii+1]-xpeaks_sorted[ii]-10,xpeaks_sorted[ii+1]-xpeaks_sorted[ii]+10); // mu2
          fitfun->SetParameter(3, xpeaks_sorted[ii+1]);                          // mu2
          //fitfun->SetParLimits(3, xpeaks_sorted[ii+1]-10,xpeaks_sorted[ii+1]+10);       // mu2
          fitfun->SetParLimits(3, xpeaks_sorted[ii+1]-exp_noise[next], xpeaks_sorted[ii+1]+exp_noise[next]);       // mu2
          
          fitfun->SetParameter(4, ypeaks_sorted[ii+1]);              // con
          fitfun->SetParLimits(4, 0.5*ypeaks_sorted[ii+1], 1.5*ypeaks_sorted[ii+1]);         //con

          fitfun->SetParameter(5, exp_noise[next]);                       // sigma  
          fitfun->SetParLimits(5, 0.5*exp_noise[next], 1.5*exp_noise[next]);                                 // sigma 
          

          hist_c4->Fit(fitfun, "RML");                                 // Fit del pico 0 y 1, que son el pico de 0e- y el pico de 1e-

          // Guardar el puntero del TF1 para redibujar en c12
          fits_by_ext[next].push_back(fitfun);
          
          meanPeak[next][ii] = fitfun->GetParameter(1);                   // Guardar el mean del pico de 0e- en el array meanPeak    
          meanPeakErr[next][ii]=fitfun->GetParError(1);
          // if(meanPeakErr[next][ii]==0){meanPeakErr[next][ii]=0.1;}
          meanPeak[next][ii+1] = fitfun->GetParameter(3);                 // Guardar el mean del pico de 1e- en el array meanPeak
          meanPeakErr[next][ii+1]=fitfun->GetParError(1);
          // if(meanPeakErr[next][ii+1]==0){meanPeakErr[next][ii+1]=0.1;}
          sigma[next][ii] = fitfun->GetParameter(2);                      // Guardar el sigma del pico de 0e- en el array sigma 
          sigma[next][ii+1] = fitfun->GetParameter(2);                    // Guardar el sigma del pico de 1e- en el array sigma
          fitGain[next][ii]= (meanPeak[next][ii+1])-meanPeak[next][ii];   // Guardar la ganancia calculada a partir del fit en el array fitGain, 
                                                                          // que es la diferencia entre el mean del pico de 1e- y el mean del pico de 0e-        ndf[ii]=fitfun->GetNDF();
          egainPeak[next][ii]=sqrt(pow(meanPeakErr[next][ii],2)+pow(meanPeakErr[next][ii+1],2));
          egainPeak[next][ii+1]=egainPeak[next][ii];
          // if(egainPeak[next][ii]==0){egainPeak[next][ii]=0.1;}

          
          chisquare[next][ii] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
          chisquare[next][ii+1] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
          ndf[next][ii]=fitfun->GetNDF();
          ndf[next][ii+1]=fitfun->GetNDF();
          cociente[next][ii]=chisquare[next][ii]/ndf[next][ii];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
          cociente[next][ii+1]=chisquare[next][ii+1]/ndf[next][ii+1];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente

          // fitfun->SetLineColor(kBlue);
          fitfun->SetLineColor(kRed);
          // legend_c4->AddEntry(fitfun, TString::Format("fit peak %d-%d, Gain peak %d-%d=%1.3f",ii,ii+1,ii+1,ii,fitGain[next][ii]).Data(), "l");
          fitfun->Draw("lsame");
        }
        else{
          TF1 *fitfun = new TF1(Form("fitfun_%d", ii), twoGaussians, xpeaks_sorted[ii]-exp_noise[next], xpeaks_sorted[ii+1]+exp_noise[next], 6);
          fitfun->SetParameter(0, ypeaks_sorted[ii]);                            // norm
          fitfun->SetParLimits(0, 0.5*ypeaks_sorted[ii], 1.5*ypeaks_sorted[ii]);        // norm
          //hacer Setparameter(de todos parametros) tambien
          fitfun->SetParameter(1, xpeaks_sorted[ii]);                          // mu     
          fitfun->SetParLimits(1, xpeaks_sorted[ii]-exp_noise[next], xpeaks_sorted[ii]+exp_noise[next]);        // mu 
          
          fitfun->SetParameter(2, exp_noise[next]);                       // sigma  
          fitfun->SetParLimits(2, 0.5*exp_noise[next], 1.5*exp_noise[next]);                                 // sigma  
          // fitfun->SetParLimits(2, 0, exp_noise[next]);                                 // sigma  

          fitfun->SetParameter(3, xpeaks_sorted[ii+1]);
          fitfun->SetParLimits(3, xpeaks_sorted[ii+1]-exp_noise[next],xpeaks_sorted[ii+1]+exp_noise[next]); // mu2
          //fitfun->SetParLimits(4, ypeaks_sorted[ii+1]-10, ypeaks_sorted[ii+1]+10);         //norm2 
          fitfun->SetParameter(4, ypeaks_sorted[ii+1]);
          fitfun->SetParLimits(4, 0.5*ypeaks_sorted[ii+1], 1.5*ypeaks_sorted[ii+1]);

          fitfun->SetParameter(5, exp_noise[next]);                       // sigma  
          fitfun->SetParLimits(5, 0.5*exp_noise[next], 1.5*exp_noise[next]);                                 // sigma 

          hist_c4->Fit(fitfun, "RML");                                       // Fit de la siguiente pareja de picos, 

          // Guardar el puntero del TF1 para redibujar en c12
          fits_by_ext[next].push_back(fitfun);

          // norm_mu = fitfun->GetParameter(0);
          norm_mu2 = fitfun->GetParameter(4);
          meanPeak[next][ii] = fitfun->GetParameter(1);                         // Guardar el mean del primer pico de la pareja en el array meanPeak
          meanPeakErr[next][ii]=fitfun->GetParError(1);
          // if(meanPeakErr[next][ii]==0){meanPeakErr[next][ii]=0.1;}
          meanPeak[next][ii+1] = fitfun->GetParameter(3);                       // Guardar el mean del segundo pico de la pareja en el array meanPeak
          meanPeakErr[next][ii+1]=fitfun->GetParError(1);
          // if(meanPeakErr[next][ii+1]==0){meanPeakErr[next][ii+1]=0.1;}
          sigma[next][ii] = fitfun->GetParameter(2);                            // Guardar el sigma del primer pico de la pareja en el array sigma
          sigma[next][ii+1] = fitfun->GetParameter(2);                          // Guardar el sigma del segundo pico de la pareja en el array sigma
          fitGain[next][ii]= abs((meanPeak[next][ii+1]))-abs(meanPeak[next][ii]);  // Guardar la ganancia calculada a partir del fit en el array fitGain, 
                                                                          // que es la diferencia entre el mean del segundo pico de la pareja y el mean del primer pico de la pareja, en valor absoluto para evitar problemas con picos negativos
          fitGain[next][ii-1]= abs(meanPeak[next][ii])-abs(meanPeak[next][ii-1]); // Guardar la ganancia calculada del primer pico de la pareja actual y el segundo de la pareja anterior
                                                                                  // a partir del fit en el array fitGain, 

          egainPeak[next][ii]=sqrt(pow(meanPeakErr[next][ii],2)+pow(meanPeakErr[next][ii+1],2));
          egainPeak[next][ii+1]=egainPeak[next][ii];
          if(egainPeak[next][ii]==0){egainPeak[next][ii]=0.1;}

          chisquare[next][ii] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
          chisquare[next][ii+1] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
          ndf[next][ii]=fitfun->GetNDF();
          ndf[next][ii+1]=fitfun->GetNDF();
          cociente[next][ii]=chisquare[next][ii]/ndf[next][ii];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
          cociente[next][ii+1]=chisquare[next][ii+1]/ndf[next][ii+1];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
                         // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
          
          fitfun->SetLineColor(kRed);
          fitfun->Draw("lsame");
          
          
          
          // legend_c4->AddEntry(fitfun, TString::Format("fit peak %d-%d, Gain peak %d-%d=%1.3f",ii,ii+1,ii+1,ii,fitGain[next][ii] ).Data(), "l");
          // legend_c4->Draw();
        } 
        printf("ultimo pico:\n norm: %1.3f\n mu:%1.3f\nsigma: %1.3f\n,gain: %1.3f\n", norm_mu2, meanPeak[next][ii], sigma[next][ii], fitGain[next][ii]);
      }
      

      for(int ii=10; ii<peak2measure-1; ii=ii+2){ // Loop de picos usando la media calculada del segundo pico de la pareja actual y aproximar la media del primer pico de la 
                                    // siguiente pareja usando la ganancia calculada del fit
        printf("nContador de picos: %d\nposicion en X de los picos: %1.3f, %1.3f\n mu= %1.3f, mu2= %1.3f\nExpected gain: %1.3f\n", ii, meanPeak[next][ii-1]+fitGain[next][ii-2]-sigma[next][ii-1], meanPeak[next][ii-1]+2*fitGain[next][ii-2]+sigma[next][ii-1], meanPeak[next][ii-1]+fitGain[next][ii-2], meanPeak[next][ii-1]+fitGain[next][ii-2]*2, exp_gain[next]);                             
        TF1 *fitfun = new TF1(Form("fitfun_%d", ii), twoGaussians, meanPeak[next][ii-1]+fitGain[next][ii-2]-0.3*exp_gain[next], meanPeak[next][ii-1]+2*fitGain[next][ii-2]+0.3*exp_gain[next], 6);
        fitfun->SetParameter(0, norm_mu2);       // norm
        //hacer Setparameter(de todos parametros) tambien
        //fitfun->SetParameter(1, meanPeak[next][ii-1]+fitGain[next][ii-2]);
        fitfun->SetParameter(1, meanPeak[next][ii-1]+fitGain[next][ii-2]);                          // mu
        fitfun->SetParLimits(1, meanPeak[next][ii-1]+fitGain[next][ii-2]-0.3*exp_gain[next], meanPeak[next][ii-1]+fitGain[next][ii-2]+0.3*exp_gain[next]);        // mu
        // fitfun->SetParameter(1, 2490.0);                          // mu     
        // fitfun->SetParLimits(1, 2490.0-10, 2490.0+10);        // mu 
        
        fitfun->SetParameter(2, exp_noise[next]);                       // sigma  
        fitfun->SetParLimits(2, 0.5*exp_noise[next], 1.5*exp_noise[next]);                                 // sigma  
        //fitfun->SetParLimits(2, 0, exp_noise[next]);                                 // sigma  

        //fitfun->SetParameter(3, meanPeak[next][ii-1]+fitGain[next][ii-2]*2); // mu2
        fitfun->SetParameter(3, meanPeak[next][ii-1]+fitGain[next][ii-2]*2); // mu2
        fitfun->SetParLimits(3, meanPeak[next][ii-1]+fitGain[next][ii-2]*2-0.3*exp_gain[next], meanPeak[next][ii-1]+fitGain[next][ii-2]*2+0.3*exp_gain[next]); // mu2
        //fitfun->SetParLimits(4, norm_mu2*0.6, norm_mu2*0.8);         //norm2
        //fitfun->SetParLimits(4, norm_mu2*0.8, norm_mu2*0.9);   //norm2 
        fitfun->SetParameter(4, norm_mu2);
        // fitfun->SetParLimits(4, ypeaks[ii+1]-1, ypeaks[ii+1]+1);

        fitfun->SetParameter(5, exp_noise[next]);                       // sigma  
        //fitfun->SetParLimits(5, 0, exp_noise[next]);                                 // sigma 
        fitfun->SetParLimits(5, 0.5*exp_noise[next], 1.5*exp_noise[next]);
        
        // hist_c4->Fit(fitfun, "RML");                                       // Fit de la siguiente pareja de picos, 
        hist_c4->Fit(fitfun, "R+");
        // Guardar el puntero del TF1 para redibujar en c12
        fits_by_ext[next].push_back(fitfun);
        meanPeak[next][ii] = fitfun->GetParameter(1);                         // Guardar el mean del primer pico de la pareja en el array meanPeak
        meanPeakErr[next][ii]=fitfun->GetParError(1);
        // if(meanPeakErr[next][ii]==0){meanPeakErr[next][ii]=0.1;}
        meanPeak[next][ii+1] = fitfun->GetParameter(3);                       // Guardar el mean del segundo pico de la pareja en el array meanPeak
        meanPeakErr[next][ii+1]=fitfun->GetParError(3);
        // if(meanPeakErr[next][ii+1]==0){meanPeakErr[next][ii+1]=0.1;}
        sigma[next][ii] = fitfun->GetParameter(2);                            // Guardar el sigma del primer pico de la pareja en el array sigma
        sigma[next][ii+1] = fitfun->GetParameter(2);                          // Guardar el sigma del segundo pico de la pareja en el array sigma
        fitGain[next][ii]= abs((meanPeak[next][ii+1]))-abs(meanPeak[next][ii]);  // Guardar la ganancia calculada a partir del fit en el array fitGain, 
                                                                        // que es la diferencia entre el mean del segundo pico de la pareja y el mean del primer pico de la pareja, en valor absoluto para evitar problemas con picos negativos
        fitGain[next][ii-1]= abs(meanPeak[next][ii])-abs(meanPeak[next][ii-1]); // Guardar la ganancia calculada del primer pico de la pareja actual y el segundo de la pareja anterior
                                                                                // a partir del fit en el array fitGain, 

        egainPeak[next][ii]=sqrt(pow(meanPeakErr[next][ii],2)+pow(meanPeakErr[next][ii+1],2));  
        egainPeak[next][ii+1]=egainPeak[next][ii];  
        if(egainPeak[next][ii]==0){egainPeak[next][ii]=0.1;}     

        norm_mu2 = fitfun->GetParameter(4);
        chisquare[next][ii] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
        chisquare[next][ii+1] = fitfun->GetChisquare();        // Guardar el chi cuadrado del fit en el array chisquare
        ndf[next][ii]=fitfun->GetNDF();
        ndf[next][ii+1]=fitfun->GetNDF();
        cociente[next][ii]=chisquare[next][ii]/ndf[next][ii];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
        cociente[next][ii+1]=chisquare[next][ii+1]/ndf[next][ii+1];           // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
                    // Guardar el cociente entre el chi cuadrado y los grados de libertad del fit en el array cociente
        
        //fitfun->SetLineColor(kRed);
        fitfun->SetLineColor(kBlue);
        fitfun->Draw("lsame");


      }

      dirPath=Form("%s/%s/C4_Fit_2gauss", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }

      c4->Update();
      // c4->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C4_Fit_2gauss/%s_ext%i_hist_fits.png", fileroot, next));
      c4->SaveAs(Form("%s/%s/C4_Fit_2gauss/%s_ext%i_hist_fits.png", path, fileroot, fileroot, next)); 
      printf("./%s_ext%i_hist.png\n", fileroot, next);

      // ---------------------------------------------
      // Plot C12, sincronizacion del Zoom

      

      // c12->Clear();
      // c12->cd();
      // c12->Divide(1, 2);
      
      // TH1F* hist_c12 = (TH1F*)hist_c4->DrawClone("hist_c12");
      // hist_arriba = CloneHistogramRange((TH1F*)hist_c12, -250, peak2measure*exp_gain[next], "hist_arriba");

      // // TH1F* hist_arriba = (TH1F*)hist_c4->DrawClone("hist_arriba");
      // // hist_arriba->GetXaxis()->SetRangeUser(-200.0, peak2measure*exp_gain[next]);
      

      // g_abajo = new TGraph();
      // g_abajo->SetMarkerStyle(kFullCircle);
      // pointIndex = 0;
      
      
      // for(int ii=1; ii<peak2measure; ii++){
      //     g_abajo->SetPoint(pointIndex, ii, (meanPeak[next][ii]/ii));
      //     printf("Pico %d: gain=%1.3f \n", ii, (meanPeak[next][ii]/ii));
      //     pointIndex++;
      // }
      // g_abajo->SetTitle(Form("ext%i", goodext[next]));
      // g_abajo->GetXaxis()->SetTitle("Peak number");
      // g_abajo->GetYaxis()->SetTitle("Gain [ADU/e-]");
      // g_abajo->GetXaxis()->SetRangeUser(0, peak2measure);
      // //g_abajo->GetYaxis()->SetRangeUser(fitGain[next][0]-15,fitGain[next][0]+15);
      
      // c12->cd(1);
      // gPad->SetGrid();
      // gPad->SetLogy();
      // hist_arriba->Draw(); //"hist_clone"
      // // Redibujar los TF1 guardados para la extensión actual en el pad superior
      // for (size_t ifit = 0; ifit < fits_by_ext[next].size(); ++ifit) {
      //   TF1 *f = fits_by_ext[next][ifit];
      //   if (!f) continue;
      //   f->SetLineWidth(2);
      //   f->Draw("lsame");
      //   hist_arriba->GetListOfFunctions()->Add(f);
      // }
      // TExec *execZoom = new TExec("execZoom", "SincronizarZoom();");
      // hist_arriba->GetListOfFunctions()->Add(execZoom);

      // c12->cd(2);
      // pad_abajo = gPad;  // --- AHORA SÍ FUNCIONA: TVirtualPad* = TVirtualPad* ---
      // gPad->SetGrid();
      // g_abajo->Draw("AL");
      // g_abajo->GetXaxis()->SetLimits(0, peak2measure);
      // c12->Update();
    





      // ---------------------------------------------
      // Fit de los picos usando una convolucion de una gaussiana con una distribucion de Poisson,
      TH1F *hist_clone = (TH1F*)hpix[next]->Clone("hist_clone");
      c5->Clear();
      c5->cd();
      gPad->SetLogy();
      hist_clone->GetXaxis()->SetRangeUser(-250, peak2measure*exp_gain[next]); // Cambia estos valores al rango que quieras mostrar
      hist_clone->Draw("hist_clone");
      for(int ii=0; ii<nfound; ii=ii+2){
        if(ii==0){ // 
          TF1 *fitfun_poisson = new TF1(Form("fitfun_poisson_%d", ii), gauss_poisson_fit, xpeaks_sorted[ii]-200,xpeaks_sorted[ii+1]+1.5*exp_noise[next], 5);

          fitfun_poisson->SetParameter(0, ypeaks_sorted[ii]);                      // norm
          
          fitfun_poisson->SetParameter(1, xpeaks_sorted[ii]);                      // mu
          fitfun_poisson->SetParLimits(1, xpeaks_sorted[ii]-10, xpeaks_sorted[ii]+10);    // mu 
          
          fitfun_poisson->SetParameter(2, exp_noise[next]);                 // sigma
          fitfun_poisson->SetParLimits(2, exp_noise[next] - 5, exp_noise[next] + 5);  // sigma SetParLimits(2, 45,80); // sigma  
          
          fitfun_poisson->SetParameter(3, exp_lambda[next]); // poisson lambda
          fitfun_poisson->SetParLimits(3, .09,0.3); // poisson lambda
          
          fitfun_poisson->SetParameter(4, exp_gain[next]);                            // gain
          fitfun_poisson->SetParLimits(4, exp_gain[next] - 5, exp_gain[next] + 5);  // gain
          //fitfun_poisson->SetParLimits(4, xpeaks[ii+1]-xpeaks[ii]-10,xpeaks[ii+1]-xpeaks[ii]+10); // gain


          hist_clone->Fit(fitfun_poisson, "RML");
          // Guardar el puntero del TF1 poisson para redibujar en c12 si se desea
          fits_by_ext[next].push_back(fitfun_poisson);
          // double meanPeak = fitfun->GetParameter(1);
          // double meanPeakErr = fitfun->GetParError(1); 
          
          lambda[next] = fitfun_poisson->GetParameter(3); 
          
          fitfun_poisson->SetLineColor(kBlue);
          fitfun_poisson->Draw("lsame");

        }
        TLegend *legend_c5 = new TLegend(0.6,0.7,0.9,0.9);
        legend_c5->Draw();
      }
      dirPath=Form("%s/%s/C5_Fit_poisson", path,fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      //c5->SaveAs(Form("%s/MCM45 _ANSAMP500_hepeaks/C5_Fit_poisson/%s_ext%i_hist_poissonFit.png", fileroot, next));
      c5->SaveAs(Form("%s/%s/C5_Fit_poisson/%s_ext%i_hist_poissonFit.png", path, fileroot, fileroot, next));
      printf("./%s_ext%i_hist_poissonFit.png\n", fileroot, next);
      c5->Update();

      // ---------------------------------------------
      // Ganancia TSpectrum
      c6->Clear();
      c6->cd();

      TGraph *g_gain = new TGraph();
      g_gain->SetMarkerStyle(kFullCircle);
      int pointIndex = 0;
      for(int ii=1; ii<peak2measure; ii++){
        g_gain->SetPoint(pointIndex, ii, xpeaks_sorted[ii]/ii);
        printf("Pico %d: xpeak gain=%1.3f \n", ii, xpeaks_sorted[ii]/ii);
        pointIndex++;
      }
      g_gain->SetTitle(Form("ext%i", goodext[next]));
      g_gain->GetXaxis()->SetTitle("Peak number");
      g_gain->GetYaxis()->SetTitle("Gain [ADU/e-]");
      g_gain->GetXaxis()->SetRangeUser(0, peak2measure);
      g_gain->GetYaxis()->SetRangeUser(xpeaks_sorted[1]-15,xpeaks_sorted[1]+15);
      g_gain->Draw("AP");
      printf("./%s_ext%i_hist_poissonFit.png\n", fileroot, next);
      dirPath=Form("%s/%s/C6_gainTspectrum", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      //c6->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C6_gainTspectrum/%s_ext%i_xpeak_Peak_num.png", fileroot, next));
      c6->SaveAs(Form("%s/%s/C6_gainTspectrum/%s_ext%i_xpeak_Peak_num.png", path, fileroot, fileroot, next));
      c6->Update();

      // ---------------------------------------------
      // Ganancia, Media del vs numero de pico, usando el fit de dos gausianas
      c7->Clear();
      c7->cd();
      printf("Plot 7: xpeak Fit vs #Pico \n");
      TGraph *g_fit = new TGraph();
      g_fit->SetMarkerStyle(kFullCircle);
      pointIndex = 0;
      for(int ii=1; ii<peak2measure-1; ii++){
            g_fit->SetPoint(pointIndex, ii, (meanPeak[next][ii])); // ii));
            printf("Pico %d: mean on histogram =%1.3f , mean TSpectrum =%1.3f\n", ii, (meanPeak[next][ii]), (xpeaks_sorted[ii]));///ii);
            pointIndex++;
      }
      g_fit->SetTitle(Form("ext%i", goodext[next]));
      g_fit->GetXaxis()->SetTitle("Peak number");
      g_fit->GetYaxis()->SetTitle("Gain [ADU/e-]");
      g_fit->GetXaxis()->SetRangeUser(0, peak2measure-1);
      // g_fit->GetYaxis()->SetRangeUser(fitGain[next][0]-15,fitGain[next][0]+15);


      g_fit->Draw("AP");
      
      TF1 *fit_7 = new TF1("fit_7", "pol1", 0, peak2measure);
      g_fit->Fit(fit_7);

      pendiente[next] = fit_7->GetParameter(1);
      error_pendiente[next] = fit_7->GetParError(1);
      TLatex *tex = new TLatex();
      tex->SetNDC(); // coordenadas normalizadas (0-1)
      tex->DrawLatex(0.15, 0.85, Form("Pendiente = %.4f +- %.4f", pendiente[next], error_pendiente[next]));


      dirPath=Form("%s/%s/C7_gainFit", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }

      //c7->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C7_gainFit/%s__ext%i_xpeakFit_Peak_num.png", fileroot, next));
      c7->SaveAs(Form("%s/%s/C7_gainFit/%s__ext%i_xpeakFit_Peak_num.png", path, fileroot, fileroot, next));
      c7->Update();


      // ---------------------------------------------
      // Sigma del fit de dos gausianas vs numero de pico
      c8->Clear();
      c8->cd();
      printf("Plot 8: sigma Fit vs #Pico \n");
      TGraph *g_sigma = new TGraph();
      g_sigma->SetMarkerStyle(kFullCircle);
      pointIndex = 0;
      for(int ii=0; ii<peak2measure-1; ii++){
            g_sigma->SetPoint(pointIndex, ii, sigma[next][ii]);
            printf("Pico %d: sigma=%1.3f \n", ii, sigma[next][ii]);
            pointIndex++;
      }
      g_sigma->SetTitle(Form("ext%i", goodext[next]));
      g_sigma->GetXaxis()->SetTitle("Peak number");
      g_sigma->GetYaxis()->SetTitle("Sigma [ADU]");
      g_sigma->GetXaxis()->SetRangeUser(0, peak2measure-1);
      //  g_sigma->GetYaxis()->SetRangeUser(0, 100);
      g_sigma->Draw("AP");
      dirPath=Form("%s/%s/C8_sigmaFit", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      // c8->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C8_sigmaFit/%s__ext%i_sigmaFit_Peak.png", fileroot, next));
      c8->SaveAs(Form("%s/%s/C8_sigmaFit/%s__ext%i_sigmaFit_Peak.png", path, fileroot, fileroot, next));
      c8->Update();


      // ---------------------------------------------
      // Chi-cuadrado del fit de dos gausianas vs numero de pico
      c9->Clear();
      c9->cd();
      printf("Plot 9: chi-square vs #Pico \n");
      TGraph *g_chi = new TGraph();
      g_chi->SetMarkerStyle(kFullCircle);
      pointIndex = 0;
      for(int ii=0; ii<peak2measure-1; ii=ii+2){
            g_chi->SetPoint(pointIndex, ii, cociente[next][ii]);
            printf("Pico %d: chi-square/ndf=%1.3f, ndf=%d\n", ii, cociente[next][ii], ndf[next][ii]);
            pointIndex++;
      }
      g_chi->SetTitle(Form("ext%i", goodext[next]));
      g_chi->GetXaxis()->SetTitle("Peak number");
      g_chi->GetYaxis()->SetTitle("Chi-square/NDF");
      g_chi->GetXaxis()->SetRangeUser(0, peak2measure-1);
      // g_chi->GetYaxis()->SetRangeUser(0, 100);
      g_chi->Draw("AP");

      dirPath=Form("%s/%s/C9_chiSquareFit", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      // c9->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C9_chiSquareFit/%s__ext%i_chiSquareFit_Peak.png", fileroot, next));
      c9->SaveAs(Form("%s/%s/C9_chiSquareFit/%s__ext%i_chiSquareFit_Peak.png", path, fileroot, fileroot, next));
      c9->Update();

      //----------------------------------------------
      //  Ajuste lineal de la ganancia calculada a partir del fit de dos gausianas,
      //  es decir, un ajuste de la media del pico de cada numero de pico vs el numero de pico, para obtener 
      // la ganancia a partir del ajuste lineal, que seria la pendiente del ajuste lineal

      c10->Clear();
      c10->cd();

      TGraphErrors *g_gainfit = new TGraphErrors();
      g_gainfit->SetMarkerStyle(kFullCircle);

      int pointIndexx = 0;
      for(int ii = 1; ii < peak2measure-1; ii++){

        g_gainfit->SetPoint(pointIndexx, ii, (meanPeak[next][ii]-meanPeak[next][0])/ii);

        if(meanPeakErr[next][ii]*sqrt(cociente[next][ii])==0){
          g_gainfit->SetPointError(pointIndexx, 0, 0.1);
        }
        else{
          g_gainfit->SetPointError(pointIndexx, 0, meanPeakErr[next][ii]*sqrt(cociente[next][ii])/ii); // cociente -> se calcula en C4. entre el chi cuadrado y los grados de libertad del fit
        }

        
        //g_gainfit->SetPointError(pointIndexx, 0, meanPeakErr[next][ii]); // cociente -> se calcula en C4. entre el chi cuadrado y los grados de libertad del fit
        

        printf("Pico %d: mean=%1.3f \n", ii, (meanPeak[next][ii]-meanPeak[next][0])/ii);

        pointIndexx++;
      }

      g_gainfit->SetTitle(Form("ext%i;Peak number;ADC value", goodext[next]));

      g_gainfit->Draw("AP");
      g_gainfit->GetXaxis()->SetRangeUser(0, peak2measure-1);
      // g_gainfit->GetYaxis()->SetRangeUser(350, 450);
      // Ajuste lineal
      TF1 *fline = new TF1("f", "pol0", 0, peak2measure);
      fline->SetParLimits(0,300,400);
      TFitResultPtr r = g_gainfit->Fit(fline, "S");

      // Parámetros//
      
      double gainnnn   = r->Parameter(0);
      //double egain  = r->ParError(1);
      //double offset = r->Parameter(0);
      Gain_Method1[next]=gainnnn;
      // Mostrar resultado
      TLatex latex;
      latex.SetNDC();
      latex.DrawLatex(0.15,0.85,
      Form("Gain = %.3f" , gainnnn)    );
      dirPath=Form("%s/%s/C10_gainFit", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      // c10->SaveAs(Form("%s/MCM15_ANSAMP300_hepeaks/C10_gainFit/%s__ext%i_xpeakFit_Peak_num.png", fileroot, next));
      c10->SaveAs(Form("%s/%s/C10_gainFit/%s__ext%i_xpeakFit_Peak_num.png", path, fileroot, fileroot, next));
      c10->Update();


      //---------------------------------------------- (Texto generado automaticamente)
      // Ajuste lineal de la ganancia calculada a partir del fit de dos gausianas, pero usando la diferencia entre el mean del pico de cada numero de pico y el mean del pico anterior, 
      // para obtener la ganancia a partir del ajuste lineal, que seria la pendiente del ajuste lineal
      c11->Clear();
      c11->cd();
      
      TGraphErrors *g_gainfit_err = new TGraphErrors();
      g_gainfit_err->SetMarkerStyle(kFullCircle);
      pointIndex = 0;
      for(int ii=0; ii<peak2measure-2; ii++){
            g_gainfit_err->SetPoint(pointIndex, ii+1, fitGain[next][ii]);
            // if(egainPeak[next][ii]*sqrt(cociente[ii])==0){
            //   g_gainfit_err->SetPointError(pointIndex, 0, 0.1);
            // }
            // else{
            g_gainfit_err->SetPointError(pointIndex, 0, egainPeak[next][ii]*sqrt(cociente[next][ii]));// egain*sqrt(chi2/ndf));
            // }
            pointIndex++;
      }
      g_gainfit_err->SetTitle(Form("ext%i;Peak number;ADC value", goodext[next]));

      g_gainfit_err->Draw("AP");
      g_gainfit_err->GetXaxis()->SetRangeUser(0, peak2measure-2);
      g_gainfit_err->GetYaxis()->SetTitle("Gain Value [ADU/e-]");

      
      printf("./%s_ext%i_hist.png\n", fileroot, next);
      // Ajuste lineal
      TF1 *ffline = new TF1("ff", "pol0", 0, 10);
      TFitResultPtr r2 = g_gainfit_err->Fit(ffline, "S");

      // Parámetros//
      double gainnnn2   = r2->Parameter(0);
      double egain2  = r2->ParError(0);
      //double offset = r->Parameter(0);
      Gain_Method2[next]=gainnnn2;
      // Mostrar resultado
      TLatex latex2;
      latex2.SetNDC();
      latex2.DrawLatex(0.15,0.85,
      Form("Gain = %.3f#pm %.3f" , gainnnn2,egain2));
      dirPath=Form("%s/%s/C11_gainFit", path, fileroot);
      if (gSystem->AccessPathName(dirPath)) {
        std::cout << "El directorio no existe. Creándolo..." << std::endl;
        gSystem->MakeDirectory(dirPath);
      }
      
      // c11->SaveAs(Form("%s/MCM15_ANSAMP300/C11_gainFit/%s_ext%i_hist.png", fileroot, next));
      c11->SaveAs(Form("%s/%s/C11_gainFit/%s_ext%i_hist.png", path, fileroot, fileroot, next));
      c11->Update();


      for (int i = 0; i < peak2measure-1; i++){
        printf("Sigma %d: %1.3f\n", i, sigma[next][i]);
      }

     
    }
  }

char mcmID[32];
int anSamp;
TString ts = fileroot;

sscanf(ts.Data(),"MCM%[0-9A-Z]_ANSAMP%d", mcmID, &anSamp);
  

  std::ofstream gain_stream(Form("%s/%s/%s_gainMethods.tsv", path, fileroot, fileroot));

  // std::stringstream gain_stream;

  gain_stream << "MCMID\tGAIN_method1\tGAIN_method2\tGAIN_method3" << std::endl;
  gain_stream << "{";

  gain_stream << "MCM" << mcmID << "}\t{";

  for(int i = 0; i < numext; i++){
      gain_stream << Gain_Method1[i];
      if(i != numext-1) gain_stream << ", ";
  }

    gain_stream << "}\t{";

  for(int i = 0; i < numext; i++){
      gain_stream << Gain_Method2[i];
      if(i != numext-1) gain_stream << ", ";
  }

  gain_stream << "}\t{";

  for(int i = 0; i < numext; i++){  
      gain_stream << pendiente[i];  // Ganing_Method3[i] = pendiente[i];
      if(i != numext-1) gain_stream << ", ";
  }

  gain_stream << "}";


 gain_stream.close();




// Crear y abrir el archivo TSV
std::ofstream archivoSalida(Form("%s/%s/%s_allParameters.tsv", path, fileroot, fileroot));
if (!archivoSalida.is_open()) {
    std::cerr << "Error: No se pudo abrir el archivo" << std::endl;
    
  }

archivoSalida << "MCMID\tNum_Ext\tANSAMP\tGAIN_method1:Media/#pico\tEGAIN_method1\tGAIN_method2:mean1-mean2\tEGAIN_method2\tGAIN_method3:slope\tEGAIN_method3\tNoise\tChi2/ndf/" << std::endl;


// Escribir los datos en UNA SOLA FILA
  for(int next=0; next<16; next++){  // numext -> 0 - 15
    archivoSalida << mcmID << "\t";   // filtrar el numero de MCM del archivo, COLUMA 1
    archivoSalida << next+1 << "\t"; // Filtar el num de Samps del archivo  COLUMNA 2
    archivoSalida << anSamp << "\t"; // Filtar el num de Samps del archivo  COLUMNA 3
    archivoSalida << "{";       // Abrir corchete para la lista de GAIN_method1  Media/#pico
    for (int i = 1; i < peak2measure; i++) {
        if(i==1){archivoSalida << meanPeak[next][i];} //COLUMNA4     Metodo 1: Media/#pico
        else{archivoSalida << (meanPeak[next][i]-meanPeak[next][0])/i;}
        if (i < peak2measure - 1) archivoSalida << ",";
    }

    archivoSalida << "}\t{";    //  Error normalizado de la media/#pico // COLUMNA 5
    for (int i = 1; i < peak2measure; i++) {
        archivoSalida << meanPeakErr[next][i]*sqrt(cociente[next][i])/i;
        if (i < peak2measure - 1) archivoSalida << ",";
    }

    archivoSalida << "}\t{";    //ganancia medida de la (media del fit pico1- media del fit pico2)
    for (int i = 0; i < peak2measure-1; i++) { //COLUMNA 6  Metodo 2: Media del fit pico1- media del fit pico2
        archivoSalida << fitGain[next][i];
        if (i < peak2measure - 2) archivoSalida << ",";
    }

    archivoSalida << "}\t{";  //error de la ganancia medida de la (media del fit pico1- media del fit pico2)
    for (int i = 0; i < peak2measure-1; i++) {  //COLUMNA 7
        archivoSalida << egainPeak[next][i];
        if (i < peak2measure - 2) archivoSalida << ",";
    }
    
    archivoSalida << "}\t{";
    archivoSalida << pendiente[next];  // COLUMNA 8  Metodo 3: pendiente del ajuste lineal de la ganancia calculada a partir del fit de dos gausianas
    archivoSalida <<"}\t{";   
    archivoSalida <<error_pendiente[next];// Filtar el num de Samps del archivo 
    
    // archivoSalida << "}" << std::endl;

    archivoSalida << "}\t{";  //error de la ganancia medida de la (media del fit pico1- media del fit pico2)
    for (int i = 0; i < peak2measure; i++) {
        archivoSalida << sigma[next][i];
        if (i < peak2measure - 1) archivoSalida << ",";
    }
    
    archivoSalida << "}\t{";  //error de la ganancia medida de la (media del fit pico1- media del fit pico2)
    for (int i = 0; i < peak2measure; i++) {
        archivoSalida << cociente[next][i]/ndf[next][i];
        if (i < peak2measure - 1) archivoSalida << ",";
    }
    archivoSalida << "}" << std::endl;
  }
archivoSalida.close();
std::cout << "Archivo TSV creado exitosamente: " << fileroot << std::endl;
}