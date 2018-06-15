#include "H2inv_ll_br_limits.C"
#define nparam 7

//this function scan over parameters and give the values for optimal parameters
void cuts_optimization_scan(Int_t print_canvas =0){
  vector<Double_t> dtheta;
  dtheta.push_back(100);
  //dtheta.push_back(110);
  //dtheta.push_back(115);

  vector<Double_t> Pt;
  //Pt.push_back(5);
  Pt.push_back(10);
  //Pt.push_back(15);

  vector<Double_t> Pl;
  //Pl.push_back(40);
  Pl.push_back(50);
  //Pl.push_back(60);

  vector<Double_t> acopl;
  //acopl.push_back(5);
  acopl.push_back(10);
  //acopl.push_back(15);

  //vector<Double_t> helicity = {0, 5, 10};
  vector<Double_t> helicity;
  helicity.push_back(0);

  vector<Double_t> mll_up;
  //mll_up.push_back(93.2);
  //mll_up.push_back(94.2);
  mll_up.push_back(95.2);
  //mll_up.push_back(96.2);
  //mll_up.push_back(100.2);

  vector<Double_t> mll_down;
  //mll_down.push_back(85.2);
  //mll_down.push_back(86.2);
  mll_down.push_back(87.2);
  //mll_down.push_back(82.2);

  std::ofstream logfile;
  logfile.open("log/cuts_optimization_scan.log");

  Int_t tot_cicle = dtheta.size()*Pt.size()*Pl.size()*acopl.size()*helicity.size()*mll_up.size()*mll_down.size();
  logfile << "This scan will last " << tot_cicle << " cicle" << std::endl << std::flush;

  vector<Double_t> br_limit;
  vector<Double_t> br_sign5;

  Double_t min_br_limit = 10.0;
  Double_t min_br_limit_param[nparam];

  Double_t min_br_sign5 = 20.0;
  Double_t min_br_sign5_param[nparam];

  Double_t results[2];

  Int_t iteration = 1;
  for (unsigned int idt = 0; idt < dtheta.size(); idt++) {

    for (unsigned int ipt = 0; ipt < Pt.size(); ipt++) {

      for (unsigned int ipl = 0; ipl < Pl.size(); ipl++) {

        for (unsigned int iac = 0; iac < acopl.size(); iac++) {

          for (unsigned int ihe = 0; ihe < helicity.size(); ihe++) {

            for (unsigned int imup = 0; imup < mll_up.size(); imup++) {

              for (unsigned int imdw = 0; imdw < mll_down.size(); imdw++) {

                const char* cuts = Form("dtheta_ll_lab>%f&&Pt_ll>%f&&Pl_ll<%f&&acoplanarity_angle>%f&&helicity_angle>%f&&M_ll<%f&&M_ll>%f", dtheta[idt], Pt[ipt], Pl[ipl], acopl[iac], helicity[ihe], mll_up[imup] , mll_down[imdw]);

                logfile << "------ Beginning iteration " << iteration << "/" << tot_cicle << std::endl;

                H2inv_ll_br_limits(&results[0], &results[1],cuts,print_canvas);

                br_limit.push_back(results[0]);
                br_sign5.push_back(results[1]);

                if(results[0] < min_br_limit){
                  min_br_limit = results[0];
                  min_br_limit_param[0]=dtheta[idt];
                  min_br_limit_param[1]=Pt[ipt];
                  min_br_limit_param[2]=Pl[ipl];
                  min_br_limit_param[3]=acopl[iac];
                  min_br_limit_param[4]=helicity[ihe];
                  min_br_limit_param[5]=mll_up[imup];
                  min_br_limit_param[6]=mll_down[imdw];
                  logfile << "New min_br_limit found : " << min_br_limit << std::endl;
                }
                if(results[1] < min_br_sign5){
                  min_br_sign5 = results[1];
                  min_br_sign5_param[0]=dtheta[idt];
                  min_br_sign5_param[1]=Pt[ipt];
                  min_br_sign5_param[2]=Pl[ipl];
                  min_br_sign5_param[3]=acopl[iac];
                  min_br_sign5_param[4]=helicity[ihe];
                  min_br_sign5_param[5]=mll_up[imup];
                  min_br_sign5_param[6]=mll_down[imdw];
                  logfile << "New min_br_sign5 found : " << min_br_sign5 << std::endl;
                }

                iteration++;

              }

            }

          }

        }

      }

    }

  }

  logfile << "Loop is over! Saving results" << std::endl;

  logfile.close();


std::ofstream myfile;
myfile.open(Form("log/cuts_optimization_scan_results_%d.txt",time(NULL)));

myfile << "Span space:" << std::endl;
myfile << "dtheta = ";
for(unsigned int i=0; i<dtheta.size(); ++i){
  myfile << dtheta[i] << '\t';
};
myfile << std::endl;

myfile << "Pt = ";
for(unsigned int i=0; i<Pt.size(); ++i){
  myfile << Pt[i] << '\t';
};
myfile << std::endl;

myfile << "Pl = ";
for(unsigned int i=0; i<Pl.size(); ++i){
  myfile << Pl[i] << '\t';
};
myfile << std::endl;

myfile << "acopl = ";
for(unsigned int i=0; i<acopl.size(); ++i){
  myfile << acopl[i] << '\t';
};
myfile << std::endl;

myfile << "helicity = ";
for(unsigned int i=0; i<helicity.size(); ++i){
  myfile << helicity[i] << '\t';
};
myfile << std::endl;

myfile << "mll_up = ";
for(unsigned int i=0; i<mll_up.size(); ++i){
  myfile << mll_up[i] << '\t';
};
myfile << std::endl;

myfile << "mll_down = ";
for(unsigned int i=0; i<mll_down.size(); ++i){
  myfile << mll_down[i] << '\t';
};
myfile << std::endl;

myfile << "Macro for cuts:" << std::endl;
myfile << "dtheta_ll_lab> x &&Pt_ll> x &&Pl_ll< x &&acoplanarity_angle> x &&helicity_angle> x &&M_ll< x &&M_ll> x " << std::endl;

myfile << std::endl << std::endl;

myfile << "min_br_limit = " << min_br_limit << std::endl;
myfile << "Parameters: " << std::endl;
for(unsigned int i=0; i< nparam ; ++i){
  myfile << min_br_limit_param[i] << std::endl;
};
myfile << std::endl << std::endl;

myfile << "min_br_sign5 = " << min_br_sign5 << std::endl;
myfile << "Parameters: " << std::endl;
for(unsigned int i=0; i< nparam ; ++i){
  myfile << min_br_sign5_param[i] << std::endl;
};
myfile << std::endl << std::endl;

myfile << "All values obtained for br_limit:" << std::endl;
for(unsigned int i=0; i<br_limit.size(); ++i){
  myfile << br_limit[i] << std::endl;
};
myfile << std::endl << std::endl;

myfile << "All values obtained for br_sign5:" << std::endl;
for(unsigned int i=0; i<br_sign5.size(); ++i){
  myfile << br_sign5[i] << std::endl;
};
myfile << std::endl << std::endl;


myfile.close();

}
