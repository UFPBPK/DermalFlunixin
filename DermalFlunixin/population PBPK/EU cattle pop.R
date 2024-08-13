library(mrgsolve)   # for PBPK model
library(dplyr)      # for dataframe manipulation
library(knitr)      # for output table manipulation
library(truncnorm)  # R package for Truncated normal distribution
library(EnvStats)   # R package for Truncated lognormal distribution
library(ggplot2)    # R package for plotting


PBPK <- '
$PARAM @annotated
//All equations related to skin submodel are from (Talreja et al. 2001; Van der Merwe et al. 2006)
  doseevap         :0.000006 :  (rate of permeant evaporation from skin surface (mg/h))
  solventevap      :0.024    :  rate of solvent evaporation from skin surface (ml/h)
  logP             :0.75     :  Log [stratum corneum]/[solvent] 
  D                :12       :  diffusivity in stratum corneum lipid matrix (mm**2/h)
  area             :4800     :  skin surface area (mm**2)
  b                :0.0001   :  minimum solvent volume on skin surface (ml)

  noimped          :1        :  resistance to transport through membrane with no impediments 
  kt               :0.000303 :  corneocyte thickness (mm) 
  kd               :0.032    :  corneocyte diameter (mm) 
  g                :0.000033222922 :  vertical gap between corneocytes (mm)
  s                :0.000033222922 :  ateral gap between corneocytes (mm)
  N                :20       :  number of stratum corneum layers 
  d1               :0.0276   :  long leg of corneocyte overlap (mm)
  d2               :0.0044   :  short leg of corneocyte overlap (mm)
  logPsw           :0.62     :  LogP [stratum corneum]/[water]
  MTf              :6        :  Mass transfer factor
  Vmax             :0.00006  :   # maximum rate of metabolism (mg/h)
  Km               :1000.0   :  [permeant] at 50 % of Vmax (mg/l)
  dermisdepth      :0.5      :  dermatome depth setting (mm)
  Rdepth           :4.8      :  receptor chamber depth  (mm)
  
  
  // Cardiac Output and Blood flow
  
  QCC        : 5.97      :                       Cardiac output (L/h/kg), Li 2019
  QLCa       : 0.405     :                       Fraction of flow to the liver, Li 2019
  QKCa       : 0.09 	   :                       Fraction of flow to the kidneys, Li 2019
  QFCa       : 0.08      :                       Fraction of flow to the fat, Li 2019
  QMCa       : 0.18      :                       Fraction of flow to the muscle, Li 2019
  QrestCa    : 0.245     : Fractional rest of body (total sum equals to 1) (1-VLC-VKC-VFC-VMC-VbloodC)
  QrestC1a   : 0.505     : Fractional rest of body for 5-OH FLU (total sum equals to 1) (1-VLC-VKC-VbloodC)
  
  // Tissue Volume
  
  VbloodCa   : 0.04      :                       Unitless, Fraction vol. of blood; Li 2019             
  VLCa       : 0.014     :                       Unitless, Fraction vol. of liver; Li 2019
  VKCa       : 0.0025    :                       Unitless, Fraction vol. of kidney; Li 2019
  VFCa       : 0.150     :                       Unitless, Fractional fat tissue, Li 2019
  VMCa       : 0.270     :                       Unitless, Fractional muscle tissue, Li 2019
  VartCa     : 0.010     :                       Unitless, Fractional arterial blood, Li 2019
  VvenCa     : 0.030     :                       Unitless, Fractional venous blood, Li 2019
  VrestCa    : 0.5235    : Fractional rest of body (total sum equals to 1) (1-VLC-VKC-VFC-VMC-VbloodC)
  VrestC1a   : 0.9435    : Fractional rest of body for 5-OH FLU (total sum equals to 1) (1-VLC-VKC-VbloodC)
  
  
  // Mass Transfer Parameters (Chemical-specific parameters) 
  // Partition coefficients(PC, tissue:plasma)

  PL        : 10.52      :                     Liver plasma PC Li et al. 2019
  PK        : 7.89       :                     Kidney plasma PC Li et al. 2019
  PM        : 0.5        :                     Muscle plasma PC Li et al. 2019
  PF        : 0.30       :                     Fatp plasma PC Li et al. 2019
  Prest     : 6.5        :                     Rest of the body plasma PC
  
  MW        :296.24      :                     mg/mmol
  MW1       :312.24      :                     mg/mmol
  
  PL1       :	9.26       :                   Liver: plasma PC Li et al. 2019
  PK1       :4           :                   Kidney:plasma PC Li et al. 2019
  Prest1    :5           :                   Rest of the body:plasma PC
  
  // Absorption and elimination parameters
  

  Kint      : 0.400    :              
  Kfeces    : 0.500    :  
  KurineC   : 0.1      :  L/(h*kg); Rate of urine elimination from urine storage
  KurineC1  : 0.2      :  L/(h*kg); Rate of urine elimination from urine storage
  

  KmC       : 0.2      :  /h/kg, metabolic rate constant
  KehcC     : 0.05     :	/h/kg, rate constant for the enterohepatic circulation



  //Biliary Elimination Rate COnstants
  KbileC    :  0.5   :  L/h/kg
  KbileC1   :  0.1   :  L/h/kg
  
  PB        :  0.95    : Plasma protein binding coefficient for flunixin
  PB1       :  0.99    : Plasma protein binding coefficient for flunixin

  
  BW        : 223     : kg, Bodyweight of swine, 
  PDOSEde   : 3.3     : TD dose, (mg/kg)
  tinterval : 24      : h

  
$MAIN
double Free = 1-PB; //Percentage of drug not bound to plasma protein
double Free1 = 1-PB1; 

double sumQ = QLCa + QKCa + QMCa + QFCa + QrestCa; // sum up cardiac output fraction
double QLC = QLCa/sumQ; // adjusted blood flow rate fraction to liver
double QKC = QKCa/sumQ; // adjusted blood flow rate fraction to kidney
double QMC = QMCa/sumQ; // adjusted blood flow rate fraction to muscle
double QFC = QFCa/sumQ; // adjusted blood flow rate fraction to fat
double QrestC = QrestCa/sumQ; // adjusted blood flow rate fraction to rest of body
double QrestC1 = 1-QLC-QKC; // adjusted blood flow rate fraction to rest of body for 5-OH FLU

double sumV = VLCa + VKCa + VMCa + VFCa + VbloodCa + VrestCa; //sum up the tissue volumes
double VLC = VLCa/sumV; //adjusted fraction of tissue volume of liver
double VKC = VKCa/sumV; //adjusted fraction of tissue volume of kidney
double VMC = VMCa/sumV; //adjusted fraction of tissue volume of muscle
double VFC = VFCa/sumV; //adjusted fraction of tissue volume of fat
double VbloodC = VbloodCa/sumV;
double VrestC = VrestCa/sumV; //adjusted fraction of tissue volume of rest of body
double VrestC1 = 1-VbloodC-VLC-VKC; // adjusted fraction of tissue volume of rest of body for 5-OH FLU

double QC = QCC*BW; //Cardiac output
double QL = QLC*QC; //Liver
double QK = QKC*QC; //Kidney
double QF = QFC*QC; //Fat
double QM = QMC*QC; //Muscle
double QR = QrestC*QC; //Rest of body
double QR1 = QrestC1*QC; //Rest of body for 5-OH FLU

double VL = VLC*BW; //Liver
double VK = VKC*BW; //Kidney
double VF = VFC*BW; //Fat
double VM = VMC*BW; //Muscle
double VR = VrestC*BW; //Rest of body
double VR1 = VrestC1*BW; //Rest of body for 5-OH FLU
double Vblood = VbloodC*BW; //Blood

double Kmet = KmC*BW;
double Kehc = KehcC*BW;
double Kurine = KurineC*BW;
double Kurine1 = KurineC1*BW;
double Kbile = KbileC*BW;
double Kbile1 = KbileC1*BW;
double tt   = solventevap*TIME;
//double solventvol = insolventvol-tt;


//tortuosity calculator
double ah = ((kt+g)*(N-1))+kt;         //Actual SC thickness (mm)
double ar = kd/kt;                    //Corneocyte aspect ratio
double br = g/kt;               //Lipid/Corneocyte thickness ratio
double ss = kt/g;                  //Slit shape term
double cvol = (kt+kd)/(kt+kd+g+s);    //Corneocyte volume fraction of SC
double w = d1/d2;     //Corneocyte offset ratio
double Tg = (d2/((N/(N-1))*kt+g))+1; //Geometric tortuosity (Talreja et al 2001)(mm)
double term = kd/(2*s);  //Corneocyte diameter/2*junction gap
double necking = log(term)*((2*g)/ah); //Resistance associated with necking down into gaps between corneocytes
double junctionresist = (N*kd*kt)/(s*ah);   //Resistance of transport through corneocyte junction gaps
double parallelresist = (kd/(1+w))*(kd/(1+w))*(w/(ah*g))*(N-1);  //Resistance of transport in the parallel gaps between corneocytes
double tortuosity = noimped + necking + junctionresist + parallelresist;  //Effective tortuosity
double th = ah*tortuosity+dermisdepth-ah;  //Effective skin thickness (mm)
double minpath = ((kt+g+d2)*(N-1))+kt;   //Mininmum pathway (mm)
      
//SC and surface concentration calculation
double P = (pow(10,logP));  //SC/Solvent partition coefficient
double Jf1 = (P*D*area)/(th/2)*0.000001;  //Fraction of solute moving into stratum corneum (l/h)  
double Jf2 = (1/P)*(D*area)/(th/2)*0.000001;  //Fraction of solute returning to surface (l/h) 
//Stratum corneum block
double Psw = (pow(10,logPsw));  //Stratum corneum/water partition coefficient
double SClipidvol = (g/kt)*area*ah*0.000001;  //Stratum corneum lipid volume (l)
double SCepipart = 1/Psw;  //SC/viable skin partition coefficient
double Jf3 = SCepipart*MTf*0.000001;  //Fraction of solute moving into viable skin (l/h)
double Jf4 = 1/SCepipart*MTf*0.000001; //Fraction of solute returning to stratum corneum (l/h)
//Viable skin block
double RVol = area*Rdepth*0.000001;  //Volume of receptor fluid chamber (l)
double VSvol = area*(dermisdepth-ah)*0.000001+ RVol;   //Volume of VS and receptor fluid chamber (l)
double lagtime = (pow(minpath,2))/D;  //Lag time
//DOSE
double DOSEde=PDOSEde*BW;//mg
double insolventvol=0.001*DOSEde/50;          //solvent volume deposited onto skin surface (l)


$CMT AA AUCCV AA1 AUCCV1 AL AUCCL AL1 AUCCL1 Aehc Amet Abile Abile1 Acolon AI Afeces 
Aurine AK AUCCK Aurine1 AK1 AUCCK1 AM AUCCM AF AUCCF AR AUCCR AR1 AUCCR1 AS ASC AVS 

$ODE

double FVSlog = tinterval< lagtime ? 0: 1;
double solventvol1 = insolventvol-tt;
double solventvol = solventvol1 <= 0.0001 ? 0.0001: solventvol1;

double CS  = AS/solventvol; //mg/l
double Csc = ASC/SClipidvol; //mg/l
double Cvs = AVS/VSvol;//mg/l

double RS  = -CS*Jf1+Csc*Jf2-doseevap;//mg/h
double RSC = (CS*Jf1-Csc*Jf2-Csc*Jf3+Cvs*Jf4)*FVSlog;//mg/h
double RVS = (Csc*Jf3) - Cvs*Jf4;//mg/h

// Dosing, transdermal

dxdt_AS  = RS;//mg
dxdt_ASC = RSC;//mg
dxdt_AVS = RVS;//mg


// Concentrations in the tissues and in the capillary blood of the tissues

double CA = AA/Vblood; // CAfree concentration of unbound drug in the arterial blood (mg/L)
double CAfree = CA*Free;

double CL = AL/VL; // Concentration of parent drug in the tissue of liver(mg/L)
double CVL = CL/PL; // Concentration of parent drug in the capillary blood of liver

double CK = AK/VK; // Concentration of parent drug in the tissue of kidney(mg/L)
double CVK = CK/PK; // Concentration of parent drug in the capillary blood of kidney 

double CM = AM/VM; // Concentration of parent drug in the tissue of muscle(mg/L)
double CVM = CM/PM; // Concentration of parent drug in the capillary blood of muscle 

double CF = AF/VF; // Concentration of parent drug in the tissue of fat(mg/L)
double CVF = CF/PF; // Concentration of parent drug in the capillary blood of  fat

double CR = AR/VR; // Crest drug concentration in the rest of the body (mg/L)
double CVR = CR/Prest; // Concentration of parent drug in the capillary blood of the rest of body

double CA1 = AA1/Vblood; // CAfree concentration of unbound drug in the arterial blood (mg/L)
double CAfree1 = CA1*Free1;

double CL1 = AL1/VL; // Concentration of parent drug in the tissue of liver(mg/L)
double CVL1 = CL1/PL1; // Concentration of parent drug in the capillary blood of liver

double CK1 = AK1/VK; // Concentration of parent drug in the tissue of kidney(mg/L)
double CVK1 = CK1/PK1; // Concentration of parent drug in the capillary blood of kidney 

double CR1 = AR1/VR1; // Crest drug concentration in the rest of the body (mg/L)
double CVR1 = CR1/Prest1; // Concentration of parent drug in the capillary blood of the rest of body

// {flunixin distribution in each compartment}
// flunixin in venous blood compartment
double CV = (QL*CVL+QK*CVK+QF*CVF+QM*CVM+QR*CVR+RVS)/QC; // CV the drug concentration in the venous blood (mmol/L)
double RA = QC*(CV-CAfree);
dxdt_AA = RA;
double CVppm = CV*MW;//(mg/L)
dxdt_AUCCV = CVppm; // (mg*h/L)

double CV1 = (QL*CVL1+QK*CVK1+QR1*CVR1)/QC; // CV1 the 5-OH FLU concentration in the venous blood (mmol/L)
double RA1 = QC*(CV1-CAfree1);
dxdt_AA1 = RA1;
double CV1ppm = CV1*MW1;//(mg/L)
dxdt_AUCCV1 = CV1ppm; // (mg*h/L)

// FLU and 5-OH FLU in liver compartment, flow-limited model
double Rmet1 = Kmet*CL*VL; // Rmet the metabolic rate in liver (mg/h)
dxdt_Amet = Rmet1; // Amet the amount of drug metabolized in liver (mg)
double Rehc = Kehc*AI;
dxdt_Aehc = Rehc;
double Rbile = Kbile*CVL;
dxdt_Abile = Rbile;
double RL = QL*(CAfree-CVL)-Rmet1-Rbile+Rehc; // RL the changing rate of the amount of drug in liver (mg/h)
dxdt_AL = RL; // AL amount of drug in liver (mg) 
double CLppm = CL*MW; //(mg/mL)
dxdt_AUCCL = CLppm; // AUCCL area under the curve of drug concentration in liver (mg*h/L)

double Rbile1 = Kbile1*CVL1;
dxdt_Abile1 = Rbile1;
double RL1 = QL*(CAfree1-CVL1)+Rmet1-Rbile1; // RL1 the changing rate of the amount of drug in liver (mg/h)
dxdt_AL1 = RL1; // AL1 amount of drug in liver (mg) 
double CL1ppm = CL1*MW1; // //(mg/mL)
dxdt_AUCCL1 = CL1ppm; // AUCCL1 area under the curve of drug concentration in liver (mg*h/L)

// Metabolism of FLU in liver compartment

// Elimination of FLU and 5-OH FLU in liver compartment


// GI compartments and Enterohepatic Circulation
double Rcolon = Kint*AI;
dxdt_Acolon = Rcolon;
double RAI = Rbile + Rbile1-Kint*AI-Kehc*AI;
dxdt_AI = RAI;
double Rfeces = Kfeces*Acolon;
dxdt_Afeces = Rfeces;

// flunixin in kidney compartment, flow-limited model
double Rurine = Kurine*CVK;
dxdt_Aurine = Rurine;
double RK = QK*(CAfree-CVK)-Rurine; // RK the changing rate of the amount of drug in kidney (mg/h)
dxdt_AK = RK; // AK amount of drug in kidney (mg)
double CKppm = CK*MW; //(mg/mL)
dxdt_AUCCK = CKppm; // AUCCK AUC of drug concentration in kidney (mg*h/L)

double Rurine1 = Kurine1*CVK1;
dxdt_Aurine1 = Rurine1;
double RK1 = QK*(CAfree1-CVK1)-Rurine1; // RK1 the changing rate of the amount of drug in kidney (mg/h)
dxdt_AK1 = RK1; // AK amount of drug in kidney (mg)
double CK1ppm = CK1*MW1; //(mg/mL)
dxdt_AUCCK1 = CK1ppm; // AUCCK AUC of drug concentration in kidney (mg*h/L)

// flunixin in muscle compartment, flow-limited model
double RM = QM*(CAfree-CVM); //  RM the changing rate of the amount of drug in muscle (mg/h) 
dxdt_AM = RM; // AM amount of the drug in muscle (mg)
double CMppm = CM*MW;//(mg/mL)
dxdt_AUCCM = CMppm;//(mg*h/L)

// flunixin in fat compartment, flow-limited model
double RF = QF*(CAfree-CVF); // RF the changing rate of the amount of drug in fat (mg/h)
dxdt_AF = RF; // AF amount of the drug in fat (mg)
double CFppm = CF*MW;//(mg/mL)
dxdt_AUCCF = CFppm; // AUCCF AUC of drug concentration in fat (mg*h/L)

// flunixin in the compartment of rest of body, flow-limited model
double RR = QR*(CAfree-CVR); // Rrest the changing rate of the amount of drug in the rest of the body (mg/h)
dxdt_AR = RR; // Arest amount of the drug in the rest of the body (mg)
double CRppm = CR*MW;//(mg/mL)
dxdt_AUCCR = CRppm; // AUCCrest AUC of drug concentration in the rest of the body (mg*h/L)

double RR1 = QR1*(CAfree1-CVR1); // Rrest the changing rate of the amount of drug in the rest of the body (mg/h)
dxdt_AR1 = RR1; // Arest amount of the drug in the rest of the body (mg)
double CR1ppm = CR1*MW1;//(mg/mL)
dxdt_AUCCR1 = CR1ppm; // AUCCrest AUC of drug concentration in the rest of the body (mg*h/L)

// {Mass balance equations}
double Qbal = QC-QM-QR-QF-QK-QL;
double Tmass = AA+AM+AR+AF+AK+AL+Aurine+Amet+Abile;
double Input = Aehc+AVS;
double Bal = Input-Tmass;

double Qbal1 = QC-QL-QK-QR1;
double Tmass1 = AA1+AR1+AK1+AL1+Aurine1+Abile1;
double Input1 = Amet;
double Bal1 = Input1-Tmass1;

$TABLE
capture Venous = CVppm;
capture Liver = CLppm;
capture Kidney = CKppm;
capture Muscle= CMppm;
capture Fat = CFppm;
capture Venous1 = CV1ppm;
capture Liver1 =CL1ppm;
capture Kidney1 = CK1ppm;
capture jf1     = Jf1;
capture rS     = RS;
capture jf2     = Jf2;
capture jf3     = Jf3;
capture jf4     = Jf4;
capture bal     =Bal;
capture solventvoll=solventvol;
capture as    = AS;
capture asc     = ASC;
capture avs     = AVS;
'



## Load Model
mod <- mcode_cache("pbpk", PBPK)

## Set Up Simulation Subjects
##There will be N simulation subjects, with their own sampled parameters that will override the default values.  The PBPK model does not use the Omega block to simulate invidual values, instead we fill the idata dataframe with these values.  Here we only change body weight, but other parameters could also be varied.

N=1000
ndoses= 1
doseinterval=24
MW = 296.24 # mg/mmol Molecular weight of flunixin
MW1 = 312.24 # mg/mmol Molecular weight of 5-OH FLU
End_time      = 15
PDOSEde=3.3

sd.log.PL = sqrt(log(1 + 2.104 ^ 2 / 10.52 ^ 2))
m.log.PL = log(10.52) - 0.5 * sd.log.PL ^ 2

sd.log.PL1 = sqrt(log(1 + 1.852 ^ 2 / 9.26 ^ 2))
m.log.PL1 = log(9.26) - 0.5 * sd.log.PL1 ^ 2

sd.log.PK = sqrt(log(1 + 1.578 ^ 2 / 7.89 ^ 2))
m.log.PK = log(7.89) - 0.5 * sd.log.PK ^ 2

sd.log.PK1 = sqrt(log(1 + 0.8 ^ 2 / 4 ^ 2))
m.log.PK1 = log(4) - 0.5 * sd.log.PK1 ^ 2

sd.log.PM = sqrt(log(1 + 0.1 ^ 2 / 0.5 ^ 2))
m.log.PM = log(0.5) - 0.5 * sd.log.PM ^ 2

sd.log.PF = sqrt(log(1 + 0.06 ^ 2 / 0.3 ^ 2))
m.log.PF = log(0.3) - 0.5 * sd.log.PF ^ 2

sd.log.Prest = sqrt(log(1 + 1.3 ^ 2 / 6.5 ^ 2))
m.log.Prest = log(6.5) - 0.5 * sd.log.Prest ^ 2

sd.log.Prest1 = sqrt(log(1 + 1 ^ 2 / 5 ^ 2))
m.log.Prest1 = log(5) - 0.5 * sd.log.Prest1 ^ 2

sd.log.KmC = sqrt(log(1 + 0.04^ 2 / 0.2 ^ 2))
m.log.KmC = log(0.2) - 0.5 * sd.log.KmC ^ 2

sd.log.KurineC = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KurineC = log(0.1) - 0.5 * sd.log.KurineC ^ 2

sd.log.KurineC1 = sqrt(log(1 + 0.06^ 2 / 0.2 ^ 2))
m.log.KurineC1 = log(0.2) - 0.5 * sd.log.KurineC1 ^ 2

sd.log.KehcC = sqrt(log(1 + 0.015^ 2 / 0.05 ^ 2))
m.log.KehcC = log(0.05) - 0.5 * sd.log.KehcC ^ 2

sd.log.KbileC = sqrt(log(1 + 0.15^ 2 / 0.5 ^ 2))
m.log.KbileC = log(0.5) - 0.5 * sd.log.KbileC ^ 2

sd.log.KbileC1 = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KbileC1 = log(0.1) - 0.5 * sd.log.KbileC1 ^ 2

sd.log.Kint = sqrt(log(1 + 0.12^ 2 / 0.4 ^ 2))
m.log.Kint = log(0.4) - 0.5 * sd.log.Kint ^ 2

sd.log.Kfeces = sqrt(log(1 + 0.15^ 2 / 0.5 ^ 2))
m.log.Kfeces = log(0.5) - 0.5 * sd.log.Kfeces ^ 2

sd.log.logPsw = sqrt(log(1 + 0.124^ 2 / 0.62 ^ 2))
m.log.logPsw = log(0.62) - 0.5 * sd.log.logPsw ^ 2

sd.log.logP = sqrt(log(1 + 0.15^ 2 / 0.75 ^ 2))
m.log.logP = log(0.75) - 0.5 * sd.log.logP ^ 2

sd.log.PB = sqrt(log(1 + 0.285^ 2 / 0.95 ^ 2))
m.log.PB = log(0.95) - 0.5 * sd.log.PB ^ 2


set.seed(11009)
idata <- 
  data_frame(ID=1:N) %>% 
  mutate(
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = 300, sd = 0.154*300), 
                    b = qnorm(0.975, mean = 300, sd = 0.154*300), 
                    mean = 300,
                    sd = 0.154*300),
    
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 5.97, sd = 5.97*0.333),
      b = qnorm(0.975, mean = 5.97, sd = 5.97*0.333),
      mean = 5.97,
      sd = 5.97*0.333),  # Cardiac output index (L/h/kg)
    # Fracion of blood flow to organs (unitless) 
    QLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.405, sd = 0.405*0.48),
      b = qnorm(0.975, mean = 0.405, sd = 0.405*0.48),
      mean = 0.405,
      sd = 0.405*0.48
    ),
    # Fraction of blood flow to the kidneys (2016 Lin)
    QKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.09, sd =  0.09*0.3),
      b = qnorm(0.975, mean = 0.09, sd =  0.09*0.3),
      mean = 0.09,
      sd =  0.09*0.3
    ),
    #
    QMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.18, sd =  0.18*0.3),
      b = qnorm(0.975, mean = 0.18, sd =  0.18*0.3),
      mean = 0.18,
      sd =  0.18*0.3
    ),
    # Fraction of blood flow to the fat (2016 Lin)
    QFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.08, sd =  0.08*0.3),
      b = qnorm(0.975, mean = 0.08, sd =  0.08*0.3),
      mean = 0.08,
      sd =  0.08*0.3
    ),
    # Fraction of blood flow to the rest of body (total sum equals to 1)
    QrestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.245, sd = 0.245*0.3),
      b = qnorm(0.975, mean = 0.245, sd = 0.245*0.3),
      mean = 0.245,
      sd = 0.245*0.3
    ),
    QrestC1a = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.505, sd = 0.505*0.3),
      b = qnorm(0.975, mean = 0.505, sd = 0.505*0.3),
      mean = 0.505,
      sd = 0.505*0.3
    ),
    
    # Fractional organ tissue volumes (unitless)
    # Fractional liver tissue (1933 Swett)
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.014, sd =  0.014*0.12),
      b = qnorm(0.975, mean = 0.014, sd =  0.014*0.12) ,
      mean = 0.014,
      sd =  0.014*0.12
    ),
    # Fractional kidney tissue (1933 Swett)
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.0025, sd =  0.0025*0.174),
      b = qnorm(0.975, mean = 0.0025, sd =  0.0025*0.174),
      mean = 0.0025,
      sd =  0.0025*0.174
    ),
    # Fractional fat tissue (2016 Lin, 2014 Leavens)
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.15, sd =  0.15*0.3),
      b = qnorm(0.975, mean = 0.15, sd =  0.15*0.3),
      mean = 0.15,
      sd = 0.15*0.3
    ),
    # Fractional muscle tissue (2016 Lin, 2014 Leavens)
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.27, sd =  0.27*0.3),
      b = qnorm(0.975, mean = 0.27, sd =  0.27*0.3),
      mean = 0.27,
      sd =  0.27*0.3
    ),
    # Venous blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VvenCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.03, sd =  0.03*0.3),
      b = qnorm(0.975, mean = 0.03, sd =  0.03*0.3),
      mean = 0.03,
      sd =  0.03*0.3
    ),
    # Arterial blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VartCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.01, sd = 0.01*0.3),
      b = qnorm(0.975, mean = 0.01, sd = 0.01*0.3),
      mean = 0.01,
      sd =  0.01*0.3
    ),
    # Fractional rest of body (total sum equals to 1)
    VrestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.5235, sd =  0.5235*0.3),
      b = qnorm(0.975, mean = 0.5235, sd =  0.5235*0.3),
      mean = 0.5235,
      sd =  0.5235*0.3
    ),
    VrestC1a = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.9435, sd =  0.9435*0.3),
      b = qnorm(0.975, mean = 0.9435, sd =  0.9435*0.3),
      mean = 0.9435,
      sd =  0.9435*0.3
    ),
    
    # PL
    
    PL <- rlnormTrunc(N, meanlog = m.log.PL, sdlog = sd.log.PL, 
                      min = qlnorm(0.025, meanlog = m.log.PL, sdlog = sd.log.PL), 
                      max = qlnorm(0.975, meanlog = m.log.PL, sdlog = sd.log.PL)),
    
    # PL1
    
    PL1 <- rlnormTrunc(N, meanlog = m.log.PL1, sdlog = sd.log.PL1, 
                       min = qlnorm(0.025, meanlog = m.log.PL1, sdlog = sd.log.PL1), 
                       max = qlnorm(0.975, meanlog = m.log.PL1, sdlog = sd.log.PL1)),
    # PK
    
    PK <- rlnormTrunc(N, meanlog = m.log.PK, sdlog = sd.log.PK, 
                      min = qlnorm(0.025, meanlog = m.log.PK, sdlog = sd.log.PK), 
                      max = qlnorm(0.975, meanlog = m.log.PK, sdlog = sd.log.PK)),
    
    # PK1
    
    PK1 <- rlnormTrunc(N, meanlog = m.log.PK1, sdlog = sd.log.PK1, 
                       min = qlnorm(0.025, meanlog = m.log.PK1, sdlog = sd.log.PK1), 
                       max = qlnorm(0.975, meanlog = m.log.PK1, sdlog = sd.log.PK1)),
    
    # PM
    
    PM <- rlnormTrunc(N, meanlog = m.log.PM, sdlog = sd.log.PM, 
                      min = qlnorm(0.025, meanlog = m.log.PM, sdlog = sd.log.PM), 
                      max = qlnorm(0.975, meanlog = m.log.PM, sdlog = sd.log.PM)),
    
    # PF
    
    PF <- rlnormTrunc(N, meanlog = m.log.PF, sdlog = sd.log.PF, 
                      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF), 
                      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)),
    
    # Prest
    
    Prest <- rlnormTrunc(N, meanlog = m.log.Prest, sdlog = sd.log.Prest, 
                         min = qlnorm(0.025, meanlog = m.log.Prest, sdlog = sd.log.Prest), 
                         max = qlnorm(0.975, meanlog = m.log.Prest, sdlog = sd.log.Prest)),
    
    # Prest1
    
    Prest1 <- rlnormTrunc(N, meanlog = m.log.Prest1, sdlog = sd.log.Prest1, 
                          min = qlnorm(0.025, meanlog = m.log.Prest1, sdlog = sd.log.Prest1), 
                          max = qlnorm(0.975, meanlog = m.log.Prest1, sdlog = sd.log.Prest1)),
    
    # KmC
    
    KmC <- rlnormTrunc(N, meanlog = m.log.KmC, sdlog = sd.log.KmC, 
                       min = qlnorm(0.025, meanlog = m.log.KmC, sdlog = sd.log.KmC), 
                       max = qlnorm(0.975, meanlog = m.log.KmC, sdlog = sd.log.KmC)),
    
    # KurineC
    
    KurineC <- rlnormTrunc(N, meanlog = m.log.KurineC, sdlog = sd.log.KurineC, 
                           min = qlnorm(0.025, meanlog = m.log.KurineC, sdlog = sd.log.KurineC), 
                           max = qlnorm(0.975, meanlog = m.log.KurineC, sdlog = sd.log.KurineC)),
    KurineC1 = rlnormTrunc(
      N,
      meanlog = m.log.KurineC1,
      sdlog = sd.log.KurineC1,
      min = qlnorm(0.025, meanlog = m.log.KurineC1, sdlog = sd.log.KurineC1),
      max = qlnorm(0.975, meanlog = m.log.KurineC1, sdlog = sd.log.KurineC1)
    ),
    # Enterohepatic circulation rate constants
    KehcC = rlnormTrunc(
      N,
      meanlog = m.log.KehcC,
      sdlog = sd.log.KehcC,
      min = qlnorm(0.025, meanlog = m.log.KehcC, sdlog = sd.log.KehcC),
      max = qlnorm(0.975, meanlog = m.log.KehcC, sdlog = sd.log.KehcC)
    ),
    
    # Biliary elimination rate constants
    KbileC = rlnormTrunc(
      N,
      meanlog = m.log.KbileC,
      sdlog = sd.log.KbileC,
      min = qlnorm(0.025, meanlog = m.log.KbileC, sdlog = sd.log.KbileC),
      max = qlnorm(0.975, meanlog = m.log.KbileC, sdlog = sd.log.KbileC)
    ),
    KbileC1 = rlnormTrunc(
      N,
      meanlog = m.log.KbileC1,
      sdlog = sd.log.KbileC1,
      min = qlnorm(0.025, meanlog = m.log.KbileC1, sdlog = sd.log.KbileC1),
      max = qlnorm(0.975, meanlog = m.log.KbileC1, sdlog = sd.log.KbileC1)
    ),
    Kint = rlnormTrunc(
      N,
      meanlog = m.log.Kint,
      sdlog = sd.log.Kint,
      min = qlnorm(0.025, meanlog = m.log.Kint, sdlog = sd.log.Kint),
      max = qlnorm(0.975, meanlog = m.log.Kint, sdlog = sd.log.Kint)
    ),
    Kfeces = rlnormTrunc(
      N,
      meanlog = m.log.Kfeces,
      sdlog = sd.log.Kfeces,
      min = qlnorm(0.025, meanlog = m.log.Kfeces, sdlog = sd.log.Kfeces),
      max = qlnorm(0.975, meanlog = m.log.Kfeces, sdlog = sd.log.Kfeces)
    ),
    PB = rlnormTrunc(
      N,
      meanlog = m.log.PB,
      sdlog = sd.log.PB,
      min = 0.5,
      max = 1
    ),
    
    logP = rlnormTrunc(
      N,
      meanlog = m.log.logP,
      sdlog = sd.log.logP,
      min = qlnorm(0.025, meanlog = m.log.logP, sdlog = sd.log.logP),
      max = qlnorm(0.975, meanlog = m.log.logP, sdlog = sd.log.logP)
    ),
    logPsw = rlnormTrunc(
      N,
      meanlog = m.log.logPsw,
      sdlog = sd.log.logPsw,
      min = qlnorm(0.025, meanlog = m.log.logPsw, sdlog = sd.log.logPsw),
      max = qlnorm(0.975, meanlog = m.log.logPsw, sdlog = sd.log.logPsw)
    ),
    
    DOSEde=PDOSEde*BW/MW)
idata %>% head %>% kable

## Set up dosing.

ev1 <- ev (ID=1:N, amt= idata$DOSEde, ii = doseinterval,addl = ndoses-1, cmt  = "AS", replicate = FALSE)
ex <- ev1 

dim(as.data.frame(ex))
#a <- as.data.frame(ex)
## set up time grid
##We'll sample every hour up to 10 days past the time of last dose

tsamp=tgrid(0,doseinterval*ndoses+10*24,1)

## Combine data and run the simulation
set.seed(11009)
system.time(
  {
    out <- 
      mod %>% 
      data_set(ex) %>%
      idata_set(idata) %>%
      mrgsim(obsonly=TRUE, tgrid=tsamp) #, hmax = 0.01, atol = 1E-20)
  }
) -> run.time

##The run time for the simulation with `r N` subjects was `r run.time['elapsed']` seconds.
## Summarize the data
out.sum = out %>% as.data.frame %>%
  mutate(Time1 = time/24 - (ndoses-1)*doseinterval/24) %>%
  dplyr::select(ID, Time1, Venous:Fat) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Venous, 0.01), CV50 = median(Venous), CV99 = quantile(Venous, 0.99), 
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL99 = quantile(Liver, 0.99),
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK99 = quantile(Kidney, 0.99),
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM99 = quantile(Muscle, 0.99),
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF99 = quantile(Fat, 0.99))

data_tissue <- as.data.frame(out.sum)

end_time = 12


## Plot tissue concentrations
# Assuming TOL is a single value
# Assuming TOL is a single value
withdrawal_time <- min(subset(data_tissue, CL99 <= 0.3 & Time1 > 0.2) %>% dplyr::select(Time1))

p1<-ggplot(data_tissue, aes(Time1)) + 
  geom_ribbon(aes(ymin = CL1,
                  ymax = CL99
  ), fill = 'red',
  show.legend = T, 
  size = 0.2,
  alpha = 0.3) +  # alpha is transparency parameter. 
  geom_line(aes(y = CL50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CL99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CL1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-(ndoses * doseinterval / 24 - 1)):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  geom_line(aes(y = 0.3),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank()) +
  geom_vline(aes(xintercept = withdrawal_time), 
             size = 1, color = "red", linetype = 2,
             show.legend = F)+
  annotate("text", x = withdrawal_time, y = 0.3, label = paste("Liver Withdrawal Interval:",round(withdrawal_time, 2),"Days"),hjust = -0.1, vjust=-0.5,angle = 0, size = 6)

p1


ggsave("p1.tiff",scale = 1.5,
       plot = p1,
       width = 25, height = 15, units = "cm", dpi=320)

withdrawal_time1 <- min(subset(data_tissue, CK99 <= 0.1 & Time1 > 0.2) %>% dplyr::select(Time1))

p2<-ggplot(data_tissue, aes(Time1)) + 
  geom_ribbon(aes(ymin = CK1,
                  ymax = CK99
  ), fill = 'red',
  show.legend = T, 
  size = 0.2,
  alpha = 0.3) +  # alpha is transparency parameter. 
  geom_line(aes(y = CK50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CK99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CK1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-(ndoses*doseinterval/24-1)):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  geom_line(aes(y = 0.1),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  geom_vline(aes(xintercept = min(subset(data_tissue, CK99<= 0.1 & Time1 > 0.2) %>% dplyr::select(Time1))), 
             size = 1, color = "red", linetype = 2,
             show.legend = F)+
  annotate("text", x = withdrawal_time1, y = 0.1, label = paste("Kidney Withdrawal Interval:",round(withdrawal_time1, 2),"Days"),hjust = -0.1, vjust = -0.5, angle = 0, size = 6)

p2


ggsave("p2.tiff",scale = 1.5,
       plot = p2,
       width = 25, height = 15, units = "cm", dpi=320)

withdrawal_time2 <- min(subset(data_tissue, CM99 <= 0.020 & Time1 > 0.2) %>% dplyr::select(Time1))

p3<-ggplot(data_tissue, aes(Time1)) + 
  geom_ribbon(aes(ymin = CM1, ymax = CM99), fill = 'red', show.legend = TRUE, size = 0.2, alpha = 0.3) +  
  geom_line(aes(y = CM50, color = '50 Percentile'), size = 1, show.legend = TRUE) + 
  geom_line(aes(y = CM99, color = '99 Percentile'), size = 1, show.legend = TRUE) +  
  geom_line(aes(y = CM1, color = '1 Percentile'), size = 1, show.legend = TRUE) + 
  scale_x_continuous(name = 'Time (Day)', breaks = c((-(ndoses * doseinterval / 24 - 1)):end_time)) + 
  ylab(expression(paste('Concentration (', mu, 'g/mL)'))) + 
  geom_line(aes(y = 0.020), color = 'black', size = 0.5, linetype = 'twodash', show.legend = FALSE) + 
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text = element_text(size = 16), legend.title = element_blank()) + 
  geom_vline(aes(xintercept = min(subset(data_tissue, CM99 <= 0.020 & Time1 > 0.2) %>% dplyr::select(Time1))), size = 1, color = "red", linetype = 2, show.legend = FALSE) + 
  annotate("text", x = withdrawal_time2, y = 0.020, label = paste("Muscle Withdrawal Interval:", round(withdrawal_time2, 2), "Days"), hjust = -0.1, vjust = -0.5, angle = 0, size = 6)
p3


ggsave("p3.tiff",scale = 1.5,
       plot = p3,
       width = 25, height = 15, units = "cm", dpi=320)

withdrawal_time3 <- min(subset(data_tissue, CF99 <= 0.03 & Time1 > 0.2) %>% dplyr::select(Time1))

p4<-ggplot(data_tissue, aes(Time1)) + 
  geom_ribbon(aes(ymin = CF1,
                  ymax = CF99
  ), fill = 'red',
  show.legend = T, 
  size = 0.2,
  alpha = 0.3) +  # alpha is transparency parameter. 
  geom_line(aes(y = CF50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CF99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CF1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-(ndoses*doseinterval/24-1)):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  geom_line(aes(y = 0.03),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  geom_vline(aes(xintercept = min(subset(data_tissue, CF99<= 0.03 & Time1 > 0.2) %>% dplyr::select(Time1))), 
             size = 1, color = "red", linetype = 2,
             show.legend = F)+
  annotate("text", x = withdrawal_time3, y = 0.03, label = paste("Fat Withdrawal Interval:",round(withdrawal_time3, 2),"Days"),hjust = -0.1, vjust = -0.5, angle = 0, size = 6)
p4


ggsave("p4.tiff",scale = 1.5,
       plot = p4,
       width = 25, height = 15, units = "cm", dpi=320)
run.time

