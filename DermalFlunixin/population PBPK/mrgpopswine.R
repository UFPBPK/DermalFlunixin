library(mrgsolve)   # for PBPK model
library(dplyr)      # for dataframe manipulation
library(knitr)      # for output table manipulation
library(truncnorm)  # R package for Truncated normal distribution
library(EnvStats)   # R package for Truncated lognormal distribution
library(ggplot2)    # R package for plotting


PBPK <- '
$PARAM @annotated
//skin parameters based on Van der Merwe et al. 2006
  doseevap         :0.000006 :  (rate of permeant evaporation from skin surface (mg/h))
  solventevap      :0.024    :  rate of solvent evaporation from skin surface (ml/h)
  logP             :1.19     :  Log [stratum corneum]/[solvent] 
  D                :19.17       :  diffusivity in stratum corneum lipid matrix (mm**2/h)
  area             :4800      :  skin surface area (mm**2)
  e                :2.718282 :  base of natural logarithm
  b                :0.0001   :  minimum solvent volume on skin surface (ml)

  noimped          :1        :  resistance to transport through membrane with no impediments 
  kt               :0.000303 :  corneocyte thickness (mm) 
  kd               :0.032    :  corneocyte diameter (mm) 
  g                :0.000033222922 :  vertical gap between corneocytes (mm)
  s                :0.000033222922 :  ateral gap between corneocytes (mm)
  N                :21.9     :  number of stratum corneum layers 
  d1               :0.0276   :  long leg of corneocyte overlap (mm)
  d2               :0.0044   :  short leg of corneocyte overlap (mm)
  logPsw           :2.13    :  LogP [stratum corneum]/[water]
  MTf              :3.62        :  Mass transfer factor
  Vmax             :0.00006  :   # maximum rate of metabolism (mg/h)
  Km               :1000.0   :  [permeant] at 50 % of Vmax (mg/l)
  dermisdepth      :0.5      :  dermatome depth setting (mm)
  Rdepth           :4.8      :  receptor chamber depth  (mm)
  
  
  
  // Cardiac Output and Blood flow
  
  QCC        : 8.543      :                       Cardiac output (L/h/kg), Li 2019
  QLCa       : 0.273     :                       Fraction of flow to the liver, Li 2019
  QKCa       : 0.116	   :                       Fraction of flow to the kidneys, Li 2019
  QFCa       : 0.128     :                       Fraction of flow to the fat, Li 2019
  QMCa       : 	0.293     :                       Fraction of flow to the muscle, Li 2019
  QrestCa    : 	0.19     : Fractional rest of body (total sum equals to 1) (1-VLC-VKC-VFC-VMC-VbloodC)
  QrestC1a   : 0.611     : Fractional rest of body for 5-OH FLU (total sum equals to 1) (1-VLC-VKC-VbloodC)
  
  // Tissue Volume
  
  VbloodCa   : 0.06      :                       Unitless, Fraction vol. of blood; Li 2019              
  VLCa       : 	0.023    :                       Unitless, Fraction vol. of liver; Li 2019
  VKCa       :0.0045    :                       Unitless, Fraction vol. of kidney; Li 2019
  VFCa       :0.235    :                       Unitless, Fractional fat tissue, Li 2019
  VMCa       : 	0.355    :                       Unitless, Fractional muscle tissue, Li 2019
  VartCa     : 0.016     :                       Unitless, Fractional arterial blood, Li 2019
  VvenCa     : 0.044     :                       Unitless, Fractional venous blood, Li 2019
  VrestCa    : 0.3225     : Fractional rest of body (total sum equals to 1) (1-VLC-VKC-VFC-VMC-VbloodC)
  VrestC1a   : 0.9125     : Fractional rest of body for 5-OH FLU (total sum equals to 1) (1-VLC-VKC-VbloodC)
  
  
  // Mass Transfer Parameters (Chemical-specific parameters) 
  // Partition coefficients(PC, tissue:plasma)

  PL        : 10.52     :                     Liver plasma PC 
  PK        : 	4       :                     Kidney plasma PC 
  PM        : 0.5       :                     Muscle plasma PC 
  PF        : 0.6       :                     Fatp plasma PC 
  Prest     :8          :                     Rest of the body plasma PC
  
  MW        :296.24    :                     mg/mmol
  MW1       :312.24    :                     mg/mmol
  
  PL1       :9.26      :                   Liver: plasma PC 
  PK1       :4         :                   Kidney:plasma PC 
  Prest1    :5         :                   Rest of the body:plasma PC
  
  // Absorption and elimination parameters
  

  Kint      : 0.400    :              
  Kfeces    : 0.5      :  
  KurineC   : 0.1      :  L/(h*kg); Rate of urine elimination from urine storage
  KurineC1  : 0.1      :  L/(h*kg); Rate of urine elimination from urine storage
  

  KmC       : 0.2    :  /h/kg, metabolic rate constant
  KehcC     : 0.15    :	/h/kg, rate constant for the enterohepatic circulation



  //Biliary Elimination Rate COnstants
  KbileC    :  0.1     :  L/h/kg
  KbileC1   :  0.1     :  L/h/kg
  
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
double Tg = (d2/((N/(N-1))*kt+g))+1; //Geometric tortuosity (Talreja et al 2001)
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
double solventvol = solventvol1 <= 0 ? 0.0001: solventvol1;

double CS  = AS/solventvol; //mg/l
double Csc = ASC/SClipidvol; 
double Cvs = AVS/VSvol;

double RS  = -CS*Jf1+Csc*Jf2-doseevap;//mg/h
double RSC = (CS*Jf1-Csc*Jf2-Csc*Jf3+Cvs*Jf4)*FVSlog;
double RVS = (Csc*Jf3) - Cvs*Jf4;

// Dosing, transdermal

dxdt_AS  = RS;
dxdt_ASC = RSC;
dxdt_AVS = RVS;


// Concentrations in the tissues and in the capillary blood of the tissues

double CA = AA/Vblood; // CAfree concentration of unbound drug in the arterial blood (mg/L)
double CAfree = CA*Free;

double CL = AL/VL; // Concentration of parent drug in the tissue of liver
double CVL = CL/PL; // Concentration of parent drug in the capillary blood of liver

double CK = AK/VK; // Concentration of parent drug in the tissue of kidney
double CVK = CK/PK; // Concentration of parent drug in the capillary blood of kidney 

double CM = AM/VM; // Concentration of parent drug in the tissue of muscle
double CVM = CM/PM; // Concentration of parent drug in the capillary blood of muscle 

double CF = AF/VF; // Concentration of parent drug in the tissue of fat
double CVF = CF/PF; // Concentration of parent drug in the capillary blood of  fat

double CR = AR/VR; // Crest drug concentration in the rest of the body (mg/L)
double CVR = CR/Prest; // Concentration of parent drug in the capillary blood of the rest of body

double CA1 = AA1/Vblood; // CAfree concentration of unbound drug in the arterial blood (mg/L)
double CAfree1 = CA1*Free1;

double CL1 = AL1/VL; // Concentration of parent drug in the tissue of liver
double CVL1 = CL1/PL1; // Concentration of parent drug in the capillary blood of liver

double CK1 = AK1/VK; // Concentration of parent drug in the tissue of kidney
double CVK1 = CK1/PK1; // Concentration of parent drug in the capillary blood of kidney 

double CR1 = AR1/VR1; // Crest drug concentration in the rest of the body (mg/L)
double CVR1 = CR1/Prest1; // Concentration of parent drug in the capillary blood of the rest of body

// {flunixin distribution in each compartment}
// flunixin in venous blood compartment
double CV = (QL*CVL+QK*CVK+QF*CVF+QM*CVM+QR*CVR+RVS)/QC; // CV the drug concentration in the venous blood (mmol/L)
double RA = QC*(CV-CAfree);
dxdt_AA = RA;
double CVppm = CV*MW;// 
dxdt_AUCCV = CV; // mg/L

double CV1 = (QL*CVL1+QK*CVK1+QR1*CVR1)/QC; // CV1 the 5-OH FLU concentration in the venous blood (mmol/L)
double RA1 = QC*(CV1-CAfree1);
dxdt_AA1 = RA1;
double CV1ppm = CV1*MW1;// 
dxdt_AUCCV1 = CV1; // mg/L

// FLU and 5-OH FLU in liver compartment, flow-limited model
double Rmet1 = Kmet*CL*VL; // Rmet the metabolic rate in liver (mg/h)
dxdt_Amet = Rmet1; // Amet the amount of drug metabolized in liver (mg)
double Rehc = Kehc*AI;
dxdt_Aehc = Rehc;
double Rbile = Kbile*CVL;
dxdt_Abile = Rbile;
double RL = QL*(CAfree-CVL)-Rmet1-Rbile+Rehc; // RL the changing rate of the amount of drug in liver (mg/h)
dxdt_AL = RL; // AL amount of drug in liver (mg) 
double CLppm = CL*MW; // 
dxdt_AUCCL = CL; // AUCCL area under the curve of drug concentration in liver (mg*h/L)

double Rbile1 = Kbile1*CVL1;
dxdt_Abile1 = Rbile1;
double RL1 = QL*(CAfree1-CVL1)+Rmet1-Rbile1; // RL1 the changing rate of the amount of drug in liver (mg/h)
dxdt_AL1 = RL1; // AL1 amount of drug in liver (mg) 
double CL1ppm = CL1*MW1; // 
dxdt_AUCCL1 = CL1; // AUCCL1 area under the curve of drug concentration in liver (mg*h/L)

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
dxdt_AUCCK = CK; // AUCCK AUC of drug concentration in kidney (mg*h/L)
double CKppm = CK*MW; // 

double Rurine1 = Kurine1*CVK1;
dxdt_Aurine1 = Rurine1;
double RK1 = QK*(CAfree1-CVK1)-Rurine1; // RK1 the changing rate of the amount of drug in kidney (mg/h)
dxdt_AK1 = RK1; // AK amount of drug in kidney (mg)
dxdt_AUCCK1 = CK1; // AUCCK AUC of drug concentration in kidney (mg*h/L)
double CK1ppm = CK1*MW1; //

// flunixin in muscle compartment, flow-limited model
double RM = QM*(CAfree-CVM); //  RM the changing rate of the amount of drug in muscle (mg/h) 
dxdt_AM = RM; // AM amount of the drug in muscle (mg)
dxdt_AUCCM = CM;
double CMppm = CM*MW;

// flunixin in fat compartment, flow-limited model
double RF = QF*(CAfree-CVF); // RF the changing rate of the amount of drug in fat (mg/h)
dxdt_AF = RF; // AF amount of the drug in fat (mg)
dxdt_AUCCF = CF; // AUCCF AUC of drug concentration in fat (mg*h/L)
double CFppm = CF*MW;

// flunixin in the compartment of rest of body, flow-limited model
double RR = QR*(CAfree-CVR); // Rrest the changing rate of the amount of drug in the rest of the body (mg/h)
dxdt_AR = RR; // Arest amount of the drug in the rest of the body (mg)
dxdt_AUCCR = CR; // AUCCrest AUC of drug concentration in the rest of the body (mg*h/L)
double CRppm = CR*MW;

double RR1 = QR1*(CAfree1-CVR1); // Rrest the changing rate of the amount of drug in the rest of the body (mg/h)
dxdt_AR1 = RR1; // Arest amount of the drug in the rest of the body (mg)
dxdt_AUCCR1 = CR1; // AUCCrest AUC of drug concentration in the rest of the body (mg*h/L)
double CR1ppm = CR1*MW1;

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
capture solventvol2=solventvol;
capture as    = AS;
capture asc     = ASC;
capture avs     = AVS;
'



## Load Model
mod <- mcode_cache("pbpk", PBPK)

## Set Up Simulation Subjects
##There will be N simulation subjects, with their own sampled parameters that will override the default values.  The PBPK model does not use the Omega block to simulate invidual values, instead we fill the idata dataframe with these values.  Here we only change body weight, but other parameters could also be varied.

N=1000
PDOSEiv = 2.2
PDOSEim = 2.2
PDOSEsc = 2.2
ndoses= 1
doseinterval=24
timeiv = 0.005
timeim = 0.01
timesc = 0.01
MW = 296.24 # mg/mmol Molecular weight of flunixin
MW1 = 312.24 # mg/mmol Molecular weight of 5-OH FLU
End_time      = 15
BW            = 223 # KG
PDOSEde       = 3.3 # mg/kg/day
DOSE          = PDOSEde*BW/MW
area          = 4800

sd.log.PL = sqrt(log(1 + 2.104 ^ 2 / 10.52 ^ 2))
m.log.PL = log(10.52) - 0.5 * sd.log.PL ^ 2

sd.log.PL1 = sqrt(log(1 + 1.852 ^ 2 / 9.26 ^ 2))
m.log.PL1 = log(9.26) - 0.5 * sd.log.PL1 ^ 2

sd.log.PK = sqrt(log(1 + 0.8 ^ 2 / 4 ^ 2))
m.log.PK = log(4) - 0.5 * sd.log.PK ^ 2

sd.log.PK1 = sqrt(log(1 + 0.8 ^ 2 / 4 ^ 2))
m.log.PK1 = log(4) - 0.5 * sd.log.PK1 ^ 2

sd.log.PM = sqrt(log(1 + 0.1 ^ 2 / 0.5 ^ 2))
m.log.PM = log(0.5) - 0.5 * sd.log.PM ^ 2

sd.log.PF = sqrt(log(1 + 0.12 ^ 2 / 0.6 ^ 2))
m.log.PF = log(0.6) - 0.5 * sd.log.PF ^ 2

sd.log.Prest = sqrt(log(1 + 1.6 ^ 2 / 8 ^ 2))
m.log.Prest = log(8) - 0.5 * sd.log.Prest ^ 2

sd.log.Prest1 = sqrt(log(1 + 1 ^ 2 / 5 ^ 2))
m.log.Prest1 = log(5) - 0.5 * sd.log.Prest1 ^ 2

sd.log.KmC = sqrt(log(1 + 0.04^ 2 / 0.2 ^ 2))
m.log.KmC = log(0.2) - 0.5 * sd.log.KmC ^ 2

sd.log.KurineC = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KurineC = log(0.1) - 0.5 * sd.log.KurineC ^ 2

sd.log.KurineC1 = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KurineC1 = log(0.1) - 0.5 * sd.log.KurineC1 ^ 2

sd.log.KehcC = sqrt(log(1 + 0.045^ 2 / 0.15 ^ 2))
m.log.KehcC = log(0.15) - 0.5 * sd.log.KehcC ^ 2

sd.log.KbileC = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KbileC = log(0.1) - 0.5 * sd.log.KbileC ^ 2

sd.log.KbileC1 = sqrt(log(1 + 0.03^ 2 / 0.1 ^ 2))
m.log.KbileC1 = log(0.1) - 0.5 * sd.log.KbileC1 ^ 2

sd.log.Kint = sqrt(log(1 + 0.12^ 2 / 0.4 ^ 2))
m.log.Kint = log(0.4) - 0.5 * sd.log.Kint ^ 2

sd.log.Kfeces = sqrt(log(1 + 0.15^ 2 / 0.5 ^ 2))
m.log.Kfeces = log(0.5) - 0.5 * sd.log.Kfeces ^ 2

sd.log.logPsw = sqrt(log(1 + 0.426^ 2 /2.13 ^ 2))
m.log.logPsw = log(2.13) - 0.5 * sd.log.logPsw ^ 2

sd.log.logP = sqrt(log(1 + 0.24^ 2 / 1.2 ^ 2))
m.log.logP = log(1.2) - 0.5 * sd.log.logP ^ 2

sd.log.PB = sqrt(log(1 + 0.285^ 2 / 0.95 ^ 2))
m.log.PB = log(0.95) - 0.5 * sd.log.PB ^ 2


set.seed(11009)
idata <- 
  data_frame(ID=1:N) %>% 
  mutate(
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = 223, sd = 0.3*223), 
                    b = qnorm(0.975, mean = 223, sd = 0.3*223), 
                    mean = 223,
                    sd = 0.3*223),
    
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 8.543, sd = 8.543*0.224),
      b = qnorm(0.975, mean = 8.543, sd = 8.543*0.224),
      mean = 8.543,
      sd = 8.543*0.224),  # Cardiac output index (L/h/kg)
    # Fracion of blood flow to organs (unitless) 
    QLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.273, sd = 0.3*0.273),
      b = qnorm(0.975, mean = 0.273, sd = 0.3*0.273),
      mean = 0.273,
      sd = 0.3*0.273
    ),
    # Fraction of blood flow to the kidneys (2016 Lin)
    QKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.116, sd =  0.149*0.116),
      b = qnorm(0.975, mean = 0.116, sd =  0.149*0.116),
      mean = 0.116,
      sd =  0.149*0.116
    ),
    #
    QMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.293, sd =  0.144*0.293),
      b = qnorm(0.975, mean = 0.293, sd =  0.144*0.293),
      mean = 0.293,
      sd =  0.144*0.293
    ),
    # Fraction of blood flow to the fat (2016 Lin)
    QFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.128, sd =  0.128*0.3),
      b = qnorm(0.975, mean = 0.128, sd =  0.128*0.3),
      mean = 0.128,
      sd =  0.128*0.3
    ),
    # Fraction of blood flow to the rest of body (total sum equals to 1)
    QrestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.19, sd = 0.19*0.3),
      b = qnorm(0.975, mean = 0.19, sd = 0.19*0.3),
      mean = 0.19,
      sd = 0.19*0.3
    ),
    QrestC1a = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.611, sd = 0.611*0.3),
      b = qnorm(0.975, mean = 0.611, sd = 0.611*0.3),
      mean = 0.611,
      sd = 0.611*0.3
    ),
    
    # Fractional organ tissue volumes (unitless)
    # Fractional liver tissue (1933 Swett)
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.023, sd =  0.023*0.015),
      b = qnorm(0.975, mean = 0.023, sd =  0.023*0.015) ,
      mean = 0.023,
      sd =  0.023*0.015
    ),
    # Fractional kidney tissue (1933 Swett)
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.0045, sd =  0.0045*0.038),
      b = qnorm(0.975, mean = 0.0045, sd =  0.0045*0.038),
      mean = 0.0045,
      sd =  0.0045*0.038
    ),
    # Fractional fat tissue (2016 Lin, 2014 Leavens)
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.235, sd =  0.235*0.077),
      b = qnorm(0.975, mean = 0.235, sd =  0.235*0.077),
      mean = 0.235,
      sd = 0.235*0.077
    ),
    # Fractional muscle tissue (2016 Lin, 2014 Leavens)
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.355, sd =  0.355*0.007),
      b = qnorm(0.975, mean = 0.355, sd =  0.355*0.007),
      mean = 0.355,
      sd =  0.355*0.007
    ),
    # Venous blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VvenCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.044, sd =  0.044*0.3),
      b = qnorm(0.975, mean = 0.044, sd =  0.044*0.3),
      mean = 0.044,
      sd =  0.044*0.3
    ),
    # Arterial blood volume, fraction of blood volume (2016 Lin# 2008 Leavens)
    VartCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.016, sd = 0.016*0.3),
      b = qnorm(0.975, mean = 0.016, sd = 0.016*0.3),
      mean = 0.016,
      sd =  0.016*0.3
    ),
    # Fractional rest of body (total sum equals to 1)
    VrestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.3225, sd =  0.3225*0.3),
      b = qnorm(0.975, mean = 0.3225, sd =  0.3225*0.3),
      mean = 0.3225,
      sd =  0.3225*0.3
    ),
    VrestC1a = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.9125, sd =  0.9125*0.3),
      b = qnorm(0.975, mean = 0.9125, sd =  0.9125*0.33),
      mean = 0.9125,
      sd =  0.9125*0.3
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
   
    
    DOSEiv = PDOSEiv*BW/MW, DOSEim = PDOSEim*BW/MW, DOSEsc = PDOSEsc*BW/MW,DOSEde=PDOSEde*BW/MW)
    idata %>% head %>% kable

## Set up dosing.
##We'll give a fixed dose of 100 mg, given QDx5 into the SC compartment.  Two dose records are used, one to the Amtsitesc and one to DOSEppgsc, and the Fracsc parameter modifies the dose in each compartment.
##I noted some trick used in the original code to apply an infusion over the integration interval.  This is not necessary if the reason was to avoid inserting a state change in LSODA.  With mrgsolve we can apply dose events and they will stop the integrator to safely change the state.

ev1 <- ev(ID=1:N, amt= idata$DOSEiv, ii=doseinterval, addl=ndoses-1, cmt="AA", rate=(idata$DOSEiv)/timeiv, replicate = FALSE)
ev2 <- ev(ID=1:N, amt= idata$DOSEim, ii=doseinterval, addl=ndoses-1, cmt="Amtsiteim", rate=(idata$DOSEim)/timeim, replicate = FALSE)
ev3 <- ev(ID=1:N, amt= idata$DOSEsc, ii=doseinterval, addl=ndoses-1, cmt="Amtsitesc", rate=(idata$DOSEsc)/timesc, replicate = FALSE)
ev4 <- ev (ID=1:N, amt= idata$DOSEde, ii = doseinterval,addl = ndoses-1, cmt  = "AS", replicate = FALSE)
ex <- ev4 #ev1 + ev2 + ev3+ 

dim(as.data.frame(ex))
#a <- as.data.frame(ex)
## set up time grid
##We'll sample every hour up to 8 days past the time of last dose

tsamp=tgrid(0,doseinterval*ndoses+15*24,1)

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
  mutate(Time1 = time/24 ) %>%
  dplyr::select(ID, Time1, Venous:Kidney1) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Venous, 0.025), CV50 = median(Venous), CV99 = quantile(Venous, 0.975), 
                   CL1 = quantile(Liver, 0.025), CL50 = median(Liver), CL99 = quantile(Liver, 0.975),
                   CK1 = quantile(Kidney, 0.025), CK50 = median(Kidney), CK99 = quantile(Kidney, 0.975),
                   CM1 = quantile(Muscle, 0.025), CM50 = median(Muscle), CM99 = quantile(Muscle, 0.975),
                   CF1 = quantile(Fat, 0.025), CF50 = median(Fat), CF99 = quantile(Fat, 0.975),
                   CV11 = quantile(Venous1, 0.025), CV150 = median(Venous1), CV199 = quantile(Venous1, 0.975), 
                   CL11 = quantile(Liver1, 0.025), CL150 = median(Liver1), CL199 = quantile(Liver1, 0.975),
                   CK11 = quantile(Kidney1, 0.025), CK150 = median(Kidney1), CK199 = quantile(Kidney1, 0.975))

data_tissue <- as.data.frame(out.sum)

end_time = 16


## Plot tissue concentrations
ggplot(data_tissue, aes(x=Time1)) + 
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
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-4):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  annotate("text", x = 8, y = 1, label = paste("Liver"), 
           hjust = -0.5, vjust = 0, angle = 0, size = 6) +
  geom_line(aes(y = 0.125),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16))

ggplot(data_tissue, aes(x=Time1)) + 
  geom_ribbon(aes(ymin = CM1,
                  ymax = CM99
  ), fill = 'red',
  show.legend = T, 
  size = 0.2,
  alpha = 0.3) +  # alpha is transparency parameter. 
  geom_line(aes(y = CM50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CM99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CM1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-4):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  annotate("text", x = 8, y = 1, label = paste("Muscle"), 
           hjust = -0.5, vjust = 0, angle = 0, size = 6) +
  geom_line(aes(y = 0.020),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16))

ggplot(data_tissue, aes(x=Time1)) + 
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
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-4):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  
  annotate("text", x = 8, y = 1, label = paste("Fat"), 
           hjust = -0.5, vjust = 0, angle = 0, size = 6) +
  geom_line(aes(y = 0.03),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16))

ggplot(data_tissue, aes(x=Time1)) + 
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
  scale_x_continuous(name = c('Time (Day)'),breaks = c((-4):end_time)) + 
  ylab(expression(paste('Concentration (',mu,'g/mL)'))) + 
  
  annotate("text", x = 8, y = 1, label = paste("Kidney"), 
           hjust = -0.5, vjust = 0, angle = 0, size = 6) +
  geom_line(aes(y = 0.1),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  scale_y_log10() + theme_bw() + theme(axis.text=element_text(size=16))

run.time