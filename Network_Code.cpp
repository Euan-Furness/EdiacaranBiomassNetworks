// This code is assocaited with "Sloppy feeding by predators destabilises early animal food webs in the Ediacaran".
// In order to run this code, please modify the filepath located at flag "A)", and modify the settings located at flags "B)",
// "C)", "D)", and "E)" as required to replicate the desired experiments.

#pragma comment(linker, "/HEAP:8000000")
#pragma comment(linker, "/STACK:8000000")

#include <iostream>
#include <fstream>
#include <math.h>
#include <C:/Eigen/Eigenvalues>
#include <time.h>

using namespace std;

#define maxSpecies 16 // 16 for Nilpena
#define maxLoops 8192 //Static maximum number of loops. Should be 2^(maxSpecies+1), or 2*(maxSpecies^3).
#define nReplicates 216000 //Static number of loops performed by main() to explore parameter space.

bool WeCareAboutLoopWeights = true;
bool WeREALLYCareAboutLoopWeights = false;
bool ThreeLengthLoopsOnly = true;
bool WriteLog = true;
bool ConsoleOutputs = false;
bool ShortCircuitsAreNotAllowed = false;
bool WeightsAreScaledByd = true; // Boolean value which controls whether or not loop weights are divided through by death rates or not
bool OutputdValues = false;
bool OutputEassValues = false;
bool OutputEprodValues = false;
bool NetLoopEffects = false; //If true, loop weights are calculated by subtracting one loop from its inverse. Otherwise, they are absolute values.

int ReplicateNumber; //Global to prevent compile errors due to saved array loading blocks at the bottom of the script

//Globals to avoid stack overflow
double StabilityArray[nReplicates];

class SpeciesNode {
    public:
    string NodeName = "Name Not Assigned";
    int ID = -1; //A unique ID for the species, must be less than maxSpecies, and should always be the same as the species' position in the interactions arrays. -1 is the flag for nonexistant species
    
    double EqBiomass = 0; //B: Equilibrium biomass in the species
    double NaturalDeathRate = 0; //d: Death rate of the species independent of consumption.
    double AssimilationEfficiency = 0; //e^ass: Efficiency with which the species consumes acquired resources.
    double MetabolicEfficiency = 0; //e^prod: Efficiency with which the species converts acquired resources into biomass.
    double FeedingRate = 0; //Fi
    double UnnaturalDeathRate = 0; //Mi: death rate from feeding. Sum of all Fij values.
    double EffectOnDetritus = 0; //aDj: the effect of the species on detritus.
    double TrophicLevel = 0; //The trophic levels of the species, as determined based on linkages.

    bool isBag = false; //Bool dictates whether or not a species is recorded as a bag.

    bool RandomiseBiomass = false; //Bool tracks whether or not the species' biomass should be randomised during batch runs.
    bool RandomiseBagness = false; //Bool tracks whether or not the species' bag (vs solid) tissue identity is randomised during batch runs.
    bool MiDefined = false; //Bool tracks whether or not Mi has been calculated for this species yet.
    bool FiDefined = false; //Bool tracks whether or not Fi has beeen calculated for this species yet.

    int CreateSpecies(std::string, double, double, double, double); //Returns the ID of the species created.
    void KillSpecies(int);
    void CalculateFi();
    void CalculateMi();
};

class TrophicInteraction {
    public: //Here, -1 is the flag for no interaction, but 0 also denotes no interaction and will come up more often
    int consumer = -1;
    int resource = -1;
    double affinity100 = 0; //wij * 100; multiplication allows for handling as an int
    double rate = 0; //Fij parameter
    bool FijDefined = false; //Bool tracks whether or not Fij has been calculated for this interaction yet.
    double InteractionStrengthAonB = 0; //The aij value for the interaction
    void AssignLinkage(double);
    int CalculateFij(int);
};

class IntList {
    public:
    int Entries[maxSpecies];
    void append(int);
    void prefill();
    int findSpecies(int); //This function compensates for the empty -1 values in the loop object by returning the nth non-empty value
};

class DoubleList {
    public:
    double Entries[maxSpecies];
    void append(double);
    void prefill();
};

class Loop {
    public:
    int LoopID = -1; //ID value of the loop itself
    int partnerID = -1; //ID value of the reverse loop.
    int Species[maxSpecies]; //ID values of species in the loop
    int length = 0; //Number of species in the loop
    double weight = 0; //Loop weight
    double net_weight = 0; //Difference between the absolute weight of this loop and the absolute weight of the same loop in he other direction.
    int type = 1; //The loop form: 0: detrital; 1: trophic
    int polarity = 0; //The polarity of the loop: -1: negative feedback, 0: positive feedback
    IntList pathway;
    void Reset(); //This function functionally removes a loop from the array
};

//SpeciesNode* SpeciesArray = (SpeciesNode*) malloc(62 * sizeof(SpeciesNode));

SpeciesNode SpeciesArray[maxSpecies];
TrophicInteraction TrophicInteractionsArray[maxSpecies][maxSpecies];
int NextID = 0;

Loop LoopArray[maxLoops]; //Global to ensure access throughout the program.
int WeightOrderedLoops[maxLoops]; //A loop to order the trophic loops by weight, for display. Global prevents stack overflow.
int NextLoopID = 0;

bool IsConnected(int SpeciesA, int SpeciesB, int DetritusID) {
    if (TrophicInteractionsArray[SpeciesA][SpeciesB].InteractionStrengthAonB != 0) {
        //cout << "Trophic connection between " << SpeciesArray[SpeciesA].NodeName << " and " << SpeciesArray[SpeciesB].NodeName << endl;
        return true;
    }
    else if (((SpeciesA == DetritusID) & (SpeciesB != DetritusID)) || ((SpeciesA != DetritusID) & (SpeciesB == DetritusID))) { //Detritus is connected to everything through the detrital pathways
        //cout << "Detrital connection between " << SpeciesArray[SpeciesA].NodeName << " and " << SpeciesArray[SpeciesB].NodeName << endl;
        return true; //XOR is used here because otherwise detritus would always short circuit the loop by being connected to itself
    }
    else {
        //cout << "No connection between " << SpeciesArray[SpeciesA].NodeName << " and " << SpeciesArray[SpeciesB].NodeName << endl;
        return false;
    }
}

IntList SearchPathways(int FirstSpecies, int CurrentSpecies, int PreviousSpecies, IntList SpeciesToFind, IntList SpeciesFound, int DetritalBox) {
    IntList BlankReturn;
    BlankReturn.prefill();
    if (SpeciesToFind.findSpecies(1) == -1) { //i.e. if there are no entries left in the species to find list
        if (IsConnected(FirstSpecies, CurrentSpecies, DetritalBox)) { //Detrital box argument is required so that connections through the detritus box can be checked for
            return SpeciesFound;
        }
        else {
            return BlankReturn;
        }
    }

    if (ShortCircuitsAreNotAllowed) {
        for (int n = 0; n < maxSpecies; n++) { //Short circuit detector
            if ((SpeciesFound.Entries[n] != -1) & (SpeciesFound.Entries[n] != PreviousSpecies) & (SpeciesFound.Entries[n] != CurrentSpecies)) {
                if (IsConnected(CurrentSpecies, SpeciesFound.Entries[n], DetritalBox)) { //Short circuit in the loop
                    return BlankReturn;
                }
            }
        }
    }

    if (SpeciesToFind.findSpecies(1) != -1) { //Otherwise, we need to try to find another species in the SpeciesToFind list
        for (int n = 0; n < maxSpecies; n++) {
            if (SpeciesToFind.Entries[n] != -1) { //We can try to find this species
                if (IsConnected(CurrentSpecies, SpeciesToFind.Entries[n], DetritalBox)) {
                    IntList NewSpeciesToFind = SpeciesToFind;
                    IntList NewSpeciesFound = SpeciesFound;
                    NewSpeciesToFind.Entries[n] = -1;
                    NewSpeciesFound.append(SpeciesToFind.Entries[n]);
                    IntList FoundRoute = SearchPathways(FirstSpecies, SpeciesToFind.Entries[n], CurrentSpecies, NewSpeciesToFind, NewSpeciesFound, DetritalBox); //Recursive
                    bool ARouteWasFound = false;
                    for (int n = 0; n < maxSpecies; n++) {
                        if (FoundRoute.Entries[n] != -1) {
                            ARouteWasFound = true;
                        }
                    }
                    if (ARouteWasFound) {
                        return FoundRoute;
                    }
                }
            }
        }
    }
    return BlankReturn;
}

void IntList::append(int NextValue) {
    for (int n = 0; n < maxSpecies; n++) {
        if (Entries[n] == -1) {
            Entries[n] = NextValue;
            break;
        }
    }
}

void IntList::prefill() {
    for (int n = 0; n < maxSpecies; n++) {
        Entries[n] = -1;
    }
}

void DoubleList::append(double NextValue) {
    for (int n = 0; n < maxSpecies; n++) {
        if (Entries[n] == -1) {
            Entries[n] = NextValue;
            break;
        }
    }
}

void DoubleList::prefill() {
    for (int n = 0; n < maxSpecies; n++) {
        Entries[n] = -1;
    }
}

int UpdateNextID() {
    for (int i = 0; i < maxSpecies; i++) {
        if (SpeciesArray[i].ID == -1) {
            return i;
        }
    }
    if (ReplicateNumber == 0) {
        std::cout << "Warning: max species reached." << endl;
    }
    return -1;
}

int IntList::findSpecies(int whichSpecies) {
    int CurrentSpecies = 0;
    for (int n = 0; n < maxSpecies; n++) {
        if (Entries[n] != -1) {
            CurrentSpecies++;
            if (CurrentSpecies == whichSpecies) {
                return n;
            }
        }
    }
    return -1;
}

void Loop::Reset() {
    for (int i = 0; i < maxSpecies; i++) {
        Species[i] = -1;
    }
    length = 0;
    LoopID = -1;
    partnerID = -1;
    weight = 0;
    type = 1;
    polarity = 0;
    pathway.prefill();
}

int SpeciesNode::CreateSpecies(string PassedNodeName, double PassedEqBiomass, double PassedNaturalDeathRate, double PassedAssimilationEfficiency, double PassedMetabolicEfficiency) {
    ID = NextID;
    NextID = UpdateNextID();

    NodeName = PassedNodeName;
    EqBiomass = PassedEqBiomass;
    NaturalDeathRate = PassedNaturalDeathRate;
    AssimilationEfficiency = PassedAssimilationEfficiency;
    MetabolicEfficiency = PassedMetabolicEfficiency;

    return ID;
}

void SpeciesNode::KillSpecies(int IDIndex) {
    ID = -1;
    NodeName = "Name Not Assigned";
    EqBiomass = 0;
    NaturalDeathRate = 0;
    AssimilationEfficiency = 0;
    MetabolicEfficiency = 0;
    FeedingRate = 0;
    UnnaturalDeathRate = 0;
    TrophicLevel= 0;
    EffectOnDetritus = 0;
    MiDefined = false;
    FiDefined = false;
    RandomiseBiomass = false;
    RandomiseBagness = false;
    isBag = false;
    NextID = UpdateNextID();
        for (int j = 0; j < maxSpecies; j++) {
            TrophicInteractionsArray[IDIndex][j].affinity100 = 0;
            TrophicInteractionsArray[IDIndex][j].consumer = -1;
            TrophicInteractionsArray[IDIndex][j].resource = -1;
            TrophicInteractionsArray[IDIndex][j].rate = 0;
            TrophicInteractionsArray[IDIndex][j].InteractionStrengthAonB = 0;
            TrophicInteractionsArray[IDIndex][j].FijDefined = false;
            TrophicInteractionsArray[j][IDIndex].affinity100 = 0;
            TrophicInteractionsArray[j][IDIndex].consumer = -1;
            TrophicInteractionsArray[j][IDIndex].resource = -1;
            TrophicInteractionsArray[j][IDIndex].rate = 0;
            TrophicInteractionsArray[j][IDIndex].InteractionStrengthAonB = 0;
            TrophicInteractionsArray[j][IDIndex].FijDefined = false;
        }
}

void TrophicInteraction::AssignLinkage(double PassedAffinity100) {
    affinity100 = PassedAffinity100;
}

void SpeciesNode::CalculateFi() {
    if (MiDefined) {
        FeedingRate = ((NaturalDeathRate * EqBiomass) + UnnaturalDeathRate) / (AssimilationEfficiency * MetabolicEfficiency);
        FiDefined = true;
    }
}

int TrophicInteraction::CalculateFij(int PriorAssignedFij) {

    if ((SpeciesArray[consumer].FiDefined) & (FijDefined == false)) {
        double FijDenominator = 0;
        for (int i = 0; i < maxSpecies; i++) {
            FijDenominator += (TrophicInteractionsArray[i][consumer].affinity100 * SpeciesArray[i].EqBiomass);
        }
        if (FijDenominator > 0.00000000001) {
            rate = TrophicInteractionsArray[resource][consumer].affinity100 * SpeciesArray[resource].EqBiomass * SpeciesArray[consumer].FeedingRate / FijDenominator;
        }
        PriorAssignedFij++;
        FijDefined = true;
    }
    
    return PriorAssignedFij;
}

void SpeciesNode::CalculateMi() {
    double Mi = 0;
    bool AllFijFound = true;
    for (int i = 0; i < maxSpecies; i++) {
        if (TrophicInteractionsArray[ID][i].affinity100 > 0.000000000001) {
            if (TrophicInteractionsArray[ID][i].FijDefined) {
                Mi += TrophicInteractionsArray[ID][i].rate;
            }
            else {
                AllFijFound = false;
            }
        }
    }
    if (AllFijFound) {
        UnnaturalDeathRate = Mi;
        MiDefined = true;
    }
}

Loop HeaviestPositive3Loop;
Loop HeaviestNegative3Loop;
Loop HeaviestPositiveLoop;
Loop HeaviestNegativeLoop;
float HeaviestPositive3LoopWeight;
float HeaviestNegative3LoopWeight;
float HeaviestPositiveLoopWeight;
float HeaviestNegativeLoopWeight;
string NamesList[maxSpecies]; //Global to allow access in the output block

int main() {

    if ((WeREALLYCareAboutLoopWeights == true) & (WeCareAboutLoopWeights == false)) {
        cout << "Do we care about loop weights or not? Inconsistent settings.\n";
        return 0;
    }

    ofstream OutputFile; // A) Modify the below filepath as required; a text file will be created there, containing the output data cube.
    if (WriteLog) {
        OutputFile.open("C:/Users/Local User/Data Outputs/Outputs/Network Sims/Experimental for manuscript/Cubes/Systematically Varied Parameter Values/Nilpena_1TF/output.txt");
    }

    srand(time(NULL));

    if ((maxLoops != pow(2, maxSpecies + 1)) & (ThreeLengthLoopsOnly == false)) {
        std::cout << "Warning! maxLoops and maxSpecies may not be aligned. ";
        std::cout << "Expected maxLoops value of " << pow(2, maxSpecies + 1) << " for a maxSpecies value of " << maxSpecies << "." << endl;
    }
    else if ((maxLoops != 2 * pow(maxSpecies, 3)) & (ThreeLengthLoopsOnly)) {
        std::cout << "Warning! maxLoops and maxSpecies may not be aligned. ";
        std::cout << "Expected maxLoops value of " << 2 * pow(maxSpecies, 3) << " for a maxSpecies value of " << maxSpecies << "." << endl;
    }

    int CubeDim1 = -1;
    int CubeDim2 = -1;
    int CubeDim3 = -1;
    ReplicateNumber = 0;
    while (ReplicateNumber < nReplicates) { //Loop is not indented. Includes all of main.

    CubeDim1++;
    if ((ReplicateNumber % 60) == 0) {
        CubeDim1 = 0;
        CubeDim2++;
    }
    if ((ReplicateNumber % 3600) == 0) {
        CubeDim2 = 0;
        CubeDim3++;
    }

    for (int i = 0; i < maxSpecies; i++) {
        SpeciesArray[i].KillSpecies(i);
    }

    for (int i = 0; i < maxSpecies; i++) { //Initial linkages are all absent
        for (int j = 0; j < maxSpecies; j++) {
            TrophicInteractionsArray[i][j].consumer = j;
            TrophicInteractionsArray[i][j].resource = i;
        }
    }
    for (int n = 0; n < maxLoops; n++) { //Initial loops are all empty (species = -1)
        LoopArray[n].Reset();
        for (int i = 0; i < maxSpecies; i++) {
            LoopArray[n].Species[i] = -1;
        }
    }

    //Begin array loading

    int propOsmo = 25;

    // B) The below bool controls the presence or absence of model predators in the network.
    bool ArtificialPredators = true;
    // C) The below float controls the biomass of predators in the network.
    float PredatorBiomass = 8;
    float dPred = exp(1 - (0.25 * log(5)));
    // D) The below bool controls whether or not predators feed on detritus in addition to live prey.
    bool FacultativePredators = false;
    // E) The below bool controls whether or not microorganisms have high death rates
    bool FastPlankton = false;

    // All solid tissue values
    float TriradialPIB = 0.22986; //Tribrachidium and Rugoconites
    float FrondPIB = 13.8226;
    float SpriggPIB = 0.1781; //Spriggina and Andiva
    float ParvanPIB = 0.084667;
    float SpongePIB = 0.79331; //Palaeopascichnus and Coronacollina
    float FunisiaPIB = 8.29455;
    float DickinPIB = 0.435424;
    float AlgaePIB = 11.517;
    float AuloPIB = 259.18134;

    //Hardcoded per square metre abundance values
    float TriradialPerSqMAbundance = 0.641025641;
    float FrondPerSqMAbundance = 1.025641026;
    float SpriggPerSqMAbundance = 1.794871795;
    float ParvanPerSqMAbundance = 0.7692307692;
    float SpongePerSqMAbundance = 0.5982905983;
    float FunisiaPerSqMAbundance = 0.2136752137;
    float DickinPerSqMAbundance = 1.965811966;
    float AlgaePerSqMAbundance = 0.7692307692;
    float AuloPerSqMAbundance = 1.068376068;
    
    float DOCBiomass = 240; //2.4 for modern, 1200 for Ediacaran? Maybe 240 is better?
    float AutoBiomass = 21;
    float HeteroBiomass = 11;
    float AmoebaeBiomass = 1.5;
    float MatBiomass = 2000; //600 is low. Appropriate for the Avalon, but should probably be an order of magnitude higher in the Nama

    //Total system biomass calculations
    float TriradialBiomass = TriradialPIB * TriradialPerSqMAbundance;
    float FrondBiomass = FrondPIB * FrondPerSqMAbundance;
    float SpriggBiomass = SpriggPIB * SpriggPerSqMAbundance;
    float ParvanBiomass = ParvanPIB * ParvanPerSqMAbundance;
    float SpongeBiomass = SpongePIB * SpongePerSqMAbundance;
    float FunisiaBiomass = FunisiaPIB * FunisiaPerSqMAbundance;
    float DickinBiomass = DickinPIB * DickinPerSqMAbundance;
    float AlgaeBiomass = AlgaePIB * AlgaePerSqMAbundance;
    float AuloBiomass = AuloPIB * AuloPerSqMAbundance;

    //Natural death rate calculations
    float dTriradial = exp(1 - (0.25 * log(TriradialPIB)));
    float dFrond = exp(1 - (0.25 * log(FrondPIB)));
    float dSprigg = exp(1 - (0.25 * log(SpriggPIB)));
    float dParvan = exp(1 - (0.25 * log(ParvanPIB)));
    float dSponge = exp(1 - (0.25 * log(SpongePIB)));
    float dFunisia = exp(1 - (0.25 * log(FunisiaPIB)));
    float dDickin = exp(1 - (0.25 * log(DickinPIB)));
    float dAlgae = exp(1 - (0.25 * log(AlgaePIB)));
    float dAulo = exp(1 - (0.25 * log(AuloPIB)));

    //Species creation
    int DOC = SpeciesArray[NextID].CreateSpecies("DOC",DOCBiomass,0,1,1);
    int Autotrophic_Plankton;
    int Heterotrophic_Plankton;
    int Amoebae;
    if (FastPlankton) {
        Autotrophic_Plankton = SpeciesArray[NextID].CreateSpecies("Autotrophic plankton",AutoBiomass,12,1,1);
        Heterotrophic_Plankton = SpeciesArray[NextID].CreateSpecies("Heterotrophic plankton",HeteroBiomass,12,1,0.4);
        Amoebae = SpeciesArray[NextID].CreateSpecies("Amoebae",AmoebaeBiomass,60,0.95,0.4);
    }
    else {
        Autotrophic_Plankton = SpeciesArray[NextID].CreateSpecies("Autotrophic plankton",AutoBiomass,1.2,1,1);
        Heterotrophic_Plankton = SpeciesArray[NextID].CreateSpecies("Heterotrophic plankton",HeteroBiomass,1.2,1,0.4);
        Amoebae = SpeciesArray[NextID].CreateSpecies("Amoebae",AmoebaeBiomass,6,0.95,0.4);
    }
    int Microbial_Mat = SpeciesArray[NextID].CreateSpecies("Microbial mat",MatBiomass,1.2,1,1);

    int Triradial = SpeciesArray[NextID].CreateSpecies("Triradial",TriradialBiomass,dTriradial,0.8,0.4);
    int Frond = SpeciesArray[NextID].CreateSpecies("Frond",FrondBiomass,dFrond,0.8,0.4);
    int Sprigg = SpeciesArray[NextID].CreateSpecies("Sprigg",SpriggBiomass,dSprigg,0.6,0.27);
    int Parvan = SpeciesArray[NextID].CreateSpecies("Parvan",ParvanBiomass,dParvan,0.8,0.4); //Parvancorina is motile, but only just. Efficiency is coded as sessile.
    int Sponge = SpeciesArray[NextID].CreateSpecies("Sponge",SpongeBiomass,dSponge,0.8,0.4);
    int Funisia = SpeciesArray[NextID].CreateSpecies("Funisia",FunisiaBiomass,dFunisia,0.8,0.4);
    int Dickin = SpeciesArray[NextID].CreateSpecies("Dickin",DickinBiomass,dDickin,0.6,0.27);
    int Algae = SpeciesArray[NextID].CreateSpecies("Algae",AlgaeBiomass,dAlgae,1,1);
    int Aulo = SpeciesArray[NextID].CreateSpecies("Aulo",AuloBiomass,dAulo,0.8,0.4);

    int Predators;
    if (ArtificialPredators) {
        Predators = SpeciesArray[NextID].CreateSpecies("Predators",PredatorBiomass,dPred,0.8,0.27);
    }

    //Bag identities
    SpeciesArray[Aulo].isBag = true;

    for (int i = 0; i < maxSpecies; i++) {
        if (SpeciesArray[i].isBag) {
            SpeciesArray[i].EqBiomass = SpeciesArray[i].EqBiomass * 0.01;
            SpeciesArray[i].NaturalDeathRate = SpeciesArray[i].NaturalDeathRate * 3.16; // 3.16 calculated as the approximate conversion factor
        }
    }

    //Cube, in this case
    float DOCValues[61];
    float AutValues[61];
    float HetValues[61];
    DOCValues[0] = log(SpeciesArray[DOC].EqBiomass * 0.01);
    DOCValues[60] = log(SpeciesArray[DOC].EqBiomass * 100);
    AutValues[0] = log(SpeciesArray[Autotrophic_Plankton].EqBiomass * 0.01);
    AutValues[60] = log(SpeciesArray[Autotrophic_Plankton].EqBiomass * 100);
    HetValues[0] = log(SpeciesArray[Heterotrophic_Plankton].EqBiomass * 0.01);
    HetValues[60] = log(SpeciesArray[Heterotrophic_Plankton].EqBiomass * 100);
    for (int i = 1; i < 60; i++) {
        DOCValues[i] = DOCValues[0] + (float(i) / 60) * (DOCValues[60] - DOCValues[0]);
        AutValues[i] = AutValues[0] + (float(i) / 60) * (AutValues[60] - AutValues[0]);
        HetValues[i] = HetValues[0] + (float(i) / 60) * (HetValues[60] - HetValues[0]);
    }
    for (int i = 0; i < 61; i++) {
        DOCValues[i] = exp(DOCValues[i]);
        AutValues[i] = exp(AutValues[i]);
        HetValues[i] = exp(HetValues[i]);
    }
    SpeciesArray[DOC].EqBiomass = (DOCValues[CubeDim1 + 1] + DOCValues[CubeDim1]) / 2;
    SpeciesArray[Autotrophic_Plankton].EqBiomass = (AutValues[CubeDim2 + 1] + AutValues[CubeDim2]) / 2;
    SpeciesArray[Heterotrophic_Plankton].EqBiomass = (HetValues[CubeDim3 + 1] + HetValues[CubeDim3]) / 2;
    
    //Interaction Assignment
    TrophicInteractionsArray[Heterotrophic_Plankton][Amoebae].AssignLinkage(33);
    TrophicInteractionsArray[Autotrophic_Plankton][Amoebae].AssignLinkage(33);
    TrophicInteractionsArray[DOC][Amoebae].AssignLinkage(33);

    TrophicInteractionsArray[DOC][Heterotrophic_Plankton].AssignLinkage(50);
    TrophicInteractionsArray[Autotrophic_Plankton][Heterotrophic_Plankton].AssignLinkage(50);

    TrophicInteractionsArray[DOC][Triradial].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Triradial].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Triradial].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Triradial].AssignLinkage((100-propOsmo)/3);

    TrophicInteractionsArray[DOC][Frond].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Frond].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Frond].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Frond].AssignLinkage((100-propOsmo)/3);

    TrophicInteractionsArray[DOC][Sprigg].AssignLinkage(100);

    TrophicInteractionsArray[DOC][Parvan].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Parvan].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Parvan].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Parvan].AssignLinkage((100-propOsmo)/3);

    TrophicInteractionsArray[DOC][Sponge].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Sponge].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Sponge].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Sponge].AssignLinkage((100-propOsmo)/3);

    TrophicInteractionsArray[DOC][Funisia].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Funisia].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Funisia].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Funisia].AssignLinkage((100-propOsmo)/3);

    TrophicInteractionsArray[Microbial_Mat][Dickin].AssignLinkage(100);

    TrophicInteractionsArray[DOC][Aulo].AssignLinkage(propOsmo);
    TrophicInteractionsArray[Autotrophic_Plankton][Aulo].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Heterotrophic_Plankton][Aulo].AssignLinkage((100-propOsmo)/3);
    TrophicInteractionsArray[Amoebae][Aulo].AssignLinkage((100-propOsmo)/3);

    if (ArtificialPredators) {
        TrophicInteractionsArray[Triradial][Predators].AssignLinkage(17);
        TrophicInteractionsArray[Frond][Predators].AssignLinkage(17);
        TrophicInteractionsArray[Sprigg][Predators].AssignLinkage(17);
        TrophicInteractionsArray[Parvan][Predators].AssignLinkage(17);
        TrophicInteractionsArray[Funisia][Predators].AssignLinkage(17);
        TrophicInteractionsArray[Dickin][Predators].AssignLinkage(17);
        if (FacultativePredators) {
            TrophicInteractionsArray[DOC][Predators].AssignLinkage(17);
        }
    }

    int DetritusBox = DOC;

    //End array loading

    double Complexity = 0; //Mean links per species.
    double Connectivity = 0; //Links per species squared (approximately links per possible link, but this ignores the matrix diagonal elements...)
    double TrophicCoherence = 0; //Degree to which the network forms discrete trophic levels. I'm not sure if this value is quite the same as is described by theory, but it does the job for now.

    //Compute complexity and connectivity
    double nSpecies = 0;
    double nLinks = 0;
    for (int i = 0; i < maxSpecies; i++) {
        if (SpeciesArray[i].ID != -1) {
            nSpecies++;
        }
        for (int j = 0; j < maxSpecies; j++) {
            if (TrophicInteractionsArray[i][j].affinity100 != 0) {
                nLinks++;
            }
        }
    }
    Complexity = nLinks/nSpecies;
    Connectivity = nLinks/(nSpecies*nSpecies);

    int AssignedFijCount = 0;
    while (AssignedFijCount < (maxSpecies * (maxSpecies - 1))) {
        for (int i = 0; i < maxSpecies; i++) {
            SpeciesArray[i].CalculateMi();
            SpeciesArray[i].CalculateFi();
            for (int j = 0; j < maxSpecies; j++) {
                if (i != j) {
                    AssignedFijCount = TrophicInteractionsArray[i][j].CalculateFij(AssignedFijCount);
                }
            }
        }
    }

    for (int i = 0; i < maxSpecies; i++) { //Calculate aij values
        for (int j = 0; j < maxSpecies; j++) {
            if (TrophicInteractionsArray[i][j].affinity100 > 0.00000000001) {
                TrophicInteractionsArray[i][j].InteractionStrengthAonB = SpeciesArray[j].AssimilationEfficiency * SpeciesArray[j].MetabolicEfficiency * TrophicInteractionsArray[i][j].rate / SpeciesArray[i].EqBiomass; //Emily's thesis is in disagreement with the De Ruiter paper about this equation. I think the paper is right.
                TrophicInteractionsArray[j][i].InteractionStrengthAonB = -1 * TrophicInteractionsArray[i][j].rate / SpeciesArray[j].EqBiomass;
            }
        }
    }

    //Compute trophic levels
    bool MissingTrophicLevels = true;
    while (MissingTrophicLevels) {
        MissingTrophicLevels = false;
        for (int i = 0; i < maxSpecies; i++) {
            double LinksStrengthToNode = 0;
            if ((SpeciesArray[i].TrophicLevel < -0.00000000001) & (SpeciesArray[i].ID != -1)) {
                MissingTrophicLevels = true;
                bool MissingUnderlyingTrophicLevels = false;
                for (int j = 0; j < maxSpecies; j++) {
                    if ((TrophicInteractionsArray[j][i].InteractionStrengthAonB > 0.00000000001) & (i != j)) {
                        if (SpeciesArray[j].TrophicLevel < -0.00000000001) {
                            MissingUnderlyingTrophicLevels = true;
                        }
                        else {
                            LinksStrengthToNode += TrophicInteractionsArray[j][i].InteractionStrengthAonB;
                            SpeciesArray[i].TrophicLevel += TrophicInteractionsArray[j][i].InteractionStrengthAonB * SpeciesArray[j].TrophicLevel;
                        }
                    }
                }
                if (LinksStrengthToNode > 0.00000000001) {
                    SpeciesArray[i].TrophicLevel /= LinksStrengthToNode;
                }
                SpeciesArray[i].TrophicLevel += 1;
                if (MissingUnderlyingTrophicLevels) {
                    SpeciesArray[i].TrophicLevel = 0;
                }
            }
        }
    }
    double LValue = 0;
    for (int i = 0; i < maxSpecies; i++) {
        for (int j = 0; j < maxSpecies; j++) {
            if (TrophicInteractionsArray[i][j].InteractionStrengthAonB > 0.00000000001) {
                LValue += TrophicInteractionsArray[i][j].InteractionStrengthAonB;
                TrophicCoherence += (TrophicInteractionsArray[i][j].InteractionStrengthAonB) * ((SpeciesArray[j].TrophicLevel - SpeciesArray[i].TrophicLevel) - 1) * ((SpeciesArray[j].TrophicLevel - SpeciesArray[i].TrophicLevel) - 1);
            }
        }
    }
    TrophicCoherence /= LValue;
    TrophicCoherence = sqrt(TrophicCoherence);
    TrophicCoherence = 1 - TrophicCoherence;

    for (int i = 0; i < maxSpecies; i++) { //Calculate detrital effect values.
        if ((SpeciesArray[i].ID != -1) & (SpeciesArray[i].ID != DetritusBox)) {
            SpeciesArray[i].EffectOnDetritus = SpeciesArray[i].NaturalDeathRate - (TrophicInteractionsArray[DetritusBox][i].rate/SpeciesArray[i].EqBiomass);
            for (int j = 0; j < maxSpecies; j++) {
                SpeciesArray[i].EffectOnDetritus += (1 - SpeciesArray[i].AssimilationEfficiency) * TrophicInteractionsArray[j][i].rate / SpeciesArray[i].EqBiomass;
                SpeciesArray[i].EffectOnDetritus += (1 - SpeciesArray[j].AssimilationEfficiency) * TrophicInteractionsArray[i][j].rate / SpeciesArray[i].EqBiomass; //Is this last term correct? This is as written in the equations, but it seems strange to flip all terms but this one.
            }
        }
        if (SpeciesArray[i].ID == DetritusBox) {
            SpeciesArray[i].EffectOnDetritus = 0;
            for (int j = 0; j < maxSpecies; j++) {
                SpeciesArray[i].EffectOnDetritus -= SpeciesArray[j].AssimilationEfficiency * TrophicInteractionsArray[DetritusBox][j].rate / SpeciesArray[DetritusBox].EqBiomass;
            }
        }
    }

    for (int i = 0; i < maxSpecies; i++) { //This block replaces the original effects of species on detritus with their new effects
        if (SpeciesArray[i].ID != -1) {
            TrophicInteractionsArray[i][DetritusBox].InteractionStrengthAonB = SpeciesArray[i].EffectOnDetritus;
        }
    }

    // Calculate net detritus flow
    float NetDetritusFlux = 0;
    for (int i = 0; i < maxSpecies; i++) {
        if (i != DetritusBox) {
            NetDetritusFlux += TrophicInteractionsArray[i][DetritusBox].InteractionStrengthAonB * SpeciesArray[i].EqBiomass;
        }
    }

    if (ConsoleOutputs) {
        std::cout << "aij Values:\n";
        for (int i = 0; i < maxSpecies; i++) {
            if (SpeciesArray[i].ID != -1) {
                for (int j = 0; j < maxSpecies; j++) {
                        if (SpeciesArray[j].ID != -1) {
                        std::cout << SpeciesArray[i].NodeName << " - " << SpeciesArray[j].NodeName << ": " << TrophicInteractionsArray[i][j].InteractionStrengthAonB << endl;
                    }
                }
            }
        }
    }

    //Calculate degree distribution: I'm not sure if these actually correct: they never allow for two-way connections, which perhaps they should
    int InDegreeDistribution[maxSpecies];
    int OutDegreeDistribution[maxSpecies];
    for (int i = 0; i < maxSpecies; i++) {
        InDegreeDistribution[i] = 0;
        OutDegreeDistribution[i] = 0;
    }
    for (int i = 0; i < maxSpecies; i++) {
        int InDegree = 0;
        int OutDegree = 0;
        if (SpeciesArray[i].ID != -1) {
            for (int j = 0; j < maxSpecies; j++) {
                if (SpeciesArray[i].ID != -1) {
                    if (TrophicInteractionsArray[j][i].InteractionStrengthAonB > 0.00000000001) {
                        InDegree++;
                    }
                    else if (TrophicInteractionsArray[j][i].InteractionStrengthAonB < -0.00000000001) {
                        OutDegree++;
                    }
                }
            }
            InDegreeDistribution[InDegree]++;
            OutDegreeDistribution[OutDegree]++;
        }
    }

    //Here, we calculate the stability metric for the matrix
    Eigen::Matrix<double, maxSpecies, maxSpecies> A; //Create eigen matrix
    A.setZero();
    
    for (int i = 0; i < maxSpecies; i++) { //Add interaction strengths to the matrix
        for (int j = 0; j < maxSpecies; j++) {
            A(i,j) = TrophicInteractionsArray[i][j].InteractionStrengthAonB;
        }
    }

    //Block to home in on the best value of s:
    double PossibleSValues[100000]; //Increasing the size of this array (and associated other values), in order to increase precision, causes an error.
    int maxSn = 99999;
    int minSn = 0;
    for (int n = 0; n < 100000; n++) {
        PossibleSValues[n] = n * -0.00005;
    }
    while ((minSn != maxSn) & (minSn != maxSn - 1)) { //Loop terminates when s is constrained between two adjacent bounds
        int testSn = (maxSn + minSn) * 0.5; //Coerces to an int, which is fine.
        double testS = PossibleSValues[testSn];
        for (int i = 0; i < maxSpecies; i++) {
            if ((SpeciesArray[i].ID != -1) & (SpeciesArray[i].ID != DetritusBox)) {
                A(i,i) = testS * SpeciesArray[i].NaturalDeathRate;
            }
        }
        Eigen::EigenSolver<Eigen::Matrix<double, maxSpecies, maxSpecies> > s(A);
        bool biggerS = false;
        for (int n = 0; n < maxSpecies; n++) {
            if (s.eigenvalues()[n].real() > 0) { //We only care about the real eigenvalue
                biggerS = true;
            }
        }
        if (biggerS) {
            minSn = testSn;
        }
        else {
            maxSn = testSn;
        }
    }

    Eigen::EigenSolver<Eigen::Matrix<double, maxSpecies, maxSpecies> > s(A);
    if (ConsoleOutputs) {
        std::cout << "Interactions matrix" << endl;
        std::cout << A << endl << endl;
        std::cout << "eigenvalues:" << endl;
        std::cout << s.eigenvalues() << endl << endl;;
        std::cout << "s lies between " << PossibleSValues[maxSn] << " and " << PossibleSValues[minSn] << "." << endl;
        std::cout << endl;

        std::cout << "Complexity: " << Complexity << endl;
        std::cout << "Connectivity: " << Connectivity << endl;
        std::cout << "Trophic Coherence: " << TrophicCoherence << endl;

        std::cout << "In Degree Distribution: ";
        for (int i = 0; i < maxSpecies; i++) {
            std::cout << InDegreeDistribution[i] << ",";
        }
        std::cout << endl;
        std::cout << "Out Degree Distribution: ";
        for (int i = 0; i < maxSpecies; i++) {
            std::cout << OutDegreeDistribution[i] << ",";
        }
        std::cout << endl;
    }

    StabilityArray[ReplicateNumber] = PossibleSValues[maxSn];

    if (WeCareAboutLoopWeights) {

        if (ThreeLengthLoopsOnly) {
            for (int i = 0; i < maxSpecies; i++) {
                for (int j = 0; j < maxSpecies; j++) {
                    for (int k = 0; k < maxSpecies; k++) {
                        if ((i != j) & (i != k) & (j != k) & (SpeciesArray[i].ID != -1) & (SpeciesArray[j].ID != -1) & (SpeciesArray[k].ID != -1)) {
                            for (int p = 0; p < maxSpecies; p++) {
                                LoopArray[NextLoopID].Species[p] = -1;
                            }
                            LoopArray[NextLoopID].Species[0] = SpeciesArray[i].ID;
                            LoopArray[NextLoopID].Species[1] = SpeciesArray[j].ID;
                            LoopArray[NextLoopID].Species[2] = SpeciesArray[k].ID;
                            LoopArray[NextLoopID].length = 3;
                            LoopArray[NextLoopID].LoopID = NextLoopID;
                            NextLoopID++;
                        }
                    }
                }
            }
        }
        else {
            for (int i = 0; i < maxSpecies; i++) { //Here, we construct the longest possible loop.
                LoopArray[0].Species[i] = SpeciesArray[i].ID;
                if (SpeciesArray[i].ID != -1) {
                    LoopArray[0].length += 1;
                }
            }
            LoopArray[0].LoopID = 0;
            NextLoopID += 1;

            //Here, we construct all other possible loops
            bool LoopMakerArray[maxSpecies];
            int MaximumLoopBinary = 1;
            int ActualSpecies = 0;
            for (int i = 0; i < maxSpecies; i++) { //Defines the maximum number of loops
                LoopMakerArray[i] = false;
                if (SpeciesArray[i].ID != -1) {
                    MaximumLoopBinary = MaximumLoopBinary * 2;
                    ActualSpecies++;
                }
            }

            int CurrentLoopBinary = 2;
            int ProcessingLoopBinary;
            while (CurrentLoopBinary < MaximumLoopBinary) {
                for (int i = 0; i < ActualSpecies; i++) { //Resets the loop-maker array
                    LoopMakerArray[i] = false;
                }
                ProcessingLoopBinary = CurrentLoopBinary;
                while (ProcessingLoopBinary != 0) { //This block will recur for each trophic loop until all binary components are detected
                    int k = 0; //Power, used as the indexer for binary digit position
                    int k2 = 1; //Two to the power of k: the value of the binary digit at k
                    while (k2 <= ProcessingLoopBinary) {
                        k2 = k2 * 2;
                        k++;
                    }
                    k2 = k2 / 2;
                    k--;
                    LoopMakerArray[k] = true;
                    ProcessingLoopBinary -= k2;
                }
                int InsertionPositon = 0;
                for (int i = 0; i < ActualSpecies; i++) {
                    if (LoopMakerArray[i] == true) {
                        LoopArray[NextLoopID].Species[InsertionPositon] = LoopArray[0].Species[i];
                        InsertionPositon++;
                    }
                }
                LoopArray[NextLoopID].length = InsertionPositon;
                LoopArray[NextLoopID].LoopID = NextLoopID;
                NextLoopID++;
                CurrentLoopBinary++;
            }
        }

        for (int n = 0; n < maxLoops; n++) { //Now, for every loop, we must determine whether or not it is biologically feasible
            if (LoopArray[n].LoopID != -1) { //If a loop exists
                IntList SpeciesInLoop;
                SpeciesInLoop.prefill();
                for (int m = 0; m < maxSpecies; m++) { //SpeciesInLoop is to be a list of all species that need to be in the loop
                    if (LoopArray[n].Species[m] != -1) {
                        SpeciesInLoop.append(LoopArray[n].Species[m]);
                    }
                }
                IntList SpeciesSearched; //Species searched is a list that will start empty, and be filled by the searching algorithm
                SpeciesSearched.prefill();
                SpeciesSearched.Entries[0] = SpeciesInLoop.Entries[0]; //Put the first species into the "found" list
                SpeciesInLoop.Entries[0] = -1;
                IntList Pathway = SearchPathways(SpeciesSearched.Entries[0], SpeciesSearched.Entries[0], -1, SpeciesInLoop, SpeciesSearched, DetritusBox); //Identifies real loops
                LoopArray[n].pathway = Pathway;

                bool ThereIsAPathway = false;
                for (int n = 0; n < maxSpecies; n++) {
                    if (Pathway.Entries[n] != -1) {
                        ThereIsAPathway = true;
                    }
                    if (Pathway.Entries[n] == DetritusBox) {
                        LoopArray[n].type = 0;
                    }
                }
                if (ThereIsAPathway == false) {
                    LoopArray[n].Reset();
                }
            }
        }

        int maxLoopIDtoUse = NextLoopID;
        for (int n = 0; n < maxLoops; n++) { //Every loop must be duplicated once, and inverted, since aij need not equal aji.
            if ((LoopArray[n].LoopID != -1) & (n < maxLoopIDtoUse)) {
                LoopArray[NextLoopID].length = LoopArray[n].length;
                LoopArray[NextLoopID].LoopID = NextLoopID;
                LoopArray[NextLoopID].partnerID = LoopArray[n].LoopID;
                LoopArray[n].partnerID = LoopArray[NextLoopID].LoopID;
                LoopArray[NextLoopID].type = LoopArray[n].type;
                IntList Pathway;
                Pathway.prefill();
                for (int i = 0; i < maxSpecies; i++) {
                    LoopArray[NextLoopID].Species[i] = LoopArray[n].Species[maxSpecies - 1 - i];
                }
                for (int i = 0; i < LoopArray[NextLoopID].length; i++) {
                    Pathway.append(LoopArray[n].pathway.Entries[LoopArray[NextLoopID].length - 1 - i]);
                }
                LoopArray[NextLoopID].pathway = Pathway;
                NextLoopID++;
            }
        }

        SpeciesArray[DetritusBox].NaturalDeathRate = TrophicInteractionsArray[DetritusBox][DetritusBox].InteractionStrengthAonB;
        //Now that we only have the biologically feasible loops, we can calculate their weights. We have invented a detrital natural death rate based on the stability metric.
        //The alternative is to leave dDetritus = 0, which results in undefined loop weights...
        for (int n = 0; n < maxLoops; n++) {
            if (LoopArray[n].LoopID != -1) {
                LoopArray[n].weight += TrophicInteractionsArray[LoopArray[n].pathway.Entries[0]][LoopArray[n].pathway.Entries[1]].InteractionStrengthAonB;
                if (WeightsAreScaledByd) {
                    LoopArray[n].weight /= SpeciesArray[LoopArray[n].pathway.Entries[0]].NaturalDeathRate;
                }
                if (TrophicInteractionsArray[LoopArray[n].pathway.Entries[0]][LoopArray[n].pathway.Entries[1]].InteractionStrengthAonB < -0.00000000001) {
                    LoopArray[n].polarity = -1;
                }
                else {
                    LoopArray[n].polarity = 1;
                }
                for (int k = 2; k < LoopArray[n].length; k++) {
                    LoopArray[n].weight *= TrophicInteractionsArray[LoopArray[n].pathway.Entries[k-1]][LoopArray[n].pathway.Entries[k]].InteractionStrengthAonB;
                    if (WeightsAreScaledByd) {
                        LoopArray[n].weight /= SpeciesArray[LoopArray[n].pathway.Entries[k-1]].NaturalDeathRate;
                    }
                    if (TrophicInteractionsArray[LoopArray[n].pathway.Entries[k-1]][LoopArray[n].pathway.Entries[k]].InteractionStrengthAonB < -0.00000000001) {
                        LoopArray[n].polarity = LoopArray[n].polarity * -1;
                    }
                }
                LoopArray[n].weight *= TrophicInteractionsArray[LoopArray[n].pathway.Entries[LoopArray[n].length - 1]][LoopArray[n].pathway.Entries[0]].InteractionStrengthAonB;
                if (WeightsAreScaledByd) {
                    LoopArray[n].weight /= SpeciesArray[LoopArray[n].pathway.Entries[LoopArray[n].length - 1]].NaturalDeathRate;
                }
                if (TrophicInteractionsArray[LoopArray[n].pathway.Entries[LoopArray[n].length - 1]][LoopArray[n].pathway.Entries[0]].InteractionStrengthAonB < -0.00000000001) {
                    LoopArray[n].polarity = LoopArray[n].polarity * -1;
                }

                if (LoopArray[n].weight < 0) {
                    LoopArray[n].weight *= -1;
                }

                double exponent = pow(LoopArray[n].length, -1); //C++ won't divide 1 by length because length is an int, so we do it this way
                LoopArray[n].weight = pow(LoopArray[n].weight, exponent);
            }
        }

        for (int n = 0; n < maxLoops; n++) {
            if (NetLoopEffects) {
                LoopArray[n].net_weight = LoopArray[n].weight - LoopArray[LoopArray[n].partnerID].weight; // Calculate weight of the loop, minus the reverse of the same loop: the net effect.
            }
            else {
                LoopArray[n].net_weight = LoopArray[n].weight; // Otherwise, just translate acorss the absolute weight.
            }
        }

        for (int i = 0; i < maxLoops; i++) {
            WeightOrderedLoops[i] = -1;
        }
        for (int i = 0; i < maxLoops; i++) {
            if (LoopArray[i].LoopID != -1) {
                int CarriedID = LoopArray[i].LoopID;
                int NewCarriedID;
                for (int j = 0; j < maxLoops; j++) {
                    if (WeightOrderedLoops[j] == -1) {
                        WeightOrderedLoops[j] = CarriedID;
                        break;
                    }
                    else if (LoopArray[WeightOrderedLoops[j]].net_weight < LoopArray[CarriedID].net_weight) {
                        NewCarriedID = WeightOrderedLoops[j];
                        WeightOrderedLoops[j] = CarriedID;
                        CarriedID = NewCarriedID;
                    }
                }
            }
        }

        if (ConsoleOutputs) {
            int n3Loops = 0;
            std::cout << endl;
            std::cout << "Loop weights and polarities:" << endl;
            int NOutputs = 0;
            for (int n = 0; n < maxLoops; n++) {
                if ((WeightOrderedLoops[n] != -1) & (LoopArray[WeightOrderedLoops[n]].net_weight != 0)) {
                    if (LoopArray[WeightOrderedLoops[n]].length == 3) {
                        n3Loops++;
                    }
                    if ((WeREALLYCareAboutLoopWeights) & (NOutputs < 32)) {
                        NOutputs++;
                        std::cout << "Loop " << LoopArray[WeightOrderedLoops[n]].LoopID << ": " << "Polarity: " << LoopArray[WeightOrderedLoops[n]].polarity << ": ";
                        for (int i = 0; i < maxSpecies; i++) {
                            if (LoopArray[WeightOrderedLoops[n]].pathway.Entries[i] != -1) {
                                std::cout << LoopArray[WeightOrderedLoops[n]].pathway.Entries[i] << " - ";
                            }
                        }
                        std::cout << "weight: " << LoopArray[WeightOrderedLoops[n]].net_weight << endl;
                    }
                }
            }
            std::cout << endl;
            std::cout << "Loops of length 3: " << n3Loops << endl;
        }
    }

    HeaviestPositive3LoopWeight = 0;
    HeaviestNegative3LoopWeight = 0;
    HeaviestPositiveLoopWeight = 0;
    HeaviestNegativeLoopWeight = 0;
    if (WeCareAboutLoopWeights) {
        //Find the heaviest loops
        for (int n = 0; n < maxLoops; n++) {
            if (LoopArray[WeightOrderedLoops[n]].length == 3) {
                if ((LoopArray[WeightOrderedLoops[n]].polarity == 1) & (HeaviestPositive3LoopWeight == 0)) {
                    HeaviestPositive3Loop = LoopArray[WeightOrderedLoops[n]];
                    HeaviestPositive3LoopWeight = LoopArray[WeightOrderedLoops[n]].net_weight;
                }
                else if ((LoopArray[WeightOrderedLoops[n]].polarity == -1) & (HeaviestNegative3LoopWeight == 0)) {
                    HeaviestNegative3Loop = LoopArray[WeightOrderedLoops[n]];
                    HeaviestNegative3LoopWeight = LoopArray[WeightOrderedLoops[n]].net_weight;
                }
            }
            if ((LoopArray[WeightOrderedLoops[n]].polarity == 1) & (HeaviestPositiveLoopWeight == 0)) {
                HeaviestPositiveLoop = LoopArray[WeightOrderedLoops[n]];
                HeaviestPositiveLoopWeight = LoopArray[WeightOrderedLoops[n]].net_weight;
            }
            else if ((LoopArray[WeightOrderedLoops[n]].polarity == -1) & (HeaviestNegativeLoopWeight == 0)) {
                HeaviestNegativeLoop = LoopArray[WeightOrderedLoops[n]];
                HeaviestNegativeLoopWeight = LoopArray[WeightOrderedLoops[n]].net_weight;
            }
        }
    }

    if (WriteLog == true) {
        if (ReplicateNumber == 0) {
            OutputFile << "Stability";
            for (int i = 0; i < maxSpecies; i++) {
                if (SpeciesArray[i].ID != -1) {
                    OutputFile << "," << SpeciesArray[i].NodeName << " biomass";
                }
            }
            if (WeCareAboutLoopWeights) {
                OutputFile << ",Heaviest positive 3 loop,Heaviest positive 3 loop weight,Heaviest negative 3 loop,Heaviest negative 3 loop weight";
                OutputFile << ",Heaviest positive loop,Heaviest positive loop weight,Heaviest negative loop,Heaviest negative loop weight";
            }
            if (WeREALLYCareAboutLoopWeights) {
                for (int k = 0; k < maxLoops; k++) {
                    if (LoopArray[k].LoopID != -1) {
                        OutputFile << ",";
                        bool FirstEntry = true;
                        for (int p = 0; p < maxSpecies; p++) {
                            if (LoopArray[k].pathway.Entries[p] != -1) {
                                if (FirstEntry) {
                                    OutputFile << LoopArray[k].pathway.Entries[p];
                                    FirstEntry = false;
                                }
                                else {
                                    OutputFile << "-" << LoopArray[k].pathway.Entries[p];
                                }
                            }
                        }
                    }
                }
            }

            OutputFile << ",Net Flux of Detritus";
            OutputFile << ",Prop Osmotrophy";

            if (OutputdValues) {
                for (int k = 0; k < maxSpecies; k++) {
                    if (SpeciesArray[k].ID != -1) {
                        OutputFile << "," << SpeciesArray[k].NodeName << " d";
                    }
                }
            }
            if (OutputEassValues) {
                for (int k = 0; k < maxSpecies; k++) {
                    if (SpeciesArray[k].ID != -1) {
                        OutputFile << "," << SpeciesArray[k].NodeName << " Eass";
                    }
                }
            }
            if (OutputEprodValues) {
                for (int k = 0; k < maxSpecies; k++) {
                    if (SpeciesArray[k].ID != -1) {
                        OutputFile << "," << SpeciesArray[k].NodeName << " Aass";
                    }
                }
            }
            OutputFile << "\n";
        }
        OutputFile << PossibleSValues[maxSn];
        for (int k = 0; k < maxSpecies; k++) {
            if (SpeciesArray[k].ID != -1) {
                OutputFile << "," << SpeciesArray[k].EqBiomass;
            }
        }
        if (WeCareAboutLoopWeights) {
            bool FirstEntry;
            FirstEntry = true;
            OutputFile << ",";
            for (int j = 0; j < maxSpecies; j++) {
                if (HeaviestPositive3Loop.pathway.Entries[j] != -1) {
                    if (FirstEntry) {
                        OutputFile << HeaviestPositive3Loop.pathway.Entries[j];
                        FirstEntry = false;
                    }
                    else {
                        OutputFile << "-" << HeaviestPositive3Loop.pathway.Entries[j];
                    }
                }
            }
            OutputFile << "," << HeaviestPositive3LoopWeight;
            FirstEntry = true;
            OutputFile << ",";
            for (int j = 0; j < maxSpecies; j++) {
                if (HeaviestNegative3Loop.pathway.Entries[j] != -1) {
                    if (FirstEntry) {
                        OutputFile << HeaviestNegative3Loop.pathway.Entries[j];
                        FirstEntry = false;
                    }
                    else {
                        OutputFile << "-" << HeaviestNegative3Loop.pathway.Entries[j];
                    }
                }
            }
            OutputFile << "," << HeaviestNegative3LoopWeight;
            FirstEntry = true;
            OutputFile << ",";
            for (int j = 0; j < maxSpecies; j++) {
                if (HeaviestPositiveLoop.pathway.Entries[j] != -1) {
                    if (FirstEntry) {
                        OutputFile << HeaviestPositiveLoop.pathway.Entries[j];
                        FirstEntry = false;
                    }
                    else {
                        OutputFile << "-" << HeaviestPositiveLoop.pathway.Entries[j];
                    }
                }
            }
            OutputFile << "," << HeaviestPositiveLoopWeight;
            FirstEntry = true;
            OutputFile << ",";
            for (int j = 0; j < maxSpecies; j++) {
                if (HeaviestNegativeLoop.pathway.Entries[j] != -1) {
                    if (FirstEntry) {
                        OutputFile << HeaviestNegativeLoop.pathway.Entries[j];
                        FirstEntry = false;
                    }
                    else {
                        OutputFile << "-" << HeaviestNegativeLoop.pathway.Entries[j];
                    }
                }
            }
            OutputFile << "," << HeaviestNegativeLoopWeight;
        }
        if (WeREALLYCareAboutLoopWeights) {
            for (int k = 0; k < maxLoops; k++) {
                if (LoopArray[k].LoopID != -1) {
                    OutputFile << "," << (LoopArray[k].weight * LoopArray[k].polarity);
                }
            }
        }
        
        OutputFile << "," << NetDetritusFlux;
        OutputFile << "," << propOsmo;

        if (OutputdValues) {
            for (int k = 0; k < maxSpecies; k++) {
                if (SpeciesArray[k].ID != -1) {
                    OutputFile << "," << (SpeciesArray[k].NaturalDeathRate);
                }
            }
        }
        if (OutputEassValues) {
            for (int k = 0; k < maxSpecies; k++) {
                if (SpeciesArray[k].ID != -1) {
                    OutputFile << "," << (SpeciesArray[k].AssimilationEfficiency);
                }
            }
        }
        if (OutputEprodValues) {
            for (int k = 0; k < maxSpecies; k++) {
                if (SpeciesArray[k].ID != -1) {
                    OutputFile << "," << (SpeciesArray[k].MetabolicEfficiency);
                }
            }
        }
        OutputFile << "\n";
    }

    ReplicateNumber++;
    if (ConsoleOutputs == false) {
        std::cout << ReplicateNumber << endl;
    }
    //Globals need cleaned for the next run, if there is one:
    for (int n = 0; n < maxSpecies; n++) {
        SpeciesArray[n].KillSpecies(n);
    }
    for (int n = 0; n < maxLoops; n++) {
        LoopArray[n].Reset();
    }
    NextID = 0;
    NextLoopID = 0;

    } //End of loop to vary parameters in main

    OutputFile.close();

    return 0;
}