// Standard library includes
#include <iostream>

// XSecAnalyzer includes
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"


SelectionBase::SelectionBase(std::string fSelectionName_) {
  fSelectionName = fSelectionName_;
  nPassedEvents = 0;

  eventNumber = 0;

  TrueFV = {BOGUS,BOGUS,BOGUS,BOGUS,BOGUS,BOGUS};
  RecoFV = {BOGUS,BOGUS,BOGUS,BOGUS,BOGUS,BOGUS};

  STVTools = STV_Tools();

}

void SelectionBase::Setup(TTree* Tree_, bool Create_) {
  SetupTree(Tree_, Create_);
  DefineCategoryMap();
  DefineConstants();
}

void SelectionBase::ApplySelection(AnalysisEvent* Event) {
 
  Reset();
  std::cout<< "finished - Reset()" << std::endl;
  
  MC_Signal = DefineSignal(Event);
  std::cout<< "finished - DefineSignal(Event)" << std::endl;
  Selected = Selection(Event);
  std::cout<< "finished - Selection(Event)" << std::endl;
  EvtCategory = CategorizeEvent(Event);
  std::cout<< "finished - CategorizeEvent(Event)" << std::endl;
  ComputeRecoObservables(Event);
   std::cout<< "finished - ComputeRecoObservables(Event)" << std::endl;
  if (Event->is_mc_) {   //Event->is_mc_ is set in CategorizeEvent
    ComputeTrueObservables(Event);
  }
  std::cout<< "finished - ComputeTrueObservables(Event)" << std::endl;
  if (Selected) {
    nPassedEvents++;
  }
  eventNumber++;
}

void SelectionBase::Summary() {
  std::cout << fSelectionName << " has " << nPassedEvents << " events which passed" << std::endl;
}

void SelectionBase::SetupTree(TTree* Tree_, bool Create_) {
  std::cout<< "SelectionBase::SetupTree" << std::endl;
  Tree = Tree_;
  Create = Create_;

  std::string BranchName;

  BranchName = "Selected";
  SetBranch(&Selected,BranchName,kBool);
  std::cout<< "Set Branch::"<<BranchName << std::endl;

  BranchName = "MC_Signal";
  SetBranch(&MC_Signal,BranchName,kBool);
  std::cout<< "Set Branch::"<<BranchName << std::endl;

  BranchName = fSelectionName+"Category";
  SetBranch(&EvtCategory,"EventCategory",kInteger);

  std::cout<< "Set Branch::"<<BranchName << std::endl;

  DefineAdditionalInputBranches();
  std::cout<<"Finished :: DefineAdditionalInputBranches()"<< std::endl;
  DefineOutputBranches();
  std::cout<<"Finished :: DefineOutputBranches()"<< std::endl;
}

void SelectionBase::SetBranch(void* Variable, std::string VariableName, VarType VariableType) {
  SaveVariablePointer(Variable,VariableType);

  VariableName = fSelectionName+"_"+VariableName;
  std::string Leaflist = VariableName;

  switch (VariableType) {
  case kBool:
    Leaflist += "/O";
    break;
  case kDouble:
    Leaflist += "/D";
    break;
  case kFloat:
    Leaflist += "/F";
    break;
  case kInteger:
    Leaflist += "/I";
    break;
  case kTVector:
    //set_object_output_branch_address< TVector >(*Tree,VariableName,Variable,Create);
    break;
  case kSTDVector:
    //set_object_output_branch_address< std::vector<double> >(*Tree,VariableName,Variable,Create);
    break;
  case kVectorVectorFloat:
    break;
  case kVectorInteger:
    break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }

  if (Leaflist!="" /*&& VariableType != kVectorInteger&& VariableType != kVectorVectorFloat*/) {
    set_output_branch_address(*Tree,VariableName,Variable,Create,Leaflist);
  } 
  /*
  else if(Leaflist!="" && VariableType == kVectorInteger){
    set_object_output_branch_address< std::vector<int> >(*Tree,
    VariableName, Variable, Create );
  }
    else if(Leaflist!="" && VariableType == kVectorVectorFloat){
   set_object_output_branch_address< std::vector<std::vector<float > > >(*Tree,
    VariableName, Variable, Create);
  }*/
  else {
    set_output_branch_address(*Tree,VariableName,Variable,Create);
  }
}

void SelectionBase::SaveVariablePointer(void* Variable, VarType VariableType) {
  switch (VariableType) {
  case kBool:
    Pointer_Bool.push_back((bool*)Variable);
    break;
  case kDouble:
    Pointer_Double.push_back((double*)Variable);
    break;
  case kFloat:
    Pointer_Float.push_back((float*)Variable);
    break;
  case kInteger:
    Pointer_Integer.push_back((int*)Variable);
    break;
  case kTVector:
    Pointer_TVector.push_back((TVector3*)Variable);
    break;
  case kSTDVector:
    Pointer_STDVector.push_back((std::vector<double>*)Variable);
    break;
  case kVectorVectorFloat:
    Pointer_VectorVectorFloat.push_back((std::vector<std::vector<float>>*)Variable);
    break;
  case kVectorInteger:
     Pointer_VectorInteger.push_back((std::vector<int>*)Variable);
     break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }
}

void SelectionBase::Reset() {
  for (size_t i=0;i<Pointer_Bool.size();i++) {
    *(Pointer_Bool[i]) = false;
  }
  for (size_t i=0;i<Pointer_Double.size();i++) {
    *(Pointer_Double[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Float.size();i++) {
    *(Pointer_Float[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Integer.size();i++) {
    *(Pointer_Integer[i]) = BOGUS_INDEX;
  }
  for (size_t i=0;i<Pointer_TVector.size();i++) {
    /*
    for (size_t j=0;j<(*(Pointer_TVector[i])).GetNrows();j++) {
      (*(Pointer_TVector[i]))[j] = 0.;
    }
    */
    *(Pointer_TVector[i]) = TVector3(0,0,0);
  }
  for (size_t i=0;i<Pointer_STDVector.size();i++) {
    (*(Pointer_STDVector[i])).clear();
  }
  
   for (size_t i=0;i<Pointer_VectorVectorFloat.size();i++) {
   (*(Pointer_VectorVectorFloat[i])).clear();
   }
    for (size_t i=0;i<Pointer_VectorInteger.size();i++) {
    (*(Pointer_VectorInteger[i])).clear();
  }
  
  
  
  
}
