// ROOT needs these dictionary definitions to be able to manipulate TTree
// branches with specific types, such as a std::vector of TVector3 objects
#include "TVector3.h"
#include "TLorentzVector.h"

#pragma link C++ class std::vector< TVector3 >+;
#pragma link C++ class std::vector< TLorentzVector >+;
#pragma link C++ class std::vector< std::vector< float > >+;

#include <vector>
#include <string>
#include <map>

#ifdef __ROOTCLING__
#pragma link C++ class std::vector< int >+;
#pragma link C++ class std::vector< float >+;
#pragma link C++ class std::vector< double >+;
#pragma link C++ class std::string+;
#pragma link C++ class std::map<std::string,std::vector<double> >+;
#endif
