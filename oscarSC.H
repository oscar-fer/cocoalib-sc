#ifndef CoCoA_oscarSC_H
#define CoCoA_oscarSC_H

//   Copyright (c)  2013 Oscar Fernandez

//   This file "could be" part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/DynamicBitset.H"
#include "CoCoA/TmpGTypes.H"
#include "CoCoA/io.H"
#include "CoCoA/time.H"

#include <cstddef>
#include <vector>
#include <list>
#include <set>
#include <bitset>

//----- seems to be not necessary any of these:
//using std::size_t;
//using std::vector;
//using std::list;
//using std::bitset;
//#include <utility>
//using std::pair;


namespace CoCoA
{

  class PolyRing; // forward declaration -- defined in PolyRing.H
  class RingElem; // forward declaration -- defined in ring.H
  class ideal; //Oscar
  typedef DynamicBitset face;
  typedef std::list<face>::const_iterator SimplicialComplexConstIter;
  typedef std::list<face>::iterator SimplicialComplexIter;

// WARNING TODO: faces are supposed to have all the same len.
// Add a static field and some check.

// Add the PolyRing (at least Monomial) to the SimplicialComplex class.
// This helps the I/O. Or, change the I/O in such a way to get the
// PR from the PL, special case list []
  

  class SimplicialComplex
  {
  public:
    //CONSTRUCTORS
    SimplicialComplex(const long& l) {myEmptySC(l);}; //built a SC on l vertices and the empty set as the only face.
    SimplicialComplex(const PolyList&);
    SimplicialComplex(const std::list<face>&);
    SimplicialComplex(const SimplicialComplex& theSimplicialComplex, const std::vector<long>& VerticesToKeep);
    SimplicialComplex(const SimplicialComplex& theSimplicialComplex, const face& VerticesToKeep);
    SimplicialComplex(){}; // do I need this? Just for declaring without initialization
    ~SimplicialComplex(){};

    //OPERATORS
    friend std::ostream& operator<<(std::ostream&, const SimplicialComplex&);
    SimplicialComplex& operator=(const SimplicialComplex&);
    
    //FUNCTIONS
   
    friend SimplicialComplex simplexSC(const long& n);
    
    void myEmptySC(const int& );
    void myInsert(const face& );
    void myErase(const face& f) { CoCoA_ASSERT(f!=face(len(f))); myElems.remove(f); }// so far I only use this in myInsert, where it is called is to erase a facet contained in a new one so I do not need to erase vertices from myVSet. In that case the CoCoA ASSERT is not necessary
    long myNumIndets() const;
    long myNumFacets() const { return len(myElems); }
    long myDim() const { return count(myElems.back())-1; }
    std::set<face> myifaces(const long& , const long&);
    std::list<face> myFacetList() const;
    face myGetVSet() const { return myVSet; }
    face myComplF(const face& f) const { return myVSet-f; } // complement wrt appearing vertices. f should be a subset of myVSet
 
    SimplicialComplex myLinkSC(const face& ) const;
    SimplicialComplex myStarSC(const face& ) const;
    SimplicialComplex myDelSC(const face& ) const;
    SimplicialComplex myFaceDelSC(const face& ) const;
    
    std::list<face> myCCSupports() const;
    std::list<SimplicialComplex> myConnCompsSC() const;
  
    bool IamEmptySC() const { return ( len(myElems)==1 and myElems.front()==face(myNumIndets())); }
    bool IamSimplexSC() const { return ( len(myElems)==1 and myElems.front()==myGetVSet()); }
    bool IamFace(const face&) const;
    bool IamFacet(const face&) const;
    bool IamPure() const;
    bool IamConnected() const { return len(myCCSupports())==1; }
    bool IamSheddingFace(const face& ) const;
    bool IamVertexDecomposable() const;
    bool IamKDecomposable(const int&) const;

   ideal myStanleyReisnerIdeal(const SparsePolyRing&) ;

  private:
    void myLabelVSet() const;
    
    //FIELDS
    std::list<face> myElems;
    face myVSet;
    mutable std::list<long> myVLabels;
    
  };// end class SimplicialComplex


  inline SimplicialComplex operator|(const SimplicialComplex&, const SimplicialComplex&); //???????
  inline SimplicialComplex operator&(const SimplicialComplex&, const SimplicialComplex&); //???????


  //--------- creating some particular Simplicial Complexes --------//

  SimplicialComplex simplexSC(const long&);
  SimplicialComplex boundarySimplexSC(const long&);

  //-------functions on faces------------//

  bool IsInList(const face& , const std::list<face>& );
  std::list<face> ifacesF(const int& , const face& );
  std::list<face> boundaryF(const  face& ); 
  std::vector<face> vertices(const face& );

  //-------functions on SimplicialComplexes------------//

  SimplicialComplex unionSC(const SimplicialComplex&, const SimplicialComplex&);
  SimplicialComplex intersectSC(const SimplicialComplex&, const SimplicialComplex&);
  SimplicialComplex operator|(const SimplicialComplex& SC1, const SimplicialComplex& SC2) {return (unionSC(SC1,SC2));}
  SimplicialComplex operator&(const SimplicialComplex& SC1, const SimplicialComplex& SC2) {return (intersectSC(SC1,SC2));}

  //--------- other functions ------------------------//

  ideal antiCycleIdeal(const long& , const SparsePolyRing& );
  ideal polarization(const ideal& );
  SimplicialComplex ADPAC(const long& , const long& );


  //  std::vector<RingElem> FaceList2PolyList(const SparsePolyRing&, const std::list<face>&);


}// end namespace cocoa#endif


#endif



/*-----------future changes:

DONE! cambiar typedef face to simplex or face

functions on vertices ends in V
functions on faces ends in S
functions on facets ends in F
functions on simplicial complexes ends in SC

 */
