//   Copyright (c)  2013 Oscar

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


#include "CoCoA/oscarSC.H"
#include "CoCoA/DynamicBitset.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/symbol.H" 
#include "CoCoA/ExternalLibs-Frobby.H"
#include "CoCoA/TmpIsTree.H"
#include "CoCoA/RingQQ.H" 

using std::bitset;
using std::ostream;
using std::endl;
using std::list;
using std::vector;


//----- seems to be not necessary any of these:
//#include <algorithm>
//using std::find; 
//#include <bitset>
//#include <functional> 
//#include <iostream>
//using std::flush;
//#include <list>
//#include <utility> 
//using std::pair;
//#include <vector>


static const bool MAX_DEBUG = false;
 
namespace CoCoA
{
  class ideal;

  class SubsetIterator
  {
  public:
    SubsetIterator(int n, int k);
    SubsetIterator& operator++();
    const vector<int>& operator*() const { return myV; }
    bool IamFinished() const { return myV[0] == myN; }
  private: // data members
    const int myN;
    const int myK;
    std::vector<int> myV;
  }; //class SubsetIterator


//--------------------------------------------------------------------------------
//------ SubsetIterator class functions ---------------------------------------
//--------------------------------------------------------------------------------

  SubsetIterator::SubsetIterator(int n, int k):
    myN(n),
    myK(k),
    myV(k)
  {
    for (int i=0; i < k; ++i) { myV[i] = i; }
  }
  
  
  SubsetIterator& SubsetIterator::operator++()
  {
    if (IamFinished()) CoCoA_ERROR("iterator has ended", "SubsetIterator::operator++");
    if (myV[myK-1] < myN-1) { ++myV[myK-1]; return *this; }
    if (myV[0] == myN-myK) { myV[0] = myN; return *this; }
    int i=1;
    while (myV[myK-i] == myN-i) { ++i; }
    ++myV[myK-i];
    while (--i >= 1)
      myV[myK-i] = myV[myK-i-1]+1;
    return *this;
  }
  

  bool IsEnded(const SubsetIterator& S)
  {
    return S.IamFinished();
  }



//--------------------------------------------------------------------------------
//------ SimplicialComplex class functions ---------------------------------------
//--------------------------------------------------------------------------------


  SimplicialComplex::SimplicialComplex( const PolyList& PL)
  { 
    CoCoA_ASSERT(!PL.empty());
    myEmptySC(NumIndets(AsSparsePolyRing(owner(PL.front()))));
    for (PolyList::const_iterator it=PL.begin(); it!=PL.end();++it)
      myInsert(face(LPP(*it))); //need to convert a PolyRing *it into a PP for apply face()
    } //SimplicialComplex(PolyList)
  
  
  SimplicialComplex::SimplicialComplex(const list<face>& l)
  {
    CoCoA_ASSERT(!l.empty());
    myEmptySC(len(l.front()));
    for (list<face>::const_iterator it=l.begin();it!=l.end();++it)
      myInsert(*it);
  } //SimplicialComplex(list<face>)


  SimplicialComplex::SimplicialComplex(const SimplicialComplex& SC, const std::vector<long>& VerticesToKeep)
  {// build the induced simp. comp. of SC on the vertices whose indexes are in VerticesToKeep
    myEmptySC(SC.myNumIndets());   //we create a new SC with the empty set as the only face (facet)
    face mask(SC.myNumIndets());  
    for (long i=0; i < len(VerticesToKeep); ++i)
      mask.mySet(VerticesToKeep[i]); //we will consider the intersection of any facet with this face
    for (list<face>::const_iterator it=SC.myElems.begin(); it != SC.myElems.end(); ++it)
      myInsert(*it & mask);
  } //SimplicialComplex(SimplicialComplex, vector<long>)


  SimplicialComplex::SimplicialComplex(const SimplicialComplex& SC, const face& f)
  {//build the induced simp. comp. of SC on the vertices of the face f
    myEmptySC(SC.myNumIndets());
    for (list<face>::const_iterator it=SC.myElems.begin(); it != SC.myElems.end(); ++it)
      myInsert(*it & f);
  } //SimplicialComplex(SimplicialComplex, face)

  std::ostream& operator<<(std::ostream& out, const SimplicialComplex& the_SimplicialComplex)
  { 
    out << endl << " dim " << the_SimplicialComplex.myDim();
    out << endl << " Vertex Set: "<<the_SimplicialComplex.myGetVSet();
    out << endl << " V Labels  : ";
    for (long i= the_SimplicialComplex.myNumIndets()-1;i>-1; --i)
      if (the_SimplicialComplex.myVSet.Iam1At(i)) out << i;  else out << "*"; // this is a bad way for more than 10 vertices
    SimplicialComplexConstIter it=the_SimplicialComplex.myElems.begin();
    out << endl << " " <<  the_SimplicialComplex.myNumFacets() << " Facets  : " << endl << "             " << *it << endl;
    for (++it; it!=the_SimplicialComplex.myElems.end();++it)
       out  << "             " << *it << endl;
    return out;
  }//operator <<

  
  SimplicialComplex& SimplicialComplex::operator=(const SimplicialComplex& the_SimplicialComplex)
  {
    myVLabels=the_SimplicialComplex.myVLabels;
    myVSet=the_SimplicialComplex.myVSet;
    myElems=the_SimplicialComplex.myElems;
    return *this;
  }//operator =


  void SimplicialComplex::myEmptySC(const int& l)
  {// simp. comp. on l vertices with one facet: the empty set
    //for cretaing SimplicialComplex(long)
    face f(l);
    myVSet=f;
    myElems.clear(); // needed???
    myElems.push_back(f);
  } //myEmptySC


   long SimplicialComplex::myNumIndets() const
  {
    return len(myElems.front());    // Anna: was return myElems.front().myVecSize()*DynamicBitset::ourNumBitsInBlock;
  }//SimplicialComplex::myNumIndets


  void SimplicialComplex::myInsert( const face& f)
  { // add the facet f to the list of facets (if it does not belong already) in the corresponding position for the graded lexicographic order 
    // we are assuming that the list of facets is already sorted.
    bool notContained=false; //for checking if there is any facet in SC containing f
    bool foundPosition=false;
    long cardf=count(f);
    SimplicialComplexIter posInsert=myElems.end();
    for(SimplicialComplexIter it=myElems.begin(); it!=myElems.end(); ++it) // we compare with every facet in the SC, which are in increasing order
      { 
	if(!foundPosition) // we have not yet found the position for f
	  {
	    if (count(*it) < cardf or (count(*it) == cardf and *it < f)) // f greater than *it
	      {if (IsSubset(*it,f))  {it=myElems.erase(it); --it; notContained=true;}} //*it will not be a facet since it is contained in f and f will not be contained in other facet of L (=myElems) 
	    else /*    f =< *it    */ 
	      if (cardf < count(*it) or (count(*it) == cardf and f < *it)) // we must insert f in the position pointed by it
		{
		  posInsert=it;
		  foundPosition=true;
		  --it; //this can be improved: check if is contained or not 
		  if (notContained == true) // we can insert already f in it: f cannot be contained in any facet of L since it does contain one
		    {
		      myElems.insert(posInsert,f);
		      myVSet|=f; 
		      return; 
		    } 
		}
	      else{return;} // f == *it
	  }
	else //we know where should be inserted if it is not contained in the rest of facets of L
	  if (IsSubset (f,*it)) {return;} // if it is contained in other facet it will not be added
      } 
    myElems.insert(posInsert,f);
    myVSet|=f;
  }//myInsert
  

  std::list<face> SimplicialComplex::myFacetList() const
  {  // returns a copy of myElems
    list<face> l;
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end();++it)
      l.push_back(*it);
    return l;
  }//myFacetList


  std::set<face> SimplicialComplex::myifaces(const long& i0, const long& i1 ) const
  {// so far I try to avoid using this function. Namely for IamKDecomposable. Now I use it since it is faster than vertices(VSet)
   // it returns a set with all the faces of with dim between i0-1 and i1-1
    std::list<face> L;
    std::set <face> S;
    for (long i=i0; i<i1+1; ++i)
      {
	for(SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it)
	  {
	    L=ifacesF(i,*it);
	    for (SimplicialComplexConstIter it2=L.begin(); it2!=L.end(); ++it2)
	      S.insert(S.begin(),*it2);
	  }
      }
    return S;
  }//myifaces

 
  ideal  SimplicialComplex::myStanleyReisnerIdeal(const SparsePolyRing& R) //WRT the VSet. I_D=((I_D)^V)^V=(m_(F^c): F is a facet of D)^V
  {
   CoCoA_ASSERT(myNumFacets()!=1);    //if (myNumFacets()==1)  return ideal(0); how?
   std::list<face> L;
   for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
       L.push_back(myComplF(*it));
   return(AlexanderDualFrobby(ideal(FacetList2PolyList(R,L))));
  }// myStanleyReisnerIdeal


  SimplicialComplex SimplicialComplex:: myLinkSC(const face& f) const
  {
    if (f==face(myNumIndets())) return *this;
    SimplicialComplex LK(myNumIndets()); // empty SC on the same vertex set as f
    LK.myElems.pop_back(); // we remove the empty face
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
    { 
      if (f.IamSubset(*it))
      {
      LK.myElems.push_back(*it-f);
      LK.myVSet|=(*it-f); // ^f could be done just once at the end before return
      }
    }
    if (LK.myElems.empty()) LK.myElems.push_back(LK.myVSet); // if we have added no facet, we add the empty set again
    return LK;
  }//myLinkSC


  SimplicialComplex SimplicialComplex::myStarSC(const face& f) const
  {
    SimplicialComplex ST(myNumIndets()); 
    ST.myElems.pop_back(); //we remove the empty face
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
    { 
      if (f.IamSubset(*it))
      {
	ST.myElems.push_back(*it); // the facets of star are the facets of SC containing f
      ST.myVSet|=(*it);
      }
    }
    if (ST.myElems.empty()) ST.myElems.push_back(ST.myVSet); // if we have added no facet, we add the empty set again
    return ST;
  }//myStarSC


   SimplicialComplex SimplicialComplex::myDelSC(const face& f) const
  {
    face cf;
    cf=myComplF(f);
    SimplicialComplex D(*this,cf); // we consider the induced subcomplex on the vertices in VSet but not in f  
    return D;
  }//myDelSC


  SimplicialComplex SimplicialComplex::myFaceDelSC(const face& f) const
  { //I should check the lenght of f and if f is the empty set 
    SimplicialComplex newSC(myNumIndets());
    newSC.myElems.pop_back(); //we remove the empty face
    std::vector<face> L=vertices(f);
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end();++it)
      {
	if (!IsSubset(f,*it)) {newSC.myInsert(*it);}//to improve this
	else // we must remove one vertex in f to the facet and add all such new facets
	  {
	    for (vector<face>::const_iterator itL=L.begin(); itL!=L.end(); ++itL)
	      newSC.myInsert(*it-*itL); //maybe we can improve also this. Use boundary?
	  }
      }
    return newSC;
  }//myFaceDelSC


 face SimplicialComplex::myAppex() const
  {
    face appex=myVSet;
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end();++it)
      appex=appex & *it;
    return appex;
  }


  std::list<face> SimplicialComplex::myCCSupports() const
  { 
    std::list<face> connCompSupports;
    if (myNumFacets()==1) {connCompSupports.push_back(myVSet); return (connCompSupports);}// if it is a simplex or the empty SC returns a list with one element: the vertex set.
    connCompSupports.push_back(myElems.front());
    SimplicialComplexIter itb;
    SimplicialComplexIter itcc;
    for (SimplicialComplexConstIter it=++myElems.begin();
         it!=myElems.end();++it)//better to start from back
    {
      itb=connCompSupports.end();// initialized not pointing to any conn comp
      for (itcc=connCompSupports.begin(); itcc!=connCompSupports.end();++itcc)
	{
	  if (!((*it).IamDisjoint(*itcc))) 
	    {
	      if (itb==connCompSupports.end()) // this is the first conn comp that *it meets
		{
		  if (((*it)|(*itcc))==myVSet) {connCompSupports.push_back(myVSet); return (connCompSupports);} 
                  itb=connCompSupports.insert(itcc,(*it)|(*itcc));//it moves to the right but itb stay in the old place, pointing tothe new entry
		  itcc=--connCompSupports.erase(itcc);//CAREFUL!! using itcc after previous action
		}
	      else // *it meets a previous conn comp pointed by itb
		{
		  connCompSupports.insert(itb,(*itb)|(*itcc));
		  itb=--connCompSupports.erase(itb);
		  itcc=--connCompSupports.erase(itcc);
		}
	     }
	}
      if (itb==connCompSupports.end()) /* the facet does not meet any conn comp so far */ connCompSupports.push_back(*it);
    }
    return connCompSupports;
  } // myCCSupports


  std::list<SimplicialComplex> SimplicialComplex::myConnCompsSC() const
  { // returns a list of SCs built as induced SC on the support of the connected components
    std::list<SimplicialComplex> induSCs;
    if (myNumFacets()==1) {induSCs.push_back(*this); return (induSCs);}// if it is a simplex or the empty SC returns a list with one element: itself
    std::list<face> connCompSupports=myCCSupports();    
    for (SimplicialComplexIter itcc=connCompSupports.begin(); itcc!=connCompSupports.end();++itcc) 
      induSCs.push_back(SimplicialComplex(*this,*(itcc)));
    return induSCs;
  } // myConnCompsSC


  bool  SimplicialComplex::IamFacet(const face& f) const //used in IamSheddinFace
  {
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end();++it)
       if (f==*it) return true;
    return false;
  }//IamFacet


  bool  SimplicialComplex::IamFace(const face& f) const
  {
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end();++it)
      if (f.IamSubset(*it)) return true;
    return false;
  }//IamFace


  bool SimplicialComplex::IamPure() const 
  {
    return count(myElems.front())==count(myElems.back()); // facets are sorted by degrees
  }//IamPure


  bool  SimplicialComplex::IamCone() const 
  {
    face empty=face(myNumIndets());
    if (myAppex()==empty) return false;
    return true;
  }//IamCone


  bool  SimplicialComplex::IamSuspension() const 
  {
    if (IamCone()) return false;
    std::set<face>V=myifaces(0,0);
    std::vector<face>W;
    SimplicialComplex D;
    for(std::set<face>::const_iterator it=V.begin(); it!=V.end(); ++it)
      {
	D=myDelSC(*it);
	if (D.IamCone()) // maybe better if D.myAppex is not the empty face...
	  {
	    W=vertices(D.myAppex());
	    for(long i=0; i<len(W); ++i)
	      if(myLinkSC(*it).myFacetList()==myLinkSC(W[i]).myFacetList()) return true;
	  }
      }
    return false;
  }//IamSuspension

 
 bool SimplicialComplex::IamSheddingVertex(const face& f) const // if every facet of del_SC(f) is a facet of SC
 { // check f is a vertex
    std::list<face> LF;
    SimplicialComplex delSC=myDelSC(f); //
    for (SimplicialComplexConstIter it=delSC.myElems.begin(); it!=delSC.myElems.end();++it)
      if (!IamFacet(*it)) return false;
    return true;
  }//IamSheddingVertex


  bool SimplicialComplex::IamSheddingFace(const face& f) const // iff no facet of fdel_(star_SC(f))(f)) is a facet of fdel_SC(f) [Wo] Remark 3.12
  {
    std::list<face> LF;
    SimplicialComplex fdel_star=myStarSC(f); //
    fdel_star=fdel_star.myFaceDelSC(f);      //is it possible to do directly SimplicialComplex fdel_star=myStarSC(f).myFaceDelSC(f) ?
    SimplicialComplex fdel=myFaceDelSC(f);
    for (SimplicialComplexConstIter it=fdel_star.myElems.begin(); it!=fdel_star.myElems.end();++it)
      if (fdel.IamFacet(*it)) return false;
    //   std::cout << "fdstar: " << fdel_star << endl;
    //   std::cout << "fd: " << fdel <<  endl;
    return true;
  }//IamSheddingFace


  bool SimplicialComplex::IamVertexDecomposable() const
  { 
    if (IamSimplexSC()) return true;
    if (IamEmptySC()) return true;
    std::set<face> V;
    V=myifaces(0,0);
    for (std::set<face>::const_iterator it=V.begin(); it!=V.end(); ++it)
      if (IamSheddingVertex(*it) and myLinkSC(*it).IamVertexDecomposable() and myDelSC(*it).IamVertexDecomposable()) return true;
    return false;
  }//IamVertexDecomposable

  
  bool SimplicialComplex::IamKDecomposable(const int& i) const  //check if it is i- decomposable
  { // check that i is not bigger than the dimension
    if (IamSimplexSC()) return true;
    if (IamEmptySC()) return true;
    if (IamVertexDecomposable()) return true; // k-decomp implies j- decomp for k =< j =< dim SC 
    // there must exists a face f of dim(f) < i such that f it is a shedding face and link(f) and fdel(f) are i-decomposoble
    std::list<face> ifaces;
    std::list<face> checked;
    for (long j=1; j<i+1; ++j)
      {
	for (SimplicialComplexConstIter it=myElems.begin();
	     it!=myElems.end(); ++it)
	  {
	    ifaces=ifacesF(j,*it); // improve this to avoid building the list from the begining: better build one and check, build next and check...
	    for(SimplicialComplexConstIter itf=ifaces.begin(); itf!=ifaces.end(); ++itf)
	      {
		if(!IsInList(*itf,checked))
		  {
		    if (IamSheddingFace(*itf) and myLinkSC(*itf).IamKDecomposable(i) and myFaceDelSC(*itf).IamKDecomposable(i)) return true;
		    checked.push_back(*itf);
		  }	     
	      }
	  }
      }
    return false;
  }//IamKDecomposable



  //-------------private----------------------------
  //------------------------------------------------


  void SimplicialComplex::myLabelVSet() const  //not used so far
  {
    myVLabels.clear();
    for (long i=0; i < myNumIndets(); ++i)
      if (myVSet.Iam1At(i)) myVLabels.push_back(i);
  }//myLabelVSet



  //------------------------------------------------------//
  //-----------functions on faces------------------------//
  //------------------------------------------------------//
  

   bool IsInList(const face& f, const list<face>& L)
   {
     for(SimplicialComplexConstIter it=L.begin(); it!=L.end(); ++it)
       if (f==*it) return true;
     return false;
   }//IsInList


  std::list<face> ifacesF(const int& i, const face& f) //faces of dim i in f
  { 
    long l=len(f);
    face theface;
    std::list<face> ifaces;
    if (!(i<count(f))) return(ifaces);
    std::vector<face> V=vertices(f);
    for (SubsetIterator subsV(len(V),i+1); !IsEnded(subsV); ++subsV)
      {
	theface=face(l);
	for(long j=0; j<i+1; ++j)
	    theface |=  V[(*subsV)[j]];
	ifaces.push_back(theface);
      }
    return ifaces;
  }//ifacesF


  std::list<face> boundaryF(const face& f)
  {
    std::list<face> l;
    face b;
    for (long i=0; i < len(f); ++i) 
      {
	b=f;
	if (f.Iam1At(i)) {l.push_back(b.mySet(i,false));}
      }
    return l; // if f is the empty set then returns an empty list
  }//boudaryF
  

  std::vector<face> vertices(const face& f)
  { 
    //    CoCoA_ASSERT(f!=face(len(f)));
    std::vector<face> vertices;
    long l=len(f);
    for (long i=0; i < len(f); ++i)
      if (f.Iam1At(i))
	vertices.push_back(face(l).mySet(i));
    return vertices;
  }//vertices


  //------------------------------------------------------//
  //-----------functions on SimplicialComplexes-----------//
  //------------------------------------------------------//

  SimplicialComplex simplexSC(const long& n)
  {
    SimplicialComplex SC(n);
    for (long i=0; i < n; ++i)
    SC.myVSet.mySet(i);
    SC.myElems.pop_back();
    SC.myElems.push_back(SC.myVSet);
    return SC;
  }//simplexSC


  SimplicialComplex boundarySimplexSC(const long& n)
  {
    SimplicialComplex S=simplexSC(n);
    face f=S.myGetVSet();
    //list<face> l;
    //l=boundaryF(f);
    return SimplicialComplex(boundaryF(f));    
  }//boundarySimplexSC

  SimplicialComplex unionSC(const SimplicialComplex& SC1, const SimplicialComplex& SC2) 
  {
    //do independently the empty case
    if (SC1.IamEmptySC()) return SC2;
    if (SC2.IamEmptySC()) return SC1;
    std::list<face> l1,l2;
    l1=SC1.myFacetList();
    l2=SC2.myFacetList();
    SimplicialComplex SC(l1);
    for (SimplicialComplexConstIter it2=l2.begin(); it2!=l2.end(); ++it2)
      SC.myInsert(*it2);
    return SC;  
  }// unionSC
  

  SimplicialComplex intersectSC(const SimplicialComplex& SC1, const SimplicialComplex& SC2)
  {
    //do independently the empty case
    if (SC1.IamEmptySC()) return SC1;
    if (SC2.IamEmptySC()) return SC2;
    std::list<face> l1,l2;
    l1=SC1.myFacetList();
    l2=SC2.myFacetList();
    SimplicialComplex SC(SC1.myNumIndets());
    for (SimplicialComplexConstIter it1=l1.begin(); it1!=l1.end(); ++it1) 
      {   
	for (SimplicialComplexConstIter it2=l2.begin(); it2!=l2.end(); ++it2) 
	  SC.myInsert((*it1) & (*it2));
      }
    return SC; 
  }// intersectSC


  
  //---------- other functions ----------------------//
  
  ideal antiCycleIdeal(const long& l, const SparsePolyRing& R)
  { // l must be smaller than or equal to NumIndets(R)
    const vector<RingElem>& x = indets(R);
    std::vector<RingElem> L;
    for (long i=2; i<l-1; ++i )
      L.push_back(x[0]*x[i]);
    for (long j=1; j<l-2; ++j )
      {
	for ( long k=j+2; k<l; ++k )
	  L.push_back(x[j]*x[k]);
      }
    return ideal(L);
  }//antiCycleIdeal


 ideal antiPathIdeal(const long& l, const SparsePolyRing& R)
  { // l must be smaller than or equal to NumIndets(R)
    const vector<RingElem>& x = indets(R);
    std::vector<RingElem> L;
    for (long i=2; i<l; ++i ) //x[0]*x[l-1] is in the ideal
      L.push_back(x[0]*x[i]);
    for (long j=1; j<l-2; ++j )
      {
	for ( long k=j+2; k<l; ++k )
	  L.push_back(x[j]*x[k]);
      }
    return ideal(L);
  }//antiCycleIdeal
  
  
  ideal polarization(const ideal& I)
  {//check AreGensMonomials(I)
    long n=NumIndets(AsSparsePolyRing(AmbientRing (I)));
    PPMonoid M=PPM(AsSparsePolyRing(AmbientRing (I)));
    std::vector<RingElem> theGens=gens(I);
    std::vector<RingElem> newGens;// we could use only theGens
    std::vector<long> theExps(n,0);
    PPMonoidElem lcmPP(M,theExps);
    PPMonoidElem ppM(M);
    for (long i=0; i<len(theGens); ++i)
      lcmPP=lcm(lcmPP,LPP(theGens[i]));
    exponents(theExps,lcmPP);
    std::vector<symbol> newVars;
    for (long j=0; j<n; ++j)
      {
	for(long k=0; k<theExps[j]; ++k)
	  newVars.push_back(symbol("z",j,k));
      }
    SparsePolyRing RZ=NewPolyRing(RingQQ(),newVars);
    // std::cout << endl << "new variables: " << endl << indets(RZ) << endl;
    // there is a lot to improve here
    long indexVar;
    for (long i2=0; i2<len(theGens); ++i2)
      {
	indexVar=0;
	RingElem newGen=one(RZ);
	ppM=LPP(theGens[i2]);
	for (long j2=0; j2<n; ++j2)
	  {
	    for (long k2=0; k2<exponent(ppM,j2); ++k2)
	      newGen=newGen*indet(RZ,indexVar+k2);
	    indexVar=indexVar+theExps[j2];
	  }
	newGens.push_back(newGen);
      }
    return ideal(newGens);
  } //polaization
  
  
  SimplicialComplex  ADPAC(const long& n, const long& k )
  {
    std::list<face> LF;
    SparsePolyRing R(NewPolyRing(RingQQ(),SymbolRange("x",0,n-1)));
    ideal I=antiCycleIdeal(n,R);
    ideal J=I;
    for (long i=1; i<k; ++i) J=J*I;
    SimplicialComplex D(gens(polarization(J)));
    std::list<face> L=D.myFacetList();
    for (SimplicialComplexConstIter it=L.begin(); it!=L.end(); ++it)
      LF.push_back(D.myComplF(*it));
    return SimplicialComplex(LF);
  }//ADPAC

SimplicialComplex  ADPAP(const long& n, const long& k )
  {
    std::list<face> LF;
    SparsePolyRing R(NewPolyRing(RingQQ(),SymbolRange("x",0,n-1)));
    ideal I=antiPathIdeal(n,R);
    ideal J=I;
    for (long i=1; i<k; ++i) J=J*I;
    SimplicialComplex D(gens(polarization(J)));
    std::list<face> L=D.myFacetList();
    for (SimplicialComplexConstIter it=L.begin(); it!=L.end(); ++it)
      LF.push_back(D.myComplF(*it));
    return SimplicialComplex(LF);
  }//ADPAP


   // OPT: don't create a new SimplicialComplex, use a non modifying alg on the old one
   // bool SimplicialComplex::AreConnected_new(const face& g1,const face& g2)const
   // {
   //   if (g1.IamEmpty() || g2.IamEmpty()) return false;
   //   if (g1==g2 || AreDirectlyConnected(g1,g2)) return true;
   //   if (this->myNumFacets()==2) return false;// they are 2, they are not directly connected, they are not connected
   //   unsigned int g1_index=myGetIndex(g1);
   //   unsigned int g2_index=myGetIndex(g2);
   //   vector<unsigned int> P,P1;
   //   for (unsigned int i=0;i!=myNumFacetns  ();++i)
   //     P.push_back(i);
   //   list<unsigned int> xj;
   //   vector<unsigned int> g;
   //   for (unsigned int j=0;j!=myNumIndets();++j)
   //   {
   //     myMakeXj(xj,j);
   //     if (xj.size()>1)
   //     {
   //       myMakeG(g,P,xj);
   //       P1.clear();
   //       for (vector<unsigned int>::const_iterator it1=P.begin();it1!=P.end();++it1)
   //         P1.push_back(g[P[*it1]]);
   //      P1.swap(P);
   //     }//if
   //   }//for
   //   return (P[g1_index]==P[g2_index]);
   //  }// AreConnected_new


}// end namespace cocoa


/*

Some future optimization:

*/
