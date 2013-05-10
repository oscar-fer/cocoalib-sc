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
#include "CoCoA/library.H"

#include <algorithm>
using std::find;
//#include <bitset>
using std::bitset;
#include <functional>
//??
#include <iostream>
using std::ostream;
using std::endl;
using std::flush;
//#include <list>
using std::list;
#include <utility>
using std::pair;
//#include <vector>
using std::vector;

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



    SimplicialComplex::SimplicialComplex(/*const PolyRing& P,*/ const PolyList& PL)
  { 
    CoCoA_ASSERT(!PL.empty());
  
    myEmptySC(NumIndets(AsSparsePolyRing(owner(PL.front()))));
    for (PolyList::const_iterator it=PL.begin(); it!=PL.end();++it)
    {
      myInsert(facet(LPP(*it)));
    }
  } //SimplicialComplex(PolyList)


    SimplicialComplex::SimplicialComplex(const list<facet>& l)
  {
    CoCoA_ASSERT(!l.empty());

    myEmptySC(len(l.front()));
    for (list<facet>::const_iterator it=l.begin();it!=l.end();++it)
    {
      myInsert(*it);
    }
  } //SimplicialComplex(list<facet>)


    SimplicialComplex::SimplicialComplex(const SimplicialComplex& SC, const std::vector<long>& VerticesToKeep)
  {
    myEmptySC(SC.myNumIndets());
    facet mask(SC.myNumIndets());
    for (long i=0; i < len(VerticesToKeep); ++i)
      mask.mySet(VerticesToKeep[i]);
    for (list<facet>::const_iterator it=SC.myElems.begin(); it != SC.myElems.end(); ++it)
      myInsert(*it & mask);
  } //SimplicialComplex(SimplicialComplex, vector<long>)

    SimplicialComplex::SimplicialComplex(const SimplicialComplex& SC, const facet& f)
  {
    myEmptySC(SC.myNumIndets());
    for (list<facet>::const_iterator it=SC.myElems.begin(); it != SC.myElems.end(); ++it)
      myInsert(*it & f);
  } //SimplicialComplex(SimplicialComplex, facet)

   std::ostream& operator<<(std::ostream& out, const SimplicialComplex& the_SimplicialComplex)
   { 
     out << "dim " << the_SimplicialComplex.myDim() << endl;
    out << endl << "Vertex Set: "<<the_SimplicialComplex.myGetVSet();
    out << endl << "V Labels  : ";
    for (long i= the_SimplicialComplex.myNumIndets()-1;i>-1; --i)
      {if (the_SimplicialComplex.myVSet.Iam1At(i)) out << i;  else out << "*";}// this would be a bad way for more than 10 vertices
    SimplicialComplexConstIter it=the_SimplicialComplex.myElems.begin();
    out << endl << the_SimplicialComplex.mySize() << " Facets  : " << *it << endl;

    for (++it; it!=the_SimplicialComplex.myElems.end();++it)
    {
       out  << "            " << *it << endl;
    }
    return out;
  }//operator <<

  // f1=f2
  SimplicialComplex& SimplicialComplex::operator=(const SimplicialComplex& the_SimplicialComplex)
  {
    myVLabels=the_SimplicialComplex.myVLabels;
    myVSet=the_SimplicialComplex.myVSet;
    myElems=the_SimplicialComplex.myElems;
    return *this;
  }//operator =


  void SimplicialComplex::myEmptySC(const int& l)
  {
    facet f(l);
    myVSet=f;
    myElems.clear();
    myElems.push_back(f);
  } //SimplicialComplex(long)


   long SimplicialComplex::myNumIndets() const
  {
    //  if (IamEmptySC()) return 0;
    return len(myElems.front());
    // Anna: was return myElems.front().myVecSize()*DynamicBitset::ourNumBitsInBlock;
  }//SimplicialComplex::myNumIndets

  //  long SimplicialComplex::myGetIndex(const facet& f) const
  // {
  //   long i=0;
  //   for(SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it,++i)
  //   {
  //     if (*it==f)
  //        return i;
  //   }
  //   return 0;
  // }//SimplicialComplex::myGetIndex


  void SimplicialComplex::myInsert( const facet& f)
  {
    bool notContained=false;
    bool foundPosition=false;
    SimplicialComplexIter posInsert;
    posInsert=myElems.end();
    for(SimplicialComplexIter it=myElems.begin(); it!=myElems.end(); ++it) // we compare with every facet in SC1
    { 
      if(!foundPosition) // we have not yet found the position for *it2
     {
       if (count(*it) < count(f) or (count(*it) == count(f) and *it < f)) 
	{
          if (IsSubset(*it,f))
            {it=myElems.erase(it); --it; notContained=true;} //*it will not be a facet in the union since it is contained in other element and f will not be contained in other facet of L (=myElems)
	  //it can be moved to myElems.begin()-1  !!!! but seems not to be problems
        }
       else /*// f =< *it //*/ if (count(f) < count(*it) or (count(*it) == count(f) and f < *it)) 
        {
          posInsert=it;
          foundPosition=true;
          --it; //this can be improved: check if is contained or not 
          if (notContained == true ) 
          {
            myElems.insert(posInsert,f);
            myVSet|=f; 
            return; // f can not be contained in any facet of L since it does  contain one facet of L
          } 
        }
	else{return;} // f == *it
      }
      else //we know where should be inserted if it is not contained in the rest of facets of L
       {
	 if (IsSubset (f,*it)) {return;} // if it is contained in other facet it will not be added
       }
    } 
    myElems.insert(posInsert,f);
    myVSet|=f;
  }//myInsert


  list<facet> SimplicialComplex::mySimplicialComplex2FacetList() const
  {
    list<facet> l;
    for (SimplicialComplexConstIter it=myElems.begin();
         it!=myElems.end();++it)
    {
	 l.push_back(*it);
     }
     return l;
  }//SimplicialComplex::mySimplicialComplex2FacetList


std::set<facet> SimplicialComplex::myifaces(const long& i0, const long& i1 )
{// try with set
  std::list<facet> L;
  std::set<facet> S;
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
}


  facet SimplicialComplex::myGetVSet() const
  {
    return myVSet;
  }

 
  ideal  SimplicialComplex::myStanleyReisnerIdeal(const SparsePolyRing& R) //WRT the VSet. I_D=((I_D)^V)^V=(m_(F^c): F is a facet of D)^V
  {
   CoCoA_ASSERT(!IamEmptySC()); // a possible answer could be the maximal graded ideal in the polynomial ring with myNumIndets variables 
    // if (mySize()==1)  return L2; // if it is a simplex we return ideal 0
   CoCoA_ASSERT(mySize()!=1);    //if (mySize()==1)  return ideal(0);
   std::list<facet> L;
   for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
       L.push_back(myComplF(*it));
   //return(ideal(FacetList2PolyList(R,L)));
    return(AlexanderDualFrobby(ideal(FacetList2PolyList(R,L))));
  }


  SimplicialComplex SimplicialComplex:: myLinkSC(const facet& f) const
  {
    SimplicialComplex LK(myNumIndets());
    LK.myElems.pop_back(); //we remove the empty face
   
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
    { 
      if (f.IamSubset(*it))
      {
      LK.myElems.push_back(*it^f);
      LK.myVSet|=(*it^f);
      }
    }
    if (LK.myElems.empty()) LK.myElems.push_back(LK.myVSet); 
   
    return LK;
  }//myLinkSC


 SimplicialComplex SimplicialComplex::myStarSC(const facet& f) const
  {
    SimplicialComplex ST(myNumIndets()); 
    ST.myElems.pop_back(); //we remove the empty face
   
    for (SimplicialComplexConstIter it=myElems.begin(); it!=myElems.end(); ++it) 
    { 
      if (f.IamSubset(*it))
      {
      ST.myElems.push_back(*it);
      ST.myVSet|=(*it);
      }
    }
    if (ST.myElems.empty()) ST.myElems.push_back(ST.myVSet); 
   
    return ST;
  }//myStarSC


   SimplicialComplex SimplicialComplex::myDelSC(const facet& f) const
  {
    facet cf;
    cf=myComplF(f);
    SimplicialComplex D(*this,cf);
    return D;
  }//myDelSC


  SimplicialComplex SimplicialComplex::myFaceDelSC(const facet& f) const
  { //I should check the lenght of f
    SimplicialComplex newSC(myNumIndets());
    newSC.myElems.pop_back(); //we remove the empty face
    std::vector<facet> L=vertices(f);
    for ( list<facet>::const_iterator it=myElems.begin();
         it!=myElems.end();++it)
      {
	if (!IsSubset(f,*it)) {newSC.myInsert(*it);}//to improve this
	else
	  {
	    for (vector<facet>::const_iterator itL=L.begin(); itL!=L.end(); ++itL)
	      newSC.myInsert(*it-*itL); //maybe we can improve also this. Use boundary?
	  }
      }
    return newSC;
  }//myFaceDelSC



  std::list<SimplicialComplex> SimplicialComplex::connCompsSC() const
  { 
    std::list<SimplicialComplex> induSCs;
    if (mySize()==1) {induSCs.push_back(*this); return (induSCs);}

    std::list<facet> connComps;
    connComps.push_back(myElems.front());
    SimplicialComplexIter itb;
    SimplicialComplexIter itcc;
    for (SimplicialComplexConstIter it=++myElems.begin();
         it!=myElems.end();++it)//better to start from back
    {
      itb=connComps.end();
      for (itcc=connComps.begin(); itcc!=connComps.end();++itcc)
	{
	  if (!((*it).IamDisjoint(*itcc))) 
	    {
	      if (itb==connComps.end()) 
		{
		  if (((*it)|(*itcc))==myVSet) {/*induSCs.clear();*/induSCs.push_back(*this); return (induSCs);}
                  itb=connComps.insert(itcc,(*it)|(*itcc));//it moves to the right but itb stay in the old place, pointing tothe new entry
		  itcc=--connComps.erase(itcc);
		}
	      else 
		{
		  connComps.insert(itb,(*itb)|(*itcc));
		  itb=--connComps.erase(itb);
		  itcc=--connComps.erase(itcc);
		}
	     }
	}

      if (itb==connComps.end()) connComps.push_back(*it);
    }
    for (itcc=connComps.begin(); itcc!=connComps.end();++itcc) 
    	   induSCs.push_back(SimplicialComplex(*this,*(itcc)));
     return induSCs;
  }

  bool  SimplicialComplex::IamFacet(const facet& f) const //used in IamSheddinFace
  {
    for (SimplicialComplexConstIter it=myElems.begin();
         it!=myElems.end();++it)
    {
      if (f==*it) return true;
    }
     return false;
  }



  bool  SimplicialComplex::IamFace(const facet& f) const
  {
    for (SimplicialComplexConstIter it=myElems.begin();
         it!=myElems.end();++it)
    {
      if (f.IamSubset(*it)) return true;
    }
     return false;
  }

  bool SimplicialComplex::IamPure() const 
  {
    return count(myElems.front())==count(myElems.back());
  }


  bool SimplicialComplex::IamConnected() const
  { 
    if (mySize()==1) return true;

    std::list<facet> connComp;
    connComp.push_back(myElems.front());
    SimplicialComplexIter itb;
    for (SimplicialComplexConstIter it=++myElems.begin();
         it!=myElems.end();++it)//better to start from back
    {
      itb=connComp.end(); 
      for (SimplicialComplexIter itcc=connComp.begin(); itcc!=connComp.end();++itcc)
	{
	  if (!((*it).IamDisjoint(*itcc))) 
	    {
	      if (itb==connComp.end()) 
		{
		  if (((*it)|(*itcc))==myVSet) return true;
		  itb=connComp.insert(itcc,(*it)|(*itcc)); //CAREFULL!!
		  itcc=--connComp.erase(itcc);
		}
	      else
		{
		  connComp.insert(itb,(*itb)|(*itcc));
		  itb=--connComp.erase(itb);
		  itcc=--connComp.erase(itcc);
		}
	    }
	}   
      if (itb==connComp.end()) connComp.push_back(*it);
    }
     return false;
  }


  bool SimplicialComplex::IamSheddingFace(const facet& f) const // iff no facet of fdel_(star_SC(f))(f)) is a facet of fdel_SC(f) [Wo] Remark 3.12
  {
    std::list<facet> LF;
    SimplicialComplex fdel_star=myStarSC(f);
    fdel_star=fdel_star.myFaceDelSC(f); 
    SimplicialComplex fdel=myFaceDelSC(f);
    for (SimplicialComplexConstIter it=fdel_star.myElems.begin();
         it!=fdel_star.myElems.end();++it)
      if (fdel.IamFacet(*it)) return false;
    //   std::cout << "fdstar: " << fdel_star << endl;
    //   std::cout << "fd: " << fdel <<  endl;
    return true;
  }


  
   bool SimplicialComplex::IamKDecomposable(const int& i) const
   { // check that i is not bigger than the dimension
     std::list<facet> ifaces;
     std::list<facet> checked;
     for (SimplicialComplexConstIter it=myElems.begin();
	  it!=myElems.end(); ++it)
       {
	 ifaces=ifacesF(i,*it);
	 //	 std::cout<< "termina" <<endl;
	 for(SimplicialComplexConstIter itf=ifaces.begin(); itf!=ifaces.end(); ++itf)
	   {
	     //  std::cout << "iface: " << *itf << endl;
	     if(!IsInList(*itf,checked))
	       {
		 if (IamSheddingFace(*itf)) return true;
		 else { checked.push_back(*itf);}
	       }	     
	   }
       }
     return false;
   }



  //-------------private----------------------------
  //----------------------------------


   void SimplicialComplex::myGetVLabels() const
  {
    myVLabels.clear();
    for (long i=0; i < myNumIndets(); ++i)
    if (myVSet.Iam1At(i)) myVLabels.push_back(i);
  }



  //------------------------------------------------------//
  //-----------functions on facets------------------------//
  //------------------------------------------------------//
  
  // facet facetFromIndexes(const long  l, const std::vector<int> V)
  // {
  //   for (long i=0; i<len(V); ++i)
  //     {
	
  //     }
  // }

   bool IsInList(const facet& f, const list<facet>& L)
   {
     for(SimplicialComplexConstIter it=L.begin(); it!=L.end(); ++it)
       if (f==*it) return true;
     return false;
   }

  std::list<facet> ifacesF(const int& i, const facet& f) //faces of dim i in f
  { 
    long l=len(f);
    facet face;
    std::list<facet> ifaces;
    if (!(i<count(f))) return(ifaces);
    std::vector<facet> V=vertices(f);
    for (SubsetIterator subsV(len(V),i+1); !IsEnded(subsV); ++subsV)
      {
	face=facet(l);
	for(long j=0; j<i+1; ++j)
	    face |=  V[(*subsV)[j]];
	ifaces.push_back(face);
      }
    return ifaces;
  }



  std::list<facet> boundaryF(facet& f)
  {
    std::list<facet> l;
    facet b;
    for (long i=0; i < len(f); ++i) 
    {
      if (f.Iam1At(i)) {b=f.mySet(i,false); l.push_back(b); f.mySet(i);}// improve this
    }
    return l;
  }//boudaryF


  std::vector<facet> vertices(const facet& f)
  { 
    //    CoCoA_ASSERT(f!=facet(len(f)));
    std::vector<facet> vertices;
    facet emptyF(len(f));
    for (long i=0; i < len(f); ++i)
    if (f.Iam1At(i))
      {
	emptyF.mySet(i);
	vertices.push_back(emptyF);
	emptyF.mySet(i,false);
      }
    return vertices;
  }

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
    facet f=S.myGetVSet();
    //list<facet> l;
    //l=boundaryF(f);
    return SimplicialComplex(boundaryF(f));    
  }//boundarySimplexSC

  SimplicialComplex unionSC(const SimplicialComplex& SC1, const SimplicialComplex& SC2) {
    //do independently the empty case
    if (SC1.IamEmptySC()) return SC2;
    if (SC2.IamEmptySC()) return SC1;

    std::list<facet> l1,l2;
    l1=SC1.mySimplicialComplex2FacetList();
    l2=SC2.mySimplicialComplex2FacetList();
    SimplicialComplex SC(l1);
    for (SimplicialComplexConstIter it2=l2.begin(); it2!=l2.end(); ++it2)
    {
      SC.myInsert(*it2);
    }
  return SC;  
  }// unionSC
  
SimplicialComplex intersectSC(const SimplicialComplex& SC1, const SimplicialComplex& SC2)
{
//do independently the empty case
    if (SC1.IamEmptySC()) return SC1;
    if (SC2.IamEmptySC()) return SC2;

   std::list<facet> l1,l2;
   l1=SC1.mySimplicialComplex2FacetList();
   l2=SC2.mySimplicialComplex2FacetList();
   SimplicialComplex SC(SC1.myNumIndets());
   for (SimplicialComplexConstIter it1=l1.begin(); it1!=l1.end(); ++it1) 
   {   
     for (SimplicialComplexConstIter it2=l2.begin(); it2!=l2.end(); ++it2) 
     {
       SC.myInsert((*it1) & (*it2));
     }
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
}

 
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
}

  
SimplicialComplex  ADPAC(const long& n, const long& k )
  {
    std::list<facet> LF;;
    SparsePolyRing R(NewPolyRing(RingQQ(),SymbolRange("x",0,n-1)));
    ideal I=antiCycleIdeal(n,R);
    ideal J=I;
    for (long i=1; i<k; ++i) J=J*I;
    SimplicialComplex D(gens(polarization(J)));
    std::list<facet> L=D.mySimplicialComplex2FacetList();
    for (SimplicialComplexConstIter it=L.begin(); it!=L.end(); ++it)
      LF.push_back(D.myComplF(*it));
    return SimplicialComplex(LF);
  }






   // OPT: don't create a new SimplicialComplex, use a non modifying alg on the old one
   // bool SimplicialComplex::AreConnected_new(const facet& g1,const facet& g2)const
   // {
   //   if (g1.IamEmpty() || g2.IamEmpty()) return false;
   //   if (g1==g2 || AreDirectlyConnected(g1,g2)) return true;
   //   if (this->mySize()==2) return false;// they are 2, they are not directly connected, they are not connected
   //   unsigned int g1_index=myGetIndex(g1);
   //   unsigned int g2_index=myGetIndex(g2);
   //   vector<unsigned int> P,P1;
   //   for (unsigned int i=0;i!=mySize();++i)
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



     //----------------------------------------------------------------------
  // functions of class facet. Oscar: not mine

  // // Not efficient, efficiency not needed at the moment
  // RingElem Facet2RingElem(const SparsePolyRing& theP,const DynamicBitset& b)
  // {
  //   return monomial(theP, 1, NewPP(PPM(theP), b));
  // }//Facet2RingElem


  // // Not efficient, efficiency not needed at the moment
  // std::vector<RingElem> FacetList2PolyList(const SparsePolyRing& theP,const std::list<DynamicBitset>& theFL)
  // {
  //   vector<RingElem> PL;
  //   for (list<DynamicBitset>::const_iterator it=theFL.begin(); it!=theFL.end(); ++it)
  //     PL.push_back(Facet2RingElem(theP,*it));
  //   return PL;
  // }//FacetList2PolyList

 // end functions of class facet
     //-------------------------------------------------------------------------
 
}// end namespace cocoa




/*

Some future optimization:

*/
