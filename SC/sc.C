// Copyright (c) 2006  John Abbott
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a template file for example programs.  \n"
  "The program itself does nothing whatsoever.    \n";

const string LongDescription =
  "Make a copy of this file (called \"foo.C\", say) and put your code \n"
  "inside the procedure \"program\".                                  \n"
  "To compile your file in the examples directory just do this:       \n"
  "  make foo                                                         \n";
//----------------------------------------------------------------------

void program()
{
  GlobalManager CoCoAFoundations;
  

  cout << endl << "------ building complexes------" << endl;

  ring myQQ = RingQQ();
  SparsePolyRing R=NewPolyRing(myQQ,5);
  const vector<RingElem>& x = indets(R);

  cout <<"------ from PolyList------" << endl;
  
  PolyList L;
  L.push_back(x[1]*x[3]);
  L.push_back(x[0]*x[1]*x[4]);
  SimplicialComplex SC(L);
  cout << "PolyList L: " << L  << endl;
  cout << "SC(L): " << SC <<endl;

  cout << "------ from std::list<face> ------" << endl;  
  
  DynamicBitset f(LPP(x[1]*x[2]));
  DynamicBitset f2(LPP(x[0]*x[2]*x[4]));
  face f3(5);  
  f3.mySet(3,1);
  f3.mySet(2,1);
  f3.mySet(4,1);
  std::list<face> l;
  l.push_back(f);
  l.push_back(f2);
  l.push_back(f3);
  SimplicialComplex SC2(l);
  cout << "std::list<face> l: " << l <<endl;
  cout << "SC2(l): " << SC2 <<endl;

  cout << "------ the empty SC on n vertices------" << endl;
 
  SimplicialComplex emptySC(5);
  cout << "the empty SC on 5 vertices: " << emptySC <<endl;
  cout << "is empty? :: " << emptySC.IamEmptySC() <<endl;
  cout << "is a simplex? : " << emptySC.IamSimplexSC() <<endl;
  cout << endl << "------ the simplex on n vertices ------" << endl;
 
  SimplicialComplex S4=simplexSC(5);
  cout << "the simplex on 5 vertices, S4: " << S4 <<endl;
  cout << "is S4 a simplex? : " << S4.IamSimplexSC() <<endl;
  cout << "is S4 empty? :: " << S4.IamEmptySC() <<endl;
  //-----------------------------------------

  cout << endl <<"------ the boundary of the simplex on n vertices ------" << endl;

  SimplicialComplex bdS4=boundarySimplexSC(5);
  cout << "the boundary of S4: " << bdS4 <<endl;
  //------------------------------------------

  cout << "------ subcomplexes ------" << endl; 

  std::vector<long> V(0);
  V.push_back(0); 
  V.push_back(1);
  V.push_back(2);

  SimplicialComplex subSC(SC,V);
  cout << "SC: " << SC <<endl;
  cout << "vector<long> V: " << V <<endl;
  cout << "subSC(SC,V): " << subSC << endl;
  cout << "subSC.myVSet: " << subSC.myGetVSet() <<endl << endl;
 
  SimplicialComplex subSC2(SC,f2);
  cout << "SC: " << SC <<endl;
  cout << "face f: " << f  << endl;
  cout << "subSC2(SC,f): " << subSC2 << endl;
  cout << "subSC2.myVSet: " << subSC2.myGetVSet() <<endl << endl;
  //-----------------------------------------

 cout << "------myNumIndets, myNumFacets, myDim, myFacetList,myGetVSet------"<< endl << endl;
  cout << "SC: " << SC <<endl;

  cout << "SC.myNumIndets(): " << SC.myNumIndets() << endl << endl;
  cout << "SC.myNumFacets(): " << SC.myNumFacets() << endl << endl;
  cout << "SC.myDim(): " << SC.myDim()  << endl << endl;
  cout << "SC.myFacetList(): " << SC.myFacetList()  << endl << endl;
  cout << "SC.myGetVSet(): " << SC.myGetVSet() <<endl  << endl;
  //------------------------------------------

  cout <<  "------ complements ------" << endl << endl;

  std::list<face> L1=SC.myFacetList();
  face f1=L1.front(); 
  cout << "SC.myVSet: " << SC.myGetVSet() << endl;
  cout << "f1: " << f1 << endl;
  cout << "complement f1: " << SC.myComplF(f1) << endl<< endl;
  //-----------------------------------------

  cout <<  "------ link, star, del, fdel ------"<< endl << endl;

  face v0(LPP(x[0]));
  cout << "vertex 0: " << v0 << endl;
  cout << "boundary of simplex of dim 4, bdS4: " << bdS4 << endl;
  SimplicialComplex LKv=bdS4.myLinkSC(v0);
  cout << "link(v0,bdS4): " << LKv << endl;
  SimplicialComplex STv=bdS4.myStarSC(v0);
  cout << "star(v0,bdS4): " << STv << endl;
  SimplicialComplex Delv=bdS4.myDelSC(v0);
  cout << "del(v0,bdS4): " << Delv << endl;

  cout << "boundary of simplex of dim 4, bdS4: " << bdS4 << endl;
  face f0(LPP(x[0]*x[1]*x[2]*x[3]));
  cout << "f0: " << f0 << endl << endl;
  SimplicialComplex LKf=bdS4.myLinkSC(f0);
  cout << "link(f0,bdS4): " << LKf << endl;
  SimplicialComplex STf=bdS4.myStarSC(f0);
  cout << "star(f0,bdS4): " << STf << endl;
  SimplicialComplex Delf=bdS4.myDelSC(f0);
  cout << "del(f0,bdS4): " << Delf << endl;
  SimplicialComplex FDelf=bdS4.myFaceDelSC(f0);
  cout << "faceDel(f0,bdS4): " << FDelf << endl;

  face g0(LPP(x[0]*x[1]*x[2]));
  cout << "g0: " << g0 << endl << endl;
  LKf=bdS4.myLinkSC(g0);
  cout << "link(g0,bdS4): " << LKf << endl;
  STf=bdS4.myStarSC(g0);
  cout << "star(g0,bdS4): " << STf << endl;
  Delf=bdS4.myDelSC(g0);
  cout << "del(g0,bdS4): " << Delf << endl;
  FDelf=bdS4.myFaceDelSC(g0);
  cout << "faceDel(g0,bdS4): " << FDelf << endl;
  //-----------------------------------------

  cout << "------connectivity------"<< endl << endl;

  cout << "is emptySc connected?: " << emptySC.IamConnected() << endl;
  cout << "is S4 connected?: " << S4.IamConnected() << endl; 
  cout << "is bdS4? connected?: " << bdS4.IamConnected() << endl;
  std::vector<long> V1(0);
  V1.push_back(0); 
  V1.push_back(1);
  V1.push_back(2);
  std::vector<long> V2(2);
  V2[0]=3; 
  V2[1]=4;
  //cout << V2;
  SimplicialComplex newSC1(SC,V2);
  SimplicialComplex newSC2(bdS4,V1);
  cout << "newSC1: " << newSC1 <<endl;
  cout << "newSC2: " << newSC2 <<endl;
  SimplicialComplex newSC=unionSC(newSC1,newSC2);
  cout << "newSC: " << newSC <<endl;
  cout << "is newSC connected?: " << newSC.IamConnected() << endl;
  cout << "connected components of newSC: ";
  
  std::list<SimplicialComplex> CC=newSC.myConnCompsSC();
  cout << len(CC);
  for (std::list<SimplicialComplex>::const_iterator it=CC.begin(); it!=CC.end();++it) cout << *it;
  PolyList NL;
  NL.push_back(x[2]*x[3]);
  NL.push_back(x[1]*x[4]);
  NL.push_back(x[3]*x[4]);
  NL.push_back(x[0]);
  SimplicialComplex NSC(NL);

  cout << endl <<"NSC: " << NSC <<endl;
  cout << "is NSC connected?: " << NSC.IamConnected() << endl;
  cout << "connected components of NSC: ";
  
  CC=NSC.myConnCompsSC();
  cout << len(CC);
  for (std::list<SimplicialComplex>::const_iterator it=CC.begin(); it!=CC.end();++it) cout << *it;
  //-----------------------------------------

  cout << endl << "------membership------" << endl << endl;

  cout << "is f0 in link(v0,bdS4): " << LKv.IamFace(f0) << endl;  
  cout << "is f0 in link(f0,bdS4): " << LKf.IamFace(f0) << endl;
  cout << "is v0 in bdS4: " << bdS4.IamFace(v0) << endl;
  cout << "is v0 a facet in bdS4: " << bdS4.IamFacet(v0) << endl;
  face F0=bdS4.myFacetList().back();
  cout << "F0: " << F0 << endl;
  cout << "is F0 a facet in bdS4: " << bdS4.IamFacet(F0) << endl; 

 //-----------------------------------------

  cout << endl <<  "------purity------" << endl << endl;

  cout << "is emptySC pure?: " << emptySC.IamPure() << endl;
  cout << "is S4 pure?: " << S4.IamPure() << endl;  
  cout << "is bdS4 pure?: " << bdS4.IamPure() << endl;
  cout << "is SC pure?: " << SC.IamPure() << endl;
  cout << "SC: " << SC <<endl;
  //-----------------------------------------

  cout << endl << "------ k-faces ------" << endl << endl;

  cout << "face: " << f0 << endl;
  cout << "1-faces: " << ifacesF(1,f0) << endl;
  //---------------------------------------------------------

  cout << endl << "------ k-decomposability ------" << endl << endl;

  face newf=face(9);
  newf.mySet(0);
  newf.mySet(2);
  newf.mySet(3);
  newf.mySet(4);
  // newf.mySet(7);
cout << "newf: " << newf << endl;
  SimplicialComplex S8=simplexSC(9);
  SimplicialComplex D=S8.myFaceDelSC(newf);
  cout << "D: " << D << endl;
  cout << "Is D 8-decomposable?: " << S8.IamKDecomposable(4) << endl;
  D=simplexSC(5);  
  for(SimplicialComplexConstIter k=l.begin(); k!=l.end(); ++k)
    D=D.myFaceDelSC(*k);
  cout << "D: " << D << endl;
  cout << "Is D 0-decomposable?: " << D.IamKDecomposable(0) << endl;
  cout << "Is D (dim D)-decomposable?: " << D.IamKDecomposable(2) << endl;
  cout << "Is D vertex-decomposable?: " << D.IamVertexDecomposable() << endl;
  SparsePolyRing R7=NewPolyRing(myQQ,7);
  const vector<RingElem>& y = indets(R7);
  PolyList L2;
  RingElem vL2[]={y[0]*y[1]*y[2]*y[6],y[0]*y[1]*y[3]*y[6],y[0]*y[2]*y[5]*y[6],y[0]*y[3]*y[5]*y[6],y[1]*y[2]*y[4]*y[6],y[1]*y[3]*y[4]*y[6],y[2]*y[4]*y[5]*y[6],y[0]*y[1]*y[3]*y[4],y[0]*y[2]*y[3]*y[5],y[1]*y[2]*y[4]*y[5]};
  for (long i=0; i<10; ++i)
    L2.push_back(vL2[i]);
  SimplicialComplex nVD(L2);
  cout << "Is nVD 0-decomposable?: " << nVD.IamKDecomposable(0) << endl << endl;
  //---------------------------------------------------------

  cout << endl << "------ Stanley-Reisner Ideal ------" << endl;
    
  //cout << endl << "emptySC: " << emptySC <<endl;
  // cout << emptySC.myStanleyReisnerIdeal(R) << endl;
  // cout << endl << "S4: " << S4 <<endl;
  //cout << S4.myStanleyReisnerIdeal(R) << endl;
  cout << endl << "bdS4: " << bdS4 <<endl;
  cout << bdS4.myStanleyReisnerIdeal(R) << endl;
  cout << endl << "NSC: " << NSC <<endl;
  //std::list<PPMonoidElem> SRI= NSC.myStanleyReisnerIdeal();
  cout << NSC.myStanleyReisnerIdeal(R) << endl;
  //RingElem(SRI.front()); ?????
  //-----------------------------------------

  cout << endl << "------ boundaries of faces------"<< endl << endl;

  cout << "f2: " << f2 << endl;
  cout << "bound of f2: " << boundaryF(f2)  << endl;
  //-----------------------------------------
 
  cout << endl << "------vertices------"<< endl << endl;

  cout << "f2: " << f2 << endl;
  cout << "vertices of f2: " << vertices(f2) <<  endl;
  //-----------------------------------------
 
  cout << endl << "------unions------" << endl;

  cout << "SC: " << SC <<endl;
  cout << "SC2: " << SC2 << endl;
  SimplicialComplex U=SC|SC2;
  cout << "SC union SC2 :" << U << endl;
  //-----------------------------------------

  cout << "------ intersections ------" << endl;

  cout << "SC: " << SC <<endl;
  cout << "SC2: " << SC2 << endl;
  SimplicialComplex I=intersectSC(SC,SC2);
  cout << "SC intersection SC2 :" << I << endl;
  //-----------------------------------------

 cout << "------ anticycle ideals ------" << endl << endl;

 // SparsePolyRing S=NewPolyRing(myQQ,8);
 // ideal I=antiCycleIdeal(NumIndets(S),S);
 ideal AC=antiCycleIdeal(NumIndets(R),R);
 cout <<"anty-cycle of length :" << NumIndets(R) << endl << AC << endl;
  //---------------------------------------------------------

 cout << endl << "------ polarization of ideals ------" << endl << endl;
 ideal AC2=AC*AC; 
 // cout << AC2 << endl;
 SparsePolyRing S=NewPolyRing(myQQ,4);
 ideal AC4=antiCycleIdeal(NumIndets(R)-1,R);
 // cout << AC4 << endl;
 ideal J=AC2*AC4;
 cout << J << endl;
 cout << endl << "polarized ideal :" << endl << polarization(J) << endl;
 //-------------------------------------------------------


 cout << endl << "------ ADPAC ------" << endl << endl;
 //ideal JJ=ADPAC(6,3);
 SimplicialComplex Delta=ADPAC(5,3);
 cout << Delta << endl;
 // cout << ADPAC(6,3) << endl;
 //--------------------------------------------



  // ideal SRI(NSC.myStanleyReisnerIdeal());
  //cout << "ideal: "<< SRI << endl;



  // cout << "------iterators------" << endl;

  // for (PolyList::const_iterator k=NL.begin(); k!=NL.end(); ++k)
  //   {
  //     for(PolyList::const_iterator j=k; j!=NL.end(); ++j)
  // 	  { 
  // 	    std::cout<< *j <<endl;
  // 	  }
  //   }


  // cout << "------lists------" << endl;

  // std::list<face> LL;
  
  // LL=SC.myFacetList();
  // cout << "L1: " << LL << endl;
  // SimplicialComplexConstIter it1=LL.begin();
  //it1=LL.remove(it1);
  //cout << "L1: " << LL << endl;
  //cout << "*it1: " << *it1 << endl;
  //--it1;
  //++it1;
  //cout << "*it1: " << *it1 << endl;
  // f=f & f2;  
  // cout << "f: " << f << endl;
  // cout << "count f: " << count(f) << endl;
 
  // std::list<face> L1,L2;
  // L1=SC.myFacetList();
  // SimplicialComplexConstIter it1=L1.begin();
  // L2=SC2.myFacetList();
  // SimplicialComplexConstIter it2=L2.begin();

  // cout << IsSubset(*it1,*it2) << endl;
  // cout << len(*it1) << endl;
  // cout << len(*it2) << endl;
  // bool menor;
  // menor=*it1 < *it1;
  // cout << "is <: " << menor << endl;
  //  SC.myUnion(SC2);
  

  //SC.mySort();
  //cout << "SC union SC2 sorted:" << SC <<endl;
}

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }
  return 1;
}

