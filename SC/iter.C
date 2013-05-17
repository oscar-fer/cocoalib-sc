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

  // cout << endl << "------ vertex-decomposable ------" << endl << endl;
  // ring myQQ = RingQQ();
  // SparsePolyRing R7=NewPolyRing(myQQ,7);
  // const vector<RingElem>& y = indets(R7);
  // PolyList L2;
  // RingElem vL2[]={y[0]*y[1]*y[2]*y[6],y[0]*y[1]*y[3]*y[6],y[0]*y[2]*y[5]*y[6],y[0]*y[3]*y[5]*y[6],y[1]*y[2]*y[4]*y[6],y[1]*y[3]*y[4]*y[6],y[2]*y[4]*y[5]*y[6],y[0]*y[1]*y[3]*y[4],y[0]*y[2]*y[3]*y[5],y[1]*y[2]*y[4]*y[5]};
  // for (long i=0; i<10; ++i)
  //   L2.push_back(vL2[i]);
  // SimplicialComplex nVD(L2);
  // cout << "Is nVD 0-decomposable?: " << nVD.IamKDecomposable(0) << endl;

 cout << endl << "------ ADPAC ------" << endl << endl;
 //ideal JJ=ADPAC(6,3);
 long n=6;
 long k=3;
 SimplicialComplex Delta=ADPAC(n,k);
 cout << Delta.IamKDecomposable(1) << endl;
 long d=Delta.myDim();
 cout << "dim :" << d << endl;
 SimplicialComplex fdelD;
 SimplicialComplex linkD;
 std::set<facet> LiF=Delta.myifaces(5,5);
 std::vector<facet> sheddingF;
 face f6; 
 for (set<facet>::const_iterator it=LiF.begin(); it!=LiF.end(); ++it)
   if(Delta.IamSheddingFace(*it)) {f6=*it; break;}
 cout << f6 << endl;
 
 long ds=n-1;
 fdelD=Delta.myFaceDelSC(f6);
 linkD= Delta.myLinkSC(f6);
 cout << "facets fdel :" << fdelD.myNumFacets() << endl;
 cout << "facets link :" << linkD.myNumFacets() << endl; 
 
 f6=face(len(f6));
 LiF=fdelD.myifaces(ds,ds);
 for(set<facet>::const_iterator it=LiF.begin(); it!=LiF.end(); ++it)
   if(fdelD.IamSheddingFace(*it))  {f6=*it; break;}

 while(linkD.myNumFacets()==1 and fdelD.myNumFacets()>1)
   {
     cout << f6 << endl;
     fdelD=fdelD.myFaceDelSC(f6);
     linkD=fdelD.myLinkSC(f6);
     cout << "facets fdel :" << fdelD.myNumFacets() << endl;
     cout << "dim fdel :" << fdelD.myDim() << endl;
     cout << fdelD.IamKDecomposable(1) << endl;
     cout << "facets link :" << linkD.myNumFacets() << endl; 
     f6=face(len(f6));
     LiF=fdelD.myifaces(ds,ds);
     for(set<facet>::const_iterator it=LiF.begin(); it!=LiF.end(); ++it)
       if(fdelD.IamSheddingFace(*it))  {f6=*it; break;}
     while (f6==face(len(f6)))
       {
	 ds=ds-1;
         cout << fdelD.IamKDecomposable(1) << endl;
	 cout<<ds<<endl;
	 LiF=fdelD.myifaces(ds,ds);
	 for(set<facet>::const_iterator it=LiF.begin(); it!=LiF.end(); ++it)
	   if(fdelD.IamSheddingFace(*it))  {f6=*it; break;}
       }
   }
 cout << " fdel :" << fdelD << endl;
 cout << " link :" << linkD << endl; 
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-empty.C,v 1.6 2012/11/30 14:04:55 abbott Exp $
// $Log: ex-empty.C,v $
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
