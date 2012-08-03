/*************************************************************************
 *
 * Copyright (C) 2010. Mladen Nikolic (nikolic@matf.bg.ac.rs)
 *
 * This program uses ALGLIB (www.alglib.net)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (www.fsf.org); either version 2 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License is available at
 * http://www.fsf.org/licensing/licenses
 *
 *************************************************************************/

#include "GehanWilcoxon.hpp"

#include <iostream>
#include <fstream>

using std::cin;
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  if(argc!=4 && argc!=5)
  {
    cout << "usage: cmpsat data_file1 data_file2 num_used [bootstrapps]" << endl;
    cout << "Note: if bootstrapps is omitted, jackknife variance estimate is used." << endl;
    exit(-1);
  }

  int m=atoi(argv[3]);
  std::ifstream f(argv[1]);
  std::ifstream g(argv[2]);
  int bootstrapps=0;
  if(argc==5)
    bootstrapps=atoi(argv[4]);

  int n, n1;
  f >> n;
  g >> n1;
  if(n!=n1)
  {
    cout << "Numbers of observations for each solver should be the same!" << endl;
    exit(-1);
  }
  
  double zsum=0;
  double zvar=0;
  double ravg=0;
  double pavg=0;
  int num=0;
  double sup=0;

  while(1)
  {
      Sample s1, s2;
      
      bool allCensored=true;
      double max=0;
      bool eof=false;
      for(int i=0; i<n; i++)
      {
	double x;
	int c;
	f >> x;
	f >> c;
	if(!f.good())
	{
	  if(i==0)
	  {
	    eof=true;
	    break;
	  }
	  else
	  {
	    cout << "Incomplete sample!" << endl;
	    exit(-1);
	  }
	    
	}
	if(i<m)
	{
	  Datum d(x,c);
	  s1.addDatum(d);
	  if(c==1)
	    allCensored=false;
	  if(x>max)
	    max=x;
	}
      }

      if(eof)
	break;

      for(int i=0; i<n; i++)
      {
	double x;
	int c;
	g >> x;
	g >> c;
	if(!g.good())
	{
	  if(i==0)
	  {
	    eof=true;
	    break;
	  }
	  else
	  {
	    cout << "Incomplete sample!" << endl;
	    exit(-1);
	  }
	    
	}
	if(i<m)
	{
	  Datum d(x,c);
	  s2.addDatum(d);
	  if(c==1)
	    allCensored=false;
	  if(x>max)
	    max=x;
	}
      }
      
      if(eof)
	break;

      if(allCensored)
	continue;
      if(max<5)
      	continue;

      num++;

      GehanWilcoxon test(s1, s2, bootstrapps);
      double rv=test.rValue();
      ravg+=rv;
      zsum+=0.5*log((1+rv)/(1-rv));
      double var=test.variance();
      var/=(1-rv*rv)*(1-rv*rv);
      zvar+=var;
      sup+=test.probabilityOfSuperiority();
  }

  f.close();
  g.close();
  ravg/=num;
  cout << "n: " << num << endl;
  cout << "z: " << zsum/sqrt(zvar) << endl;
  cout << "r: " << ravg << endl;
  cout << "p: " << 2-2*normaldistribution(fabs(zsum/sqrt(zvar))) << endl;
  cout << "pi: " << sup/num << endl;
  return 0;
}
