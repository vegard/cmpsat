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

#ifndef _GEHAN_WILCOXON_
#define _GEHAN_WILCOXON_

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "normaldistr.h"

using std::cin;
using std::cout;
using std::endl;

class Datum
{
public:
  Datum(double value, bool known)
    : _value(value), _known(known)
  {
  }
  
  double getValue() const
  {
    return _value;
  }

  bool valueKnown() const
  {
    return _known;
  }

private:
  double _value;
  bool _known;
};

class Sample
{
public:
  Sample()
  {
  }
  
  void addDatum(Datum &datum)
  {
    _data.push_back(datum);
  }

  size_t size()
  {
    return _data.size();
  }

  int diedAtTime(double time)
  {
    size_t i;
    for(i=0; i<_data.size() && 
	  _data[i].valueKnown() && 
	  _data[i].getValue()<time; i++)
      ;
      
    if(i==_data.size() || !_data[i].valueKnown() || _data[i].getValue()!=time)
      return 0;

    int died=1;
    while((i+died<_data.size()) && 
	  _data[i+died].valueKnown() && 
	  (_data[i+died].getValue()==time))
      died++;
    return died;    
  }

  int atRisk(double time)
  {
    int numAtRisk=0;
    size_t i;
    for(i=0; i<_data.size(); i++)
      if(!_data[i].valueKnown() || _data[i].getValue()>=time)
	break;
    return size()-i;
    
  }

  double hazard(double time)
  {
    return diedAtTime(time)/(double)atRisk(time);
  }

  Datum &operator[](size_t i)
  {
    return _data[i];
  }

  std::vector<Datum>::iterator begin()
  {
    return _data.begin();
  }

  std::vector<Datum>::iterator end()
  {
    return _data.end();
  }
  
  void insert(std::vector<Datum>::iterator dest,
	      std::vector<Datum>::iterator beg,
	      std::vector<Datum>::iterator en)
  {
    return _data.insert(dest,beg,en);
  }
  

private:
  std::vector<Datum> _data;
};


class GehanWilcoxon
{
public:
  GehanWilcoxon(Sample &firstSample, Sample &secondSample, int bootstrapps)
    : _firstSample(firstSample), _secondSample(secondSample)
  {
    formPooledSample();
    _testStatistic=testStatistic();
    if(bootstrapps>0)
      _variance=bootstrapVariance(bootstrapps);
    else
      _variance=jackKnifeVariance();
    _zScore=_testStatistic/sqrt(_variance);
    _pValue=2*(1-normaldistribution(fabs(_zScore)));
    _rValue=rValue(_ranks,_indicators);
  }

  double gehanStatistic()
  {
    return _testStatistic;
  }

  double variance()
  {
    return _variance;
  }

  double zScore()
  {
    return _zScore;
  }

  double pValue()
  {
    return _pValue;
  }

  double rValue()
  {
    return _rValue;
  }

  double rCorr()
  {
    return rCorr(_ranks,_indicators);
  }

  double probabilityOfSuperiority()
  {
    return -0.5*_testStatistic+0.5;
  }


private:
  
  double rValue(std::vector<double> r1, std::vector<int> r2)
  {
    double firstSampleMean=0;
    double secondSampleMean=0;
    double r=0;
    double firstSampleVariance=0;
    double secondSampleVariance=0;
    for(int i=0; i<r1.size(); i++)
    {
      firstSampleMean+=r1[i];
      secondSampleMean+=r2[i];
    }
    firstSampleMean/=r1.size();
    secondSampleMean/=r2.size();
    for(int i=0; i<r1.size(); i++)
    {
      r+=(r1[i]-firstSampleMean)*(r2[i]-secondSampleMean);
      firstSampleVariance+=(r1[i]-firstSampleMean)*(r1[i]-firstSampleMean);    
      secondSampleVariance+=(r2[i]-secondSampleMean)*(r2[i]-secondSampleMean);
    }
   
    if(r==0)
      return 0;

    r/=sqrt(firstSampleVariance)*sqrt(secondSampleVariance);

    return r;
  }

  double rCorr(std::vector<double> r1, std::vector<int> r2)
  {
    double firstSampleMean=0;
    double secondSampleMean=0;
    double r=0;
    double rv=0;
    double firstSampleVariance=0;
    double secondSampleVariance=0;
    for(int i=0; i<r1.size(); i++)
    {
      firstSampleMean+=r1[i];
      secondSampleMean+=r2[i];
    }
    firstSampleMean/=r1.size();
    secondSampleMean/=r2.size();
    for(int i=0; i<r1.size(); i++)
    {
      r+=(r1[i]-firstSampleMean)*(r2[i]-secondSampleMean);
      firstSampleVariance+=(r1[i]-firstSampleMean)*(r1[i]-firstSampleMean);    
      secondSampleVariance+=(r2[i]-secondSampleMean)*(r2[i]-secondSampleMean);
    }
   
    if(r==0)
      return 0;

    rv=r;

    r/=sqrt(firstSampleVariance)*sqrt(secondSampleVariance);

    if(fabs(r-1)<0.00001)
      return (rv-1)/(sqrt(firstSampleVariance)*sqrt(secondSampleVariance));

    if(fabs(r+1)<0.00001)
      return (rv+1)/(sqrt(firstSampleVariance)*sqrt(secondSampleVariance));

    return r;   
  }

  void formPooledSample()
  {
    _pooledSample.insert(_pooledSample.begin(), 
			 _firstSample.begin(), _firstSample.end());
    _indicators.assign(_firstSample.size(), 1);
    _pooledSample.insert(_pooledSample.end(), 
			 _secondSample.begin(), _secondSample.end());
    _indicators.insert(_indicators.end(),_secondSample.size(), -1);

    for(size_t i=0; i<_pooledSample.size(); i++)
    {
      size_t min_ind=i;
      for(size_t j=i+1; j<_indicators.size(); j++)
	if(_pooledSample[j].getValue()<_pooledSample[min_ind].getValue())
	  min_ind=j;

      Datum tmpDat=_pooledSample[i];
      _pooledSample[i]=_pooledSample[min_ind];
      _pooledSample[min_ind]=tmpDat;

      double tmp=_indicators[i];
      _indicators[i]=_indicators[min_ind];
      _indicators[min_ind]=tmp;

    }

    for(size_t i=0; i<_pooledSample.size(); )
    {
      size_t j=i+1;
      for( ; (j<_pooledSample.size()) && 
	     (_pooledSample[i].getValue()==_pooledSample[j].getValue()); j++)
	;

      for(size_t k=i; k<j; k++)
	_ranks.push_back((i+1+j)/2.0);
      
      i+=j-i;
    }
  }

  class CompareData
  {
  public:
    bool operator() (const Datum &i, const Datum &j)
    {
      if(i.valueKnown() && !j.valueKnown())
	return true;
      if(!i.valueKnown() && j.valueKnown())
	return false;
      return i.getValue()<j.getValue();
    }
  };

  double testStatistic()
  {
    double W=0;
    size_t n=_firstSample.size();
    size_t m=_secondSample.size();
    for(size_t i=0; i<n; i++)
      for(size_t j=0; j<m; j++)
	if(_firstSample[i].getValue()<_secondSample[j].getValue())
	  W--;
	else if(_firstSample[i].getValue()>_secondSample[j].getValue())
	  W++;
    return W/(m*n);
  }

  double mantelNullVariance()
  {
    double v=0;
    size_t n1=_firstSample.size();
    size_t n2=_secondSample.size();
    for(size_t i=0; i<_pooledSample.size(); i++)
    {
      int u=0;
      for(size_t j=0; j<_pooledSample.size(); j++)
	if(_pooledSample[j].getValue()<_pooledSample[i].getValue())
	  u++;
	else if(_pooledSample[j].getValue()>_pooledSample[i].getValue())
	  u--;
      v+=u*u;
    }
    return v/(n1*n2*(n1+n2)*(n1+n2+1));
  }

  double gehanNullVariance()
  {
    std::vector<double> deathTimes;
    for(size_t i=0; i<_firstSample.size(); i++)
      if(_firstSample[i].valueKnown())
	deathTimes.push_back(_firstSample[i].getValue());
    for(size_t i=0; i<_secondSample.size(); i++)
      if(_secondSample[i].valueKnown())
	deathTimes.push_back(_secondSample[i].getValue());
    sort(deathTimes.begin(), deathTimes.end());
    std::vector<double>::iterator newEnd=unique(deathTimes.begin(), deathTimes.end());
    deathTimes.erase(newEnd, deathTimes.end());    
    
    int n=_firstSample.size();
    int m=_secondSample.size();

    int M=0;
    double var;

    for(size_t i=0; i<deathTimes.size(); i++)
    {
      int mi=_pooledSample.diedAtTime(deathTimes[i]);
      var+=mi*M*(M+1)+mi*(m+n-(M+mi))*(m+n-3*M-mi-1);
      M+=mi;      
    }
    return var/(m*n*(m+n)*(m+n-1));
  }

  double bootstrapVariance(int num)
  {
    double mean=0;
    double var=0;
    std::vector<double> first;
    std::vector<int> second;
    std::vector<double> rhos;
    for(int i=0; i<num; i++)
    {
      first.clear();
      second.clear();
      bool same=true;
      int last;
      for(int j=0; j<_ranks.size(); j++)
      {
	int t=(rand()/(1+(double)RAND_MAX))*_ranks.size();
	if(j==0)
	  last=t;
	else if(t!=last)
	  same=false;
	first.push_back(_ranks[t]);
	second.push_back(_indicators[t]);
      }
      if(same)
	rhos.push_back(0);
      else
	rhos.push_back(rCorr(first,second));
    }
    for(int i=0; i<rhos.size(); i++)
      mean+=rhos[i];
    mean/=rhos.size();
    for(int i=0; i<rhos.size(); i++)
      var+=(rhos[i]-mean)*(rhos[i]-mean);    
    var/=rhos.size()-1;
    return var;
  }

  double jackKnifeVariance()
  {
    double mean=0;
    double var=0;
    std::vector<double> rhos;
    std::vector<double> first;
    std::vector<int> second;
    for(int i=0; i<_ranks.size(); i++)
    {
      first.clear();
      second.clear();
      for(int j=0; j<_ranks.size(); j++)
      {
	if(i==j)
	  continue;
	first.push_back(_ranks[j]);
	second.push_back(_indicators[j]);
      }
      rhos.push_back(rValue(first,second));
    }
    for(int i=0; i<rhos.size(); i++)
      mean+=rhos[i];
    mean/=rhos.size();
    for(int i=0; i<rhos.size(); i++)
      var+=(rhos[i]-mean)*(rhos[i]-mean);
    var*=rhos.size()-1;
    var/=rhos.size();
    return var;
  }
  


private:
  Sample &_firstSample;
  Sample &_secondSample;
  Sample _pooledSample;
  std::vector<int> _indicators;
  std::vector<double> _ranks;
  double _testStatistic;
  double _zScore;
  double _rValue;
  double _variance;
  double _pValue;
};


#endif
