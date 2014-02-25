/*****************************************************************************
  test.cpp
  Last-modified: 19 Oct 2013 05:12:31 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[])
{
	istream *in;
	ifstream fin;
	ostream *out;
	ofstream fout;
	ostringstream sout;
	string outstring;
	int count;
	istringstream ss("test1\ntest2");
	fin.open(argv[1]);
	if(fin.is_open())
		in = &fin;
	else
		in = &ss;
	std::string line;    
	*in >> count;
	sout << count << endl;
	while (std::getline(*in, line))
		sout << line << std::endl;
	if(fin.is_open())
		fin.close();
	outstring = sout.str();
	ostream *mycout = &cout;
	*mycout << "test" << endl;
	return 0;
}

