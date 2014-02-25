/*****************************************************************************
  test2.cpp
  Last-modified: 19 Oct 2013 03:48:22 PM

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

using namespace std;
#define ctheaderlength 250

int main(int argc, char* argv[])
{
	char base[2],header[ctheaderlength],temp[1000];
	ifstream fin;
	fin.open(argv[1]);
	int count;
	fin >> count;
	cout << count << endl;
	fin.getline(header, ctheaderlength -1);
	cout << header << endl;
	fin.close();
    return 0;
}

