/*****************************************************************************
  myDraw.cpp
  Last-modified: 19 Oct 2013 05:51:43 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include "myDraw.h"

using namespace std;

std::string get_file_contents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
  return (string(""));
}





int main(int argc, char* argv[])
{
    string infile, outfile;
	string contents, result;
	infile = argv[1];
	outfile = argv[2];
	contents = get_file_contents(argv[1]);
	StructureImageHandler* handler = new StructureImageHandler();
	//handler->setNucleotidesCircled( true );
	handler->readRadial( contents, 1); 
	result = handler->writePostscript();
	cout << result << endl;
	//handler->writeSVG( outfile );
//	handler->writePostscript( outputFile);
	delete handler;
	return 0;
}

