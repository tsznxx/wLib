/*
 * An implementation file for a program that parses Unix-like command lines to
 * determine the data being given for input.
 *
 * (c) 2008  Mathews Lab, University of Rochester Medical Center
 * Redone in 2012.
 * Written by Jessica S. Reuter
 */

#include "ParseCommandLine.h"
#include <ctype.h>

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
ParseCommandLine::ParseCommandLine( string name ) {
	//software version
  	version = 5.6;

	// Initialize the error flag and the specialized usage flags.
	error = false;
	specializedUsage = false;
	specializedUsageSet = false;

	// Initialize the interface name and usage string.
	interfaceName = name;
	usageString = "USAGE: " + name + " ";

	// Add the help flag, which is found in all text interfaces.
	standardHelpFlags.push_back( "-h" );
	standardHelpFlags.push_back( "-H" );
	standardHelpFlags.push_back( "--help" );
	addOptionFlagsNoParameters( standardHelpFlags, "Display the usage details message." );

	// Add the version flag, which is found in all text interfaces.
	standardVersionFlags.push_back( "-v" );
	standardVersionFlags.push_back( "-V" );
	standardVersionFlags.push_back( "--version" );
	addOptionFlagsNoParameters( standardVersionFlags, "Display version and copyright information for this interface." );
}

///////////////////////////////////////////////////////////////////////////////
// Add a group of option flags that don't have parameters.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addOptionFlagsNoParameters(
	vector<string> list, string text ) {

	// Build the vector list into a string.
	string flags = "";
	for( unsigned int i = 1; i <= list.size(); i++ ) {
		string flag( list[i-1] );
		flags += ( flag + " " );
	}

	// Add the full option description.
	descriptionsOfOptionsNoFlags[flags] = text;

	// Create a transformed copy of the flag list in lower case, then set
	// that transformed copy as the list.
	// While creating the transformed copy, add the transformed flags to the
	// flag set.
	vector<string> copy;
	for( unsigned int i = 1; i <= list.size(); i++ ) {
		string flag( list[i-1] );
		transform( flag.begin(), flag.end(), flag.begin(), ::tolower );
		copy.push_back( flag );
		lowerOptionsNoFlags.insert( flag );
	}
	list = copy;
}

///////////////////////////////////////////////////////////////////////////////
// Add a group of option flags that have parameters.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addOptionFlagsWithParameters(
	vector<string> list, string text ) {

	// Build the vector list into a string.
	string flags = "";
	for( unsigned int i = 1; i <= list.size(); i++ ) {
		string flag( list[i-1] );
		flags += ( flag + " " );
	}

	// Add the full option description.
	descriptionsOfOptionsWithFlags[flags] = text;

	// Create a transformed copy of the flag list in lower case, then set
	// that transformed copy as the list.
	// While creating the transformed copy, add the transformed flags to the
	// flag set.
	vector<string> copy;
	for( unsigned int i = 1; i <= list.size(); i++ ) {
		string flag( list[i-1] );
		transform( flag.begin(), flag.end(), flag.begin(), ::tolower );
		copy.push_back( flag );
		lowerOptionsWithFlags.insert( flag );
	}
	list = copy;
}

///////////////////////////////////////////////////////////////////////////////
// Add a parameter description.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::addParameterDescription(
	string id, string description ) {

	// Create the parameter string.
	string paramString = "<" + id + ">";

	// Update the usage string.
	usageString += ( paramString + " " );

	// Add the full parameter description.
	descriptionsOfParameters.push_back( make_pair( paramString, description ) );
}

///////////////////////////////////////////////////////////////////////////////
// Check whether the parser contains an option in a list.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::contains( vector<string> list ) {

	// Go through the possible flags that denote the option.
	// If one is found, return true.
	for( int i = 1; i <= list.size(); i++ ) {
		for (int j = 0; j < list[i-1].length();++j) {
			list[i-1][j] = tolower(list[i-1][j]);
		}
		if( parsedData.find( list[i-1] ) != parsedData.end() ) { 
			return true; 
		}
	}

	// Return false if no option was found.
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Get an optional value as a string.
///////////////////////////////////////////////////////////////////////////////
string ParseCommandLine::getOptionString( vector<string> list, bool exists ) {

	// Go through the possible flags that denote the option.
	for( int i = 1; i <= list.size(); i++ ) {
		for (int j = 0; j < list[i-1].length();++j) {
			list[i-1][j] = tolower(list[i-1][j]);
		}
		// Get the next option if it exists.
		string option = list[i-1];
		if( parsedData.find( option ) != parsedData.end() ) {

			// If string is meant to be an existing file name, but the file
			// doesn't exist, set an error and return an empty string.
			string file = parsedData[option];
			if( ( exists == true ) && ( ifstream( file.c_str() ) == NULL ) ) {
				cerr << "File required for option " << option
				     << " does not exist." << endl;
				setError();
				return "";
			}

			// Return the file.
			return file;
		}
	}

	// Return an empty string if no option was found.
	return "";
}

///////////////////////////////////////////////////////////////////////////////
// Get a parameter.
///////////////////////////////////////////////////////////////////////////////
string ParseCommandLine::getParameter( int number ) {

	// Create a string that would be this parameter's key in the data map.
	stringstream stream( stringstream::in | stringstream::out );
	stream << "param" << number;

	// Return the parameter's value.
	return parsedData[stream.str()];
}

///////////////////////////////////////////////////////////////////////////////
// Get whether the parser encountered an error.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::isError() {

	return error;
}

///////////////////////////////////////////////////////////////////////////////
// Get whether the class doesn't need to handle usage.
///////////////////////////////////////////////////////////////////////////////
bool ParseCommandLine::isSpecializedUsage() {

	return specializedUsageSet;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line into its pieces.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::parseLine( int argc, char* argv[] ) {

	// Add the options piece to the usage string, because there always are at
	// least two types of options (help and version information).
	usageString += "[options]";

	// Convert the command line into strings.
	// During this process, transform all flags into lower case.
	vector<string> cmdLine;
	for( int i = 1; i < argc; i++ ) {
		string next( argv[i] );
		if( next[0] == '-' ) {
			transform( next.begin(), next.end(), next.begin(), ::tolower );
		}
		cmdLine.push_back( next );
	}
	unsigned int length = cmdLine.size();

	// Check if a help flag exists.
	// If one does, print out a usage message, set an error, and return.
	// This isn't a real error, it's just to stop parsing and running later.
	vector<string>::iterator helpIt;
	for( unsigned int i = 1; i <= standardHelpFlags.size(); i++ ) {
		string next = standardHelpFlags[i-1];
		helpIt = find( cmdLine.begin(), cmdLine.end(), next );
		if( helpIt != cmdLine.end() ) {
			usage();
			setError();
			return;
		}
	}

	// Check if a version flag exists.
	// If one does, print out a version message, set an error, and return.
	// This isn't a real error, it's just to stop parsing and running later.
	vector<string>::iterator versionIt;
	for( unsigned int i = 1; i <= standardVersionFlags.size(); i++ ) {
		string next = standardVersionFlags[i-1];
		versionIt = find( cmdLine.begin(), cmdLine.end(), next );
		if( versionIt != cmdLine.end() ) {
			cout << interfaceName << ": Version " << version << "." << endl
			     << "Copyright Mathews Lab, University of Rochester." << endl;
			setError();
			return;
		}
	}

	// Populate the parsed data map.
	unsigned int numParams = 0;
	for( unsigned int i = 1; i <= cmdLine.size(); i++ ) {

		// Get the next command line argument.
		string next = cmdLine[i-1];

		// If the argument doesn't start with "-", call it a parameter.
		if( next[0] != '-' ) {
			numParams++;
			stringstream paramStream( stringstream::in | stringstream::out );
			paramStream << "param" << numParams;
			parsedData[paramStream.str()] = next;
		}

		// Otherwise, process the flag.
		else {
			// Create a variable to track whether the flag exists.
			bool exists = false;

			// Check the flag against the possible flags with no parameters.
			// Set the flag's value if found.
			set<string>::iterator it1;
			set<string> s1 = lowerOptionsNoFlags;
			for( it1 = s1.begin(); it1 != s1.end(); it1++ ) {
				string flag = *it1;
				if( next == flag ) {
					parsedData[next] = "";
					exists = true;
				}
			}

			// Check the flag against the possible options with parameters.
			// Set the flag's value if found.
			// Only do this if a flag wasn't previously found.
			if( exists == false ) {
				set<string>::iterator it2;
				set<string> s2 = lowerOptionsWithFlags;
				for( it2 = s2.begin(); it2 != s2.end(); it2++ ) {
					string flag = *it2;
					if( next == flag ) {
						bool noParam =
							( i == cmdLine.size() ) ||
							( cmdLine[i][0] == '-' );
						if( noParam == true ) {
							double test = 0.0;
							stringstream testStream( cmdLine[i] );
							if( testStream >> test ) { noParam = false; }
						}

						if( noParam ) {
							cerr << "Option missing for flag: " << next << endl;
							setError();
							return;
						} else {
							parsedData[next] = cmdLine[i];
							exists = true;
							i++;
						}
					}
				}
			}

			// If the flag doesn't exist, show an error and return.
			if( exists == false ) {
				cerr << "Flag " << next << " does not exist." << endl;
				setError();
				return;
			}
		}
	}

	// If the number of parameters in the command line isn't what it should
	// be, show an error and return.
	if( numParams != descriptionsOfParameters.size() ) {
		cerr << "Incorrect number of required parameters given." << endl
		     << usageString << endl
		     << "Use any of the following options to get a help message: ";
		for( unsigned int i = 1; i <= standardHelpFlags.size(); i++ ) {
			cerr << standardHelpFlags[i-1] << " ";
		}
		cerr << endl;
		setError();
		return;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setError() {

	error = true;
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setError( string type ) {

	setError();
	cerr << "Invalid " << type << " given." << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Set that an error occurred in parsing.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setErrorSpecialized( string errString ) {

	setError();
	cerr << errString << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as a double.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setOptionDouble(
	vector<string> list, double& defaultValue ) {

	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to a double and return that value.
	for( int i = 1; i <= list.size(); i++ ) {

		for (int j = 0; j < list[i-1].length();++j) {
			list[i-1][j] = tolower(list[i-1][j]);
		}
		string option = list[i-1];
		if( parsedData.find( option ) != parsedData.end() ) {

			// Try to read the given value as a double.
			// If that can't be done, show an error and set the default value
			// to infinity.
			stringstream stream( parsedData[option] );
			if( !( stream >> defaultValue ) ) {
				cerr << "Non-numeric input given for flag " << option << endl;
				setError();
				defaultValue = numeric_limits<double>::infinity();
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as a float.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setOptionFloat(
	vector<string> list, float& defaultValue ) {

	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to a float and return that value.
	for( int i = 1; i <= list.size(); i++ ) {

		for (int j = 0; j < list[i-1].length();++j) {
			list[i-1][j] = tolower(list[i-1][j]);
		}
		string option = list[i-1];
		if( parsedData.find( option ) != parsedData.end() ) {

			// Try to read the given value as a float.
			// If that can't be done, show an error and set the default value
			// to infinity.
			stringstream stream( parsedData[option] );
			if( !( stream >> defaultValue ) ) {
				cerr << "Non-numeric input given for flag " << option << endl;
				setError();
				defaultValue = numeric_limits<float>::infinity();
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Set an optional value as an integer.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setOptionInteger(
	vector<string> list, int& defaultValue ) {

	// Go through the possible flags that denote the option.
	// If one is found, attempt to convert it to an int and return that value.
	for( int i = 1; i <= list.size(); i++ ) {

		for (int j = 0; j < list[i-1].length();++j) {
			list[i-1][j] = tolower(list[i-1][j]);
		}
		string option = list[i-1];
		if( parsedData.find( option ) != parsedData.end() ) {

			// Try to read the given value as an int.
			// If that can't be done, show an error and set the default value
			// to infinity.
			stringstream stream( parsedData[option] );
			if( !( stream >> defaultValue ) ) {
				cerr << "Non-numeric input given for flag " << option << endl;
				setError();
				defaultValue = numeric_limits<int>::infinity();
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Tell the class that it doesn't need to handle usage.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::setSpecializedUsage() {

	specializedUsage = true;
}

///////////////////////////////////////////////////////////////////////////////
// Print out a detailed usage message.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::usage() {

	// If the standard usage message should not be printed out, return.
	if( specializedUsage == true ) {
		specializedUsageSet = true;
		return;
	}

	// Print out a short usage string and general flag information.
	cout << usageString << endl
	     << "All flags are case-insensitive, "
	     << "and grouping of flags is not allowed." << endl << endl;

	// Print out the required parameters.
	cout << "=============================" << endl
	     << "==== Required Parameters ====" << endl
	     << "=============================" << endl;
	unsigned int numParams = descriptionsOfParameters.size();
	for( unsigned int i = 1; i <= numParams; i++ ) {
		pair<string, string> nextPair = descriptionsOfParameters[i-1];
		cout << nextPair.first << endl;
		wrapString( nextPair.second );
	}

	// Print out option flags without parameters, if they exist.
	unsigned int numOptionsNoParams = descriptionsOfOptionsNoFlags.size();
	if( numOptionsNoParams >= 1 ) {
		cout << "=========================================" << endl
		     << "==== Option Flags Without Parameters ====" << endl
		     << "=========================================" << endl;
		map<string, string, compareOptions>::iterator it2;
		map<string, string, compareOptions> m2 = descriptionsOfOptionsNoFlags;
		for( it2 = m2.begin(); it2 != m2.end(); it2++ ) {
			cout << it2->first << endl;
			wrapString( it2->second );
		}
	}

	// Print out option flags with parameters, if they exist.
	unsigned int numOptionsWithParams = descriptionsOfOptionsWithFlags.size();
	if( numOptionsWithParams >= 1 ) {
		cout << "======================================" << endl
		     << "==== Option Flags With Parameters ====" << endl
		     << "======================================" << endl
		     << "All parameters must follow their associated flag directly."
		     << endl
		     << "Failure to do so may result in aberrant program behavior."
		     << endl << endl;
		map<string, string, compareOptions>::iterator it1;
		map<string, string, compareOptions> m1 = descriptionsOfOptionsWithFlags;
		for( it1 = m1.begin(); it1 != m1.end(); it1++ ) {
			cout << it1->first << endl;
			wrapString( it1->second );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Word wrap a string.
///////////////////////////////////////////////////////////////////////////////
void ParseCommandLine::wrapString( string text ) {

	// Determine the maximum width of the string to wrap.
	unsigned int bufferSize = 78;

	// Declare a variable that handles words in strings and a variable that
	// handles the description spacer.
	string word;
	string spacer = "    ";

	// Split the string into sentences.
	// Sentences are terminated by any word ending in a period.
	vector<string> sentences;
	string piece = "";
	stringstream sentenceStream( text );
	while( sentenceStream >> word ) {
		piece += word;
		if( piece.find_last_of( '.' ) != ( piece.length() - 1 ) ) {
			piece += " ";
		} else {
			sentences.push_back( spacer + piece );
			piece = "";
		}
	}
	unsigned int numSentences = sentences.size();

	// For each sentence, wrap it as necessary and print it out.
	for( unsigned int i = 1; i <= numSentences; i++ ) {

		// Get the next sentence.
		string sentence = sentences[i-1];

		// If the sentence is greater than the buffer size, add spacing to
		// wrap it properly.
		string fullSentence = "";
		unsigned int index = 0;
		int partSize = bufferSize - spacer.length();
		int pieceEnd = 0;
		if( sentence.length() > bufferSize ) {

			// While the current index is less than the sentence length, break
			// the sentence into pieces.
			while( index < sentence.length() ) {
				string piece = "";
				pieceEnd = index + partSize;

				// If the end of the next possible sentence piece is longer
				// than the sentence, just make the next piece the rest of
				// the sentence.
				if( pieceEnd > sentence.length() ) {
					piece = sentence.substr( index );
					index = sentence.length();
				}

				// Otherwise, pull out the next sentence piece and chop it off
				// at its last space.
				else {
					string fullPiece = sentence.substr( index, partSize );
					int lastSpace = fullPiece.find_last_of( ' ' ) + 1;
					piece = sentence.substr( index, lastSpace );
					index += lastSpace;
				}

				// Add spacers or new lines as necessary.
				if( piece[0] == ' ' ) {
					fullSentence += ( piece + '\n' );
				} else if( index == sentence.length() ) {
					fullSentence += ( spacer + piece );
				} else {
					fullSentence += ( spacer + piece + '\n' );
				}
			}

			// Set the newly spaced sentence as the real sentence.
			sentence = fullSentence;
		}

		// Print out the string.
		cout << sentence << endl;
	}
	cout << endl;
}
