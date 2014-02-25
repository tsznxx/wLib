/*
 * A class that holds structure image data and can write those images to files if necessary.
 *
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#ifndef STRUCTURE_IMAGE_HANDLER_H
#define STRUCTURE_IMAGE_HANDLER_H

#include <limits>
#include <map>
#include <vector>

#include "RNA.h"
#include "DrawingConstants.h"
#include "ErrorChecker.h"

using namespace std;

class StructureImageHandler {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes private variables.
         */
	StructureImageHandler();

	/*
	 * Name:        addAnnotationProbability
	 * Description: Add probability annotation to the structure.
	 * Arguments:
	 *     1. file
	 *        The file that holds probability annotation information.
	 *     2. text
	 *        True if the annotation file is a dot plot text file, false if not.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string addAnnotationProbability( string file, bool text = false );

	/*
	 * Name:        addAnnotationSHAPE
	 * Description: Add SHAPE annotation to the structure.
	 * Arguments:
	 *     1. file
	 *        The file that holds SHAPE annotation information.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string addAnnotationSHAPE( string file );

	/*
	 * Name:        flipHorizontally
	 * Description: Flip the structure image horizontally.
	 */ 
	void flipHorizontally();

	/*
	 * Name:        readCircular
	 * Description: Read a structure and determine its coordinates for a
	 *              circular layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readCircular( string file, int number );

	/*
	 * Name:        readLinear
	 * Description: Read a structure and determine its coordinates for a
	 *              linear layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readLinear( string file, int number );

	/*
	 * Name:        readRadial
	 * Description: Read a structure and determine its coordinates for a
	 *              radial layout.
	 * Arguments:
	 *     1. file
	 *        The file that holds structure information.
	 *     2. number
	 *        The structure to get information about.
	 * Returns:
	 *     A string showing the completion status.
	 */ 
	string readRadial( string file, int number );

	/*
	 * Name:        removeAnnotation
	 * Description: Remove annotation from nucleotides.
	 */
	void removeAnnotation();

	/*
	 * Name:        setNucleotidesCircled
	 * Description: Set nucleotides to be surrounded by circles when drawn.
	 * Arguments:
	 *     1. encircle
	 *        True if nucleotides should be circled, false if not.
	 */
	void setNucleotidesCircled( bool encircle );

	/*
	 * Name:        toString
	 * Description: Return a string representation of the structure handled by
	 *              this object. This is usually a Java convention, not a C++
	 *              one, but it is useful to get a snapshot of this structure
	 *              and access all information about it with a single easy call.
	 */
	string toString();

	/*
	 * Name:        writePostscript
	 * Description: Write a structure as a Postscript image.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 *     2. append
	 *        True if structures should be appended on to this file, false
	 *        if a new structure file should be created for writing.
	 */
	void writePostscript( string file);
	string writePostscript();
	/*
	 * Name:        writeSVG
	 * Description: Write a structure as an SVG image.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 */
	void writeSVG( string file );
	string writeSVG();

 private:
	// Private image writing workhorse function.

	/*
	 * Name:        writeImageFile
	 * Description: Write a structure as an image.
	 *              This method encapsulates the common strategy for
	 *              writing structure images in different formats.
	 * Arguments:
	 *     1. file
	 *        The file name to write.
	 *     2. append
	 *        True if structures should be appended on to this file, false
	 *        if a new structure file should be created for writing.
	 *     3. isSVG
	 *        True if the file should be SVG, false if not.
	 */
	void writeImageFile( string file, bool append, bool isSVG );
	// return a string of the svg or ps file.
	string writeImageString ( bool isSVG);

 protected:
	// Protected data structures.

	// Annotation colors.
	vector<string> annotations;

	// Backbone, nucleotide, and label coordinates.
	vector<string> coordinates;

	// Legend text.
	vector<string> legend;

	// Legend text colors.
	vector<string> legendColors;

	// Any extra data that is necessary to place on the image.
	// Note that this data is assumed to be formatted properly for a given
	// image type.
	vector<string> extras;

	// Pairing curve coordinates.
	vector<string> pairs;

 protected:
	// Protected variables.

	// A boolean flag specifying if the structure is bimolecular (true) or
	// not (false).
	bool bimolecular;

	// A boolean flag specifying if the nucleotides are circled (true) or
	// not (false).
	bool circleNucs;

	// The description of the image.
	string description;

	// The maximum X bound of the image.
	double maxX;

	// The maximum Y bound of the image.
	double maxY;
};

#endif /* STRUCTURE_IMAGE_HANDLER_H */
