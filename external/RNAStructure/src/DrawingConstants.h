#ifndef DRAWCONSTANTS_H
#define DRAWCONSTANTS_H

#include <sstream>
using namespace std;

/******************************************/
/* Text sizing.                           */
/******************************************/

const int GLYPHSIZE = 15;
const int TEXTSIZE = 24;
const string TEXTSIZE_STRING = "24";
const int TEXTSIZE_LEGEND = 16;
const string TEXTSIZE_LEGEND_STRING = "16";
const int TEXTSIZE_BIGGER = 26;
const string TEXTSIZE_BIGGER_STRING = "26";
const int TEXTSIZE_LEGEND_BIGGER = 18;
const string TEXTSIZE_LEGEND_BIGGER_STRING = "18";

/******************************************/
/* Miscellaneous layout properties.       */
/******************************************/

const int BORDER = 36;
const string BORDER_STRING = "36";
const int CIRCLE_RADIUS = 12;
const string CIRCLE_RADIUS_STRING = "12";
const double NUC_LABEL_HEIGHT = 24;
const string NUC_LABEL_HEIGHT_STRING = "24";
const int NUC_POS_ADJUSTMENT = 4;
const double PI = 3.14159;

/******************************************/
/* Color identification strings.          */
/******************************************/

// Generic color strings.
const string BLACK = "Black";
const string WHITE = "White";
const string GRAY = "Gray";
const string RED = "Red";
const string PINK = "Pink";
const string ORANGE = "Orange";
const string YELLOW = "Yellow";
const string LIGHT_GREEN = "Light Green";
const string GREEN = "Green";
const string LIGHT_BLUE = "Light Blue";
const string BLUE = "Blue";
const string PURPLE = "Purple";

// Inline function to get a particular type of color string.
inline string getColorString( string type, bool isSVG ) {
	if( type == WHITE ) {
		if( !isSVG ) { return "1.00 1.00 1.00"; }
		else { return "\"rgb(255,255,255)\""; }
	} else if( type == GRAY ) {
		if( !isSVG ) { return "0.67 0.67 0.67"; }
		else { return "\"rgb(171,171,171)\""; }
	} else if( type == RED ) {
		if( !isSVG ) { return "1.00 0.00 0.00"; }
		else { return "\"rgb(255,0,0)\""; }
	} else if( type == PINK ) {
		if( !isSVG ) { return "1.00 0.50 1.00"; }
		else { return "\"rgb(255,128,255)\""; }
	} else if( type == ORANGE ) {
		if( !isSVG ) { return "1.00 0.50 0.00"; }
		else { return "\"rgb(255,171,0)\""; }
	} else if( type == YELLOW ) {
		if( !isSVG ) { return "0.83 0.83 0.17"; }
		else { return "\"rgb(212,212,44)\""; }
	} else if( type == LIGHT_GREEN ) {
		if( !isSVG ) { return "0.00 1.00 0.00"; }
		else { return "\"rgb(0,255,0)\""; }
	} else if( type == GREEN ) {
		if( !isSVG ) { return "0.00 0.50 0.00"; }
		else { return "\"rgb(0,128,0)\""; }
	} else if( type == LIGHT_BLUE ) {
		if( !isSVG ) { return "0.00 0.67 1.00"; }
		else { return "\"rgb(0,171,255)\""; }
	} else if( type == BLUE ) {
		if( !isSVG ) { return "0.00 0.00 1.00"; }
		else { return "\"rgb(0,0,255)\""; }
	} else if( type == PURPLE ) {
		if( !isSVG ) { return "0.50 0.00 0.50"; }
		else { return "\"rgb(128,0,128)\""; }
	} else {
		if( !isSVG ) { return "0.00 0.00 0.00"; }
		else { return "\"rgb(0,0,0)\""; }
	}
}

/******************************************/
/* Image size constants.                  */
/******************************************/

// Maximum Postscript description length and image bounds.
const unsigned int DESC_PS = 39;
const int XBOUND_PS = 612;
const int YBOUND_PS = 792;

// Maximum SVG description length and image bounds.
const unsigned int DESC_SVG = 45;
const int XBOUND_SVG = 790;
const int YBOUND_SVG = 905;

/******************************************/
/* Drawing object template pieces.        */
/******************************************/

// Definitions of standardized variables used in structure element templates.
const string BACKGROUND = "BACKGROUND";
const string COLOR = "COLOR";
const string CONTROLX = "CONTROLX";
const string CONTROLY = "CONTROLY";
const string CURVEWEIGHT = "CURVEWEIGHT";
const string ENDX = "ENDX";
const string ENDY = "ENDY";
const string HEIGHT = "HEIGHT";
const string LINEWEIGHT = "LINEWEIGHT";
const string LOCX = "LOCX";
const string LOCY = "LOCY";
const string OUTLINE = "OUTLINE";
const string RADIUS = "RADIUS";
const string SCALEFACTOR = "SCALEFACTOR";
const string STARTX = "STARTX";
const string STARTY = "STARTY";
const string TEXTSTRING = "TEXTSTRING";
const string WIDTH = "WIDTH";
const string X1 = "X1";
const string X2 = "X2";
const string Y1 = "Y1";
const string Y2 = "Y2";

// Postscript color template.
const string COLOR_TEMPLATE_PS = "RED GREEN BLUE";

// Postscript scaling region markers. 
const string SCALE_OPEN_PS = "gsave " + SCALEFACTOR + " " + SCALEFACTOR + " scale";
const string SCALE_CLOSE_PS = "grestore";

// Postscript legend resizing markers.
const string LEGEND_RESIZE_START_PS = "[" + TEXTSIZE_LEGEND_STRING + " 0 0 -" + TEXTSIZE_LEGEND_STRING + " 0 0] /Courier-Bold sfm";
const string LEGEND_RESIZE_END_PS = "";

// Postscript syntax strings to draw specific structural elements.
const string CIRCLE_PS = LINEWEIGHT + " setlinewidth newpath " + OUTLINE + " setrgbcolor " + LOCX + " " + LOCY + " " + RADIUS + " 0 360 arc closepath gsave " + BACKGROUND + " setrgbcolor fill grestore stroke";
const string CURVE_PS = COLOR + " setrgbcolor " + CURVEWEIGHT + " setlinewidth " + X1 + " " + Y1 + " moveto " + X1 + " " + Y1 + " " + CONTROLX + " " + CONTROLY + " " + X2 + " " + Y2 + " curveto stroke";
const string LINE_PS = COLOR + " setrgbcolor " + LINEWEIGHT + " setlinewidth newpath " + STARTX + " " + STARTY + " moveto " + ENDX + " " + ENDY + " lineto closepath stroke";
const string RECTANGLE_PS = COLOR + " setrgbcolor newpath " + LOCX + " " + LOCY + " moveto 0 " + HEIGHT + " rlineto " + WIDTH + " 0 rlineto 0 -" + HEIGHT + " rlineto closepath fill";
const string TEXT_PS = LOCX + " " + LOCY + " moveto " + COLOR + " setrgbcolor (" + TEXTSTRING + ") show";

// SVG color template.
const string COLOR_TEMPLATE_SVG = "\"rgb(RED,GREEN,BLUE)\"";

// SVG scaling region markers.
const string SCALE_OPEN_SVG = "<g transform=\"scale(" + SCALEFACTOR + ")\">";
const string SCALE_CLOSE_SVG = "</g>";

// SVG legend resizing markers.
const string LEGEND_RESIZE_START_SVG = "<g font-size=\"" + TEXTSIZE_LEGEND_STRING + "\">";
const string LEGEND_RESIZE_END_SVG = "</g>";

// SVG syntax strings to draw specific structural elements.
const string CIRCLE_SVG = "<circle style=\"stroke-width:" + LINEWEIGHT + "\" cx=\"" + LOCX + "\" cy=\"" + LOCY + "\" r=\"" + RADIUS + "\" fill=" + BACKGROUND + " stroke=" + OUTLINE + "/>";
const string CURVE_SVG = "<path style=\"fill:none;stroke-width:" + CURVEWEIGHT + "\" stroke=" + COLOR + " d=\"M" + X1 + "," + Y1 + " Q" + CONTROLX + "," + CONTROLY + " " + X2 + "," + Y2 + "\"/>";
const string LINE_SVG = "<line style=\"fill:none;stroke-width:" + LINEWEIGHT + "\" stroke=" + COLOR + " x1=\"" + STARTX + "\" y1=\"" + STARTY + "\" x2=\"" + ENDX + "\" y2=\"" + ENDY + "\"/>";
const string RECTANGLE_SVG = "<rect x=\"" + LOCX + "\" y=\"" + LOCY + "\" width=\"" + WIDTH + "\" height=\"" + HEIGHT + "\" fill=" + COLOR + " stroke=" + COLOR + "/>";
const string TEXT_SVG = "<text x=\"" + LOCX + "\" y=\"" + LOCY + "\" fill=" + COLOR + " stroke=" + COLOR + ">" + TEXTSTRING + "</text>";

/******************************************/
/* File boundary marker creation.         */
/******************************************/

// Postscript file markers.
inline string createStartPS() {
	stringstream startMarkerStream( stringstream::in | stringstream::out );
	startMarkerStream << "%!" << endl
	                  << "0 " << YBOUND_PS << " translate 1 -1 scale" << endl
	                  << "/sfm { findfont exch makefont setfont } bind def" << endl
	                  << "[" << TEXTSIZE << " 0 0 " << -TEXTSIZE << " 0 0] /Courier-Bold sfm";
	return startMarkerStream.str();
}
const string END_MARKER_PS = "showpage";

// SVG file markers.
inline string createStartSVG() {
	stringstream startMarkerStream( stringstream::in | stringstream::out );
	startMarkerStream << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl
	                  << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
	                  << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
	                  << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
	                  << "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
	                  << "xml:space=\"preserve\" font-family=\"monospace\" font-size=\"" << TEXTSIZE << "\" "
	                  << "fill=" << getColorString( WHITE, true ) << " stroke=" << getColorString( BLACK, true ) << " "
	                  << "viewBox=\"0 0 " << XBOUND_SVG << " " << YBOUND_SVG << "\">";
	return startMarkerStream.str();
}
const string END_MARKER_SVG = "</svg>";

#endif /* DRAWINGCONSTANTS_H */
