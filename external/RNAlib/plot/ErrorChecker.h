/*
 * A header file for a template that checks errors in a variety of objects used
 * throughout RNAstructure.
 * Note that the isErrorStatus methods can be treated as returning either an
 * integer or a boolean, because a status of no error is always 0.
 * Therefore, an error can be treated as true, while no error is treated as
 * false. Hence, an error code not equal to 0 is a (true) error status.
 *
 * (c) 2009  Mathews Lab, University of Rochester Medical Center
 * Written by Jessica Reuter
 */

#ifndef ERRORCHECKER_H
#define ERRORCHECKER_H

#include <iostream>
using namespace std;

template <typename T>
class ErrorChecker {
 public:
  // Public constructor and methods

  /*
   * Name:        Constructor
   * Description: Initalizes an error checker with an object.
   * Arguments:
   *     1. newObject
   *        The object to check
   */
  ErrorChecker( T* newObject );

  /*
   * Name:        isErrorStatus
   * Description: One of two overloaded methods that checks the viability of a
   *              strand of nucleic acids. This particular variant checks the
   *              result of the object method GetErrorCode.
   * Arguments:
   *     1. print
   *        True if error message should be printed to stdout, false if not.
   *        Default is true (printing enabled).
   * Returns:
   *        Error code that was checked: zero if no error, nonzero if error.
   */
  int isErrorStatus( bool print = true );

  /*
   * Name:        isErrorStatus
   * Description: One of two overloaded methods that checks the viability of a
   *              strand of nucleic acids. This particular variant checks an
   *              error code explicitly given to it.
   * Arguments:
   *     1. code
   *        The explicit error code.
   *     2. print
   *        True if error message should be printed to stdout, false if not.
   *        Default is true (printing enabled).
   * Returns:
   *        Error code that was checked: zero if no error, nonzero if error.
   */
  int isErrorStatus( int code, bool print = true );

  /*
   * Name:        returnError
   * Description: One of two overloaded methods that checks the viability of a
   *              strand of nucleic acids. Unlike the isErrorStatus methods,
   *              this method returns the error string, rather than printing it
   *              out. This particular variant checks the result of the object
   *              method GetErrorCode.
   * Returns:
   *        The error message.
   */
  string returnError();

  /*
   * Name:        returnError
   * Description: One of two overloaded methods that checks the viability of a
   *              strand of nucleic acids. Unlike the isErrorStatus methods,
   *              this method returns the error string, rather than printing it
   *              out. This particular variant checks an error code explicitly
   *              given to it.
   * Arguments:
   *     1. code
   *        The explicit error code.
   * Returns:
   *        The error message.
   */
  string returnError( int code );

 private:
  // The object which is checked for errors.
  T* object;
};

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
ErrorChecker<T>::ErrorChecker( T* newObject ) {

  object = newObject;

}

///////////////////////////////////////////////////////////////////////////////
// Check the error status with the general GetErrorCode method.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
int ErrorChecker<T>::isErrorStatus( bool print ) {

  // If object has been initialized, get its error message.
  if( object != 0 ) {
    return isErrorStatus( object->GetErrorCode(), print );
  }

  // Otherwise, if the object is not initialized in the checker, show an error
  // and return -1.
  else {
    if( print ) {
      cerr << endl
	   << "Object is uninitialized; cannot be checked for errors." << endl;
    }
    return -1;
  }

}

///////////////////////////////////////////////////////////////////////////////
// Check the error status of an explicit error code.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
int ErrorChecker<T>::isErrorStatus( int code, bool print ) {

  // If error code isn't equal to 0, prepare to get a specific error message.
  if( code != 0 && print ) {

    // If object has been initialized, get its error message.
    if( object != 0 ) {
      cerr << endl << object->GetErrorMessage( code ) << endl;
    }

    // Otherwise, if the object is initialized in the checker, show an error
    // and set the error code to -1.
    else {
      cerr << endl
	   << "Object is uninitialized; cannot be checked for errors." << endl;
      code = -1;
    }
  }

  // Return the error code.
  return code;

}

///////////////////////////////////////////////////////////////////////////////
// Return an error from the general GetErrorCode method.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
string ErrorChecker<T>::returnError() {

  // If object has been initialized, get its error message.
  if( object != 0 ) {
    return returnError( object->GetErrorCode() );
  }

  // Otherwise, if the object is not initialized in the checker, show an error.
  else {
    return "Object is uninitialized; cannot be checked for errors.";
  }

}

///////////////////////////////////////////////////////////////////////////////
// Return an error from an explicit error code.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
string ErrorChecker<T>::returnError( int code ) {

  // Initialize the error string.
  string errorString = "";

  // If error code isn't equal to 0, prepare to get a specific error message.
  if( code != 0 ) {

    // If object has been initialized, get its error message.
    if( object != 0 ) {
      errorString = object->GetErrorMessage( code );
    }

    // Otherwise, if the object is initialized in the checker, show an error.
    else {
      errorString = "Object is uninitialized; cannot be checked for errors.";
    }
  }

  // Return the error string.
  return errorString;

}

#endif /* ERRORCHECKER_H */
