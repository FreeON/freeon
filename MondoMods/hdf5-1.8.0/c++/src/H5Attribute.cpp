/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#include "H5Include.h"
#include "H5Exception.h"
#include "H5IdComponent.h"
#include "H5PropList.h"
#include "H5Object.h"
#include "H5AbstractDs.h"
#include "H5Attribute.h"
#include "H5DcreatProp.h"
#include "H5CommonFG.h"
#include "H5DataType.h"
#include "H5DataSpace.h"
#include "H5private.h"

#ifndef H5_NO_NAMESPACE
namespace H5 {
#ifndef H5_NO_STD
    using std::cerr;
    using std::endl;
#endif  // H5_NO_STD
#endif

//--------------------------------------------------------------------------
// Function:	Attribute default constructor
///\brief	Default constructor: Creates a stub attribute
// Programmer	Binh-Minh Ribler - May, 2004
//--------------------------------------------------------------------------
Attribute::Attribute() : AbstractDs() {}

//--------------------------------------------------------------------------
// Function:	Attribute copy constructor
///\brief	Copy constructor: makes a copy of the original Attribute object.
///\param	original  - IN: Original Attribute object to copy
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
Attribute::Attribute( const Attribute& original ) : AbstractDs( original ) {}

//--------------------------------------------------------------------------
// Function:	Attribute overloaded constructor
///\brief	Creates an Attribute object using the id of an existing
///		attribute.
///\param	existing_id - IN: Id of an existing attribute
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
Attribute::Attribute(const hid_t existing_id) : AbstractDs(existing_id) {}

//--------------------------------------------------------------------------
// Function:	Attribute::write
///\brief	Writes data to this attribute.
///\param	mem_type  - IN: Attribute datatype (in memory)
///\param	buf       - IN: Data to be written
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
void Attribute::write( const DataType& mem_type, const void *buf ) const
{
   herr_t ret_value = H5Awrite( id, mem_type.getId(), buf );
   if( ret_value < 0 )
   {
      throw AttributeIException("Attribute::write", "H5Awrite failed");
   }
}

//--------------------------------------------------------------------------
// Function:	Attribute::write
///\brief	This is an overloaded member function, provided for convenience.
///		It writes a \a std::string to this attribute.
///\param	mem_type  - IN: Attribute datatype (in memory)
///\param	strg      - IN: Data to be written
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - Apr, 2003
//--------------------------------------------------------------------------
void Attribute::write( const DataType& mem_type, const H5std_string& strg ) const
{
   // Convert string to C-string
   const char* strg_C;
   strg_C = strg.c_str();  // strg_C refers to the contents of strg as a C-str

   herr_t ret_value = H5Awrite( id, mem_type.getId(), &strg_C );
   if( ret_value < 0 )
   {
      throw AttributeIException("Attribute::write", "H5Awrite failed");
   }
}

//--------------------------------------------------------------------------
// Function:	Attribute::read
///\brief	Reads data from this attribute.
///\param	mem_type -  IN: Attribute datatype (in memory)
///\param	buf      - OUT: Buffer for read data
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
void Attribute::read( const DataType& mem_type, void *buf ) const
{
   herr_t ret_value = H5Aread( id, mem_type.getId(), buf );
   if( ret_value < 0 )
   {
      throw AttributeIException("Attribute::read", "H5Aread  failed");
   }
}

//--------------------------------------------------------------------------
// Function:	Attribute::read
///\brief	This is an overloaded member function, provided for convenience.
///		It reads a \a std::string from this attribute.
///\param	mem_type  - IN: Attribute datatype (in memory)
///\param	strg      - IN: Buffer for read string
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - Apr, 2003
// Modification
//		2006/12/9 - H5Aread allocates memory for character string
//			    buffer with malloc, therefore, no allocation here,
//			    but HDfree is needed. - BMR
//--------------------------------------------------------------------------
void Attribute::read( const DataType& mem_type, H5std_string& strg ) const
{
   char* strg_C;  // temporary C-string for C API

   // call C API to get the attribute string of chars
   herr_t ret_value = H5Aread( id, mem_type.getId(), &strg_C);
   if( ret_value < 0 )
   {
      throw AttributeIException("Attribute::read", "H5Aread failed");
   }
   strg = strg_C;	// get 'string' from the C char*
   HDfree(strg_C);
}

//--------------------------------------------------------------------------
// Function:	Attribute::getSpace
///\brief	Gets a copy of the dataspace for this attribute.
///\return	Dataspace instance
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
DataSpace Attribute::getSpace() const
{
   // Calls C function H5Aget_space to get the id of the dataspace
   hid_t dataspace_id = H5Aget_space( id );

   // If the dataspace id is valid, create and return the DataSpace object
   if( dataspace_id > 0 )
   {
      DataSpace dataspace( dataspace_id );
      return( dataspace );
   }
   else
   {
      throw AttributeIException("Attribute::getSpace", "H5Aget_space failed");
   }
}

//--------------------------------------------------------------------------
// Function:	Attribute::p_get_type (private)
// Purpose	Gets the datatype of this attribute.
// Return	Id of the datatype
// Exception	H5::AttributeIException
// Description
// 		This private function is used in AbstractDs.
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
hid_t Attribute::p_get_type() const
{
   hid_t type_id = H5Aget_type( id );
   if( type_id > 0 )
      return( type_id );
   else
   {
      throw AttributeIException("", "H5Aget_type failed");
   }
}

//--------------------------------------------------------------------------
// Function:	Attribute::getName
///\brief	Gets the name of this attribute, returning its length.
///\param	buf_size  -  IN: Desired length of the name
///\param	attr_name - OUT: Buffer for the name string
///\return	Length of the attribute name
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - Nov, 2001
//--------------------------------------------------------------------------
ssize_t Attribute::getName( size_t buf_size, H5std_string& attr_name ) const
{
   char* name_C = new char[buf_size+1];  // temporary C-string for C API

   // Calls C routine H5Aget_name to get the name of the attribute
   ssize_t name_size = H5Aget_name( id, buf_size, name_C );

   // If H5Aget_name returns a negative value, raise an exception,
   if( name_size < 0 )
   {
      throw AttributeIException("Attribute::getName", "H5Aget_name failed");
   }
   // otherwise, convert the C attribute name and return
   attr_name = name_C;
   delete []name_C;
   return( name_size );
}

//--------------------------------------------------------------------------
// Function:	Attribute::getName
///\brief	This is an overloaded member function, provided for convenience.
///		It differs from the above function in that it returns the
///		attribute's name, not the length.
///\return	Name of the attribute
///\param	buf_size  -  IN: Desired length of the name
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - 2000
//--------------------------------------------------------------------------
H5std_string Attribute::getName( size_t buf_size ) const
{
   H5std_string attr_name;
   ssize_t name_size = getName( buf_size, attr_name );
   return( attr_name );
   // let caller catch exception if any
}

//--------------------------------------------------------------------------
// Function:	Attribute::getName
///\brief	This is an overloaded member function, provided for convenience.
///		It differs from the above functions in that it doesn't take
///		any arguments and returns the attribute's name.
///\return	Name of the attribute
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - May, 2004
//--------------------------------------------------------------------------
H5std_string Attribute::getName() const
{
   // Try with 256 characters for the name first, if the name's length
   // returned is more than that then, read the name again with the
   // appropriate space allocation
   char* name_C = new char[256];  // temporary C-string for C API
   ssize_t name_size = H5Aget_name(id, 255, name_C);

   H5std_string attr_name;
   if (name_size >= 256)
      name_size = getName(name_size, attr_name);

   // otherwise, convert the C attribute name and return
   else
      attr_name = name_C;

   delete []name_C;
   return( attr_name );
}

//--------------------------------------------------------------------------
// Function:	Attribute::getStorageSize
///\brief	Returns the amount of storage size required for this attribute.
///\return	Size of the storage or 0, for no data
///\exception	H5::AttributeIException
// Note:	H5Dget_storage_size returns 0 when there is no data.  This
//		function should have no failure. (from SLU)
// Programmer	Binh-Minh Ribler - Mar, 2005
//--------------------------------------------------------------------------
hsize_t Attribute::getStorageSize() const
{
   hsize_t storage_size = H5Aget_storage_size(id);
   return (storage_size);
}

//--------------------------------------------------------------------------
// Function:	Attribute::close
///\brief	Closes this attribute.
///
///\exception	H5::AttributeIException
// Programmer	Binh-Minh Ribler - Mar 9, 2005
//--------------------------------------------------------------------------
void Attribute::close()
{
    if (p_valid_id(id))
    {
	herr_t ret_value = H5Aclose(id);
	if( ret_value < 0 )
	{
	    throw AttributeIException("Attribute::close", "H5Aclose failed");
	}
	// reset the id because the attribute that it represents is now closed
	id = 0;
    }
}

//--------------------------------------------------------------------------
// Function:	Attribute destructor
///\brief	Properly terminates access to this attribute.
// Programmer	Binh-Minh Ribler - 2000
// Modification
//		- Replaced resetIdComponent() with decRefCount() to use C
//		library ID reference counting mechanism - BMR, Jun 1, 2004
//		- Replaced decRefCount with close() to let the C library
//		handle the reference counting - BMR, Jun 1, 2006
//--------------------------------------------------------------------------
Attribute::~Attribute()
{
    try {
	close();
    }
    catch (Exception close_error) {
	cerr << "Attribute::~Attribute - " << close_error.getDetailMsg() << endl;
    }
}

#ifndef H5_NO_NAMESPACE
} // end namespace
#endif
