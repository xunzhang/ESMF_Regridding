// $Id: Sacado_CacheFad_DFad.hpp,v 1.1 2007/08/07 20:45:56 dneckels Exp $ 
// $Source: /cvsroot/esmf/esmf/src/Infrastructure/Mesh/include/sacado/Sacado_CacheFad_DFad.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_DFAD_HPP
#define SACADO_CACHEFAD_DFAD_HPP

#include "Sacado_CacheFad_GeneralFadExpr.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"
#include "Sacado_CacheFad_DFadTraits.hpp"

namespace Sacado {

  namespace CacheFad {

    /*! 
     * \brief Forward-mode AD class using dynamic memory allocation and
     * caching expression templates.
     */
    /*!
     * This is a user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::CacheFad::GeneralFad.  It is similar
     * to Sacado::Fad::DFad, except it uses the caching expression templates
     * that cache the results of val() calculations for later dx() 
     * calculations.  Preliminary testing indicates this is more efficient for 
     * deep expressions with large numbers of derivative components, but is
     * less efficient otherwise.
     *
     * The class is templated on two types, \c ValueT and \c ScalarT.  Type
     * \c ValueT is the type for values the derivative class holds, while
     * type \c ScalarT is the type of basic scalars in the code being
     * differentiated (usually \c doubles).  When computing first derivatives, 
     * these two types are generally the same,  However when computing
     * higher derivatives, \c ValueT may be DFad<double> while \c ScalarT will
     * still be \c double.  Usually \c ScalarT does not need to be explicitly
     * specified since it can be deduced from \c ValueT through the template
     * metafunction ScalarValueType.
     */
    template <typename ValueT, 
	      typename ScalarT = typename ScalarValueType<ValueT>::type >
    class DFad : public Expr< GeneralFad<ValueT,
					 Fad::DynamicStorage<ValueT> > > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DFad() : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DFad(const ValueT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty
       */
      DFad(const ScalarT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DFad(const int sz, const ValueT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DFad(const int sz, const int i, const ValueT & x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(sz,i,x) {}

      //! Copy constructor
      DFad(const DFad& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DFad(const Expr<S>& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //@}

      //! Destructor
      ~DFad() {}

      //! Assignment operator with constant right-hand-side
      DFad& operator=(const ValueT& val) {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(val);
	return *this;
      }

      //! Assignment operator with constant right-hand-side
      DFad& operator=(const ScalarT& val) {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(ValueT(val));
	return *this;
      }

      //! Assignment operator with DFad right-hand-side
      DFad& operator=(const DFad& x) {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(static_cast<const GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> DFad& operator=(const Expr<S>& x) 
      {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(x);
	return *this;
      }
	
    }; // class DFad<ScalarT,ValueT>

    /*! 
     * \brief Forward-mode AD class using dynamic memory allocation and
     * caching expression templates.
     */
    /*!
     * This is the specialization of DFad<ValueT,ScalarT> for when
     * \c ValueT and \c ScalarT are the same type.  It removes an extra
     * constructor that would be duplicated in this case.
     */
    template <typename ValueT>
    class DFad<ValueT,ValueT> : 
      public Expr< GeneralFad<ValueT,Fad::DynamicStorage<ValueT> > > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DFad() : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DFad(const ValueT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DFad(const int sz, const ValueT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DFad(const int sz, const int i, const ValueT& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(sz,i,x) {}

      //! Copy constructor
      DFad(const DFad& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DFad(const Expr<S>& x) : 
	Expr< GeneralFad< ValueT,Fad::DynamicStorage<ValueT> > >(x) {}

      //@}

      //! Destructor
      ~DFad() {}

      //! Assignment operator with constant right-hand-side
      DFad& operator=(const ValueT& val) {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(val);
	return *this;
      }

      //! Assignment operator with DFad right-hand-side
      DFad& operator=(const DFad& x) {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(static_cast<const GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> DFad& operator=(const Expr<S>& x) 
      {
	GeneralFad< ValueT,Fad::DynamicStorage<ValueT> >::operator=(x);
	return *this;
      }

    }; // class DFad<ValueT>

  } // namespace CacheFad

} // namespace Sacado

#endif // SACADO_CACHEFAD_DFAD_HPP
