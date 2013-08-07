/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

// Header: osra_stl.h
//
// STL Helpers
//

#include <vector> // std::vector
#include <iostream> // std::ostream

#include "osra_labels.h"
#include "osra_fragments.h"

// Function: operator<<()
//
// Helper template method to print vectors.
namespace std
{
std::ostream& operator<<(std::ostream &os, const letters_t &letter);

std::ostream& operator<<(std::ostream &os, const label_t &label);

std::ostream& operator<<(std::ostream &os, const atom_t &atom);

std::ostream& operator<<(std::ostream &os, const bond_t &bond);

std::ostream& operator<<(std::ostream &os, const fragment_t &fragment);

template <class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
{
  os << '[';
  if (!v.empty())
    {
      typedef typename std::vector<T>::const_iterator const_iterator;

      const_iterator last = v.end();
      std::copy(v.begin(), --last, std::ostream_iterator<T>(os, ", "));
      os << *last;
    }
  os << ']';

  return os;
}
}
