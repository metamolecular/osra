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

#include "osra_stl.h"

// Function: operator<<()
//
// Helper template method to print various object types.
namespace std
{
std::ostream& operator<<(std::ostream &os, const letters_t &letter)
{
  os << "{letter char:" << letter.a << " x:" << letter.x << " y:" << letter.y << " r:" << letter.r << " free:" << letter.free << '}';

  return os;
}

std::ostream& operator<<(std::ostream &os, const label_t &label)
{
  os << "{label s:" << label.a << " box:" << label.x1 << "x" << label.y1 << "-" << label.x2 << "x" << label.y2 << '}';

  return os;
}

std::ostream& operator<<(std::ostream &os, const atom_t &atom)
{
  os << "{atom label:" << atom.label << " x:" << atom.x << " y:" << atom.y << " n:" << atom.n << " anum:" << atom.anum << '}';

  return os;
}

std::ostream& operator<<(std::ostream &os, const bond_t &bond)
{
  os << "{bond a:" << bond.a << " b:" << bond.b << " type:" << bond.type << " exists:" << bond.exists
     << " arom:" << bond.arom << " hash:" << bond.hash << " wedge:" << bond.wedge << " up:" << bond.up << " down:" << bond.down<< '}';

  return os;
}

std::ostream& operator<<(std::ostream &os, const fragment_t &fragment)
{
  os << "{fragment " << fragment.x1 << "x" << fragment.y1 << "-" << fragment.x2 << "x" << fragment.y2 << " atoms.size:" << fragment.atom.size() << '}';

  return os;
}
}
