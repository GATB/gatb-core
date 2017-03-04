/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "Ukl_loader.hpp"
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <fstream>

using namespace std;



/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

	UKL_loader::~UKL_loader()
	{
		if(_ukl_set!=0)
			delete _ukl_set;
	}
	
	simpleBitArray * UKL_loader::loadSetFromFile(int kmerSize, int mmerSize)
	{
		const char* datadir = 0;

		datadir = getenv ("DATA_DSK_UKL");
		if(datadir==0)
		{
			printf("Environment variable DATA_DSK_UKL not set\n");
			printf("This should point to the gatb-core/data_ukl folder\n");
			_ukl_set=0;
			return _ukl_set;
		}
		std::string filename(datadir);
		filename += tools::misc::impl::Stringify::format ("ukl_%i_%i.bin", mmerSize,kmerSize);
		std::ifstream is(filename, std::ios::binary);
		_ukl_set = new simpleBitArray();
		_ukl_set->load(is);

		return _ukl_set;
	}
	
	simpleBitArray * UKL_loader::getUKL()
	{
		return _ukl_set;
	}


	
	
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/


