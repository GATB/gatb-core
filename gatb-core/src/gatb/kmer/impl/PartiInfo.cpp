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

#include <gatb/kmer/impl/PartiInfo.hpp>
#include <algorithm>

// We use the required packages
using namespace std;

#define DEBUG(a) // printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Repartitor::computeDistrib (const PartiInfo<5>& extern_pInfo)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    std::vector<ipair> bin_size_vec;
    std::priority_queue< itriple, std::vector<itriple>,compSpaceTriple > pq;

    //sum total bins size
    u_int64_t sumsizes =0;
    for (int ii=0; ii< _nb_minims; ii++)
    {
        // sumsizes +=   extern_pInfo.getNbSuperKmer_per_minim(ii); // _binsize[ii];
        // bin_size_vec.push_back(ipair( extern_pInfo.getNbSuperKmer_per_minim(ii) ,ii));
        sumsizes +=   extern_pInfo.getNbKxmer_per_minim(ii); // _binsize[ii];
        bin_size_vec.push_back(ipair( extern_pInfo.getNbKxmer_per_minim(ii) ,ii));
    }
    u_int64_t mean_size =  sumsizes /  _nbpart;

    DEBUG (("Repartitor : mean size per parti should be :  %lli  (total %lli )\n",mean_size,sumsizes));

    //init space left
    for (int jj = 0; jj < _nbpart; jj++)  {  pq.push (itriple(jj,0,0));  }

    //sort minim bins per size
    std::sort (bin_size_vec.begin (), bin_size_vec.end (), comp_bins);

    DEBUG (("Repartitor : 20 largest estimated bin sizes \n"));
    for (size_t ii=0; ii<20 &&  ii< bin_size_vec.size(); ii++ )
    {
        DEBUG (("binsize [%llu] = %llu \n",bin_size_vec[ii].second,bin_size_vec[ii].first));
    }

    //GC suggestion : put the largest in the emptiest (la plus grosse dans la plus grosse)

    itriple smallest_parti;

    int cur_minim = 0;
    while (cur_minim < _nb_minims)
    {
        //get emptiest parti
        smallest_parti = pq.top(); pq.pop();

        //put largest bin in it
        _repart_table[bin_size_vec[cur_minim].second] = smallest_parti.first;

        //update space used in this bin, push it back in the pq
        smallest_parti.second += bin_size_vec[cur_minim].first;
        smallest_parti.third ++; // how many minimizers are in this bin (just for info)

        pq.push (smallest_parti);

        DEBUG (("Repartitor : affected minim %llu to part %llu  space used %llu  (msize %llu) \n",
            bin_size_vec[cur_minim].second,smallest_parti.first,
            smallest_parti.second , bin_size_vec[cur_minim].first
        ));

        cur_minim++;
    }
}

	
void Repartitor::NaiveMegak (const PartiInfo<5>& extern_pInfo)
{
	
	u_int64_t _num_mm_bins = _nb_minims ;
	
	
	/** We allocate a table whose size is the number of possible minimizers. */
	_repart_table.resize (_nb_minims);
	
	std::vector<ipair> bin_size_vec;
	std::priority_queue< itriple, std::vector<itriple>,compSpaceTriple > pq;
	
	
	//sum total bins size
	u_int64_t sumsizes =0;
	for (int ii=0; ii< _nb_minims; ii++)
	{
		sumsizes +=   extern_pInfo.getNbKxmer_per_minim(ii); // _binsize[ii];
	}
	u_int64_t mean_size =  sumsizes /  _nbpart;
	
	printf ("Repartitor  NaiveMegak : mean size per parti should be :  %lli  (total %lli )\n",mean_size,sumsizes);

	
	
	/////////////////////////////////////////test
	//search for largest superk link
	
	double max_link_val = 0;
	int max_from, max_to;
	
	for (int ii=0; ii<_num_mm_bins; ii++) {
		for (int jj=0; jj<_num_mm_bins; jj++) {
			if(extern_pInfo.getMinimMatrix(ii,jj) > max_link_val)
			{
				max_link_val = extern_pInfo.getMinimMatrix(ii,jj);
				max_from = ii;
				max_to = jj;
			}
		}
	}
	
	printf("frist strongest link from %i to %i : %f \n",max_from,max_to,max_link_val );
	std::set<u_int64_t> dejavu;

	dejavu.insert(max_from);
	dejavu.insert(max_to);

	//init the part with these two
	_repart_table[max_from]= 0;
	_repart_table[max_to]= 0;
//update part size
	u_int64_t acc_size = extern_pInfo.getNbKxmer_per_minim(max_from) + extern_pInfo.getNbKxmer_per_minim(max_to);
	
	
	//remplacer getMinimMatrix par getRelMinimMatrix pour avoir par % des liens totaux sortants
	//put into pq all links of these two
	std::priority_queue< itriple, std::vector<itriple_float>, compLinkfloat > pqm;
	int from = max_from;
	for (int jj=0; jj<_num_mm_bins; jj++) {
		if(extern_pInfo.getMinimMatrix(from,jj) !=0)
		{
 				pqm.push (itriple_float( from  , extern_pInfo.getMinimMatrix(from,jj) ,    jj)); // from, val, to
		}
	}
	from = max_to;
	for (int jj=0; jj<_num_mm_bins; jj++) {
		if(extern_pInfo.getMinimMatrix(from,jj) !=0)
		{
			pqm.push (itriple_float( from  , extern_pInfo.getMinimMatrix(from,jj) ,    jj)); // from, val, to
		}
	}
	
	printf("first acc size %llu \n",acc_size);
	
	//tres approx
	while (acc_size < mean_size)
	{
		itriple_float strongest_link    ;
		strongest_link = pqm.top(); pqm.pop();
		int from = strongest_link.first;
		int to = strongest_link.third;

		int newminim = -1 ;
		if(dejavu.find(from) == dejavu.end()) //si pas deja dedans
		{
			newminim = from;
		}
		else if (dejavu.find(to) == dejavu.end())
		{
			newminim = to;
		}
		else
		{
			//remove this link fro mthe pq
			;
		}
		
		if(newminim!=-1)
		{
			dejavu.insert(newminim);
			_repart_table[newminim] = 0;
			acc_size += extern_pInfo.getNbKxmer_per_minim(newminim) ;
			//add new links
			for (int jj=0; jj<_num_mm_bins; jj++) {
				if(extern_pInfo.getMinimMatrix(newminim,jj) !=0)
				{
 					pqm.push (itriple_float( newminim  , extern_pInfo.getMinimMatrix(newminim,jj) ,    jj)); // from, val, to
				}
			}
			
			printf("added %i  new size %llu link strength %f \n",newminim, acc_size, strongest_link.second);

		}
	}
	
	//on a fini une 'bonne' parti, juste pour test
	//on termine le reste avec ancienne methode
	
	u_int64_t  nb_minims_left = 0 ;

	for (int ii=0; ii< _nb_minims; ii++)
	{
			if (dejavu.find(ii) == dejavu.end()) //que ceux qui restent
			{
				bin_size_vec.push_back(ipair( extern_pInfo.getNbKxmer_per_minim(ii) ,ii));
				nb_minims_left++;
			}
	}
	
	//init space left
	for (int jj = 1; jj < _nbpart; jj++)  {  pq.push (itriple(jj,0,0));  } //tout sauf part 0
	
	//sort minim bins per size
	std::sort (bin_size_vec.begin (), bin_size_vec.end (), comp_bins);

	////////
	
	
	itriple smallest_parti;
	
	int cur_minim = 0;
	while (cur_minim < nb_minims_left)
	{
	//	printf("cur_minim %i \n",cur_minim);
		//get emptiest parti
		smallest_parti = pq.top(); pq.pop();
		
		//put largest bin in it
		_repart_table[bin_size_vec[cur_minim].second] = smallest_parti.first;
		
		//update space used in this bin, push it back in the pq
		smallest_parti.second += bin_size_vec[cur_minim].first;
		smallest_parti.third ++; // how many minimizers are in this bin (just for info)
		
		pq.push (smallest_parti);
		
	/*	printf ("Repartitor : affected minim %llu to part %llu  space used %llu  (msize %llu) \n",
				bin_size_vec[cur_minim].second,smallest_parti.first,
				smallest_parti.second , bin_size_vec[cur_minim].first
				);*/
		
		cur_minim++;
	}

	

	
	
}
	
	
	
// simple version of the code above in the case where we use frequency-based minimizers, and we just want to group minimizers according to their ordering
void Repartitor::justGroupNaive (const PartiInfo<5>& extern_pInfo, std::vector <std::pair<int,int> > &counts)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (int ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = 0; // important to have a consistent repartition for unseen (in the sample approximation) minimizers
    }

    int step = counts.size() / _nbpart;
    
    for (unsigned int i = 0; i < counts.size(); i++)
    {
        _repart_table[counts[i].second] = std::min((int)(i / step), _nbpart-1);
    }
    
}



	

// much more effective version of the function above, using estimation of number of kmers per bucket
void Repartitor::justGroup (const PartiInfo<5>& extern_pInfo, std::vector <std::pair<int,int> > &counts)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (int ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = 0; // important to have a consistent repartition for unseen (in the sample approximation) minimizers
    }

    //sum total count size
    u_int64_t total_counts =0;
    for (unsigned int i = 0; i < counts.size(); i++)
        total_counts += counts[i].first;

    u_int64_t sumsizes =0;
    for (int ii=0; ii< _nb_minims; ii++)
        sumsizes += extern_pInfo.getNbKmer_per_minim(ii);
 
    u_int64_t mean_size = sumsizes / _nbpart;
    
    u_int64_t acc = 0, j = 0;
    for (unsigned int i = 0; i < counts.size(); i++)
    {
	//	printf("mmer freq %i \n",counts[i].first);
        _repart_table[counts[i].second] = j;

        acc += extern_pInfo.getNbKmer_per_minim(counts[i].second);
        if (acc > mean_size)
        {
            acc = 0;
            if (j < _nbpart)
                j++;
        }
    }
}


//cet a cause de cette fonc approx quil y a des parti vides au bout je pense
// lexi case
void Repartitor::justGroupLexi (const PartiInfo<5>& extern_pInfo)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (int ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = 0; // important to have a consistent repartition for unseen (in the sample approximation) minimizers
    }

    u_int64_t sumsizes =0;
    for (int ii=0; ii< _nb_minims; ii++)
        sumsizes += extern_pInfo.getNbKmer_per_minim(ii);
 
    u_int64_t mean_size = sumsizes / _nbpart;
    
    u_int64_t acc = 0, j = 0;
    for (unsigned int i = 0; i < _nb_minims; i++)
    {
		//printf("mmer count %i \n",extern_pInfo.getNbKmer_per_minim(i));
        _repart_table[i] = j;
        acc += extern_pInfo.getNbKmer_per_minim(i);
        if (acc > mean_size)
        {
            acc = 0;
            if (j < _nbpart)
                j++;
        }

    }
    
}




/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Repartitor::load (tools::storage::impl::Group& group)
{
    tools::storage::impl::Storage::istream is (group, "minimRepart");
    is.read ((char*)&_nbpart,     sizeof(_nbpart));
    is.read ((char*)&_mm,         sizeof(_mm));
    is.read ((char*)&_nb_minims,  sizeof(_nb_minims));

    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    is.read ((char*)_repart_table.data(), sizeof(Value) * _nb_minims);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Repartitor::save (tools::storage::impl::Group& group)
{
    tools::storage::impl::Storage::ostream os (group, "minimRepart");
    os.write ((const char*)&_nbpart,                sizeof(_nbpart));
    os.write ((const char*)&_mm,                    sizeof(_mm));
    os.write ((const char*)&_nb_minims,             sizeof(_nb_minims));
    os.write ((const char*)_repart_table.data(),    sizeof(Value) * _nb_minims);
    os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Repartitor::printInfo ()
{
    size_t nbMinimizers = 1 << (_mm*2);
    printf("Repartitor : nbMinimizers=%ld\n", nbMinimizers);
    for(int ii=0; ii<nbMinimizers; ii++ )  {  printf("   table[%i] = %i \n",ii,_repart_table[ii]); }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
