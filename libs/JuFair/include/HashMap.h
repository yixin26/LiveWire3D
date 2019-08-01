#ifndef HASHMAP_H
#define HASHMAP_H

#include <iostream>
using namespace std;
#include <stdio.h>
#include <windows.h>

/* A hash table hashed by two int integers
 *
 */

const int MAX_HASH = 1<<16;
const int HASH_LENGTH = 16;
const int HASH_BIT_1 = 8;
const int HASH_BIT_2 = 8;

struct HashElement
{
/// Key of hash element
	int key[2];
/// Actually content of hash element
	int index ;
/// Link when collision
	HashElement  *nextHash;
};

class HashMap
{
	/// Hash table
	HashElement *table[MAX_HASH];

	/// Create hash key
	int createKey( int k1, int k2 )
	{
		int ind = ((( k1 & ( (1<<HASH_BIT_1) - 1 )) << ( HASH_BIT_2  )) |
					( k2 & ( (1<<HASH_BIT_2) - 1)) ) & ((1<<HASH_LENGTH) - 1);

		return ind;
	}

public:

	/// Constructor
	HashMap ( )
	{
		for (int i = 0; i < MAX_HASH; i ++)
		{
			table[i] = NULL;
		}
	};

	/// Lookup Method
	int findInsert( int k1, int k2, int index )
	{
		/// Create hash key
		int ind = createKey ( k1, k2 );

		/// Find it in the table
		HashElement *p = table[ind];
		while (p)
		{
			if ((p->key[0] == k1) && (p->key[1] == k2))
			{
				return p->index ;
			}
			p = p->nextHash;
		}

		// Not found
		p = new HashElement ;
		p->key[0] = k1;
		p->key[1] = k2;
		p->index = index ;
		p->nextHash = table[ind];
		table[ind] = p;

		return index ;
	};


	int findInsertSort( int k1, int k2, int index )
	{
		/// Create hash key
		if ( k1 > k2 )
		{
			int temp = k1 ;
			k1 = k2 ; 
			k2 = temp ;
		}

		int ind = createKey ( k1, k2 );

		/// Find it in the table
		HashElement *p = table[ind];
		while (p)
		{
			if ((p->key[0] == k1) && (p->key[1] == k2))
			{
				return p->index ;
			}
			p = p->nextHash;
		}

		// Not found
		p = new HashElement ;
		p->key[0] = k1;
		p->key[1] = k2;
		p->index = index ;
		p->nextHash = table[ind];
		table[ind] = p;

		return index ;
	};

	int findInsertSortReplace( int k1, int k2, int index )
	{
		/// Create hash key
		if ( k1 > k2 )
		{
			int temp = k1 ;
			k1 = k2 ; 
			k2 = temp ;
		}

		int ind = createKey ( k1, k2 );

		/// Find it in the table
		HashElement *p = table[ind];
		while (p)
		{
			if ((p->key[0] == k1) && (p->key[1] == k2))
			{
				p->index = index;	//replace it!
				return index;
			}
			p = p->nextHash;
		}

		// Not found
		p = new HashElement ;
		p->key[0] = k1;
		p->key[1] = k2;
		p->index = index ;
		p->nextHash = table[ind];
		table[ind] = p;

		return index ;
	};
	/*void copy(const HashMap& hashmap)
	{
		
	}*/
	void clear()
	{
		HashElement *p, *pp;
		for (int i = 0; i < MAX_HASH; i ++)
		{
			p = table[i];
			while (p)
			{
				pp = p->nextHash;
				delete p;
				p = pp;
			}
		}
		for (int i = 0; i < MAX_HASH; i ++)
		{
			table[i] = NULL;
		}
	}
	// Destruction method
	~HashMap()
	{
		HashElement *p, *pp;

		for (int i = 0; i < MAX_HASH; i ++)
		{
			p = table[i];

			while (p)
			{
				pp = p->nextHash;
				delete p;
				p = pp;
			}
		}
	};

};


#endif