#ifndef __binaryh__
#define __binaryh__

#define FORMAT_CHECK 1337

struct binaryHeader
{
	int formatCheck;
	int nind; // the number of indicators in the file.
	int nframes;
};

struct elementDescriptor
{
	char resname[256];
	char segid[256];
	int res;
	int chain;
};

#define MAX_BORDERS 15

struct binaryElement
{	
	int type; // 
	int nborders;
	float op;
	int cluster; // which cluster is this assigned to?
	float borderLen[MAX_BORDERS];
	float borderTheta[MAX_BORDERS];
	short borderID[MAX_BORDERS];
};



#endif
