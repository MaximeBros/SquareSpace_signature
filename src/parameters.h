#ifndef PARAMETERS_H
#define PARAMETERS_H


/*
#define P 1451
#define M 5
#define R 2
#define D1 0
#define D2 1
#define D3 5
#define C1 1449
#define C2 20
#define C3 1
*/

/*
#define P 11
#define M 5
#define R 2
#define D1 0
#define D2 2
#define D3 5
#define C1 9
#define C2 10
#define C3 1
*/

/*
#define P 2097169
#define M 167
#define R 3
#define D1 
#define D2 
#define D3 
#define C1 
#define C2 
#define C3
*/

/*
#define P 65521
#define M 131
#define R 4
#define D1 
#define D2 
#define D3 
#define C1 
#define C2 
#define C3
*/

/*
#define P 79
#define M 251
#define R 5
#define D1 
#define D2 
#define D3 
#define C1 
#define C2 
#define C3
*/


#define P 1451
#define M 131
#define R 4
#define D1 0
#define D2 1
#define D3 131
#define C1 60
#define C2 2
#define C3 1


#define R2 R * R
#define T (R * (R+1)) / 2
#define T2 ((R * R + 1) * R * R) / 2
#define COMMIT_SIZE 16 * 11 * M * T
#define PUB_KEY_SIZE (11 * M * T) / 8 + 1
#define MESSAGE_SIZE 16
#define HASH_SIZE COMMIT_SIZE + PUB_KEY_SIZE + MESSAGE_SIZE

#endif
