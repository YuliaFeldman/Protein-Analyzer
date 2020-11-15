#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>


#define MIN_ARGUMENT 2
#define NUM_OF_COORDINATES 3
#define MAX_ATOMS 20000
#define MAX_LENGTH_OF_INPUT_LINE 80
#define MIN_LENGTH_OF_ATOM_INPUT_LINE 61
#define ATOM_LINE_PREFIX "ATOM  "
#define X_COORDINATE 0
#define Y_COORDINATE 1
#define Z_COORDINATE 2
#define COORD_STR_SIZE 8
#define X_COL_NUM 30
#define Y_COL_NUM 38
#define Z_COL_NUM 46


/* Error messages */
#define NO_INPUT "Usage: AnalyzeProtein <pdb1> <pdb2> ...\n"
#define OPEN_FILE_FAILED "Error opening file: %s\n"
#define N0_ATOMS_WERE_FOUND "Error - 0 atoms were found in the file %s\n"
#define CONVERT_FAILED "Error in coordinate conversion %s!\n"
#define ATOM_LINE_TOO_SHORT "ATOM line is too short %lu characters\n"
#define CLOSE_FILE_FAILED "Error closing file: %s\n"


/* Output messages */
#define OUTPUT_MSG "PDB file %s, %d atoms were read\n"
#define CG_MSG "Cg = %.3f %.3f %.3f\n"
#define RG_MSG "Rg = %.3f\n"
#define DMAX_MSG "Dmax = %.3f\n"


/* Return values of functions */
#define FINISHED_SUCCESSFULLY 0
#define OPEN_INPUT_FILE_FAILED -1
#define INVALID_INPUT -2
#define CLOSE_FAILURE -3
#define FAILURE_READING_LINE -4
#define CONVERSION_FAILURE -5

/* Exit values */
#define INPUT_FAILURE -1
#define PROGRAM_FAILURE -2


/**
 * Function calculates distance between two given atoms
 * @param atom1 first atom
 * @param atom2 second atom
 * @return the distance between two atoms
 */
float calculateDistance(const float atom1[NUM_OF_COORDINATES], const float atom2[NUM_OF_COORDINATES])
{
	return sqrt(pow(atom1[X_COORDINATE] - atom2[X_COORDINATE], 2)
	             + pow(atom1[Y_COORDINATE] - atom2[Y_COORDINATE], 2)
	             + pow(atom1[Z_COORDINATE] - atom2[Z_COORDINATE], 2));
}


/**
 * Function calculates the maximal distance of all the atoms given in the input file.
 * @param numOfAtoms number of atoms read from input file
 * @param coordinatesOfAtoms the coordinates of all the atoms read from file
 * @return maximal distance
 */
float calculateDmax(const int numOfAtoms, const float coordinatesOfAtoms[NUM_OF_COORDINATES][MAX_ATOMS])
{
	float maxDistance = 0.0, currDistance = 0.0;
	for(int i = 0; i < numOfAtoms; i++)
	{
		float atom1[NUM_OF_COORDINATES];
		atom1[X_COORDINATE] = coordinatesOfAtoms[X_COORDINATE][i];
		atom1[Y_COORDINATE] = coordinatesOfAtoms[Y_COORDINATE][i];
		atom1[Z_COORDINATE] = coordinatesOfAtoms[Z_COORDINATE][i];

		for(int j = 0; j < numOfAtoms; j++)
		{
			float atom2[NUM_OF_COORDINATES];
			atom2[X_COORDINATE] = coordinatesOfAtoms[X_COORDINATE][j];
			atom2[Y_COORDINATE] = coordinatesOfAtoms[Y_COORDINATE][j];
			atom2[Z_COORDINATE] = coordinatesOfAtoms[Z_COORDINATE][j];

			currDistance = calculateDistance(atom1, atom2);
			if(currDistance > maxDistance)
			{
				maxDistance = currDistance;
			}

		}
	}

	return maxDistance;
}



/**
 * Function calculates Rg of the atoms given in the input file.
 * @param numOfAtoms number of atoms read from input file
 * @param coordinatesOfAtoms the coordinates of all the atoms read from file
 * @param Cg Cg of the atoms
 * @return Rg value
 */
float calculateRg(const int numOfAtoms, const float coordinatesOfAtoms[NUM_OF_COORDINATES][MAX_ATOMS],
				  const float Cg[NUM_OF_COORDINATES])
{
	float sum = 0;
	for(int i = 0; i < numOfAtoms; i++)
	{
		float atom[NUM_OF_COORDINATES];
		atom[X_COORDINATE] = coordinatesOfAtoms[X_COORDINATE][i];
		atom[Y_COORDINATE] = coordinatesOfAtoms[Y_COORDINATE][i];
		atom[Z_COORDINATE] = coordinatesOfAtoms[Z_COORDINATE][i];
		sum = sum + pow(calculateDistance(Cg, atom), 2);
	}
	return sqrt(sum / numOfAtoms);
}

/**
 * Function calculates Cg of the atoms given in the input file.
 * @param numOfAtoms number of atoms read from input file
 * @param coordinatesOfAtoms the coordinates of all the atoms read from file
 * @param Cg Cg of atoms will be saved in this parameter
 * @return FINISHED_SUCCESSFULLY when calculation finished
 */
float calculateCg(const int numOfAtoms, const float coordinatesOfAtoms[NUM_OF_COORDINATES][MAX_ATOMS],
		          float Cg[NUM_OF_COORDINATES])
{
	for(int i = 0; i < NUM_OF_COORDINATES; i++)
	{
		float sumOfCoordinates = 0;
		for (int j = 0; j < numOfAtoms; j++)
		{
			sumOfCoordinates += coordinatesOfAtoms[i][j];
		}
		Cg[i] = sumOfCoordinates / numOfAtoms;
	}
	return (float)FINISHED_SUCCESSFULLY;
}

/**
 * Function receives a string representation of one of the atoms' values and converts it into a float number.
 * @param str a string representation of one of the atoms' values
 * @param num a numeric representation of the value will be saved in num
 * @return
 */
int convertStrToNumber(char str[COORD_STR_SIZE + 1], float *num)
{

	char *endp;
	errno = 0;
	*num = strtod(str, &endp);
	if( (*num == 0 && (errno != 0 || endp == str)) || strlen(endp) != 0 )
	{
		fprintf(stderr, CONVERT_FAILED, str);
		return CONVERSION_FAILURE;
	}
	return FINISHED_SUCCESSFULLY;

}

/**
 * Function receives an input line taken from the input file and extracts from it the values of a single atom.
 * @param inputLine an input line taken from the input file
 * @param coordinatesOfAtoms a table in which the values of the the new atom will be saved
 * @param numOfAtoms the number of atoms that have been read from file so far
 * @return FAILURE_READING_LINE if failed to receive values from input line.
 * 	        FINISHED_SUCCESSFULLY, otherwise.
 */
int getCoordinatesFromLine(char inputLine[], float coordinatesOfAtoms[NUM_OF_COORDINATES][MAX_ATOMS],
		                   int numOfAtoms)
{
	char strCoordinateX[COORD_STR_SIZE + 1], strCoordinateY[COORD_STR_SIZE + 1], strCoordinateZ[COORD_STR_SIZE + 1];
	float xCoordinate, yCoordinate, zCoordinate;

	memcpy(strCoordinateX, &inputLine[X_COL_NUM], COORD_STR_SIZE);
	memcpy(strCoordinateY, &inputLine[Y_COL_NUM], COORD_STR_SIZE);
	memcpy(strCoordinateZ, &inputLine[Z_COL_NUM], COORD_STR_SIZE);
	strCoordinateX[COORD_STR_SIZE] = '\0';
	strCoordinateY[COORD_STR_SIZE] = '\0';
	strCoordinateZ[COORD_STR_SIZE] = '\0';

	if(convertStrToNumber(strCoordinateX, &xCoordinate) == CONVERSION_FAILURE ||
	   convertStrToNumber(strCoordinateY, &yCoordinate) == CONVERSION_FAILURE ||
	   convertStrToNumber(strCoordinateZ, &zCoordinate) == CONVERSION_FAILURE)
	{
	 	return FAILURE_READING_LINE;
	}

	coordinatesOfAtoms[X_COORDINATE][numOfAtoms] = xCoordinate;
	coordinatesOfAtoms[Y_COORDINATE][numOfAtoms] = yCoordinate;
	coordinatesOfAtoms[Z_COORDINATE][numOfAtoms] = zCoordinate;

	return FINISHED_SUCCESSFULLY;
}

/**
 * Function receives an input file path, opens the file, and reads from it the values of atoms.
 * The values of the atoms will be saved in a table - coordinatesOfAtoms.
 * If a failure will occur in the process of opening a file and reading from it,
 * than a proper massage will be printed and the function will return an error value.
 * @param filePath the path of the input file which contains values of atoms
 * @param coordinatesOfAtoms a table in which the atom values will be saved.
 * @return OPEN_INPUT_FILE_FAILED if failed to oped the input file,
 * 			INVALID_INPUT if failed to read from input file the values of atoms,
 * 			CLOSE_FAILURE if failed to close input file,
 * 			number of atoms read from file, in no failure occurred.
 */
int readFromFile(const char filePath[], float coordinatesOfAtoms[NUM_OF_COORDINATES][MAX_ATOMS])
{
	FILE *inputFile;
	inputFile = fopen(filePath, "r");

	if(inputFile == NULL)
	{
		fprintf(stderr, OPEN_FILE_FAILED, filePath);
		return OPEN_INPUT_FILE_FAILED;
	}

	int numOfAtoms = 0;
	char inputLine[MAX_LENGTH_OF_INPUT_LINE];

	while(fgets(inputLine, MAX_LENGTH_OF_INPUT_LINE, inputFile) != NULL && numOfAtoms <= MAX_ATOMS)
	{
		if(strncmp(inputLine, ATOM_LINE_PREFIX, strlen(ATOM_LINE_PREFIX)) == 0)
		{
			if(strlen(inputLine) < MIN_LENGTH_OF_ATOM_INPUT_LINE)
			{
				fprintf(stderr, ATOM_LINE_TOO_SHORT, strlen(inputLine));
				return INVALID_INPUT;
			}

			if(getCoordinatesFromLine(inputLine, coordinatesOfAtoms, numOfAtoms) == FAILURE_READING_LINE)
			{
				return INVALID_INPUT;
			}
			numOfAtoms ++;
		}
	}

	if(numOfAtoms == 0)
	{
		fprintf(stderr, N0_ATOMS_WERE_FOUND, filePath);
		return INVALID_INPUT;
	}

	if(fclose(inputFile))
	{
		fprintf(stderr, CLOSE_FILE_FAILED, filePath);
		return CLOSE_FAILURE;
	}
	return numOfAtoms;
}

/**
 * The main function of this program receives input files with values of atoms, and prints for each input
 * file its Cg, Rg and Dmax.
 * @param argc number of arguments
 * @param argv arguments
 * @return if program finish running successfully, it will print the calulations and return 0.
 * Otherwise, it will exit the program with the proper error code.
 */
int main(const int argc, const char *argv[])
{

	if(argc < MIN_ARGUMENT)
	{
		fprintf(stderr, NO_INPUT);
		exit(PROGRAM_FAILURE);
	}

	int numOfFiles = 1;

	while(numOfFiles < argc )
	{

		float coordinatesOfAtoms [NUM_OF_COORDINATES][MAX_ATOMS];

		int numOfAtoms = readFromFile(argv[numOfFiles], coordinatesOfAtoms);
		if( numOfAtoms ==  OPEN_INPUT_FILE_FAILED || numOfAtoms == CLOSE_FAILURE)
		{
			exit(PROGRAM_FAILURE);
		}

		if(numOfAtoms == INVALID_INPUT)
		{
			exit(INPUT_FAILURE);
		}

		float Cg[NUM_OF_COORDINATES];
		calculateCg(numOfAtoms, coordinatesOfAtoms, Cg);

		float Rg = calculateRg(numOfAtoms, coordinatesOfAtoms, Cg);
		float Dmax = calculateDmax(numOfAtoms, coordinatesOfAtoms);


		printf(OUTPUT_MSG, argv[numOfFiles], numOfAtoms);
		printf(CG_MSG, Cg[X_COORDINATE], Cg[Y_COORDINATE], Cg[Z_COORDINATE]);
		printf(RG_MSG, Rg);
		printf(DMAX_MSG, Dmax);

		numOfFiles = numOfFiles + 1;

	}

	return FINISHED_SUCCESSFULLY;
}