// SIMPLE: Calculates simplicity score from a DNA or protein sequence.
// It uses a moving window and gives scores to mono-, di-, tri- and tetra-
// elements which are found repeated within the window.
// The results are compared to those obtained using randomized sequences.

// Altered Sep 2008 by Simon Greenaway (S.Greenaway@har.mrc.ac.uk) for repeats upto 10 sequences in length

// In this version arguments need to be entered:
// path
// executable_program
// file_name
// seq_nature (y/n/p)
// file_type (1 GenBank,2 EMBL Swissprot,3 crude sequence)
// score_mono (int)
// score_di (int)
// score_tri (int)
// score_tetra (int)
// score_penta (int)    Sep 08 SJG
// score_hexa (int)    Sep 08 SJG
// score_hepta (int)    Sep 08 SJG
// score_octa (int)    Sep 08 SJG
// score_nona (int)    Sep 08 SJG
// score_deca (int)    Sep 08 SJG
// window (int)
// num_random (int)
// randomization method (1,2,3,4)
// stringency (0.9, 0.99, 0.999)
// graphical display (y/n)
// directory (job's name)   Aug 99


#include <fstream.h> // file manipulation
#include <string.h> // string manipulation
#include <stdlib.h> // C++ functions...
#include <math.h> // math operations
#include <iomanip.h> // classes, struct
#include <ctype.h> // toupper (GenBank sequences)
#include <stdio.h> // C library: fprintf..
#include <iostream.h> // cout..


struct sequence{
	long length; // sequence length
   char *seq_id1; // sequence identifications
   char *seq_id2;
   char *seq_id3;
   char *seq_id4;
   char *seq; // sequence
   char* elements; // specific elements to be considered
   int** di_count; // count of di-elements
   int* count; // count of elements
   int** count_three_frames; // count elements in the three different frames
   };

struct param_process{
	int score_mono; // score for each repeated element
   int score_di; // score for each repeated di-element
   int score_tri; // score for each repeated tri-element
   int score_tetra; // score for each repeated tetra-element
   int score_penta; // score for each repeated penta-element
   int score_hexa; // score for each repeated hexa-element
   int score_hepta; // score for each repeated hepta-element
   int score_octa; // score for each repeated octa-element
   int score_nona; // score for each repeated nona-element
   int score_deca; // score for each repeated deca-element
   int random_type; // type of randomization process
   int window; // length of moving window
   long max; // the last reference nucleotide for which to calculate simplicity scores
   int num_random; // number of random sequences
   int FileType; // type of file: GenBank, EMBL, crude sequence
	int NbrOfElements; // AGCT for DNA sequences, etc.
	char seq_nature[2]; // DNA/RNA(n), purine/pyrimidine(y) or protein(p)
	char graphical_display[1]; // option to have a postcript file
   float stringency; // to determine which motifs are associated with high simplicity
   };


struct simple_scores{
	long *ss_test;  // simplicity scores for each position
   float *score_frequency; // frequency of scores
   float SimplicityFactor; // simplicity of sequence
   int MaxScore; // maximum score in sequence
   };


#include "include/files_net.cpp"
#include "include/general_net.cpp"

#define LINE_LEN 81
#define MAX 1000000
#define MAX_DIR 61  // Aug 99


//===================================================================
//
//           Main program
//
//===================================================================

int main (int argc, char *argv[])
{

// check that the correct number of parameters (=20) has been entered  // Aug. 99. Mod Sep08 SJG

if (argc!=20)
	{
   cerr << "ERROR: You have to enter 19 parameters after the name of the executable,\n" <<
           "instead you have entered " << (argc-1) << " !";
   exit(1);
   }


sequence *seq1;          // allocate pointers to struct in heap
seq1 = new sequence;
param_process *userp;
userp = new param_process;
simple_scores *test;
test = new simple_scores;
simple_scores *random;
random = new simple_scores;
if ((!seq1)||(!userp)||(!test)||(!random))
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }

char OutputFile[LINE_LEN];  // Output files
char RandomFile[LINE_LEN];
char OutputMotifs[LINE_LEN];
char OutputScores[LINE_LEN];
char GraphicalDisplay[LINE_LEN];

// Output file SIMPDATA for ROMPLOT, Aug 99
char SimpScores[LINE_LEN];


// directory, defined using the user ID, Aug 99
char directory[MAX_DIR];

char InputFile[LINE_LEN]; // input file


// ========= FUNCTION CALLS =========================

Get_command_arguments(argv, argc -1, InputFile, userp, directory);

/*// Define user settable parameters
User_settable_parameters (userp, InputFile);*/

// Allocate names to the Output Files
Allocate_names (OutputFile, RandomFile, OutputMotifs,  // Aug 99
					 OutputScores, GraphicalDisplay, SimpScores, directory);

// Assign elements to the sequence

   cerr << '\n' << "before assign_elements";
Assign_elements (seq1,userp);
   cerr << '\n' << "after assign_elements";

// Get sequence from file
Get_sequence (seq1, userp, InputFile);

   cerr << '\n' << "after reading from file";

// Determine sequence characteristics (define max, count elements frequency,etc.)
Determine_sequence_characteristics (seq1, userp->NbrOfElements, userp);

// Output for test sequence
Output_test_sequence (seq1, userp, OutputFile, InputFile);

// Calculate simplicity of test sequence
Simplicity_test_sequence (seq1->seq, userp, test, SimpScores);

// Calculations on random sequences, outputs significance of simplicity
Random_sequence_calculations (seq1, userp, random, RandomFile, OutputFile,InputFile,
										test->MaxScore,test->SimplicityFactor);

// Determine significant scores,motifs and segments
Parameters_simplicity (seq1,userp,test,random,OutputScores,OutputMotifs,InputFile);

// Final screen output
cout << "\n\n" << "Process terminated!" <<
		"\n\n" << "Two files have been generated which contain information on the" <<
      " simplicity" <<
      "\nof the sequence:" <<
      '\n' << OutputFile << " -- Sequence characteristics, relative simplicity score, and" <<
      "\nsignificance" <<
      '\n' << OutputMotifs << " -- Motifs and segments involved in the generation of simplicity" <<
      '\n' << "SIMPDATA -- Coordinates for graphical display" << '\n' << OutputScores << " -- Significant scores in test sequence\n";
      //	<< '\n' << RandomFile << " -- Sequence composition of random sequences" <<
      //"\n\n" << "Press a key to leave the screen";


exit(0);
}


