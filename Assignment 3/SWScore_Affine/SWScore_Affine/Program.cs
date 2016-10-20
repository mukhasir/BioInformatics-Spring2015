/***********************************************************************************************
 * Date: 18 Feb 2015.  Author : Mukhasir Shah Syed
 * Using the below code we will calculat ooptimum local alignment score by implementing 
   Smith-Waterman algorithm using affine gap scoring tehnique.
 * Protein data is being used for this segment of code implementation.
 * In this code we will read one query file (EColi-query1.fa) which contains protein data for
   one sequence .
 * We will also read database file (swissprot.fa) which contains multiple protein sequences.
 * Once reading from the query file is completed we implement the algorithm for getting
 * optimum local alignment score using Smith-Waterman Affine gap scoring algorithm.
 * We take 2 rows in this implementation where the first row will be defaulted to zero and in
   the second row we we check the score and place it based on alignment.
 * Here as we are dealing with protein data we will consider GAP_OPEN as 11 and GAP_EXTEND as 1.
 * To calculate score we will run loop for database field data and check it with Single protein
 * query sequence, where SmithWatermanScore_Affine() method will be called.
 * We will used BLOSUM62 matrix(for protein) and Affine gap opening and extension penality
   which are from NCBI BLASTp implementation. 
 * Input : We will provide protein query sequence (EColi-query1.fa) and database fasta file
   (swissprot.fa). Note: swissprot.fa is relatively small protein sequence database which
   contains about 460,670 known protein sequences.
 * Output: We will list file showing optimum local alignment score (Smith-Waterman score) for
   each of the reference protein sequence.
 * Output would be displayed as mentioned below:
 * --- query len = 378 
    >gi|1684788|gb|AAB36530.1| 4-phosphoerythronate dehydrogenase [Escherichia coli str. K-12 substr. W3110]
    MKILVDENMPYARDLFSRLGEVTAVPGRPIPVAQLADADALMVRSVTKVNESLLAGKPIKFVGTATAGTDHVDEAWLKQAGIGFSAAPGCNAIAVVEYVFSSLL
    MLAERDGFSLYDRTVGIVGVGNVGRRLQARLEALGIKTLLCDPPRADRGDEGDFRSLDELVQRADILTFHTPLFKDGPYKTLHLADEKLIRSLKPGAILINACR
    GAVVDNTALLTCLNEGQKLSVVLDVWEGEPELNVELLKKVDIGTSHIAGYTLEGKARGTTQVFEAYSKFIGHEQHVALDTLLPAPEFGRITLHGPLDQPTLKRL
    VHLVYDVRRDDAPLRKVAGIPGEFDKLRKNYLERREWSSLYVICDDASAASLLCKLGFNAVHHPAR

    >gi|745997998|sp|P0DKH9.1|AREP1_ARATH RecName: Ful	(len=40)
    -- SW_score = 20(Query Index = 47, Database Index = 13)
    >gi|745997997|sp|P0DKH8.1|AMP2_FAGES RecName: Full	(len=40)
    -- SW_score = 17(Query Index = 88, Database Index = 13)....and so on.
 * To complie/run the code we have to create a Console Application using Visual Studio 2013 
   and compile/build this application which will create an executable file then run the 
   ".exe" application to execute the code.
 **********************************************************************************************/
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SWScore_Affine
{
    class Program
    {
        /// <summary>
        /// This is the method that gets hit first when this application is run.
        /// In this Main() method we will read single protein sequence from EColi-query1.fa file
        /// We will loop through the database data and place them in key value pair.
        /// We will initialize BLOSUM62 scoring matrix and then loop through the database Key-Value pairs
        /// We will call SmithWatermanScore_Affine() method to get the optimum local alignment score
        /// and store output in a file.
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            //Variables used to store label and data of query protein sequence
            string querylabel = string.Empty;
            string querydata= string.Empty;
            //Read label and data from Protein sequence file and store in above mentioned variables respectiely.
            foreach(string lines in File.ReadAllLines(@"E:\Bio Informatics\Assignment 3\EColi-query1.fa"))
            {
                if(lines.Contains(">"))
                {
                    querylabel = lines;//stores label name
                }
                else
                {
                    querydata += lines;//stores data/bases of protein.
                }
            }
            //Variables used to store label and data of database file.
            string Database_label = string.Empty;
            string Database_Seq = string.Empty;
            //Initialize KeyValuePair List for storing database sequences.
            List<KeyValuePair<string, string>> Database_Data = new List<KeyValuePair<string, string>>();
            //Loop throug database file to store the label and bases respectivel in key-value pair format.
            foreach (string Database in File.ReadLines(@"E:\Bio Informatics\Assignment 3\swissprot.fa"))
            {
                if (Database.Contains(">"))
                {
                    if (Database_Seq != string.Empty)
                    {
                        Database_Data.Add(new KeyValuePair<string, string>(Database_label, Database_Seq));
                        Database_Seq = string.Empty;
                    }
                    Database_label = Database;
                }
                else
                {
                    Database_Seq += Database;
                }
            }
            Database_Data.Add(new KeyValuePair<string, string>(Database_label, Database_Seq));
            //Initializing and declaration of ScoreMatrix (BLOSUM62)
            int[,] ScoreMatrix = new int[25, 25]{
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4}, 
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4}, 
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4},
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4},
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4},
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4},
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4},
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4},
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4},
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4},
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4},
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4},
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4},
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4},
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4},
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4},
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4},
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4},
{-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4},
{-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4},
{-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4},
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1} 
};
            //Store the output format text into variable.
            string FinalText = "--- query len = " + querydata.Length +"  \n" + querylabel + "\n" + querydata + "\n\n";
            //Loop through key-value pair of database sequence and then call SmithWatermanScore_Affine() to get the score and indexes 
            //of database and query.
            foreach (KeyValuePair<string, string> DatabaseValue in Database_Data)
            {
                //As we are processing protein data the value for gao opening is '11' and gap extension is '1'
                int gap_open = 11;
                int gap_ext = 1;
                //Forming string to show the output in required format.
                //Store the database label and bases length into string for displaying it in output.
                string header = DatabaseValue.Key.Substring(0, 50) + "\t" + "(len=" + DatabaseValue.Value.Length + ")\n";
                //Call SmithWatermanScore_Affine() method to get the score and indexes
                string vallue = SmithWatermanScore_Affine(DatabaseValue.Value, DatabaseValue.Value.Length, querydata, querydata.Length, gap_open, gap_ext, ScoreMatrix);
                string ScoreText = "-- SW_score = " + vallue;
                FinalText += header + ScoreText;
            }
            //Store the output/complete result for sequences into Result.txt file.
            File.WriteAllText(@"E:\Bio Informatics\Assignment 3\Result.txt", FinalText);
        }
        /// <summary>
        /// We pass database sequence, query sequence, lengths of database sequence and query sequence, gap open, gap extend
        /// and score matrix values.
        /// In this we will initialize two seperate 1D arrays. Where first will be set value as 'zero'.
        /// and set would be set -(gap_open) as value in it.
        /// Loop through the query sequence character and database characters and calculate optimum local alignment score
        /// using Smith-Waterman affine scoring algorithm.
        /// </summary>
        /// <param name="dbSeq"></param>
        /// <param name="dbSeq_len"></param>
        /// <param name="querySeq"></param>
        /// <param name="querySeq_len"></param>
        /// <param name="gap_open"></param>
        /// <param name="gap_ext"></param>
        /// <param name="ScoreMatrix"></param>
        /// <returns></returns>
        public static string SmithWatermanScore_Affine(string dbSeq, int dbSeq_len,string querySeq, int querySeq_len, int gap_open, int gap_ext, int[,] ScoreMatrix)
        {
            //intialize two 1 dimensional arrays of size each database sequence+ 1  that is passed.
            int[] nogap = new int[dbSeq.Length + 1];
            int[] b_gap = new int[dbSeq.Length + 1];
            //Variables used in the code.
            int i =0 , j=0;
            int last_nogap, prev_nogap;
            int score = 0;
            int GlobalScore = 0;
            int queryIndex = 0;
            int dbIndex = 0;
            //loop through b_gap array and set the value as -(gap_open)
            for (int counter = 0; counter <= dbSeq.Length; counter++)
            {
                b_gap[counter] = -(gap_open);
            }
            //loop through length of query sequence
            for(i=0;i<querySeq_len;++i)
            {
                //declare value for last_nogap, prev_nogap as 0
                int a_gap;
                last_nogap = prev_nogap = 0;
                a_gap = -(gap_open);
                //loop through length of database sequence
                for(j=0;j<dbSeq_len;++j)
                {
                    //Smith-Waterman Affine Scoring algorithm
                    a_gap = Math.Max((last_nogap-gap_open-gap_ext),(a_gap-gap_ext));
                    b_gap[j] = Math.Max((nogap[j]-gap_open-gap_ext),(b_gap[j]-gap_ext));
                    //We will call BLOSUM62Value() method to get the gap penality value from scorematrix.
                    last_nogap = Math.Max((prev_nogap + BLOSUM62Value(ScoreMatrix,dbSeq[j],querySeq[i])), 0);
                    last_nogap = Math.Max(last_nogap,a_gap);
                    last_nogap = Math.Max(last_nogap,b_gap[j]);
                    prev_nogap = nogap[j];
                    nogap[j] = last_nogap;
                    score = Math.Max(score,last_nogap);
                    //Check if the new score is greater than or equal to older one then place it in GlobalScore and store the indexes.
                    if (score >= GlobalScore)
                    {
                        GlobalScore = score;
                        queryIndex = i;
                        dbIndex = j;
                    }
                }
            }
            //Send bacl score , query index and database index.
            string ScoreTextValue = GlobalScore.ToString() + "(Query Index = " + queryIndex + ", Database Index = " + dbIndex + ")\n";
            return ScoreTextValue;
        }
        /// <summary>
        /// This method takes scorematrix, database sequence character and query sequence character
        /// and get penality value and return it for optimum local alignment score calculation.
        /// </summary>
        /// <param name="ScoreMatrix"></param>
        /// <param name="dbSeqChar"></param>
        /// <param name="queryChar"></param>
        /// <returns></returns>
        public static int BLOSUM62Value(int[,] ScoreMatrix, char dbSeqChar, char queryChar)
        {
            int BLOSUMscore = 0;
            //Intialize and declare the character set for protein sequences.(BLOSUM62)
            char[] Alphabet = new char[25] { 'A', 'R','N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*'};
            int dbIndex =- 1 , queryIndex = -1;
            //Loop through the array and get the index of character in DB and query
            for (int k = 0; k < Alphabet.Length; k++)
            {
                if (Alphabet[k] == dbSeqChar)
                {
                    dbIndex = k;
                }
                if (Alphabet[k] == queryChar)
                {
                    queryIndex = k;
                }
            }
            //Check if the character is not presnt in Alphabet array then assigh index= length - 1
            if (dbIndex == -1)
            {
                dbIndex = Alphabet.Length - 1;
            }
            if (queryIndex==-1)
            {
                queryIndex = Alphabet.Length - 1;
            }
            //Get the penality value.
            BLOSUMscore = ScoreMatrix[dbIndex, queryIndex];
            //Return the penality score
            return BLOSUMscore;
        }
    }
}
