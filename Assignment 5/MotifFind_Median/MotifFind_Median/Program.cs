/*****************************************************************************************************
 * Date : Mar-14-2015         Author: Mukhasir Shah Syed
 * Program written will be used to find median string and motif from given set of DNA sequences and 
   the length of the motif. We would report the best 5 median strings and their corresponding motif
   consensus strings and positions. 
 * INPUT:
 *      A set of DNA sequences (arbitary length each) from "HMP-617.fa" file and motif length (L-mer)
 *      of size 6.
 * OUTPUT:
 *      We will display 5 best median strings(with total distance), corresonding motif consensus 
 *      strings(with conseensus scores), and motif positions (for each motif string).
 * Firstly, we would make list which contains all different sequences of median strings(4^L[L=6]).6-mer's
 * Secondly, we will read all DNA sequences from "HMP-617.fa" and place them into a list of key value pair.
 * As now we have got all median strings and all sequences from file, we will now loop through DNA 
   sequences for each Median string. 
 * We will loop through the each substring of length as 6 in each DNA sequence and find out hamming 
   score and place median string value and score into a Dictonary.
 * After placing the value into dictionary we will sort disctionary based on the score which is placed
   in value of key value pair in ascending order.
 * After sorting is completed, we can find best 5 median strings for which motif consensus string, 
   motif consensus score and motif position in each DNA sequences with less hamming score.
 * TotalDistance() method is used to get all median strings and their hamming score for each substring
   of size "6" in each DNA sequence provided. We will pass median value and DNA sequence and dictionary 
   to store Median Values and their corresponding score. Score would be update with addition 
   of new score to old score. We check for duplicacy of Median Value and keep only unique value.
 * FindMotifPositions() method will find one motif position in each DNA sequence with low hamming score
   and place that value and its position in sequence into list. We would pass median value from top 5
 * median values selected, DNA sequences and List of key value pair.
 * FindConsensus() method is used to find motif consensus score and motif consensus score.
 * To complie/run the code we have to create a Console Application using Visual Studio 2013 
   and compile/build this application which will create an executable file then run the 
   ".exe" application to execute the code.
 ****************************************************************************************************/
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MotifFind_Median
{
    class Program
    {
        static void Main(string[] args)
        {
            //Variables usd in the program further.
            string motifCon_String = string.Empty;
            string FinalValue = string.Empty;
            string HMP_Label = string.Empty;
            string HMP_Seq = string.Empty;
            //As we consider 6-mers to use as median, we will select length of string as 6.
            int length = Convert.ToInt32(Math.Pow(4,6));
            //List to place 4096 median  values.
            List<string> medians = new List<string>();
            //Linq method which will create random strings.
            Random random = new Random();
            bool add = true;
            int Size = 6;
            string input = "ACGT";
            StringBuilder builder = null;
            char ch;
            int i =0;
            int j = 0;
            //Loop through and create 6-Mer median values.
            for (j = 0; j < length; j++)
            {
                add = true;
                builder = new StringBuilder();
                for (i = 0; i < Size; i++)
                {
                    ch = input[random.Next(0, input.Length)];
                    builder.Append(ch);
                }
                //We will check for duplicacy, if present then we will not add.
                foreach(string v in medians)
                {
                    if(v == builder.ToString())
                    {
                        j = j - 1;
                        add = false;
                        break;
                    }
                }
                //If not duplicate then we wil add the value to list.
                if(add)
                {
                    medians.Add(builder.ToString());
                }
            }
            //Sort the list such that it starts from AAAAAA
            medians.Sort();
            //Initialize KeyValuePair List for storing database sequences.
            List<KeyValuePair<string, string>> HMP_Data = new List<KeyValuePair<string, string>>();
            //Loop throug database file to store the label and bases respectively in key-value pair format.
            foreach (string Database in File.ReadLines(@"E:\Bio Informatics\Assignment 5\HMP-617.fa"))
            {
                if (Database.Contains(">"))
                {
                    if (HMP_Seq != string.Empty)
                    {
                        HMP_Data.Add(new KeyValuePair<string, string>(HMP_Label, HMP_Seq));
                        HMP_Seq = string.Empty;
                    }
                    HMP_Label = Database;
                }
                else
                {
                    HMP_Seq += Database;
                }
            }
            HMP_Data.Add(new KeyValuePair<string, string>(HMP_Label, HMP_Seq));
            //Initialize Dictionary for storing Median Values along with their respective hamming score.
            Dictionary<string, int> MedianDistanceData = new Dictionary<string, int>();
            //Loop thru each DNA sequence.
            foreach(KeyValuePair<string,string> HMPMOtif_Seq in HMP_Data)
            {
                //if (d == 1)
                //{
                //    Console.WriteLine(HMPMOtif_Seq.Value.Substring(2034, 6));
                //    break;
                //}
                //d++;
                //if (HMPMOtif_Seq.Key.Contains(">ECSE_P6-0003;") || HMPMOtif_Seq.Key.Contains(">ECSE_P5-0006;") || HMPMOtif_Seq.Key.Contains(">ECSE_P5-0007;")
                //            || HMPMOtif_Seq.Key.Contains(">ECSE_P3-0027;") || HMPMOtif_Seq.Key.Contains(">ECSE_P3-0030;"))
                //{
                    //Loop thru each median value of size 6.
                    foreach (string V in medians)
                    {
                        //This method is used to traverse through each median with each sequence and each substring of
                        //size "6" in the sequence. This will provide us a list with median value and their corresponding score.
                        //string V = "TTGCT";//"CGTAA";
                        TotalDistance(V, HMPMOtif_Seq.Value, ref MedianDistanceData);
                    }
                //}
            }
            //Linq query to sort data present in disctionary based on value(which is score) in ascending order.
            var items = from pair in MedianDistanceData
                        orderby pair.Value ascending
                        select pair;
            //Initialize KeyValuePair List for storing median values with their score.
            List<KeyValuePair<string, int>> med = new List<KeyValuePair<string, int>>(items);
            //Variables used to store and display values.
            string median_string = string.Empty;
            string motifps = string.Empty;
            string cValueFind = string.Empty;
            //Initialize KeyValuePair List for storing motif data and their positions.
            List<KeyValuePair<string, int>> MotifPositionData = null;
            //Loop through 5 best median values.
            for (int counter = 0; counter < 5; counter++)
            {
                motifps = "";
                MotifPositionData =  new List<KeyValuePair<string, int>>();
                //Loop through DNA Sequences
                foreach (KeyValuePair<string, string> b in HMP_Data)
                {
                    //Build a string value with median string and its total distance
                    median_string = "median string: " + med[counter].Key + " (total_distance = " + med[counter].Value + ")";
                    //This method is used to find motif string and its position in DNA sequence which have low 
                    //hamming distance.
                    FindMotifPositions(med[counter].Key, b.Value, ref MotifPositionData);
                }
                cValueFind = string.Empty;
                //Loop through all motif strings for a particular median string and form a string to display.
                foreach (KeyValuePair<string, int> motifSeq in MotifPositionData)
                {
                    //Concatinate motif strings into one string so that it could be passed to FindConsensus to find 
                    //motif consensus string and motif consensus score.
                    cValueFind += motifSeq.Key + ";";
                    //Concatinate position and string to diplay in output.
                    motifps += motifSeq.Value + "(" + motifSeq.Key + "), ";
                    motifps = motifps.Trim();
                }
                //Build format and concatinate final value into the string to be displayed.
                motifps = "motif positions/string s=(s1,s2,..,st-1,st): " + "\n\t" + motifps.Substring(0,motifps.Length-1);
                //This method is used to find out motif consensus string and motif consensus score by passing 
                //all motif strings found above. 
                motifCon_String = FindConsensusScore(cValueFind);
                //Concatinate all the string to display the final outpout format.
                FinalValue += median_string + "\n" + motifCon_String + "\n" + motifps + "\n\n"; 
            }
            //Write final result into the text file.
            //File.WriteAllText(@"E:\Bio Informatics\Assignment 5\Result.txt", FinalValue);
        }
        /// <summary>
        /// This method is used to find the Motif Consensus Score and Motif Consensus String from all
        /// motif strings retrieved for each median string.
        /// </summary>
        /// <param name="sequences"></param>
        /// <returns></returns>
        public static string FindConsensusScore(string sequences)
        {
            string modifyString = string.Empty;
            int score = 0;
            string consensus_string= string.Empty;
            //We will split string based on splitter ";".
            string[] motifstrings = sequences.Split(';');
            //Loop through each character present in each string.
            for (int j = 0; j < 6; j++)
            {
                //Loop through each motif string to get character at "i" position to form a string.
                for (int i = 0; i < motifstrings.Length; i++)
                {
                    if (motifstrings[i] != "")
                    {
                        modifyString += motifstrings[i][j].ToString();
                    }
                }
                //We will add occurence of character in each column and get numeric number and add the sum of all value to get consensus score.
                score += modifyString.GroupBy(x => x).OrderByDescending(x => x.Count()).First().Count();
                //We will get maximun occured character and concatinate all values to get consensus string.
                consensus_string += modifyString.GroupBy(x => x).OrderByDescending(x => x.Count()).First().Key.ToString();
                modifyString = "";
            }
            //Form a value in string to be displayed in output.
            string mcs = "motif consensus string: " + consensus_string + " (consensus_score = " + score.ToString() + ")";
            //return final string.
            return mcs;
        }
        /// <summary>
        /// This method is used to find motif positions and motif strings in for the 5 best median strings.
        /// </summary>
        /// <param name="medianValue"></param>
        /// <param name="SequenceValue"></param>
        /// <param name="MotifPositions"></param>
        public static void FindMotifPositions(string medianValue, string SequenceValue, ref List<KeyValuePair<string,int>> MotifPositions)
        {
            //Variables used in this method.
            bool setValue = true;
            string v1 = string.Empty;
            int v2 = 0;
            int k = 0;
            int score = 0;
            int substring_position = 0;
            //Loop through DNA sequence and get substring of sequence size "6" and find out hamming score.
            for (k = 0; k <= (SequenceValue.Length - medianValue.Length); k++)
            {
                score = 0;
                //Get substring of size of l(L=6)
                string tempValue = SequenceValue.Substring(k, medianValue.Length);
                //Convert median and substring into character array for comparing each character at same position.
                char[] median = medianValue.ToArray();
                char[] tmpvalue = tempValue.ToArray();
                //Loop through character array to find score of mismatch.
                for (int l = 0; l < medianValue.Length; l++)
                {
                    if (median[l] != tmpvalue[l])
                    {
                        score = score + 1;
                    }
                }
                int maxscore = score;
                //Below written if statement will set minimun hamming score.
                if (maxscore != substring_position)
                {
                    if (setValue)
                    {
                        substring_position = maxscore;
                        setValue = false;
                    }
                }
                //We will check if the new hamming score value is less than the old hamming score value then set new value
                //into variables so that that can be update into list.
                if (maxscore < substring_position)
                {
                    substring_position = maxscore;
                    v1 = tempValue;
                    v2 = k;
                }
                //If hamming score value is zero/same value then we use the below code set last string as value to below
                //variables.
                if(maxscore==substring_position)
                {
                    v1 = tempValue;
                    v2 = k;
                }
            }
            //Add motif strings and their corresponding score into list.
            MotifPositions.Add(new KeyValuePair<string, int>(v1, v2));
        }
        /// <summary>
        /// This method is used to find out median string and its score. We will add new score to the existing score if that
        /// median has already been placed into dictionary.
        /// </summary>
        /// <param name="medianValue"></param>
        /// <param name="SequenceValue"></param>
        /// <param name="searchValue"></param>
        public static void TotalDistance(string medianValue, string SequenceValue, ref Dictionary<string, int> searchValue)
        {
            //Variables to be used in this method
            bool setValue = true;
            int k = 0;
            int score = 0;
            int substring_position = 0;
            //Loop through DNA sequence till the last substring string has size "6" and find out hamming score.
            for (k = 0; k <= (SequenceValue.Length - medianValue.Length); k++)
            {
                score = 0;
                //Get substring of size of l(L=6)
                string tempValue = SequenceValue.Substring(k, medianValue.Length);
                //Convert median and substring into character array for comparing each character at same position.
                char[] median = medianValue.ToArray();
                char[] tmpvalue = tempValue.ToArray();
                //Loop through character array to find score of mismatch.
                for (int l = 0; l < medianValue.Length; l++)
                {
                    if (median[l] != tmpvalue[l])
                    {
                        score = score + 1;
                    }
                }
                int maxscore = score;
                //We will check if the new hamming score value is less than the old hamming score value then set new value
                //into variables so that that can be update into list.
                if (maxscore != substring_position)
                {
                    if (setValue)
                    {
                        substring_position = maxscore;
                        setValue = false;
                    }
                }
                //get the score for each substring.
                if (maxscore <= substring_position)
                {
                    substring_position = maxscore;
                    setValue = false;
                }
            }
            //We will check for duplicacy and if that string is not present then we will add to that list.
            //if we find that then we will add that new score to existing score.
            if (searchValue.ContainsKey(medianValue))
            {
                //Add new score to existing score for particular sustring from sequence.
                searchValue[medianValue] = searchValue[medianValue] + substring_position;
            }
            else
            {
                //Add substring from sequences and its score.
                searchValue.Add(medianValue, substring_position);
            }
        }
    }
}
