/*****************************************************************************************************
 * Date : Mar-25-2015         Author: Mukhasir Shah Syed
 * Program written will be used to find median string and motif from given set of DNA sequences and 
   the length of the motif by implementing L-mer tree way. We would report the best 5 median strings 
   and their corresponding motif consensus strings and positions. 
 * INPUT:
 *      A set of DNA sequences (arbitary length each) from "HMP-617.fa" file and motif length (L-mer)
 *      of size 6.
 * OUTPUT:
 *      We will display 5 best median strings(with total distance), corresonding motif consensus 
 *      strings(with conseensus scores), and motif positions (for each motif string).
 * Firstly, we will read all DNA sequences from "HMP-617.fa" and place them into a list of key value pair.
 * Secondly, we would traverse through L-mer tree structure by starting at "AAAAAA".
 * We will loop through Vertexes in L-mer tree in each DNA sequence and find out hamming 
   score and based on conditions if total distance for the vertex is greater than bestdistance of 6-mer 
 * then ByPass the vertex, else go to next vertex. 
 * If length of vertex is 6 then get the calculation to be done for that 6-mer leaf and then go to next vertex.
 * Place 6-mers and their distance into a dictionary and sort them based on Total distance in ascending order.
 * After sorting is completed, we can find best 5 median strings for which motif consensus string, 
   motif consensus score and motif position in each DNA sequences with less hamming score.
 * TotalDistance() method is used to get all median strings and their hamming score for each substring
   of size "6" in each DNA sequence provided. We will pass Vertex value and DNA sequence, Score would 
   be update with addition of new score to old score.
 * FindMotifPositions() method will find one motif position in each DNA sequence with low hamming score
   and place that value and its position in sequence into list. We would pass median value from top 5
 * median values selected, DNA sequences and List of key value pair.
 * FindConsensus() method is used to find motif consensus score and motif consensus score.
 * ByPass() method is used to traverse through next vertex at same level.
 * NextVertex() method is used to traverse through next vertex starting at next level. If level is at 6
 * then go to next vertex at same level. If vertex is last child then go to next vertex a level up.
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
using System.Windows.Forms;

namespace MotifFind_Tree
{
    class Program
    {
        static void Main(string[] args)
        {
            string HMP_Label = string.Empty;
            string HMP_Seq = string.Empty;
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
            
            //variables used further
            string st = string.Empty;
            int bestDistance = 0;
            int optimisticDistance = 0;
            int TotalDist = 0;
            string prefix = "AAAAAA";
            int x = 6;
            bool Flag = true;
            st = prefix;
            int TreeValue_Level = prefix.Length;
            string bestWord = string.Empty;
            //Initialize a List to store visited vertexes of 6-mers
            List<string> VisitedNodes = new List<string>();
            //Initialize a dictionary to store Median Values processed.
            Dictionary<string, int> MedianDistanceData = new Dictionary<string, int>();
            //Loop through to check each vertex applicable based on conditions
            while (x > 0)
            {
                //if level of vertex is less than 6 then go through else go into else block.
                if (x < 6)
                {
                    //Loop through each DNA sequence by passing vertex from L-mer tree and find out total distance.
                    foreach (KeyValuePair<string, string> HMPMOtif_Seq in HMP_Data)
                    {
                        optimisticDistance += TotalDistance(prefix, HMPMOtif_Seq.Value);
                    }
                    st = prefix;
                    //Check if total distance is greater than best distance then go through and bypass vertex.
                    //else goto next vertex method to find next vertex in L-mer tree.
                    if (optimisticDistance > bestDistance)
                    {
                        //This method will by pass through next vertex at same level or other level based on vertex position and level.
                        ByPass(ref st, ref TreeValue_Level, 6, 4);
                        prefix = st;
                        x = TreeValue_Level;
                    }
                    else
                    {
                        //This method will go to next vertex in L-mer tree structure.
                        NextVertex(ref st, ref TreeValue_Level, 6, 4);
                        prefix = st;
                        x = TreeValue_Level;
                    }
                    optimisticDistance = 0;
                }
                else
                {
                    //Loop through each DNA sequence by passing vertex from L-mer tree and find out total distance.
                    foreach (KeyValuePair<string, string> HMPMOtif_Seq in HMP_Data)
                    {
                        TotalDist += TotalDistance(prefix, HMPMOtif_Seq.Value);
                    }
                    VisitedNodes.Add(prefix);
                    //Add median string and its total distance for respective median string into dictionary.
                    MedianDistanceData.Add(prefix, TotalDist);
                    //Check if Total distance is less than best distance then set new total distance as best distance.
                    //set the median as best word
                    if (TotalDist < bestDistance)
                    {
                        bestDistance = TotalDist;
                        bestWord = st;
                    }
                    if (Flag)
                    {
                        bestDistance = TotalDist;
                        bestWord = st;
                        Flag = false;                        
                    }
                    //After setting best distance and best word now move to next vertex in L-mer tree structure.
                    NextVertex(ref st, ref TreeValue_Level, 6, 4);
                    prefix = st;
                    x = TreeValue_Level;
                    TotalDist = 0;
                }
            }
            string FinalVisitedNodes = string.Join("\n", VisitedNodes.ToArray());
            File.WriteAllText(@"E:\Bio Informatics\Assignment 6\TotalVisited_6mers.txt", FinalVisitedNodes);
            Console.WriteLine(MedianDistanceData.Count);
            Console.WriteLine(VisitedNodes.Count);
            //Linq query to sort data present in disctionary based on value(which is score) in ascending order.
            var items = from pair in MedianDistanceData
                        orderby pair.Value ascending
                        select pair;
            ////Initialize KeyValuePair List for storing median values with their score.
            List<KeyValuePair<string, int>> med = new List<KeyValuePair<string, int>>(items);
            string median_string = string.Empty;
            string motifps = string.Empty;
            string cValueFind = string.Empty;
            string motifCon_String = string.Empty;
            string FinalValue = string.Empty;
            //Initialize KeyValuePair List for storing motif data and their positions.
            List<KeyValuePair<string, int>> MotifPositionData = null;
            //Loop through 5 best median values.
            for (int counter = 0; counter < 5; counter++)
            {
                motifps = "";
                MotifPositionData = new List<KeyValuePair<string, int>>();
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
                motifps = "motif positions/string s=(s1,s2,..,st-1,st): " + "\n\t" + motifps.Substring(0, motifps.Length - 1);
                //This method is used to find out motif consensus string and motif consensus score by passing 
                //all motif strings found above. 
                motifCon_String = FindConsensusScore(cValueFind);
                //Concatinate all the string to display the final outpout format.
                FinalValue += median_string + "\n" + motifCon_String + "\n" + motifps + "\n\n";
            }
            File.WriteAllText(@"E:\Bio Informatics\Assignment 6\ConsensusResult.txt", FinalValue);
        }
        /// <summary>
        /// This method is used to ByPass through Vertexes in L-mer tree structure.
        /// </summary>
        /// <param name="CurrentLeaf"></param>
        /// <param name="level"></param>
        /// <param name="length"></param>
        /// <param name="maxValue"></param>
        public static void ByPass(ref string CurrentLeaf, ref int level, int length, int maxValue)
        {
            //Convert string into character array
            char[] CurLeaf_int = StringNumber(CurrentLeaf, "N").ToArray();
            bool FinalValue = true;
            int a = 0;
            //Loop through character array and increment the position to next number.
            for (a = (level-1); a >= 0; a--)
            {
                //Check if Value at particular location in character array less than max value (4) then 
                //increment the value to next value.
                if (Convert.ToInt16(CurLeaf_int[a].ToString()) < maxValue)
                {
                    CurLeaf_int[a] = Convert.ToChar((Convert.ToInt16(CurLeaf_int[a].ToString()) + 1).ToString());
                    FinalValue = false;
                    break;
                }
            }
            //Change character array into string.
            string Final = String.Join("", CurLeaf_int.ToArray());
            //Convert string in number to string in alphabets.
            Final = StringNumber(Final, "S");
            CurrentLeaf = Final.Substring(0, a + 1);
            level = CurrentLeaf.Length;
            if (FinalValue)
            {
                a = Final.Length-1;
                CurrentLeaf = Final.Substring(0, a + 1);
                level = 0;
            }
        }
        /// <summary>
        /// This method is used to traverse through next vertex in L-mer tree structurel
        /// </summary>
        /// <param name="s"></param>
        /// <param name="level"></param>
        /// <param name="length"></param>
        /// <param name="maxvalue"></param>
        public static void NextVertex(ref string s,ref int level, int length, int maxvalue)
        {
            //Variables to be used in this method
            string NextV = string.Empty;
            int counter = 0;
            bool FinalValue = true;
            //Convert string with aplhabets into number character array.
            char[] CurLeaf_int = StringNumber(s, "N").ToArray();
            //Check if level is less than length then go to next vertex starting at next level..
            if(level<length)
            {
                //Loop through character array 
                for(int i=0;i<CurLeaf_int.Length;i++)
                {
                    NextV += CurLeaf_int[i];
                }
                NextV += "1";
                level = NextV.Length;
                s = StringNumber(NextV, "S");
            }
            else
            {
                //loop through character array.
                for (counter = length-1; counter >= 0; counter--)
                {
                    //Check if value at certain index is less than max value (4)
                    //then add "1" to value.
                    if (Convert.ToInt16(CurLeaf_int[counter].ToString()) < maxvalue)
                    {
                        CurLeaf_int[counter] = Convert.ToChar((Convert.ToInt16(CurLeaf_int[counter].ToString()) + 1).ToString());
                        FinalValue = false;
                        level = counter + 1;
                        break;
                    }
                }
                //Change character array to String.
                string Final = String.Join("", CurLeaf_int.ToArray());
                //Convert string in number to string in alphabets.
                Final = StringNumber(Final, "S");
                //Set the values in s and level to be returned.
                s = Final.Substring(0, counter + 1);
                level = s.Length;
                if (FinalValue)
                {
                    counter = Final.Length-1;
                    s = Final.Substring(0, counter + 1);
                    level = 0;
                }
            }
        }
        /// <summary>
        /// This method is used to convert string into number format and viceversa.
        /// </summary>
        /// <param name="value"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public static string StringNumber(string value, string type)
        {
            string IntValue = string.Empty;
            //Check if type is "N" or "S". If "N" then we will prepare a string with alphabets with numbers 
            //for string with corresponding alphhabet values.
            //If "S" then convert Number string into alphabet format.
            if (type == "N")
            {
                //Convert string to character array
                char[] value_array = value.ToArray();
                foreach (char c in value_array)
                {
                    switch (c)
                    {
                        case 'A':
                            IntValue += "1";
                            break;
                        case 'C':
                            IntValue += "2";
                            break;
                        case 'G':
                            IntValue += "3";
                            break;
                        case 'T':
                            IntValue += "4";
                            break;
                    }
                }
            }
            else if(type=="S")
            {
                //Convert string to character array
                char[] value_array = value.ToArray();
                foreach (char c in value_array)
                {
                    switch (c)
                    {
                        case '1':
                            IntValue += "A";
                            break;
                        case '2':
                            IntValue += "C";
                            break;
                        case '3':
                            IntValue += "G";
                            break;
                        case '4':
                            IntValue += "T";
                            break;
                    }
                }
            }
            //return the respective value requested.
            return IntValue;
        }
        /// <summary>
        /// This method is used to find out median string and its score. We will add new score to the existing score if that
        /// median has already been placed into dictionary.
        /// </summary>
        /// <param name="medianValue"></param>
        /// <param name="SequenceValue"></param>
        public static int TotalDistance(string medianValue, string SequenceValue)
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
            //Return score for each sequence for corresponding vertex passed.
            return substring_position;
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
            string consensus_string = string.Empty;
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
        public static void FindMotifPositions(string medianValue, string SequenceValue, ref List<KeyValuePair<string, int>> MotifPositions)
        {
            //Variables used in this method.
            bool setValue = true;
            string v1 = string.Empty;
            int v2 = 0;
            int k = 0;
            int score = 0;
            int substring_position = 0;
            int z = 0;
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
                //If hamming score value is zero/same value then we use the below code set first string as value to below
                //variables.
                if (maxscore == substring_position)
                {
                    if (z == 0)
                    {
                        v1 = tempValue;
                        v2 = k;
                    }
                    z++;
                }
            }
            //Add motif strings and their corresponding score into list.
            MotifPositions.Add(new KeyValuePair<string, int>(v1, v2));
        }        
    }
}
