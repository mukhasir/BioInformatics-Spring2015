/*****************************************************************************************************************
 * Date : April 14 2015                           Author : Mukhasir Shah Syed
 * Below program is used for checking the mappability of the 250 gene sequences within the scope of HG19 chr1.
 * In this program we used 3 different read lengths, which are 50, 70 and 100. As reads50.fa and reads70.fa are 
 * provided , we have prepared reads100.fa by creating program for that.
 * Now we have used Bowtie to get the mappability reads for each gene by using reads50.fa, reads70.fa and 
 * reads100.fa files individually.
 * We would create BowTie output (BED format) for each readlength.
 * In this program we will use BowTie output and get the Unique GeneID's and
 * mapped read count and GeneLength to calculate the Mappability Percentile value for each GeneID.
 * Mappability Percentile = ((# of mapped reads)/(GeneLength-readLength+1))*100 formula used for calculation.
 * Above step would be repeated for Bowtie Output of read length 50,70 and 100.
 * After getting mappability values for all 3 read length 50, 70 and 100, we will combile all values into one file.
 * Final file would contains all values for read length 50, 70 & 100 for respective geneId's
 * Output file format below:
 * geneID	readLen50	readLen70	readLen100
   NR_046018	100.0	100.0	100.0
   NR_026818	94.66	99.79	100.0
   NM_001005484	100.0	100.0	100.0
   NR_039983	0.57	1.87	6.46
 * To complie/run the code we have to create a Console Application using Visual Studio 2013 
   and compile/build this application which will create an executable file then run the 
   ".exe" application to execute the code.
 ****************************************************************************************************************/
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneMap_Bowtie_R
{
    class Program
    {
        static void Main(string[] args)
        {
            string Seq = string.Empty;
            string ID = string.Empty;
            //Initialize KeyValuePair List for storing query sequences.
            List<KeyValuePair<string, string>> CHRData = new List<KeyValuePair<string, string>>();
            //Loop through chr1-250.fa file to store all query sequences into list.
            foreach (string Chr250Data in File.ReadLines(@"E:\Bio Informatics\Assignment 7\chr1-250.fa"))
            {
                if (Chr250Data.Contains(">"))
                {
                    if (Seq != string.Empty)
                    {
                        CHRData.Add(new KeyValuePair<string, string>(ID, Seq));
                        Seq = string.Empty;
                    }
                    ID = Chr250Data;
                }
                else
                {
                    Seq += Chr250Data;
                }
            }
            CHRData.Add(new KeyValuePair<string, string>(ID, Seq));
            //Initialize dictionary to store unique gene ID's from BowTie output (BED format) file.
            Dictionary<string, int> GeneData = new Dictionary<string, int>();
            //Initialize dictionary to store Gene ID's and its gene length.
            Dictionary<string, int> GeneDataLength = new Dictionary<string, int>();
            //Variable used to store header of combined file.
            string Header = "geneID" + "\t" + "readLen50" + "\t" + "readLen70" + "\t" + "readLen100" + "\n";
            string Values_Check = string.Empty;
            string Values50 = string.Empty;
            string Values70 = string.Empty;
            string Values100 = string.Empty;
            //Call MappabilityFunction meenthod to get the mappability percentages for respective GeneID's with read length as '50'
            MappabilityFunction(@"E:\Bio Informatics\Assignment 7\BTout-BED-50", GeneData, GeneDataLength, CHRData, ref Values50, 50);
            //Call MappabilityFunction meenthod to get the mappability percentages for respective GeneID's with read length as '70'
            MappabilityFunction(@"E:\Bio Informatics\Assignment 7\BTout-BED-70", GeneData, GeneDataLength, CHRData, ref Values70, 70);
            //Call MappabilityFunction meenthod to get the mappability percentages for respective GeneID's with read length as '100'
            MappabilityFunction(@"E:\Bio Informatics\Assignment 7\BTout-BED-100", GeneData, GeneDataLength, CHRData, ref Values100, 100);
            //Variable and arrays used to combile data for read lengths '50','70' & '100' and display them in one file.
            int i=0;
            string[] Arr50_GeneId = Values50.Split(';');
            string[] GeneIDs = new string[Arr50_GeneId.Length-1];
            string[] Arr50 = new string[Arr50_GeneId.Length - 1];
            for(int j=0;j < (Arr50_GeneId.Length-1);j++)
            {
                GeneIDs[j] = Arr50_GeneId[j].Split('$')[0].Split('.')[3];
                Arr50[j] = Arr50_GeneId[j].Split('$')[1];
            }
            string[] Arr70_GeneId = Values70.Split(';');
            string[] Arr70 = new string[Arr70_GeneId.Length - 1];//Values70.Split(';');
            for (int j = 0; j < (Arr70_GeneId.Length-1); j++)
            {
                Arr70[j] = Arr70_GeneId[j].Split('$')[1];
            }
            string[] Arr100_GeneId = Values100.Split(';');
            string[] Arr100 = new string[Arr100_GeneId.Length - 1];
            for (int j = 0; j < (Arr100_GeneId.Length-1); j++)
            {
                Arr100[j] = Arr100_GeneId[j].Split('$')[1];
            }
            //Loop through GeneId's for which mappability value has been calculated.
            for(int k=0; k<(GeneIDs.Length);k++)
            {
                string GeneName_Value = GeneIDs[k];
                if(i<Arr50.Length)
                    Values_Check += GeneName_Value + "\t" + Arr50[i] + "\t" + Arr70[i] + "\t" + Arr100[i] + "\n";
                i++;
            }
            //Write the GeneID and its mappability value for each read length into a file.
            File.WriteAllText(@"E:\Bio Informatics\Assignment 7\Reult1.txt", Header + Values_Check);
            File.WriteAllText(@"E:\Bio Informatics\Assignment 7\mappability-output-250", Header + Values_Check);
            Console.WriteLine("End");
        }
        /// <summary>
        /// In this method we will take one dictionary to populate GeneId's and mapped read count and another dictionary
        /// to populate GeneID's and their length. Take a list CHRData with data from chr1-250.fa file.
        /// Pass read length value and a string variable which will return the GeneID and their corresponding mappability value.
        /// </summary>
        /// <param name="FilePath"></param>
        /// <param name="GeneData"></param>
        /// <param name="GeneDataLength"></param>
        /// <param name="CHRData"></param>
        /// <param name="Values"></param>
        /// <param name="readLength"></param>
        public static void MappabilityFunction(string FilePath ,Dictionary<string, int> GeneData, Dictionary<string, int> GeneDataLength, List<KeyValuePair<string, string>> CHRData, ref string Values, int readLength)
        {
            //Clear the list as there might be values for previous read length Bowtie file.
            GeneData.Clear();
            //Loop through the Bowtie output (BED format) to get the unique GeneId's and there mapped read count            
            foreach (string GeneID in File.ReadAllLines(FilePath))
            {
                string[] BToutDetails = GeneID.Split('\t');
                string[] GeneNameDetails = BToutDetails[3].Split('.');
                //Build the gene name
                string GeneName = ">" +GeneNameDetails[0] + "." + GeneNameDetails[1] + "." + GeneNameDetails[2] + "." + GeneNameDetails[3] + "." + GeneNameDetails[4];// GeneNameDetails[3];
                //Check if there exists Gene Name in dictionary and add 1 to the value.
                if (GeneData.ContainsKey(GeneName))
                {
                    //Add 1 to existing value for particular GeneID. If Map start and end are in range of Gene ID.
                    if ((Convert.ToInt32(BToutDetails[1]) >= Convert.ToInt32(GeneNameDetails[1])) &&
                    (Convert.ToInt32(BToutDetails[2]) <= Convert.ToInt32(GeneNameDetails[2])))
                    {
                        GeneData[GeneName] = GeneData[GeneName] + 1;
                    }
                }
                else
                {
                    //Add GeneID and 1 as initial value. If Map start and end are in range of Gene ID.
                    if ((Convert.ToInt32(BToutDetails[1]) >= Convert.ToInt32(GeneNameDetails[1])) &&
                    (Convert.ToInt32(BToutDetails[2]) <= Convert.ToInt32(GeneNameDetails[2])))
                    {
                        GeneData.Add(GeneName, 1);
                    }
                }                
            }
            //Clear the dictionary before we put some information into it.
            GeneDataLength.Clear();
            //Loop through CHRData which contains information of Chr1-250.fa file and populate the unique GeneID's and their length.
            foreach (KeyValuePair<string, string> GeneLengthData in CHRData)
            {
                string[] GeneAnnotData_Array = GeneLengthData.Key.Split('.');
                string GeneName = GeneLengthData.Key;//GeneAnnotData_Array[3];
                int GeneLength = 0;
                //If we have duplicate GeneID then we will consider only the last GeneID found and its corresponding length.
                if (GeneDataLength.ContainsKey(GeneName))
                {
                    //Replace new length with the existing Value for particular GeneID which is duplicate.                    
                    GeneDataLength[GeneName] = Convert.ToInt32(GeneAnnotData_Array[2]) - Convert.ToInt32(GeneAnnotData_Array[1]);
                }
                else
                {
                    //Add GeneID and its length.
                    GeneLength = Convert.ToInt32(GeneAnnotData_Array[2]) - Convert.ToInt32(GeneAnnotData_Array[1]);
                    GeneDataLength.Add(GeneName, GeneLength);
                }
            }
            //Linq query to get data present in disctionary.
            var items = from pair in GeneData
                        select pair;
            //Initialize KeyValuePair List for storing data from GeneData dictionary.
            List<KeyValuePair<string, int>> med1 = new List<KeyValuePair<string, int>>(items);
            //Linq query to get data present in disctionary.
            var items1 = from pair in GeneDataLength
                         select pair;
            //Initialize KeyValuePair List for storing data from GeneDataLength dictionary.
            List<KeyValuePair<string, int>> med2 = new List<KeyValuePair<string, int>>(items1);
            //Set read_Length value.
            int read_length = readLength;
            double FinalValue = 0.0;
            bool present = false;
            //Loop through the list which contains GeneID and its length values.
            foreach (KeyValuePair<string, int> GeneLengthValue  in med2)
            {
                //Loop through the list which contains GeneID and its mapped reads.
                foreach (KeyValuePair<string, int> GeneDataValue in med1)
                {
                    //Check if GeneID exists or not
                    if (GeneDataValue.Key == GeneLengthValue.Key)//&& GeneLengthValue.Key == "NR_026818")
                    {
                        //Calculate the mappability value in percentile.
                        double totaltiles = 0.0;
                        //if gene length less than or equal to current read length then set total tiles as 1 other wise calculate.
                        if (GeneLengthValue.Value <= read_length)
                        {
                            totaltiles = 1.0;
                        }
                        else
                        {
                            totaltiles = GeneLengthValue.Value - read_length + 1;
                        }
                        FinalValue = ((Convert.ToDouble(GeneDataValue.Value) / totaltiles)) * 100;                        
                        string vl = Math.Round(FinalValue, 2).ToString();
                        if (vl.IndexOf('.') != -1)
                        {
                            Values += GeneDataValue.Key + "$" + vl + ";";
                        }
                        else
                        {
                            Values += GeneDataValue.Key + "$" + vl + ".0" + ";";
                        }
                        present = true;                        
                    }
                }
                //Based on flag value above, if GeneID doesnt exists then we will set value as '0.0'
                if(!present)
                {
                    Values += GeneLengthValue.Key + "$" + "0.0;";
                    present = false;
                }
                else
                {
                    present = false;
                }
            }
        }
    }
}
