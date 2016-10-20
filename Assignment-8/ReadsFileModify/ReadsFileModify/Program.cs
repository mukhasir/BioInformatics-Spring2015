using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ReadsFileModify
{
    class Program
    {
        static void Main(string[] args)
        {
            int number = 1;
            bool readStat = true;
            TextWriter tsw = new StreamWriter(@"E:\Bio Informatics\Assignment 8\ERR030893-1_Modified.fa");
            foreach (string Read in File.ReadAllLines(@"E:\Bio Informatics\Assignment 8\ERR030893-1.fq"))
            {
                if(readStat)
                {
                    tsw.Write(Read + "\n");
                    if (number == 2)
                    {
                        readStat = false;
                        number = 0;
                    }
                    number++;
                }
                else
                {
                    if (number == 2)
                    {
                        readStat = true;
                        number = 0;
                    }
                    number++;
                }
                //Console.WriteLine(Read);
                //if (i == 9)
                //    break;
                //i++;
            }
            tsw.Close();
        }
    }
}
