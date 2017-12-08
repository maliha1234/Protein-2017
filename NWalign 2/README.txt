//-----------------------------------------------------------------
// 
// Description: Implementation of Needleman-Wunsch global alignment.
// This code is written by Ren-Xiang Yan in China Agricultural University and is originally based on 
// the fortran implementation from Dr. Yang Zhang (http://zhanglab.ccmb.med.umich.edu/NW-align/).
// Last update is in 2010/08/14. 
//
//  Usage:
//      java -jar NWAlign.jar F1.fasta F2.fasta  (align two sequences in fasta file)
//		java -jar NWAlign.jar F1.pdb F2.pdb    1 (align two sequences in PDB file)
//		java -jar NWAlign.jar F.fasta F.pdb  2 (align sequences 1 in fasta and 1 in pdb)
//		java -jar NWAlign.jar GKDGL EVADELVSE    3 (align two sequences in plain text)
//		java -jar NWAlign.jar GKDGL F.fasta  4 (align sequences 1 in text and 1 in fasta)
//		java -jar NWAlign.jar GKDGL F.pdb    5 (align sequences 1 in text and 1 in pdb)
//  
//   Note: You also could complied the code by yourself. 
//         Decompress the NWAlign.jar file and you can get the source code in the NWAlign folder.
//   The program can be compiled by 
//              javac NWAlign.java
//   Then you could use the program by the following commands:
//		java NWAlign F1.fasta F2.fasta  (align two sequences in fasta file)
//		java NWAlign F1.pdb F2.pdb    1 (align two sequences in PDB file)
//		java NWAlign file1.fasta file2.pdb  2 (align sequences 1 in fasta and 1 in pdb)
//		java NWAlign GKDGL EVADELVSE    3 (align two sequences in plain text)
//		java NWAlign GKDGL F.fasta  4 (align sequences 1 in text and 1 in fasta)
//		java NWAlign GKDGL F.pdb    5 (align sequences 1 in text and 1 in pdb)            
//-----------------x-------------------x-------------------x----------