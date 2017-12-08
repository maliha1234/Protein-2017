
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author Maliha
 */
public class Neddle_wnch {

    TemplateMatching.Template t;
    List<String> sequence1;
    List<String> sequence2;
    List<String> mSeqA;
    List<String> mSeqB;
    int[][] mD;
    int mScore;
    String mAlignmentSeqA = "";
    String mAlignmentSeqB = "";
    String file_name;

    public Neddle_wnch(String file_name) {

        // this.t=t;
        this.file_name = file_name;
    }

    public Neddle_wnch() {

        // this.t=t;

    }

    int Matched(String queryFile, String pairFile) throws IOException {
        System.out.println("\n \n \n *********** ");
        String sCurrentLine;


        //String file_path="E:\\4-2\\Thesis 4-2\\TemplateMatching_1\\";
        String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/";

        int strLen, flag = 0, total_match = 0;
        sequence2 = new ArrayList<>();
        sequence1 = new ArrayList<>();
        BufferedReader br = null;
        String final_path;
        try {
            final_path = file_path + queryFile + "/" + pairFile + ".txt";
            System.out.println(final_path);
            // br = new BufferedReader(new FileReader(file_path+this.file_index+".txt"));
            br = new BufferedReader(new FileReader(final_path));

            int j = 0;

            while ((sCurrentLine = br.readLine()) != null) {


                int i, k;
                String tmp = sCurrentLine;


                if (!sCurrentLine.startsWith(">") && flag == 1) {


                    for (i = j, k = 0; k < tmp.length(); i++, k++) {
                        if (tmp.charAt(k) != '\u0000') {
                            sequence1.add(new String(String.valueOf(tmp.charAt(k))));
                            //             System.out.print(tmp.charAt(k));
                        }
                    }

                    j = i - 1;

                } else if (!sCurrentLine.startsWith(">") && flag == 0) {

                    for (i = j, k = 0; k < tmp.length(); i++, k++) {
                        if (tmp.charAt(k) != '\u0000') {
                            sequence2.add(new String(String.valueOf(tmp.charAt(k))));
                            //           System.out.print(tmp.charAt(k));
                        }
                    }

                    j = i - 1;

                } else if (sCurrentLine.startsWith(">") && flag == 0) {

                    flag = 1;
                    continue;


                } else if (sCurrentLine.startsWith(">") && flag == 1) {

                    flag = 0;
                    j = 0;
                    continue;


                }


            }
            //   System.out.print("wrong");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TemplateMatching.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println("\nfinal seq 1  " + sequence1.size());

        for (int i = 0; i < sequence1.size(); i++) {

            System.out.print(sequence1.get(i));
        }

        System.out.println("\nfinal seq 2  " + sequence2.size());
        for (int i = 0; i < sequence2.size(); i++) {
            System.out.print(sequence2.get(i));
        }

        init(sequence1, sequence2);
        process();
        backtrack();
        //   printMatrix();
        System.out.println();
        total_match = printScoreAndAlignments();

        return total_match;

    }


    void init(List<String> seqA, List<String> seqB) {
        mSeqA = seqA;
        mSeqB = seqB;
        // System.out.print("helllooooo"+mSeqA.length);
        mD = new int[mSeqA.size() + 1][mSeqB.size() + 1];
        for (int i = 0; i <= mSeqA.size(); i++) {
            for (int j = 0; j <= mSeqB.size(); j++) {
                if (i == 0) {
                    mD[i][j] = -j;
                } else if (j == 0) {
                    mD[i][j] = -i;
                } else {
                    mD[i][j] = 0;
                }
            }
        }
    }

    void process() {
        for (int i = 1; i <= mSeqA.size(); i++) {
            for (int j = 1; j <= mSeqB.size(); j++) {
                int scoreDiag = mD[i - 1][j - 1] + weight(i, j);
                int scoreLeft = mD[i][j - 1] - 1;
                int scoreUp = mD[i - 1][j] - 1;
                mD[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
            }
        }
    }

    void backtrack() {
        int i = mSeqA.size();
        int j = mSeqB.size();
        mScore = mD[i][j];
        while (i > 0 && j > 0) {
            if (mD[i][j] == mD[i - 1][j - 1] + weight(i, j)) {
                mAlignmentSeqA += mSeqA.get(i - 1);
                mAlignmentSeqB += mSeqB.get(j - 1);
                i--;
                j--;
                continue;
            } else if (mD[i][j] == mD[i][j - 1] - 1) {
                mAlignmentSeqA += "-";
                mAlignmentSeqB += mSeqB.get(j - 1);
                j--;
                continue;
            } else {
                mAlignmentSeqA += mSeqA.get(i - 1);
                mAlignmentSeqB += "-";
                i--;
                continue;
            }
        }
        mAlignmentSeqA = new StringBuffer(mAlignmentSeqA).reverse().toString();
        mAlignmentSeqB = new StringBuffer(mAlignmentSeqB).reverse().toString();
    }

    private int weight(int i, int j) {
        if (sequence1.get(i - 1).equalsIgnoreCase(sequence2.get(j - 1))) {
            return 1;
        } else {
            return -1;
        }
    }

    void printMatrix() {
        System.out.println("D =");
        for (int i = 0; i < mSeqA.size() + 1; i++) {
            for (int j = 0; j < mSeqB.size() + 1; j++) {
                System.out.print(String.format("%4d ", mD[i][j]));
            }
            System.out.println();
        }
        System.out.println();
    }

    int printScoreAndAlignments() {
        System.out.println("Score: " + mScore);
        System.out.println("Sequence A: " + mAlignmentSeqA);
        System.out.println("Sequence B: " + mAlignmentSeqB);
        System.out.println();
        int c = newSeq(mAlignmentSeqA, mAlignmentSeqB);
        System.out.println("\n %%%%done print" + c);

        return c;
    }

    public int newSeq(String A, String B) {

        int l = mAlignmentSeqA.length();
        char[] array = new char[l + 1];
        int i, j = 0, count;
        char[] sequenceA = A.toCharArray();
        char[] sequenceB = B.toCharArray();

        for (i = 0; i < l; i++)

        {

            if (sequenceA[i] != '-' && sequenceB[i] != '-') {
                array[j] = sequenceA[i];
                j++;

            }

        }

        for (i = 0; i < j; i++) {
            System.out.print(array[i]);
        }


        count = j;
        System.out.println("\ndone" + count);
        return count;

    }


}
