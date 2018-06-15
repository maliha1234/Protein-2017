import java.io.*;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

public class MAKE_SEQ_FROM_PDB {
    public static void main(String[] args) throws IOException {
        readPdb(null);
    }

    static AVG_PDB.Template readPdb(String pdb_file) {


        BufferedReader br = null;
        AVG_PDB.Template t = null;

        try {

            String sCurrentLine;
            String firstLetterArray = "";


            String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/Data_comparison_2018/Score_comparision_May_2018/T0866/T0866.pdb.txt";

            br = new BufferedReader(new FileReader(file_path));
            while ((sCurrentLine = br.readLine()) != null) {
                StringTokenizer st = new StringTokenizer(sCurrentLine);
                ArrayList l = new ArrayList();
                while (st.hasMoreElements()) {

                    l.add(st.nextElement());
                }

                System.out.print(l.get(0));
                try {
                    firstLetterArray += l.get(0).toString().subSequence(0,1);
                    System.out.println(firstLetterArray);


                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            writePDBCordinates(firstLetterArray);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) br.close();

            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        return t;

    }


    static void writePDBCordinates(String firstLetterArray) throws IOException {

        File fout = new File("/Users/maliha.sarwat/Desktop/NE/Thesis_2017/Data_comparison_2018/Score_comparision_May_2018/T0866/SEQUENCE FROM PDB CA/T0866_seq.txt");
        FileOutputStream fos = null;
        try {
            fos = new FileOutputStream(fout);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
        bw.write(firstLetterArray);
        bw.close();

    }

    public static String[] fileRead(String queryFolder) throws IOException {

        System.out.print(queryFolder);

        String[] pairs = new String[10];
        String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/";
        BufferedReader br = null;
        int i = 0;

        try {

            br = new BufferedReader(new FileReader(file_path + queryFolder + "/homologs" + ".txt"));
            String temp = null;

            while ((temp = br.readLine()) != null && !temp.startsWith(">")) {


                pairs[i++] = temp;
                System.out.print(pairs[i]);


            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TemplateMatching.class.getName()).log(Level.SEVERE, null, ex);

        }

        return pairs;
    }
}
