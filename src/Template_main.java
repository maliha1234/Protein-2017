import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by maliha.sarwat on 11/9/17.
 */
public class Template_main {

    public static void main(String[] args) throws IOException {


        //*3zqc, 2k9n, 3osf, 1h88, 2dim, 2d9n*//
        String[] fileName = {"3ey7+3ghj", "3ey7+3r4q", "3ghj+3rri", "3rri+4rt5"};
        int k, i;

          /*   for(int i=0;i<4;i++)
             {
              Neddle_wnch nw = new Neddle_wnch(fileName[i]);
              nw.Matched(fileName[i]);
             }*/


        String[] queryProteins = {"T0669", "T0680", "T0696"};
        String p = null;
        Neddle_wnch nw = new Neddle_wnch();

        String[] pairProteins = fileRead(queryProteins[0]);

        for (i = 0, k = 0; i < 5; i++) {

            System.out.println(pairProteins.length + ":" + queryProteins[k] + pairProteins[i] + "test");
            nw.Matched(queryProteins[k], pairProteins[i]);
            if (pairProteins[i] == "null") {
                break;
                    /* k++;
                     i=0;

                     pairProteins=fileRead(queryProteins[k]);*/
            }
        }
    }


    public static String[] fileRead(String queryFolder) throws IOException {


        String[] pairs = new String[10];
        String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/";
        BufferedReader br = null;
        int i = 0;

        try {

            br = new BufferedReader(new FileReader(file_path + queryFolder + "/pairs" + ".txt"));
            String temp = null;

            while ((temp = br.readLine()) != null && !temp.startsWith(">")) {


                pairs[i++] = temp;


            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TemplateMatching.class.getName()).log(Level.SEVERE, null, ex);

        }

        return pairs;
    }


}

