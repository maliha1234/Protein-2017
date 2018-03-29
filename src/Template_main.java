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

    public static String[] queryProteins = {"T0669", "T0680", "T0696", "T0859","T0921", "T0866", "T0868", "T0877", "T0888", "T0891", "T0895", "T0908", "T0915", "T0935", "T0948"};
    public static int k = 14;
    public static void main(String[] args) throws IOException {

        int i;

          /*   for(int i=0;i<4;i++)
             {
              Neddle_wnch nw = new Neddle_wnch(fileName[i]);
              nw.Matched(fileName[i]);
             }*/

        String p = null;

        String[] pairProteins = fileRead(queryProteins[k]);

        for (i = 0; i < 5; i++) {

            System.out.println(pairProteins.length + ":" + queryProteins[k] + pairProteins[i] + "test");
            Neddle_wnch nw = new Neddle_wnch();

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
        String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/";
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

