import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by maliha.sarwat on 12/1/17.
 */
public class HomologsRead {
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
