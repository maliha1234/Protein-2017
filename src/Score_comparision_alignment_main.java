import java.io.IOException;

/**
 * Created by maliha.sarwat on 5/27/18.
 */
public class Score_comparision_alignment_main {
    public static void main(String[] args) throws IOException {

        int i;

          /*   for(int i=0;i<4;i++)
             {
              Neddle_wnch nw = new Neddle_wnch(fileName[i]);
              nw.Matched(fileName[i]);
             }*/

        String p = null;

        //String[] pairProteins = fileRead(queryProteins[k]);

       // for (i = 0; i < 5; i++) {

          //  System.out.println(pairProteins.length + ":" + queryProteins[k] + pairProteins[i] + "test");
            Neddle_wnch_for_cmparision nw = new Neddle_wnch_for_cmparision();

            nw.Matched("T0866", "Target+Homolog");
           // if (pairProteins[i] == "null") {
              //  break;
                    /* k++;
                     i=0;

                     pairProteins=fileRead(queryProteins[k]);*/
        //    }
       // }
    }

}
