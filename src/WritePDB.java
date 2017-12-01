import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by maliha.sarwat on 11/18/17.
 */
public class WritePDB {
    public static void main(String[] args) throws IOException {

        //  String[] TO680_proteins = new String[]{"1YO7.pdb", "3U8V.pdb", "4FM3.pdb", "4k0D.pdb"};
        //  String[] TO669_proteins = new String[]{"1PQX.pdb", "2K1H.pdb", "2LTL.pdb", "2LTM.pdb"};
//        String[] TO669_proteins = new String[]{"1PQX.pdb", "2K1H.pdb", "2LTL.pdb", "2LTM.pdb"};


        //   for(int i=0; i<4; i++)
        // readPdb(TO669_proteins[i]);

        //readPdb("2K1H.pdb");
        String[] homologs = fileRead(Template_main.queryProteins[Template_main.k]);

        for (int i = 0; i < 4; i++) {
            System.out.print(homologs[i]);
            readPdb(homologs[i]);
        }

    }

    static AVG_PDB.Template readPdb(String pdb_file) {


        BufferedReader br = null;
        AVG_PDB.Template t = null;
        //  t = new Template();

        try {

            String sCurrentLine, lines;

            int seqid, i = 0, length;
            int[] matched_align;
            float x, y, z;
            String type, symbol;
            // String file_path="/Users/maliha.sarwat/Desktop/Thesis_2017/T0669/";
            //      String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/T0669/";
            String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/" + Template_main.queryProteins[Template_main.k] + "/";

            br = new BufferedReader(new FileReader(file_path + pdb_file));
            List<AVG_PDB.Atom> atom = new ArrayList<AVG_PDB.Atom>();
            //atom = new List();
            List<AVG_PDB.AminoAcid> aminoacid = new ArrayList<AVG_PDB.AminoAcid>();
            while ((sCurrentLine = br.readLine()) != null) {

                if (sCurrentLine.startsWith("TER")) {

                    break;
                }

                if (sCurrentLine.startsWith("ATOM") || sCurrentLine.startsWith("HETATM")) {


                    StringTokenizer st = new StringTokenizer(sCurrentLine);
                    ArrayList l = new ArrayList();
                    //System.out.println("---- Split by space ------");
                    while (st.hasMoreElements()) {

                        l.add(st.nextElement());

                    }

                    if (l.get(2).equals("CA")) {
                        //   System.out.println(sCurrentLine);
                        //System.out.println(i);
                        x = Float.valueOf(l.get(6).toString());
                        y = Float.valueOf(l.get(7).toString());
                        z = Float.valueOf(l.get(8).toString());
                        type = "CA";
                        symbol = l.get(3).toString();
                        seqid = Integer.parseInt(l.get(5).toString());
                        atom.add(i, new AVG_PDB.Atom(x, y, z, type));
                        //  atom.get(i).printAtoms(atom.get(i));
                        aminoacid.add(i, new AVG_PDB.AminoAcid(atom.get(i), symbol, seqid));
                        //      aminoacid.get(i).printAminoAcids(aminoacid.get(i));
                        //      System.out.println(i);
                        i++;
                        //  System.out.print(x+"____"+y+"_____"+z+"\n");
                    }
                }
            }
            length = i;
            //    Template t[file_index]=new Template(length);


            t = new AVG_PDB.Template(pdb_file, aminoacid);
            writePDBCordinates(t);
            double d[][] = AVG_PDB.Template.distance(t.amino_acid, t.amino_acid.size());
            t = new AVG_PDB.Template(t.amino_acid, d);

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


    static void writePDBCordinates(AVG_PDB.Template template) throws IOException {

        //  String fileName = template.template_name.split(".")[0];
        File fout = new File("/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/" + Template_main.queryProteins[Template_main.k] + "/Coordinates/" + template.template_name + ".txt");
        FileOutputStream fos = null;
        try {
            fos = new FileOutputStream(fout);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));

        for (int i = 0; i < template.amino_acid.size(); i++) {
            bw.write(template.amino_acid.get(i).symbol
                    + " " + template.amino_acid.get(i).atom.x + " " + template.amino_acid.get(i).atom.y + " " + template.amino_acid.get(i).atom.z);
            bw.newLine();
        }

        bw.close();




      /*  BufferedWriter bw = null;
        FileWriter fw = null;

        try {

            String content = "This is the content to write into file\n";

            fw = new FileWriter(FILENAME);
            bw = new BufferedWriter(fw);
            bw.write(content);

            System.out.println("Done");

        } catch (IOException e) {

            e.printStackTrace();

        } finally {

            try {

                if (bw != null)
                    bw.close();

                if (fw != null)
                    fw.close();

            } catch (IOException ex) {

                ex.printStackTrace();

            }

        }*/
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
