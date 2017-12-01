
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.lang.Math;
import java.util.List;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author hp user
 */
public class AVG_PDB {


    public static Template input_Template;


    static class Atom {

        double x, y, z;
        String type;

        Atom(double x, double y, double z, String type) {

            this.x = x;
            this.y = y;
            this.z = z;
            this.type = type;

        }

    }

    static class AminoAcid {

        Atom atom;
        String symbol;
        int seqid;

        void printAminoAcids(AminoAcid aa) {

            System.out.println("Symbol :" + aa.symbol + " seqid :" + aa.seqid);
        }

        AminoAcid(Atom atom, String symbol, int seqid) {

            this.atom = atom;
            this.symbol = symbol;
            this.seqid = seqid;
        }


    }

    static class Structure {

        AminoAcid amino_acid;

        Structure(AminoAcid amino_acid) {

            this.amino_acid = amino_acid;

        }


    }

    static class Template {

        float score;
        String template_name;
        int index;
        int[] matched_align;
        float weight;
        char[] sequence;
        double[][] d;
        // Structure structure;
        List<AminoAcid> amino_acid;

        /**
         * this method creates intra-distance matrix for a template
         ***/
        static double[][] distance(List<AminoAcid> a, int length) {

            int i = 0, j = 0;
            double[][] d = new double[length][length];

            for (i = 0; i < length; i++)
                for (j = i + 1; j < length; j++) {

                    double temp1, temp2, temp3, sum = 0.00;
                    temp1 = (a.get(i).atom.x - a.get(j).atom.x);
                    temp2 = Math.pow(temp1, 2);
                    sum += temp2;
                    temp1 = (a.get(i).atom.y - a.get(j).atom.y);
                    temp2 = Math.pow(temp1, 2);
                    sum += temp2;
                    temp1 = (a.get(i).atom.z - a.get(j).atom.z);
                    temp2 = Math.pow(temp1, 2);
                    sum += temp2;
                    d[i][j] = Math.sqrt(sum);


                }
         
    /*   for(i=0; i<length; i++)
        for(j=i+1; j<length; j++){
            
            System.out.println(" "+d[i][j]+" ");
        }*/
            return d;
        }


        Template() {

            this.amino_acid = null;
            this.d = null;
            this.index = 0;
            this.matched_align = null;
            this.score = 0;
            this.sequence = null;
            this.template_name = null;
            this.weight = 0;
        }

        /***
         * this constructor is for building input template i.e. target protein
         ***/
        Template(float score, String template_name, double[][] d, List<AminoAcid> amino_acid) {

            this.score = score;
            this.template_name = template_name;
            //this.matched_align=align;
            this.d = d;
            //  this.weight=weight;
            this.sequence = sequence;
            //this.structure=structure;
            this.amino_acid = amino_acid;

        }

        Template(List<AminoAcid> amino_acid, double[][] d) {

            this.amino_acid = amino_acid;
            this.d = d;

        }

        /***
         * this constructor is for building 10 matched templates
         ***/
        Template(String template_name, List<AminoAcid> amino_acid) {

            //  this.score=score;
            this.template_name = template_name;
            // this.matched_align=align;
            this.d = d;
            //   this.weight=weight;
            this.sequence = sequence;
            //this.structure=structure;
            this.amino_acid = amino_acid;  //readPdb theke aminoacid array ekhane ashbe, for one template, we get an array of CA aminoacids


        }

    }

    static Template readPdb(String pdb_file) {


        BufferedReader br = null;
        Template t = null;
        //  t = new Template();

        try {

            String sCurrentLine, lines;

            int seqid, i = 0, length;
            int[] matched_align;
            float x, y, z;
            String type, symbol;
            //     String file_path="/Users/maliha.sarwat/Desktop/NE/Thesis_2017/T0669/";
            ///   String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/T0680/";

            String file_path = "/Users/maliha.sarwat/Desktop/NE/Thesis_2017/NewPDB/" + Template_main.queryProteins[Template_main.k] + "/";

            br = new BufferedReader(new FileReader(file_path + pdb_file));
            List<Atom> atom = new ArrayList<Atom>();
            //atom = new List();
            List<AminoAcid> aminoacid = new ArrayList<AminoAcid>();
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
                        atom.add(i, new Atom(x, y, z, type));
                        //  atom.get(i).printAtoms(atom.get(i));
                        aminoacid.add(i, new AminoAcid(atom.get(i), symbol, seqid));
                        //      aminoacid.get(i).printAminoAcids(aminoacid.get(i));
                        //      System.out.println(i);
                        i++;
                        //  System.out.print(x+"____"+y+"_____"+z+"\n");
                    }
                }
            }
            length = i;
            //    Template t[file_index]=new Template(length);


            t = new Template(pdb_file, aminoacid);
            double d[][] = Template.distance(t.amino_acid, t.amino_acid.size());
            t = new Template(t.amino_acid, d);

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


    static void writePDBCordinates(Template template) throws IOException {

        File fout = new File("/Users/maliha.sarwat/Desktop/NE/Thesis_2017/" + template.template_name + ".txt");
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

    public static void main(String[] args) throws IOException {
        
       
     /* Template template1=new Template();
        Template template2=new Template();
        
        template1=readPdb("2K1H.pdb");
        System.out.print("***********\n");  
        
        template2=readPdb("1PQX.pdb");
        System.out.print("***********\n"); */

        //*3zqc, 2k9n, 3osf, 1h88, 2dim, 2d9n*//
        Template[] templates = new Template[5];

      /*  templates[0]=readPdb("2LTL.pdb"); // 669
        System.out.print("***********\n");  
        
        templates[1]=readPdb("2LTM.pdb");
        System.out.print("***********\n");
        
        templates[2]=readPdb("1PQX.pdb");
        System.out.print("***********\n");  
        
        templates[3]=readPdb("2K1H.pdb");
        System.out.print("***********\n"); // 669*/

        templates[0] = readPdb("3U8V.pdb"); // 680
        System.out.print("***********\n");

        templates[1] = readPdb("4FM3.pdb");
        System.out.print("***********\n");

        templates[2] = readPdb("4K0D.pdb");
        System.out.print("***********\n");

        templates[3] = readPdb("1YO7.pdb");
        System.out.print("***********\n");
        Template concensus = new Template();
        // concensus=readPdb("3RRI.pdb");

        concensus = templates[0];
        Template c2 = new Template();
        c2 = find_avg(templates[0], templates[1], concensus);

        double d[][] = Template.distance(c2.amino_acid, c2.amino_acid.size());
        c2 = new Template(c2.amino_acid, d);

        double seq_drmsd, struct_drmsd;


        Neddle_wnch nw = new Neddle_wnch();
        int matched_count = nw.Matched("T0680", "3u8v+4fm3");
        // int matched_count= nw.Matched("T0669","2k1h+1pqx");
        //    seq_drmsd=mg.calc_initial_drmsd(templates, concensus, matched_count);
        System.out.println("initial drmsd:");
        System.out.println("matched count:" + matched_count);
        ModelGeneration mg = new ModelGeneration();
        seq_drmsd = mg.calc_initial_drmsd(templates, concensus, matched_count);

        int count = 0;
        System.out.println("drmsd after structural avg:");

        c2 = concensus;
        for (int i = 0; i < 10; i++) {

            c2 = find_avg(templates[0], templates[1], c2);


            double d2[][] = Template.distance(c2.amino_acid, c2.amino_acid.size());
            c2 = new Template(c2.amino_acid, d2);
            struct_drmsd = mg.calc_initial_drmsd(templates, c2, matched_count);

            if ((seq_drmsd - struct_drmsd) == 0.00) break;
            else count++;
            // System.out.println(struct_drmsd);

        }
        System.out.println(count);
    /*   for(int i=0;i<concensus.amino_acid.size();i++)
        {
            System.out.println("original");
            printAtoms(concensus.amino_acid.get(i).atom);
            
            System.out.println("avg");
            printAtoms(c2.amino_acid.get(i).atom);
       
        }*/

    }

    public static void printAtoms(Atom atom) {

        //  System.out.println("X :"+atom.x+" Y :"+atom.y+" Z :"+atom.z+" type :"+atom.type);
    }


    private static Template rotation_maty(char axis, Template template, Double angle) {

        Template tc = new Template();
        Template[] t = new Template[3];

        t[0] = template;
//            t[1]=template2;
//            t[2]=tc1;

        int length = t[0].amino_acid.size();
//            for(int j=1; j<3; j++){
//
//            if(length > t[j].amino_acid.size()){
//
//                length = t[j].amino_acid.size();
//            }
//
//           }

        Double xnew, ynew, znew;

        for (int i = 0; i < length; i++) {
            //if(template1.amino_acid.get(i)==null || template2.amino_acid.get(i)==null)
            // break;

            if (angle > 0.00) {
                template.amino_acid.get(i).atom.x = ((template.amino_acid.get(i).atom.x) * Math.cos(angle) + (template.amino_acid.get(i).atom.z) * Math.sin(angle));

                template.amino_acid.get(i).atom.z = ((template.amino_acid.get(i).atom.z) * Math.cos(angle) - (template.amino_acid.get(i).atom.x) * Math.sin(angle));

            } else {

                template.amino_acid.get(i).atom.x = ((template.amino_acid.get(i).atom.x) * Math.cos(angle) + (template.amino_acid.get(i).atom.y) * Math.sin(angle));

                template.amino_acid.get(i).atom.y = ((template.amino_acid.get(i).atom.y) * Math.cos(angle) - (template.amino_acid.get(i).atom.x) * Math.sin(angle));


            }

//            tc1.amino_acid.get(i).atom.y=(template1.amino_acid.get(i).atom.y+template2.amino_acid.get(i).atom.y)/2;
//            tc1.amino_acid.get(i).atom.z=(template1.amino_acid.get(i).atom.z+template2.amino_acid.get(i).atom.z)/2;
//        

        }

        return template;
    }


    private static Template rotation_matx(char axis, Template template, Double angle) {

        Template tc = new Template();
        Template[] t = new Template[3];

        t[0] = template;
//            t[1]=template2;
//            t[2]=tc1;

        int length = t[0].amino_acid.size();
//            for(int j=1; j<3; j++){
//
//            if(length > t[j].amino_acid.size()){
//
//                length = t[j].amino_acid.size();
//            }
//
//           }

        Double xnew, ynew, znew;

        for (int i = 0; i < length; i++) {
            //if(template1.amino_acid.get(i)==null || template2.amino_acid.get(i)==null)
            // break;

            if (angle > 0.00) {
                template.amino_acid.get(i).atom.x = ((template.amino_acid.get(i).atom.x) * Math.cos(angle) - (template.amino_acid.get(i).atom.y) * Math.sin(angle));

                template.amino_acid.get(i).atom.y = ((template.amino_acid.get(i).atom.x) * Math.sin(angle) + (template.amino_acid.get(i).atom.y) * Math.cos(angle));

            } else {

                template.amino_acid.get(i).atom.x = ((template.amino_acid.get(i).atom.x) * Math.cos(angle) + (template.amino_acid.get(i).atom.y) * Math.sin(angle));

                template.amino_acid.get(i).atom.y = ((template.amino_acid.get(i).atom.y) * Math.cos(angle) - (template.amino_acid.get(i).atom.x) * Math.sin(angle));


            }

//            tc1.amino_acid.get(i).atom.y=(template1.amino_acid.get(i).atom.y+template2.amino_acid.get(i).atom.y)/2;
//            tc1.amino_acid.get(i).atom.z=(template1.amino_acid.get(i).atom.z+template2.amino_acid.get(i).atom.z)/2;
//        

        }

        return template;
    }

    private static Template find_avg(Template template1, Template template2, Template tc1) {


        Template[] t = new Template[3];
//        t[0]=template1;
//        t[1]=template2;
        t[2] = tc1;

        char axis = 'x';

        t[0] = template1;
        t[1] = template2;
//        t[0]= rotation_matx(axis,template1, 270.00);
//        t[1]= rotation_matx(axis,template2, -270.00);

        int length = t[0].amino_acid.size();
        for (int j = 0; j < 2; j++) {

            if (length > t[j].amino_acid.size()) {

                length = t[j].amino_acid.size();
            }

        }


        for (int i = 0; i < length; i++) {
            //if(template1.amino_acid.get(i)==null || template2.amino_acid.get(i)==null)
            // break;
            tc1.amino_acid.get(i).atom.x = (t[0].amino_acid.get(i).atom.x + t[1].amino_acid.get(i).atom.x) / 2;
            tc1.amino_acid.get(i).atom.y = (t[0].amino_acid.get(i).atom.y + t[1].amino_acid.get(i).atom.y) / 2;
            tc1.amino_acid.get(i).atom.z = (t[0].amino_acid.get(i).atom.z + t[1].amino_acid.get(i).atom.z) / 2;
        }
        return tc1;
    }

}
