
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
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
 *
 * @author hp user
 */
public class TestClass{

  
   public static Template input_Template;

 
   static class Atom{
       
        double x,y,z;
        String type;
       
        
        Atom(double x,double y,double z,String type){
            
            this.x=x;
            this.y=y;
            this.z=z;
            this.type=type;
            
        }
        
    }
    
    static class AminoAcid{
        
        Atom atom;
        String symbol;
        int seqid;
        void printAminoAcids(AminoAcid aa){
            
            System.out.println("Symbol :"+aa.symbol+" seqid :"+aa.seqid);
        }
        
        AminoAcid(Atom atom, String symbol,int seqid){
            
            this.atom=atom;
            this.symbol=symbol;
            this.seqid=seqid;
        }
        
        
        
        
    }
    
    static class Structure{
    
        AminoAcid amino_acid;
        
        Structure(AminoAcid amino_acid){
            
            this.amino_acid=amino_acid;
            
        }
    
    
}
    
    static class Template{
        
        float score;
        String template_name;
        int index;
        int[] matched_align;
        float weight;
        char[] sequence;
        double[][] d;
       // Structure structure;
        List<AminoAcid> amino_acid;
      
        /** this method creates intra-distance matrix for a template ***/
      static double[][] distance(List<AminoAcid> a,int length){
        
        int i=0,j=0;
        double[][] d = new double[length][length];
        
        for(i=0; i<length; i++)
        for(j=i+1; j<length; j++){
        
            double temp1,temp2,temp3,sum=0.00;
            temp1=(a.get(i).atom.x-a.get(j).atom.x);
            temp2=Math.pow(temp1,2);
            sum+=temp2;
            temp1=(a.get(i).atom.y-a.get(j).atom.y);
            temp2=Math.pow(temp1,2);
            sum+=temp2;
            temp1=(a.get(i).atom.z-a.get(j).atom.z);
            temp2=Math.pow(temp1, 2);
            sum+=temp2;
            d[i][j]=Math.sqrt(sum);
       
        
        }
         
     /*  for(i=0; i<length; i++)
        for(j=i+1; j<length; j++){
            
            System.out.println(" "+d[i][j]+" ");
        }*/
       return d;
    }
  

      Template(){
          
          this.amino_acid=null;
          this.d=null;
          this.index=0;
          this.matched_align=null;
          this.score=0;
          this.sequence=null;
          this.template_name=null;
          this.weight=0;
      }
      
      /***this constructor is for building input template i.e. target protein ***/
      Template(float score,String template_name,double[][] d,List<AminoAcid> amino_acid){
          
            this.score=score;
            this.template_name=template_name;
          //this.matched_align=align;
            this.d=d;
        //  this.weight=weight;
            this.sequence=sequence;
           //this.structure=structure;
            this.amino_acid=amino_acid; 
          
      }
      
      Template(List<AminoAcid> amino_acid){
          
          this.amino_acid=amino_acid;
         // this.d=d;
          
      }
        
       /***this constructor is for building 10 matched templates***/
        Template(String template_name,List<AminoAcid> amino_acid){
           
          //  this.score=score;
            this.template_name=template_name;
           // this.matched_align=align;
            this.d=d;
        //   this.weight=weight;
            this.sequence=sequence;
           //this.structure=structure;
            this.amino_acid=amino_acid;  //readPdb theke aminoacid array ekhane ashbe, for one template, we get an array of CA aminoacids
        
        
        }
        
        
        
        
    }
    
    static Template readPdb(String pdb_file){
            
            
        BufferedReader br = null;
        Template t = null;
      //  t = new Template();
 
		try {
 
			String sCurrentLine,lines;
                        
                        int seqid,i = 0,length;
                        int[] matched_align;
                        float x,y,z;
                        String type,symbol;
                        String file_path="E:\\4-2\\Thesis 4-2\\TemplateMatching_1\\";
                   //   while(pdb_files!=null){   
			
                        br = new BufferedReader(new FileReader(file_path+pdb_file));
                        List<Atom> atom =new ArrayList<Atom>();
                        //atom = new List();
                        List<AminoAcid> aminoacid = new ArrayList<AminoAcid>();
			while ((sCurrentLine = br.readLine()) != null) {
                            
                            if(sCurrentLine.startsWith("TER")){
                                
                                break;
                            }
                            
                            if(sCurrentLine.startsWith("ATOM") || sCurrentLine.startsWith("HETATM")){
				
                                
                            StringTokenizer st = new StringTokenizer(sCurrentLine);
                            ArrayList l=new ArrayList();
                            //System.out.println("---- Split by space ------");
                            while (st.hasMoreElements()) {

                                l.add(st.nextElement());

                            }
                           
                            if(l.get(2).equals("CA"))
                            {   
                             //   System.out.println(sCurrentLine);
                               //System.out.println(i);
                                x=Float.valueOf(l.get(6).toString());
                                y=Float.valueOf(l.get(7).toString());
                                z=Float.valueOf(l.get(8).toString());
                                type="CA";
                                symbol=l.get(3).toString();
                                seqid=Integer.parseInt(l.get(5).toString());
                                atom.add(i, new Atom(x,y,z,type));
                           //  atom.get(i).printAtoms(atom.get(i));
                                aminoacid.add(i,new AminoAcid(atom.get(i),symbol,seqid));
                          //      aminoacid.get(i).printAminoAcids(aminoacid.get(i));
                          //      System.out.println(i);
                                i++;
                              //  System.out.print(x+"____"+y+"_____"+z+"\n");
                            }
                            }
			}
                         length=i;
                    //    Template t[file_index]=new Template(length);
                       
      
                       
                        t=new Template(pdb_file,aminoacid);
                    
                      
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
                              
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
            return t;  
                
   }
   
  
    
   public static void main(String[] args) throws IOException{
        
   
        Template template1=new Template();
        Template template2=new Template();
        Template concensus=new Template();
        
        template1=readPdb("2ltm.pdb");
        System.out.print("***********\n");  
        
        template2=readPdb("1pqx.pdb");
        System.out.print("***********\n");
        
        
        concensus=readPdb("2M8W.pdb");
       
        Template c2=new Template();
        c2=find_avg(template1,template2,concensus);
        
       
        
       for(int i=0;i<concensus.amino_acid.size();i++)
        {
            System.out.print("original:");
           printAtoms(concensus.amino_acid.get(i).atom);
            
          //  write_pdb(template2.amino_acid.get(i).atom);
          //  System.out.println(template2.amino_acid.get(i).symbol+" "+template2.amino_acid.get(i).seqid);
          //System.out.print("avg:");
            printAtoms(c2.amino_acid.get(i).atom);
        }
        
      
    //  write_pdb();
      //   write_pdb(template2);
        
    }
    
    
    public static void write_pdb(Template t){
        
        
        BufferedReader br = null;
        
 
	try {
 
			String sCurrentLine,lines;
                        
                        int seqid,i = 0,length;
                        int[] matched_align;
                        float x,y,z;
                        String type,symbol;
                        
                        String file_path="E:\\4-2\\Thesis 4-2\\TemplateMatching_1\\";
                     
			
                        br=new BufferedReader(new FileReader(file_path+"2M8W.pdb"));
                        List<Atom> atom =new ArrayList<Atom>();
                        
                        List<AminoAcid> aminoacid = new ArrayList<AminoAcid>();
			while ((sCurrentLine = br.readLine()) != null) {
                            
                        File file= new File("E:\\4-2\\Thesis 4-2\\TemplateMatching_1\\newfile.txt");
			
			if (!file.exists()) {
				file.createNewFile();
			}

			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
                        
                        bw.write(sCurrentLine);
                        bw.newLine();
                        bw.close();
                            
                        if(sCurrentLine.startsWith("TER")){
                                
                                break;
                               // bw.close();
                        }
                            
                         /*   if(sCurrentLine.startsWith("ATOM") || sCurrentLine.startsWith("HETATM")){
				
                                
                            StringTokenizer st = new StringTokenizer(sCurrentLine);
                            ArrayList l=new ArrayList();
                            //System.out.println("---- Split by space ------");
                            while (st.hasMoreElements()) {

                                l.add(st.nextElement());

                            }
                           
                            if(l.get(2).equals("CA"))
                            {   
                                x=Float.valueOf(l.get(6).toString());
                                y=Float.valueOf(l.get(7).toString());
                                z=Float.valueOf(l.get(8).toString());
                                type="CA";
                                symbol=l.get(3).toString();
                                seqid=Integer.parseInt(l.get(5).toString());
                                atom.add(i, new Atom(x,y,z,type));
                           
                                aminoacid.add(i,new AminoAcid(atom.get(i),symbol,seqid));
                                i++;
                            }
                            }*/
			}
                         length=i;
                      
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
                              
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
           
        
        
 /* try {

			String content=null;

			//File file = new File("/users/mkyong/filename.txt");
                        File file= new File("E:\\4-2\\Thesis 4-2\\TemplateMatching_1\\filename.txt");
			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}

			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
                        
                        for(int i=0; i<t.amino_acid.size(); i++){
                        
                            
                        content = t.amino_acid.get(i).symbol+"  "+t.amino_acid.get(i).atom.x+"   "+t.amino_acid.get(i).atom.y+"   "+t.amino_acid.get(i).atom.z+"\n";
			bw.write(content);
                        bw.newLine();
                        
                        }
			bw.close();

			System.out.println("Done");

		} catch (IOException e) {
			e.printStackTrace();
		}
    */    
        
  }
    
    
    
 public static void printAtoms(Atom atom){
            
            System.out.println("X :"+atom.x+" Y :"+atom.y+" Z :"+atom.z+" type :"+atom.type);
 }
    
    
       private static Template find_avg(Template template1, Template template2, Template tc1) {
           
        int length=0;
        if(template1.amino_acid.size()>template2.amino_acid.size())
            length=template2.amino_acid.size();
        
        else
            length=template1.amino_acid.size();
        
        double x,y,z;
        for(int i=0;i<length;i++)
        {
            //if(template1.amino_acid.get(i)==null || template2.amino_acid.get(i)==null)
               // break;
            tc1.amino_acid.get(i).atom.x=(template1.amino_acid.get(i).atom.x+template2.amino_acid.get(i).atom.x)/2;
            tc1.amino_acid.get(i).atom.y=(template1.amino_acid.get(i).atom.y+template2.amino_acid.get(i).atom.y)/2;
            tc1.amino_acid.get(i).atom.z=(template1.amino_acid.get(i).atom.z+template2.amino_acid.get(i).atom.z)/2;

        }
         return tc1;
    }
      
}
