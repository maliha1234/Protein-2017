
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;



public class ModelGeneration {
    
        int matched_count;
      //  int[] align;
        double[][] in_distance;
        double[][] t_distance;
        TemplateMatching.Template input, Tc;
        TemplateMatching.Template[] t;
    
    public ModelGeneration(TemplateMatching.Template[] t,TemplateMatching.Template input){
        
        this.input=input;
        this.t=t;
     //   this.align=t.matched_align;
       // this.in_distance=input.d;
     //   this.t_distance=t.d;
       
  }
    
      public ModelGeneration(){
        
       
       
  }
   public double calc_initial_drmsd(AVG_PDB.Template[] templates,AVG_PDB.Template Tc, int count) throws IOException{
       
        double drmsd = 0,difference,sum_drmsd=0, sum=0, root_drmsd=0;
        int matched_count=count,min,i;
         
       
          
        min=Tc.amino_acid.size();
          
        for(i=0; i<4; i++){
        
            if(min>templates[i].amino_acid.size()) 
                min=templates[i].amino_acid.size();
        }
          
          for(int k=0; k<4; k++){
          
              for(i=0; i<min; i++){
              
             
                  for(int j=i+1; j<min; j++){
                     
                          
                          difference=templates[k].d[i][j]-Tc.d[i][j];
                          sum_drmsd+=Math.pow(difference, 2);
                          
                  }  
                 
                  
                }
             
     
                   drmsd=sum_drmsd/(matched_count*(matched_count-1));
          
                
                   root_drmsd=Math.sqrt(drmsd);
                 // System.out.println("drmsd "+ " "+ root_drmsd); 
                  sum+=root_drmsd;
                 
     
          }
          
           System.out.println("sum "+ sum);
           return sum;
   }
   
   
 /*  void WritePdb(TemplateMatching.Template Tc) throws IOException{
       
       File file = new File("E/4-1/Thesis/filename.txt");
 
			// if file doesnt exists, then create it
			if (!file.exists()) {
				//file.createNewFile();
			}
 
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			bw.write(" this is the file ");
			bw.close();
 
       
       
   }*/
   
   
   public TemplateMatching.Template generateConsensus(double[] weights){
       
       int i,j,k,min_length;
       double sumx=0,sumy=0,sumz=0,newx,newy,newz,score=0;
    
       TemplateMatching.Atom t_atom;
       TemplateMatching.Template Tc = null;
       
       min_length=t[0].amino_acid.size();
       for(i=1; i<2; i++){
       
           // System.out.println("min :"+ min_length + "t%d :"+ t[i].amino_acid.size());
            if(min_length>t[i].amino_acid.size()) {
                min_length=t[i].amino_acid.size();
                System.out.println("min :"+ min_length );
            }
       }
   
         System.out.println("final min :"+ min_length );
          double[] matched_x = new double[min_length+1];
          double[] matched_y = new double[min_length+1];
          double[] matched_z = new double[min_length+1];
  
       List<TemplateMatching.Atom> atoms=new ArrayList<TemplateMatching.Atom>();
       List<TemplateMatching.AminoAcid> aminoacids=new ArrayList<TemplateMatching.AminoAcid>();
       for(i=0; i<2; i++)
           
       {
         // int[] align_indices = new int[t[i].matched_align.length];
          
            
          for(j=0; j<min_length; j++){
              
            
             t_atom = t[i].amino_acid.get(j).atom; 
           
          
         
         /*  if(t[i].matched_align[j]==1){
               //System.out.println(j);
             
           //  System.out.println("t_atom x:"+t_atom.x);
             score=t_atom.x*weights[i];
           //  System.out.println("weight :" + weights[i]);
             matched_x[j]+=score;
             matched_y[j]+=score;
             matched_z[j]+=score;
          //    System.out.println("scored matched_x: "+matched_x[j]);
             if(j==min_length-1)  { System.out.println(i+"---------------- ");;}
          
           }*/
         //  else {
               matched_x[j]+=t_atom.x; matched_y[j]+=t_atom.y; matched_z[j]+=t_atom.z;
            //   System.out.println("matched_x: "+matched_x[j]);
           //}
           
          
          }
       }
       
      int seqid=0; String symbol;
      for(k=0; k<min_length; k++){
       newx=matched_x[k]/2;               // 55.86 is the sum of all the Z-Scores
      // System.out.println(" "+ k + "  newx: "+ newx);
       
       newy=matched_y[k]/2; 
       newz=matched_z[k]/2;
       
       atoms.add(k, new TemplateMatching.Atom(newx,newy,newz,"CA"));
       symbol=input.amino_acid.get(k).symbol;
       seqid=input.amino_acid.get(k).seqid;
       aminoacids.add(k, new TemplateMatching.AminoAcid(atoms.get(k),symbol,seqid));
      }
     /*  for(int k=0; k<Tc.d.length; k++)
                  for(int p=k+1; p<Tc.d.length; p++)
                      System.out.println(Tc.d[k][p]+":");*/
       //System.out.println("Concensus Model"+Tc.d);
       //System.out.println(null);
      
      Tc=new TemplateMatching.Template(aminoacids);
      Tc.d=Tc.distance(Tc.amino_acid, min_length);
       return Tc;
   }
    
   
    
    public TemplateMatching.Template changeConsensus(TemplateMatching.Template T){
        
        int i;
        double x1,y1,z1,x,y,z;
        TemplateMatching.Template Tc1=new TemplateMatching.Template(30,"consensus",T.d, T.amino_acid);
        
        for(i=0; i<T.amino_acid.size(); i++){
            
            x1=T.amino_acid.get(i).atom.x;
            x=(x1+TemplateMatching.input_Template.amino_acid.get(i).atom.x)/2;
            
            y1=T.amino_acid.get(i).atom.y;
            y=(y1+TemplateMatching.input_Template.amino_acid.get(i).atom.y)/2;
            
            z1=T.amino_acid.get(i).atom.z;
            z=(z1+TemplateMatching.input_Template.amino_acid.get(i).atom.z)/2;
            
                        
          /*  Tc1.amino_acid.get(i).atom.x=x;
            Tc1.amino_acid.get(i).atom.y=y;
            Tc1.amino_acid.get(i).atom.z=z;
            */
            // System.out.println("T er x1 :"+ x1+"input x :"+TemplateMatching.input_Template.amino_acid.get(i).atom.x);
            Tc1.amino_acid.get(i).atom=new TemplateMatching.Atom(x,y,z,"CA");
      //   System.out.println("x :"+Tc1.amino_acid.get(i).atom.x + " y "+ Tc1.amino_acid.get(i).atom.y+" z:"+Tc1.amino_acid.get(i).atom.z);
      //  System.out.println("T er x1 :"+ x1+"input x :"+TemplateMatching.input_Template.amino_acid.get(i).atom.x);
          
        }
      
        
        
        return Tc1;
    }
    
}
