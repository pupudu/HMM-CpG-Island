import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
****Notes****
    * plus implies a CpG island
    * minus implies a non CpG island
    * Two strings sP and sM are calculated at each iteration to do Viterbi decoding easily 
    * If a letter other than ACG or T is appeared, it will be considered as a T type neucleotide
    * Emmission matrix is not defined here, since the probability of emmiting N from N+ and N- are 1 while from M+ and M- are 0(M be any neucleotide except N) 
*/

public class CpGIsland {
    
    public static void main(String[] args) {
        String fileName = "input.txt";
        if(args.length>0){
            fileName = args[0];
        }
        String test = getString(fileName);
        
        double trans[][] = {
            {0.0000000000, 0.0725193101, 0.1637630296, 0.1788242720, 0.0754545682, 0.1322050994, 0.1267006624, 0.1226380452, 0.1278950131},
            {0.0010000000, 0.1762237762, 0.2682517483, 0.4170629371, 0.1174825175, 0.0035964036, 0.0054745255, 0.0085104895, 0.0023976024},
            {0.0010000000, 0.1672435130, 0.3599201597, 0.2679840319, 0.1838722555, 0.0034131737, 0.0073453094, 0.0054690619, 0.0037524950},
            {0.0010000000, 0.1576223776, 0.3318881119, 0.3671328671, 0.1223776224, 0.0032167832, 0.0067732268, 0.0074915085, 0.0024975025},
            {0.0010000000, 0.0773426573, 0.3475514486, 0.3759440559, 0.1781818182, 0.0015784216, 0.0070929071, 0.0076723277, 0.0036363636},
            {0.0010000000, 0.0002997003, 0.0002047952, 0.0002837163, 0.0002097902, 0.2994005994, 0.2045904096, 0.2844305694, 0.2095804196},
            {0.0010000000, 0.0003216783, 0.0002977023, 0.0000769231, 0.0003016983, 0.3213566434, 0.2974045954, 0.0778441558, 0.3013966034},
            {0.0010000000, 0.0002477522, 0.0002457542, 0.0002977023, 0.0002077922, 0.2475044955, 0.2455084915, 0.2974035964, 0.2075844156},
            {0.0010000000, 0.0001768232, 0.0002387612, 0.0002917083, 0.0002917083, 0.1766463536, 0.2385224775, 0.2914165834, 0.2914155844}
        };
        
        double aP;  // probability of current level plus state
        double aM;  //probability of current level Minus state
        double bP; //temp storage for aP
        double bM;  //temp storage for aM
        
        String sP="+"; // For doing the viterbi decoding
        String sM="-"; // For doing the viterbi decoding
        
        int n1 = getPos(test, 0);   //initial plus state index(ACGT)
        int n2 = n1+4;              //initial Minus state index(ACGT)
        
        aP = trans[0][n1]; //initial plus state probabilty
        aM = trans[0][n2];  //initial Minus state probability 
        
        int c1 = n1;    //previous level plus state index
        int c2 = n2;    //previous level minus state index
        
        for(int i=1;i<test.length();i++){   // Go neucleotide by neucleotide for whole genome
            n1 = getPos(test,i);
            n2= n1+4;
            if(aP*trans[c1][n1]> aM*trans[c2][n1]){
                bP = aP*trans[c1][n1];
                sP = sP+"+"; // to a Plus state from a Plus state
            }else{
                bP = aM*trans[c2][n1];
                sP = sM+"+"; // to a Plus state from a Minus state
            }

            if(aP*trans[c1][n2]> aM*trans[c2][n2]){
                bM = aP*trans[c1][n2];
                sM = sP + "-"; // to a Minus state from a Plus state
            }else{
                bM = aM*trans[c2][n2];
                sM = sM+ "-"; // to a Minus state from a Minus state
            }
            
            c1=n1; //define previous plus state
            c2=n2;  //define previous minus state
            aP = bP;  // update plus state probability
            aM = bM;    //update minus state probability 
            
        }   
        if(aP>aM){  //Viterbi decoding, see how sP and sM are calculated in each iteration
            System.out.println("States using Viterbi Decoding: "+sP); 
            writeOutput(sP);
        }else{
            System.out.println("States using Viterbi Decoding: "+sM);
            writeOutput(sM);
        }
    }
    public static int getPos(String s, int i){  // returns the neucleotide as an index to refer transition matrix
        s=s.toUpperCase();
        if(s.charAt(i)=='A'){
            return 1;
        }else if(s.charAt(i)=='C'){
            return 2;
        }else if(s.charAt(i)=='G'){
            return 3;
        }else{
            if(s.charAt(i)!='T'){
                System.out.println("Letter '"+ s.charAt(i)+"' was replaced with a letter 'T'");
            }
            return 4;
        }
    }
    
    public static String getString(String fileName){
        File file = new File(fileName);
        FileReader fr = null;
        BufferedReader br = null;
        try {
            fr = new FileReader(file);
            br = new BufferedReader(fr);
            return br.readLine();
        } catch (Exception ex) {
            System.out.println(ex);
        }
        finally{
            try {
                if(fr!=null){
                    fr.close();
                }
                if(br!=null){
                    br.close();
                }
            } catch (IOException ex) {
                System.out.println(ex);
            }
        }
        System.out.println("Error reading file. So let us use the sample string for this run");
        return "ATCGCGCGCGGCGCAAAAAAAAAAAAAAAAAAGATATATCGCGCGCGCGCGAGCTA";
    }
    
    public static void writeOutput(String out){
        File file = new File("output.txt");
        FileWriter fw = null;
        BufferedWriter bw=null;
        try {
            fw = new FileWriter(file);
            bw = new BufferedWriter(fw);
            bw.write(out);
            System.out.println("Output.txt updated");
        } catch (IOException ex) {
            Logger.getLogger(CpGIsland.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            try{
            if(bw!=null){
                bw.close();
            }
            if(fw!=null){
                fw.close();
            }
            }catch(IOException e){
                System.out.println(e);
            }
        }
    }
}

