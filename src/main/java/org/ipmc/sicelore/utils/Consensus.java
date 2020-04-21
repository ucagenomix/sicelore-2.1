package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 *  
 */
import java.util.*;
import org.ipmc.common.utils.ExecuteCmd;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.concurrent.*;

public class Consensus implements Callable<String>
{
    private List<SubConsensus> sequences;
    private String name;
    private String cons;
    
    protected static String TMPDIR;
    protected static String POAPATH;
    protected static String RACONPATH;
    protected static String MINIMAP2PATH;
    
    private final static HashMap<Character, Integer> encode;
    private final static char[] decode;

    public Consensus() {}
    
    public Consensus(String name, List<SubConsensus> lst) {
        this.name = name;
        this.sequences = lst;
    }
    
    static {
        encode = new HashMap<Character, Integer>();
        encode.put('-', 0);
        encode.put('A', 1);
        encode.put('T', 2);
        encode.put('C', 3);
        encode.put('G', 4);
        encode.put('a', 1);
        encode.put('t', 2);
        encode.put('c', 3);
        encode.put('g', 4);

        decode = new char[5];
        decode[0] = '-';
        decode[1] = 'A';
        decode[2] = 'T';
        decode[3] = 'C';
        decode[4] = 'G';
    }

    public void setStaticParams(String tmp, String poa, String racon, String minimap2){
        this.TMPDIR = tmp;
        this.POAPATH = poa;
        this.RACONPATH = racon;
        this.MINIMAP2PATH = minimap2;
    }

    public static int getIndexForChar(char base) { return (int) encode.get(base); }
    public static char getCharForIndex(int index) { return (char) decode[index]; }

    public List<SubConsensus> getSequences() { return sequences; }
    public String getName() { return this.name; }
    public String getCons() { return this.cons; }
    
    public String toString(){ return "name:"+this.name+",nbSeq="+ this.sequences.size(); }
    public String toFasta(){ return ">"+this.name+"\n"+ this.cons+"\n"; }
    
    public String call() throws Exception
    {
        String line;
        DataOutputStream os=null;
        
        try {
            os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+this.name+"_reads.fa"));
            Iterator<SubConsensus> iterator2 = this.sequences.iterator();
            while(iterator2.hasNext()){
                SubConsensus sc = (SubConsensus) iterator2.next();
                os.writeBytes(">"+sc.getName()+"\n"+new String(sc.getSequence())+"\n");
            }
            os.close();
            
            // compute consensus with POA
            String[] commande = {"bash", "-c" , ""};
            String poadir = POAPATH.replaceAll("/poa$","");
            commande[2] = POAPATH + " -read_fasta "+TMPDIR+"/"+this.name+"_reads.fa -pir "+TMPDIR+"/"+this.name+".pir " + poadir + "/blosum80.mat -hb -best";
            //System.out.println(commande[2]);
            ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
            executeCmd.run();

            // get CONSENS0 in fasta file .pir
            this.cons = "";
            BufferedReader file = new BufferedReader(new FileReader(TMPDIR+"/"+ this.name + ".pir"));
            line = file.readLine();
            while(! ">CONSENS0".equals(line)){ line = file.readLine(); }
            while(line != null && !"".equals(line)){ line = file.readLine(); this.cons += line; }                
            file.close();
            
            this.cons = this.cons.replaceAll("-", "");
            this.cons = this.cons.replaceAll("\\n", "");
            this.cons = this.cons.replaceAll("null", "");
                
            // do we still need to do that after poa consensus ?
            os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+this.name+"_consensus.fa"));
            os.writeBytes(">"+this.name+"\n"+this.cons+"\n");
            os.close();

            commande[2] = MINIMAP2PATH + " --secondary=no -ax map-ont "+TMPDIR+"/"+this.name+"_consensus.fa "+TMPDIR+"/"+this.name+"_reads.fa > "+TMPDIR+"/"+this.name+"_overlap.sam";
            executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
            executeCmd.run();

            commande[2] = RACONPATH + " "+TMPDIR+"/"+this.name+"_reads.fa "+TMPDIR+"/"+this.name+"_overlap.sam "+TMPDIR+"/"+this.name+"_consensus.fa > "+TMPDIR+"/"+this.name+"_corrected_consensus.fa";
            executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
            executeCmd.run();
                
            BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+ this.name + "_corrected_consensus.fa"));
            line = fichier.readLine();
            this.cons = fichier.readLine();
            fichier.close();
                
            commande[2] = "rm "+TMPDIR+"/"+this.name+"_reads.fa "+TMPDIR+"/"+this.name+".pir "+TMPDIR+"/"+this.name+"_overlap.sam "+TMPDIR+"/"+this.name+"_consensus.fa "+TMPDIR+"/"+this.name+"_corrected_consensus.fa";
            executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
            executeCmd.run();
            
        }catch(Exception e){ e.printStackTrace(); }
        
        return this.toFasta();
    }
}
