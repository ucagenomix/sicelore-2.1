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

public class Molecule implements Callable<String>
{
    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private HashSet<int[]> junctionSet;
    private String barcode;
    private String umi;
    private int rn;
    private Float pctId;
    private String consensus = "";
    private String geneId = "undef";
    private String transcriptId = "undef";
    private int supporting_reads=0;
    
    protected static String TMPDIR;
    protected static String POAPATH;
    protected static String RACONPATH;
    protected static String MINIMAP2PATH;
    
    private final static HashMap<Character, Integer> encode;
    private final static char[] decode;

    public Molecule() {}
    
    public Molecule(String barcode, String umi, int rn) {
        this.longreads = new ArrayList<Longread>();
        this.geneIds = new HashSet<String>();
        this.junctionSet = new HashSet<int[]>();
        this.barcode = barcode;
        this.umi = umi;
        this.rn = rn;
    }
    
    public Molecule(String barcode, String umi, String consensus, int rn) {
        //this.longreads = new ArrayList<Longread>();
        //this.geneIds = new HashSet<String>();
        this.junctionSet = new HashSet<int[]>();
        this.barcode = barcode;
        this.consensus = consensus;
        this.umi = umi;
        this.rn=rn;
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
        //System.out.println(TMPDIR);
    }

    public static int getIndexForChar(char base) {
        return (int) encode.get(base);
    }

    public static char getCharForIndex(int index) {
        return (char) decode[index];
    }

    public List<Longread> getLongreads() {
        return longreads;
    }

    public HashSet<String> getGeneIds() {
        return geneIds;
    }
    
    public String[] getGeneIdsArray() {
        String[] array = new String[this.geneIds.size()];
        this.geneIds.toArray(array);
        return array;
    }

    public HashSet<int[]> getJunctionSet() {
        return this.junctionSet;
    }

    public String getBarcode() {
        return this.barcode;
    }

    public String getUmi() {
        return this.umi;
    }

    public int getRn() {
        return this.rn;
    }

    public String getConsensus() {
        return this.consensus;
    }

    public String getGeneId() {
        return this.geneId;
    }
    public Float getPctId() {
        return this.pctId;
    }
    
    public int getSupporting_reads() {
        return this.supporting_reads;
    }

    public void setRn(int rn) {
        this.rn = rn;
    }

    public String getTranscriptId() {
        return this.transcriptId;
    }
    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId=transcriptId;
    }
    
    public void setSupporting_reads(int r) {
        this.supporting_reads=r;
    }
    
    public String getLabel() {
        return this.barcode + "-" + this.umi + "-" + this.longreads.size();
    }
    
    public int getNumberOfReads()
    {
        if(this.rn > 1)
            return this.rn;
        else
            return this.longreads.size();
    }
    
    public String toString()
    {
        String str = "cell="+this.barcode+",umi="+ this.umi + ", reads=" + this.longreads.size() + ", " + this.geneIds.toString();
        Iterator<Longread> it = this.longreads.iterator();
        while(it.hasNext()){
            Longread lr = (Longread) it.next();
            List<LongreadRecord> lrr = lr.getLongreadrecords();
            str += "\n\t" + lr.getName() + ", SAMrecords="+lrr.size();
        }
        
        return str;
    }
    
    public void addLongread(Longread lr)
    {
        this.longreads.add(lr);
        this.pctId = new Float(1) - lr.getLongreadrecords().get(0).getDe();
        
        Iterator<String> iterator = lr.getGeneIds().iterator();
        while (iterator.hasNext())
            this.geneIds.add((String)iterator.next());
    }
    
    public void addJunction(int[] j)
    {
        this.junctionSet.add(j);
    }
    
    public String call() throws Exception
    {
        String line;
        DataOutputStream os=null;
        //DataOutputStream os10=null;
        String prefix = this.barcode + "" + this.umi;
        
        // if only one read -> get best record as consensus
        if(this.longreads.size() == 1){
            LongreadRecord lrr = this.longreads.get(0).getBestRecord();
            this.consensus = new String(lrr.getCdna());
        }
        // if 2 reads, take the one with lowest de
        else if(this.longreads.size() == 2){
            LongreadRecord lrr1 = this.longreads.get(0).getBestRecord();
            LongreadRecord lrr2 = this.longreads.get(1).getBestRecord();
            this.consensus = (lrr1.getDe() < lrr2.getDe())?new String(lrr1.getCdna()):new String(lrr2.getCdna());
        }
        // if more do consensus calling
        else{            
            try {
                //int count=0;
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads.fa"));
                //os10 = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads10.fa"));
                Iterator<Longread> iterator2 = this.longreads.iterator();
                while(iterator2.hasNext()){
                    Longread lr = (Longread) iterator2.next();
                    LongreadRecord lrr = lr.getBestRecord();
                    os.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                    //if(count++<11)
                    //    os10.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                }
                os.close();
                //os10.close();
                
                // compute consensus with POA
                String[] commande = {"bash", "-c" , ""};
                String poadir = POAPATH.replaceAll("/poa$","");
                commande[2] = POAPATH + " -read_fasta "+TMPDIR+"/"+prefix+"_reads.fa -pir "+TMPDIR+"/"+prefix+".pir " + poadir + "/blosum80.mat -hb -best";
                //System.out.println(commande[2]);
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                // get CONSENS0 in fasta file .pir
                this.consensus = "";
                BufferedReader file = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + ".pir"));
                line = file.readLine();
                while(! ">CONSENS0".equals(line)){ line = file.readLine(); }
                while(line != null && !"".equals(line)){ line = file.readLine(); this.consensus += line; }                
                file.close();
                this.consensus = this.consensus.replaceAll("-", "");
                this.consensus = this.consensus.replaceAll("\\n", "");
                this.consensus = this.consensus.replaceAll("null", "");
                
                // do we still need to do that after poa consensus ?
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_consensus.fa"));
                os.writeBytes(">"+this.barcode + this.umi+"\n"+this.consensus+"\n");
                os.close();

                commande[2] = MINIMAP2PATH + " --secondary=no -ax map-ont "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_reads.fa > "+TMPDIR+"/"+prefix+"_overlap.sam";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                commande[2] = RACONPATH + " "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa > "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                
                BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + "_corrected_consensus.fa"));
                line = fichier.readLine();
                this.consensus = fichier.readLine();
                fichier.close();
                
                commande[2] = "rm "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+".pir "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
           }
            catch(Exception e){ e.printStackTrace(); }
        }
        
        return ">" + this.getLabel() + "\n" + this.consensus + "\n";
    }
}
