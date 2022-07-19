package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */  
import java.util.*;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.annotation.Strand;
import java.util.regex.Pattern;

public class LongreadRecord implements Comparable<LongreadRecord>
{
    private String name;
    private String barcode;
    private String umi;
    private int rn;
    private String chrom;
    private Strand strand;
    private int txStart;
    private int txEnd;
    private String geneId;
    private List<int[]> exons;
    private List<Junction> junctions;
    private byte[] cdna;
    private Float de;
    private int mapqv;
    private boolean isSecondaryOrSupplementary = false;
    private boolean isChimeria = false;
    
    protected static String CELLTAG = "BC"; // BC
    protected static String UMITAG = "U8"; // U8
    protected static String RNTAG = "RN"; // RN
    protected static String GENETAG = "GE"; //GE
    protected static String TSOENDTAG = "TE"; //TE
    protected static String POLYASTARTTAG = "PS"; //PS
    protected static String USTAG = "US"; // US
    protected static String CDNATAG = "CS"; // CS
    protected static int MAXCLIP = 150; // 150
    
    //FastqLoader fastq;
    
    public LongreadRecord() { }

    public void setStaticParams(String celltag, String umitag, String genetag, String tsoendtag, String polyastart, String ustag, String cdnatag, int maxclip, String rntag){
	this.CELLTAG = celltag;
        this.UMITAG = umitag;
        this.GENETAG = genetag;
        this.TSOENDTAG = tsoendtag;
        this.POLYASTARTTAG = polyastart;
        this.USTAG = ustag;
        this.CDNATAG = cdnatag;
        this.MAXCLIP = maxclip;
        this.RNTAG = rntag;
    }
    
    //public void setFastqLoader(FastqLoader f){
	//this.fastq = f;
    //}
    
    public int compareTo(LongreadRecord lr){
        Float obj1 = new Float(((LongreadRecord) lr).getDe());
        Float obj2 = new Float(this.de);
        int retval = obj2.compareTo(obj1);
        return retval;
    }

    public static LongreadRecord fromSAMRecord(SAMRecord r, boolean load_sequence) throws LongreadParseException
    {
        LongreadRecord record = new LongreadRecord();
        
        record.geneId = (String) r.getAttribute(GENETAG);
        record.barcode = (String) r.getAttribute(CELLTAG);
        record.umi = (String) r.getAttribute(UMITAG);
        record.mapqv = r.getMappingQuality();
        
        if (record.barcode == null || r.getReadUnmappedFlag())
             return null;
        
        record.barcode=record.barcode.replace("-1","");
        try {
            record.name = r.getReadName();
            record.chrom = r.getReferenceName();
            record.txStart = r.getAlignmentStart();
            record.txEnd = r.getAlignmentEnd();
            record.strand = Strand.toStrand((r.getReadNegativeStrandFlag()) ? "-" : "+");
            
            record.isSecondaryOrSupplementary = r.isSecondaryOrSupplementary();
            record.de = ((Float) r.getAttribute("de") != null) ? (Float) r.getAttribute("de") : (Float) r.getAttribute("df"); // minimap 2.17 (de) versus 2.10 (df)
            if(record.de == null)
                record.de = new Float(1);
            record.rn = ((Integer) r.getAttribute(RNTAG) != null) ? (Integer) r.getAttribute(RNTAG) : 1;
         
            //boolean isSoftOrHardClipped = false;
            int sizeStartToClip = 0;
            int sizeEndToClip = 0;
            String exon_starts = "";
            String exon_ends = "";
            List<AlignmentBlock> blocks = r.getAlignmentBlocks();
            String cigar = r.getCigarString();
            String[] cigartype = cigar.split("[0-9]+");
            String[] cigarsize = cigar.split("[A-Z]");
            
            // detect softclipping starting or ending reads
            if ("H".equals(cigartype[1]) || "S".equals(cigartype[1])){ sizeStartToClip = new Integer(cigarsize[0]).intValue(); }
            if ("H".equals(cigartype[cigartype.length - 1]) || "S".equals(cigartype[cigartype.length - 1])) { sizeEndToClip = new Integer(cigarsize[cigarsize.length - 1]).intValue(); }
            
            if(sizeStartToClip > MAXCLIP || sizeEndToClip > MAXCLIP)
                record.isChimeria = true;
            
            // added 01/09/2021 -> cDNA tag added after barcodes assignment
            // using AddBamReadSequenceTag SEQTAG=CS
            if(load_sequence && !record.isChimeria){
                
                // directly get the cDNA sequence
                if((String)r.getAttribute(CDNATAG) != null)
                    record.cdna = ((String)r.getAttribute(CDNATAG)).getBytes();
                else{ // else compute it by substring
                    
                    // "TSO---------cdna------------AAAA-UMI-BC-ADAPTOR"
                    String readSequence = (String)r.getAttribute(USTAG);
                    int tsoEnd = ((Integer) r.getAttribute(TSOENDTAG) != null) ? (Integer) r.getAttribute(TSOENDTAG) : 0;
                    int polyAStart = ((Integer) r.getAttribute(POLYASTARTTAG) != null) ? (Integer) r.getAttribute(POLYASTARTTAG) : 0;
                    int cdna_start = (tsoEnd != 0) ? tsoEnd : 0;
                    int cdna_end = (polyAStart != 0 && polyAStart < readSequence.length()-1) ? polyAStart : readSequence.length()-1;
                    String str = readSequence.substring(cdna_start, cdna_end);
                    //System.out.println(record.name+ "\t" + readSequence.length()+ "\t" + tsoEnd + "\t" + cdna_start + "\t" + polyAEnd + "\t" + cdna_end);
                    
                    record.cdna = str.getBytes();
                    if("".equals(str))
                        record.cdna = readSequence.getBytes();
                }
            }
                
            cigar = cigar.replaceAll("[0-9]+[IS]","");
            cigartype = cigar.split("[0-9]+");
            cigarsize = cigar.split("[A-Z]");
            
            // init exons
            int i;
            AlignmentBlock currBlock;
            int block_index = 0;
            int s = blocks.get(0).getReferenceStart();
            int e = blocks.get(0).getReferenceStart();
            for(i = 0; i < cigarsize.length; i++) {
                currBlock = blocks.get(block_index);
                
                if ("M".equals(cigartype[i])) { block_index++; }
                
                if ("N".equals(cigartype[i])){
                    //System.out.println(cigarsize[i-1]+""+cigartype[i]+ " at pos "+currBlock.getReferenceStart());
                    
                    exon_starts += s + ",";
                    exon_ends += e + ",";
                    // init the start of the next exon
                    s = currBlock.getReferenceStart();
                }
                else if ("D".equals(cigartype[i]) && new Integer(cigarsize[i-1]).intValue() > 20) { // case of short intron minimap2 consider as deletion
                    exon_starts += s + ",";
                    exon_ends += e + ",";
                    // init the start of the next exon
                    s = currBlock.getReferenceStart();
                }
                
                // keep track of exon end
                if(! "D".equals(cigartype[i]))
                    e = currBlock.getReferenceStart() + currBlock.getLength() - 1;
            }
            
            exon_starts += s + ",";
            exon_ends += e + ",";
            
            //System.out.println("exon_starts="+exon_starts);
            //System.out.println("exon_ends  ="+exon_ends);
            
            int[] exonStarts = LongreadRecord.toIntArray(exon_starts);
            int[] exonEnds = LongreadRecord.toIntArray(exon_ends);
            record.exons = new ArrayList<int[]>();
            for (i = 0; i < exonStarts.length; i++)
                record.exons.add(new int[]{exonStarts[i], exonEnds[i]});
            
            record.junctions = new ArrayList();
            for (i = 1; i < record.exons.size(); i++) {
                int j = ((int[]) record.exons.get(i - 1))[1];
                int k = ((int[]) record.exons.get(i))[0];
                record.junctions.add(new Junction(j, k));
            }
            
        } catch (Exception e) { e.printStackTrace(); throw new LongreadParseException("Invalid Bam file. " + record.name + ", Can't parse: "); }

        return record;
    }
    
    
    static int polyAStartFromReadName(String name)
    {
        int PAst = 0;
        // a7f6c0db-f0e0-498a-8037-aea6986ffbf6_REV_PAst=1341_PAen=1369_AEnd=1396_TS=46_bc=CCGCGGGTACGAAGAA_ed=0_ed_sec=6_bcStart=1395_bcEnd=1380_sq=AAAAAAAAAATACCCGCAGATATTCTTCGTACCCGCGGAGA_qv=23.9_
        String[] tab = name.split("_");
        for(int i=0;i<tab.length;i++){
            if(Pattern.matches("^PAst=.*", tab[i])){
                String[] t = tab[i].split("=");
                PAst = new Integer(t[1]).intValue();
            }
        }
        return PAst;
    }
    
    /*
    public static String complementWC(String dna) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < dna.length(); i++) {
            char c = dna.charAt(i);
            if (dna.charAt(i) == 'T') {
                builder.append('A');
            }
            if (dna.charAt(i) == 'A') {
                builder.append('T');
            }
            if (dna.charAt(i) == 'C') {
                builder.append('G');
            }
            if (dna.charAt(i) == 'G') {
                builder.append('C');
            }
        }
        return builder.reverse().toString();
    }
    */
    private static int[] toIntArray(String str) throws NumberFormatException {
        str = StringUtils.stripEnd(str, ",");
        String[] vals = str.split(",");
        int[] numbers = new int[vals.length];

        for (int i = 0; i < vals.length; i++) {
            numbers[i] = Integer.valueOf(vals[i]);
        }
        return numbers;
    }

    public String toString()
    {
        //String str = "[" + chrom + ':' + txStart + '-' + txEnd + " --> " + geneId + ", "+exons.size()+" ex ]\n"+ new String(cdna) +"\n";
        String str = "lrr[" + geneId + ", " + exons.size() + " ex]\n";
       // str += "Exon Starts: " + ArrayUtils.toString(exonStarts) + "\n";
       // str += "Exon Ends: " + ArrayUtils.toString(exonEnds) + "\n";
        return str;
    }

    public String getGeneId() { return geneId; }
    public byte[] getCdna() { return this.cdna; }
    public void setGeneId(String geneId) { this.geneId = geneId; }
    public int getMapqv() { return this.mapqv; }
    public float getDe() { return this.de; }
    public boolean getIsChimeria() { return isChimeria; }
    public boolean getIsSecondaryOrSupplementary() { return isSecondaryOrSupplementary; }
    public List<int[]> getExons() { return this.exons; }   
    public List<Junction> getJunctions() { return this.junctions; }

    public String getName() { return name; }
    public void setName(String name) { this.name=name; }

    public String getBarcode() { return barcode; }
    public void setBarcode(String barcode) { this.barcode=barcode; }

    public String getUmi() {return umi; }
    public void setUmi(String umi) { this.umi=umi; }

    public int getRn() {return rn; }
    public void setRn(int rn) { this.rn=rn; }

    public String getChrom() { return chrom; }
    public Strand getStrand() { return strand; }
    public int getTxStart() { return txStart; }
    public int getTxEnd() { return txEnd; }
    
    public String printFas() {
        return ">" + name + "\n" + new String(cdna) +"\n";
    }
}
