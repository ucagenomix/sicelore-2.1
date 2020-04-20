package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.*;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SNP/editing events detection/quantification (cellBC/UMI count matrix).", oneLineSummary = "SNP/editing event detection/quantification (cellBC/UMI count matrix).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SNPMatrix extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input molecule SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "SNP", doc = "The SNP/editing events .csv file \n#-----\nname,gene,chromosome,start,end,strand,pos,ref,alt\nRG,Gria2,3,80681450,80802835,-,80692286,T,C\nQR,Gria2,3,80681450,80802835,-,80706912,T,C\n#-----")
    public File SNP;
    //@Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    //public File REFFLAT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=snp)")
    public String PREFIX = "snp";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "RN_min", doc = "Minimum number of reads for a molecule to be taken into account for SNP calling (default=0, means all molecules)")
    public int RN_min = 0;

    public CellList cellList;
    
    public SNPMatrix() {
        log = Log.getInstance(SNPMatrix.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsReadable(SNP);
        //IOUtil.assertFileIsReadable(REFFLAT);

        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t[" + this.cellList.size() + "]"});
        Matrix matrix = new Matrix(this.cellList);
        
        //UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        //String[] gria2 = {"Gria2"};
        //List<TranscriptRecord> transcripts = model.select(gria2);
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary() ;
        
        //editing.csv event list to quantify description file
        //chromosome,position,strand,name
        //11,75300577,-,Rpa1
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(SNP));
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                
                //System.out.println(tok[0] + "\t" + dictionnary.getSequence(tok[0]));
                
                if(dictionnary.getSequence(tok[0]) != null){
                    String chromosome = tok[0];
                    boolean strand = ("-".equals(tok[2]))?true:false;
                    String gene = tok[3];
                    int total=0;
                    nb++;
                    
                    if(java.util.regex.Pattern.matches(".*\\|.*", tok[1])){
                        
                        //System.out.println("2 positions:" + line);
                        
                        String[] pos = tok[1].split("\\|");
                        int pos1 = Integer.parseInt(pos[0]);
                        int pos2 = Integer.parseInt(pos[1]);
                        
                        //System.out.println(chromosome + "," + pos1 + "," + pos2);
                        
                        if(pos1 > pos2){
                            int tmp=pos2;
                            pos2=pos1;
                            pos1=tmp;
                        }
                        
                        //System.out.println(chromosome + "," + pos1 + "," + pos2);
                        
                        SAMRecordIterator iter = samReader.query(chromosome, pos1, pos2, false);
                        while(iter.hasNext()){
                            SAMRecord r = iter.next();
                            
                            //System.out.println("record:" + r);
                            
                            if(strand == r.getReadNegativeStrandFlag()){
                                pl.record(r);
                                LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, false);
                                
                                //System.out.println("good strand:" + lrr.getRn());
                                
                                if(lrr != null && lrr.getRn() >= RN_min){
                                    Longread lr = new Longread(r.getReadName());
                                    lr.addRecord(lrr);
                                    String readString = r.getReadString();
                                    int posInt1 = r.getReadPositionAtReferencePosition(pos1);
                                    int posInt2 = r.getReadPositionAtReferencePosition(pos2);
                                    
                                    //System.out.println(posInt1 + "\t" + posInt2 + "("+readString.length()+")");
                                    
                                    if(posInt1 > 0 && readString.length() > posInt1 && posInt2 > 0 && readString.length() > posInt2){
                                        String nuc1 = readString.substring(posInt1 - 1, posInt1);
                                        String nuc2 = readString.substring(posInt2 - 1, posInt2);
                                        
                                        if(r.getReadNegativeStrandFlag())
                                            nuc1 = complementBase(nuc1.charAt(0)); nuc2 = complementBase(nuc2.charAt(0));
                                        
                                        total++;
                                        Molecule molecule = new Molecule(lrr.getBarcode(), lrr.getUmi(), lrr.getRn());
                                        molecule.addLongread(lr);
                                        molecule.setGeneId(gene);
                                        molecule.setTranscriptId(chromosome + ":" + pos1 + "-" + pos2 + ".." + nuc1 + nuc2);
                                        matrix.addMolecule(molecule);
                                    }
                                }
                            }
                        }
                        iter.close();
                    }
                    else{
                        
                        //System.out.println("1 position:" + line);
                        
                        int position = Integer.parseInt(tok[1]);

                        SAMRecordIterator iter = samReader.query(chromosome, position, position, false);
                        while(iter.hasNext()){
                            SAMRecord r = iter.next();
                            
                            //System.out.println("record:" + r);
                            
                            if(strand == r.getReadNegativeStrandFlag()){
                                pl.record(r);
                                LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, false);
                                
                                //System.out.println("good strand:" + lrr.getRn());
                                
                                if(lrr != null && lrr.getRn() >= RN_min){
                                    Longread lr = new Longread(r.getReadName());
                                    lr.addRecord(lrr);
                                    String readString = r.getReadString();
                                    int posInt = r.getReadPositionAtReferencePosition(position);

                                    if(posInt > 0 && readString.length() > posInt){
                                        String nuc = readString.substring(posInt - 1, posInt);
                                        
                                        if(r.getReadNegativeStrandFlag())
                                            nuc = complementBase(nuc.charAt(0));
                                        
                                        total++;
                                        Molecule molecule = new Molecule(lrr.getBarcode(), lrr.getUmi(), lrr.getRn());
                                        molecule.addLongread(lr);
                                        molecule.setGeneId(gene);
                                        molecule.setTranscriptId(chromosome + ":" + position + ".." + nuc);
                                        matrix.addMolecule(molecule);
                                    }
                                }
                            }
                        }
                        iter.close();
                    }
                    log.info(new Object[]{"processing...\t\t" + line + "\t[" + total + " hits]\t"});
                }
                line = fichier.readLine();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }

        if(matrix.getMatrice().size() > 0){
            File MATRIX = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_matrix.txt");
            File METRICS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_metrics.txt");
            File molinfos = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_molinfos.txt");
            matrix.writeIsoformMatrix(MATRIX, METRICS, molinfos, null);
        }
        else
            log.info(new Object[]{"end of processing...\t\tnothing has been detected, check your input parameters, no output files generated\t"});
        
        return 0;
    }
    
    protected static String complementBase(char base)
    {
        String cc = "";
        if (base == 'A') cc = "T";
        if (base == 'C') cc = "G";
        if (base == 'G') cc = "C";
        if (base == 'T') cc = "A";
	
        return cc;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SNPMatrix().instanceMain(paramArrayOfString));
    }
}
