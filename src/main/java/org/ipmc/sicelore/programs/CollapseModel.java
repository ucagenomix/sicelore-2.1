package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */ 
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Collapse molecules Bam file to model Bam File", oneLineSummary = "Collapse molecules Bam file to model Bam File", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class CollapseModel extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "MINEVIDENCE", doc = "Minimum evidence for Novel isoforms (default=5 molecules)")
    public int MINEVIDENCE = 5;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=CollapseModel)")
    public String PREFIX = "CollapseModel";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "ISOFORMTAG", doc = "Isoform tag (default=IT)", optional=true)
    public String ISOFORMTAG = "IT";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    
    @Argument(shortName = "SHORT", doc = "The short read SAM or BAM file fot junction validation")
    public File SHORT;
    @Argument(shortName = "CAGE", doc = "CAGE peaks file (.bed)")
    public File CAGE;
    
    // ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.polyAs.gtf.gz
    // more +6 gencode.vM24.polyAs.gtf  | grep polyA_signal | awk '{ print $1 "\t" $4 "\t" $5 "\t" $1 ":" $4 ".." $5 "," $7 "\t1\t" $7 "\t" $4 "\t" $4+1 }' > gencode.vM24.polyAs.bed
    // sed -i -e "s/chr//g" gencode.vM24.polyAs.bed
    @Argument(shortName = "POLYA", doc = "POLYA sites file (.bed)")
    public File POLYA;

    @Argument(shortName = "cageCo", doc = "CAGE validation cutoff (default=50 bases)")
    public int cageCo = 50;
    @Argument(shortName = "polyaCo", doc = "PolyA validation cutoff (default=50 bases)")
    public int polyaCo = 50;
    @Argument(shortName = "juncCo", doc = "Junction validation cutoff (default=1 read)")
    public int juncCo = 1;
    
    public HashSet<String> allcmp = new HashSet<String>();
    BEDParser cage;
    BEDParser polyA;
    
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;
    
    UCSCRefFlatParser model;
    UCSCRefFlatParser mymodel;

    public CollapseModel()
    {
        log = Log.getInstance(CollapseModel.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        process();

        return 0;
    }

    protected void process()
    {
        this.cellList = new CellList(CSV);
        log.info(new Object[]{"\tCells detected\t\t[" + this.cellList.size() + "]"});
       
        // initialize models
        this.model = new UCSCRefFlatParser(REFFLAT);
        this.mymodel = new UCSCRefFlatParser(DELTA, MINEVIDENCE);
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary();
        
        try{
            for(SAMSequenceRecord x : dictionnary.getSequences()){
                SAMRecordIterator iter = samReader.query(x.getSequenceName(), 1, x.getSequenceLength(), false);
                //log.info(new Object[]{"\tProcessing ref. " + x.getSequenceName() + "\t[" + mymodel.getMapGenesTranscripts().size() +" genes]"});
                
                while(iter.hasNext()){
                    SAMRecord r = iter.next();
                    pl.record(r);
                    
                    String BC = (String)r.getAttribute(CELLTAG);
                    String U8 = (String)r.getAttribute(UMITAG);
                    String IG = (String)r.getAttribute(GENETAG);
                    String IT = (String)r.getAttribute(ISOFORMTAG);
                    int    RN = ((Integer) r.getAttribute(RNTAG) != null) ? (Integer) r.getAttribute(RNTAG) : 0;
                    
                    LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, true);
                    
                    // never null case if umifound or isobam from isoformMatrix pipeline used
                    // but we filter out some not reliable reads just as in isoformMatrix
                    if(lrr != null && lrr.getMapqv() > 0 && !lrr.getIsChimeria() && !lrr.getIsReversed()){
                        
                        // we have a geneId
                        if(! "undef".equals(IG)){
                            TranscriptRecord tr = new TranscriptRecord(IG, IT);
                                
                            // never seen this gene, create a new List<TranscriptRecord> including IT of lrr
                            if(! mymodel.getMapGenesTranscripts().containsKey(IG)) {
                                if(! "undef".equals(IT))
                                    tr = model.select(IG, IT);
                                
                                List<TranscriptRecord> lst = new ArrayList<TranscriptRecord>();
                                lst.add(tr);
                                mymodel.getMapGenesTranscripts().put(IG, lst);
                            }
                            // we know this gene
                            else {
                                // but we don't know this transcript
                                if(!(mymodel.getMapGenesTranscripts().get(IG)).contains(tr)){
                                    if(! "undef".equals(IT))
                                        tr = model.select(IG, IT);
                                    
                                    ((ArrayList<TranscriptRecord>) mymodel.getMapGenesTranscripts().get(IG)).add(tr);
                                }
                            }
                            
                            // finally add the molecule to the TranscriptRecord
                            mymodel.getMapGenesTranscripts().get(IG).get(mymodel.getMapGenesTranscripts().get(IG).indexOf(tr)).add(lrr);
                        }
                        // intergenic region, new genes ?
                        else{ }
                    }
                }
                iter.close();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        mymodel.collapser();
        mymodel.initialize();
        mymodel.filter(model);
        mymodel.classifier(model);
        
        if(CAGE.exists() && POLYA.exists() && SHORT.exists()){
            log.info(new Object[]{"\tPerform validation using provided CAGE bed, POLYA bed and SHORT read bam files"});
            cage = new BEDParser(CAGE);
            polyA = new BEDParser(POLYA);
            mymodel.validator(cage, polyA, SHORT, cageCo, polyaCo, juncCo);
        }
        else
            log.info(new Object[]{"\tWon't perform validation (please provide CAGE bed, POLYA bed and SHORT read bam files"});
        
        mymodel.statistics();
        
        File TXT = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".d" + DELTA + ".e" + MINEVIDENCE + ".txt");
        File GFF = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".d" + DELTA + ".e" + MINEVIDENCE + ".gff");
        File GFFVALID = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".d" + DELTA + ".e" + MINEVIDENCE + ".final.gff");
        File FAS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".d" + DELTA + ".e" + MINEVIDENCE + ".fas");
        mymodel.exportFiles(TXT,GFF,GFFVALID,FAS);
    }
                    
    public static void main(String[] args) {
        System.exit(new CollapseModel().instanceMain(args));
    }
}
