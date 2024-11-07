package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand - june 2024
 * 
 */
import gnu.trove.THashMap;
import gnu.trove.THashSet;
import htsjdk.samtools.SAMRecord;
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Statistics about parse biosciences priming polyA versus Random hexamer", oneLineSummary = "Statistics about parse biosciences priming polyA versus Random hexamer", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ParseTR extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The parse biosciences cell barcodes .csv file (ex: bc_data_n192_v4.csv) giving type and relationship for BC1")
    public File CSV;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "CELLTAG_BC", doc = "Cell tag BC1_BC2_BC3 (default=CR)", optional=true)
    public String CELLTAG_BC = "CR";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=CB)", optional=true)
    public String CELLTAG = "CB";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=pN)", optional=true)
    public String UMITAG = "pN";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=GN)", optional=true)
    public String GENETAG = "GN";
    @Argument(shortName = "XF", doc = "Gene name location, for instance INTRONIC OR INTERGENIC(default=XF)", optional=true)
    public String XF = "XF";
    @Argument(shortName = "SAMPLE", doc = "Sample name tag (default=pS)", optional=true)
    public String SAMPLE = "pS";

    public boolean DEBUG = false;
    
    private ProgressLogger pl;
    private final Log log;

    public ParseTR() {
        
        log = Log.getInstance(ParseTR.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log, 1000000, "\tProcessed\t", "Records");
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        System.out.println(htsjdk.samtools.Defaults.allDefaults());
        process();
        return 0;
    }

    protected void process()
    {
        // bci,sequence,uid,well,type
        // 1,CATTCCTA,pbs_1239,A1,T
        // 2,CTTCATCA,pbs_1205,A2,T
        // 3,CCTATATC,pbs_1247,A3,T
        
        THashMap<String, THashMap<String, Integer>> gene_matrix =  new THashMap<String, THashMap<String, Integer>>();
        THashMap<String, THashMap<String, Integer>> cell_matrix =  new THashMap<String, THashMap<String, Integer>>();
        
        Hashtable bc2type = new Hashtable();
        Hashtable bc2cond = new Hashtable();
        
        THashMap mapCell = null;
        THashMap mapGene = null;
        THashSet umiSet = null;
        
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while(line != null) {
                String[] infos = line.split(",");
                bc2type.put(infos[1], infos[4]);
                
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        String[] keys = {"CODING_T", "CODING_R", "UTR_T", "UTR_R", "INTRONIC_T", "INTRONIC_R", "INTERGENIC_T", "INTERGENIC_R"};
        
        // parse BAM
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        try{
            for(SAMRecord r : inputSam){
                pl.record(r);
                String name = r.getReadName();
                String bc1_bc2_bc3 = (String) r.getAttribute(CELLTAG_BC);
                String cell_barcode = (String) r.getAttribute(CELLTAG);
                String umi = (String) r.getAttribute(UMITAG);
                String gene = (String) r.getAttribute(GENETAG);
                String sample = (String) r.getAttribute(SAMPLE);
                String where = (String) r.getAttribute(XF);
                
                bc2cond.put(cell_barcode, sample);
                
                // do something please !
                String[] barcodes = bc1_bc2_bc3.split("_");
                String priming = (String)bc2type.get(barcodes[0]);
                
                String key = new String(where + "_" + priming);
                
                // gene information
                if ((mapGene = (THashMap) gene_matrix.get(gene)) != null) {
                    if ((umiSet = (THashSet) mapGene.get(key)) != null) {
                        umiSet.add(umi);
                    }
                    else {
                        umiSet = new THashSet();
                        umiSet.add(umi);
                        mapGene.put(key, umiSet);
                    }
                }
                else {
                    gene_matrix.put(gene, new THashMap());
                    umiSet = new THashSet();
                    umiSet.add(umi);
                    ((THashMap) gene_matrix.get(gene)).put(key, umiSet);
                }
                // cell information
                if ((mapCell = (THashMap) cell_matrix.get(cell_barcode)) != null) {
                    if ((umiSet = (THashSet) mapCell.get(key)) != null) {
                        umiSet.add(umi);
                    }
                    else {
                        umiSet = new THashSet();
                        umiSet.add(umi);
                        mapCell.put(key, umiSet);
                    }
                }
                else {
                    cell_matrix.put(cell_barcode, new THashMap());
                    umiSet = new THashSet();
                    umiSet.add(umi);
                    ((THashMap) cell_matrix.get(cell_barcode)).put(key, umiSet);
                }
            }
            inputSam.close();
            
            BufferedOutputStream os = null;
            File GENESTATS = new File(OUTDIR.getAbsolutePath() + "/gene_stats.txt");
            os = new BufferedOutputStream(new java.io.FileOutputStream(GENESTATS));
            os.write(new String("gene").getBytes());
            
            for(int i=0; i<keys.length; i++)
                os.write(new String("\t"+ keys[i]).getBytes());
            
            os.write(new String("\n").getBytes());
            
            for(String gene : gene_matrix.keySet()){
                os.write(new String(gene).getBytes());
                
                for(int i=0; i<keys.length; i++){
                    if ((umiSet = (THashSet) ((THashMap) gene_matrix.get(gene)).get(keys[i])) != null)
                        os.write(new String("\t" + umiSet.size()).getBytes());
                    else
                        os.write(new String("\t0").getBytes());
                }
                os.write(new String("\n").getBytes());
            }
            os.close();
            
            File CELLSTATS = new File(OUTDIR.getAbsolutePath() + "/cell_stats.txt");
            os = new BufferedOutputStream(new java.io.FileOutputStream(CELLSTATS));
            os.write(new String("cell\tcondition").getBytes());
            
            for(int i=0; i<keys.length; i++)
                os.write(new String("\t"+ keys[i]).getBytes());
            
            os.write(new String("\n").getBytes());
            
            for(String cell : cell_matrix.keySet()){
                os.write(new String(cell + "\t" + (String)bc2cond.get(cell)).getBytes());
                
                for(int i=0; i<keys.length; i++){
                    if ((umiSet = (THashSet) ((THashMap) cell_matrix.get(cell)).get(keys[i])) != null)
                        os.write(new String("\t" + umiSet.size()).getBytes());
                    else
                        os.write(new String("\t0").getBytes());
                }
                os.write(new String("\n").getBytes());
            }
            os.close();
            
        } catch (Exception e) { e.printStackTrace(); }
    }

    public static void main(String[] args) {
        System.exit(new ParseTR().instanceMain(args));
    }
}
