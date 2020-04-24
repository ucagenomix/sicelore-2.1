package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.util.*;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.common.utils.ExecuteCmd;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute consensus sequence per molecule.", oneLineSummary = "Compute consensus sequence per molecule.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class ComputeConsensus extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;

    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output .fasta file of consensus sequences")
    public File OUTPUT;
    @Argument(shortName = "T", doc = "The number of threads (default 20)")
    public int nThreads = 20;
    @Argument(shortName = "TMPDIR", doc = "TMPDIR")
    public String TMPDIR = "/share/data/scratch/sicelore/";
    @Argument(shortName = "CELLTAG", doc = "Cell barcode tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "UMI end tag (default=UE)", optional=true)
    public String UMIENDTAG = "UE";
    @Argument(shortName = "POLYAENDTAG", doc = "PolyA end tag (default=PE)", optional=true)
    public String POLYAENDTAG = "PE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;

    public ComputeConsensus() {
        log = Log.getInstance(ComputeConsensus.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        
        String POAPATH = this.findExecutableOnPath("poa");
        String RACONPATH = this.findExecutableOnPath("racon");
        String MINIMAP2PATH = this.findExecutableOnPath("minimap2");
        
        if(POAPATH == null)
            log.info(new Object[]{"Unable to find poa, please add it to your PATH"});
        else if(RACONPATH == null)
            log.info(new Object[]{"Unable to find racon, please add it to your PATH"});
        else if(MINIMAP2PATH == null)
            log.info(new Object[]{"Unable to find minimap2, please add it to your PATH"});
        else{
            Consensus c = new Consensus();
            c.setStaticParams(TMPDIR,POAPATH,RACONPATH,MINIMAP2PATH);
            
            LongreadRecord lrr = new LongreadRecord();
            lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,POLYAENDTAG,USTAG,MAXCLIP,RNTAG);
        
            // load_sequence = true
            // is_gene_mandatory = false
            // is_umi_mandatory = true
            LongreadParser bam = new LongreadParser(INPUT, true, false, true);
            MoleculeDataset dataset = new MoleculeDataset(bam);
            dataset.callConsensus(OUTPUT, nThreads);
        }
        
        return 0;
    }
    
    public static String findExecutableOnPath(String name)
    {
        for (String dirname : System.getenv("PATH").split(File.pathSeparator)) {
            File file = new File(dirname, name);
            if (file.isFile() && file.canExecute()) {
                return file.getAbsolutePath();
            }
        }
        return null;
    }
    
    public static void main(String[] args) {
        System.exit(new ComputeConsensus().instanceMain(args));
    }
}
