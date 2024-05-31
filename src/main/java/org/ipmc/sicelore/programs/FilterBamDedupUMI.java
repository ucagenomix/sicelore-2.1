package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashSet;
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Bam filtering to keep one read per molecule", oneLineSummary = "Bam filtering to keep one read per molecule", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterBamDedupUMI extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;

    
    public FilterBamDedupUMI() {
        log = Log.getInstance(FilterBamMF.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        THashSet<String> map = new THashSet<String>();
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader SAMFileHeader = samReader.getFileHeader();
        SAMFileWriter SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(SAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : samReader) {
                pl.record(r);
                String U8 = (String)r.getAttribute("U8");
                String BC = (String)r.getAttribute("BC");
                
                if(! map.contains(BC+"_"+U8)){
                    map.add(BC+"_"+U8);
                    SAMFileWriter.addAlignment(r);
                }
            }
            samReader.close();
            SAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new FilterBamDedupUMI().instanceMain(paramArrayOfString));
    }
}
