package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Export fastq file from bam file", oneLineSummary = "Export fastq file from bam file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class Bam2Fastq extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output FASTQ file")
    public File OUTPUT;
    
    protected static String USTAG = "US";
    protected static String QSTAG = "QS";
    
    public Bam2Fastq() {
        log = Log.getInstance(Bam2Fastq.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        BufferedOutputStream os = null;
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String seq = (String)r.getAttribute(USTAG);
                String qual = (String)r.getAttribute(QSTAG);
                String[] tmp = name.split("_");
                os.write(new String("@" + tmp[0] + "\n").getBytes());
                os.write(new String(seq + "\n+\n").getBytes());
                os.write(new String(qual + "\n").getBytes());
            }
            localSamReader.close();
            os.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new Bam2Fastq().instanceMain(paramArrayOfString));
    }
}
