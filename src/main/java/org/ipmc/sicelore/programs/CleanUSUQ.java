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

@CommandLineProgramProperties(summary = "Remove US/UQ SAM tags", oneLineSummary = "Remove US/UQ SAM tags", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class CleanUSUQ extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "USTAG", doc = "read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "UQTAG", doc = "read QV tag (default=UQ)", optional=true)
    public String UQTAG = "UQ";

    public CleanUSUQ() {
        log = Log.getInstance(CleanUSUQ.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                r.setAttribute(USTAG, "");
                r.setAttribute(UQTAG, "");
                localSAMFileWriter.addAlignment(r);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new CleanUSUQ().instanceMain(paramArrayOfString));
    }
}
