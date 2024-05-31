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
import java.util.Random;

@CommandLineProgramProperties(summary = "Bulk sample to fake single-cell experiment, adding same cell barcode and random UMI", oneLineSummary = "Bulk sample to fake single-cell experiment, adding same cell barcode and random UMI", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class Bulk2FakeSingleCell extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "The barcode cell tag (default=BC)")
    public String CELLTAG = "BC";
    @Argument(shortName = "CELL_BC", doc = "The barcode cell balu (default=AAAAAAAAAAAAAAAA)")
    public String CELL_BC = "AAAAAAAAAAAAAAAA";
    @Argument(shortName = "UMITAG", doc = "The UMI tag (default=U8)")
    public String UMITAG = "U8";

    public Bulk2FakeSingleCell() {
        log = Log.getInstance(Bulk2FakeSingleCell.class);
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
                
                byte[] array = new byte[12]; // length is bounded by 12
                new Random().nextBytes(array);
                String random_umi = this.randomDNAString(12);
                
                r.setAttribute(CELLTAG, CELL_BC);
                r.setAttribute(UMITAG, random_umi);
                
                localSAMFileWriter.addAlignment(r);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }
    public static String randomDNAString(int dnaLength) {
        Random rand = new Random();
        StringBuilder dna = new StringBuilder(dnaLength);

        for (int i = 0; i < dnaLength; i++) {
            dna.append("ACGT".charAt(rand.nextInt(4)));
        }

        return dna.toString();
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new Bulk2FakeSingleCell().instanceMain(paramArrayOfString));
    }
}