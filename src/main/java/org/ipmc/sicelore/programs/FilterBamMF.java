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
import org.ipmc.sicelore.utils.CellList;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Bam filtering for mapqv=0", oneLineSummary = "Bam filtering for mapqv=0", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterBamMF extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;

    public CellList cellList;
    
    public FilterBamMF() {
        log = Log.getInstance(FilterBamMF.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        this.cellList = new CellList(CSV);
        log.info(new Object[]{"\tCells detected\t\t[" + this.cellList.size() + "]"});
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader SAMFileHeader = samReader.getFileHeader();
        SAMFileWriter SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(SAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : samReader) {
                pl.record(r);
                
                String cell_id = (String)r.getAttribute("BC");
                
                if(this.cellList.contains(cell_id)){
                    String umi = (String)r.getAttribute("U8");
                    cell_id = cell_id + "-1";
                    r.setAttribute("CB", cell_id);
                    r.setAttribute("UB", umi);
                    
                    String name = r.getReadName();
                    String[] tab = name.split("=");
                    r.setReadName(tab[0]);
                    
                    SAMFileWriter.addAlignment(r);
                }
            }
            samReader.close();
            SAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new FilterBamMF().instanceMain(paramArrayOfString));
    }
}
