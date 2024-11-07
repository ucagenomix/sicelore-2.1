package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Bam file cell-by-cell splitter.", oneLineSummary = "Bam file cell-by-cell splitter.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SplitBamPerCell extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "CSV", doc = "The .csv cell barcode file (barcodes.tsv) no header !")
    public File CSV;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";

    public SplitBamPerCell() {
        log = Log.getInstance(SplitBamPerCell.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        HashMap localHashMap = new HashMap();

        SamReader localSamReader = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        SAMFileHeader localSAMFileHeader1 = localSamReader.getFileHeader();
        SAMFileHeader localSAMFileHeader2 = localSAMFileHeader1.clone();
        localSAMFileHeader2.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);

        try {
            BufferedReader br = new BufferedReader(new java.io.FileReader(CSV));
            String str = br.readLine();
            while (str != null) {
                String[] tmp = str.split(",");

                //if (java.util.regex.Pattern.matches(".*-1", tmp[0]))
                //    tmp[0]= tmp[0].replace("-1", "");
                
                System.out.println("-" + tmp[0] + "-");
                
                if (!localHashMap.containsKey(tmp[0]))
                    localHashMap.put(tmp[0], new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader2, true, new File(OUTPUT.getAbsolutePath() + "/" + tmp[0] + ".bam")));
 
                str = br.readLine();
            }
            br.close();
            
            for (Iterator localIterator = localSamReader.iterator(); localIterator.hasNext();) {
                SAMRecord r = (SAMRecord) localIterator.next();
                pl.record(r);
                String name = r.getReadName();
                String cell_barcode = (String) r.getAttribute(CELLTAG);
                
                if(cell_barcode != null){
                    //System.out.println(name);

                    //if (java.util.regex.Pattern.matches(".*-1", cell_barcode))
                    //    cell_barcode = cell_barcode.replace("-1", "");

                    if ((SAMFileWriter) localHashMap.get(cell_barcode) != null)
                        ((SAMFileWriter) localHashMap.get(cell_barcode)).addAlignment(r);
                }
            }
            localSamReader.close();

            Object localObject = localHashMap.keySet();
            Iterator localIterator = ((Set) localObject).iterator();
            while (localIterator.hasNext())
                ((SAMFileWriter) localHashMap.get((String) localIterator.next())).close();
            
        } catch (Exception localException) { localException.printStackTrace(); }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new SplitBamPerCell().instanceMain(paramArrayOfString));
    }
}
