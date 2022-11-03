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

@CommandLineProgramProperties(summary = "Histogram of read length distribution.", oneLineSummary = "Histogram of read length distribution.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature 
public class HistoReadLength extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output file")
    public File OUTPUT;
    @Argument(shortName = "BIN", doc = "Bin window (default=50)")
    public int BIN=50;
    @Argument(shortName = "MIN", doc = "Minimum length (default=0)")
    public int MIN=0;
    @Argument(shortName = "MAX", doc = "Minimum length (default=5000)")
    public int MAX=5000;
    
    public HashMap<Integer, Integer> histo;
    
    public HistoReadLength() {
        log = Log.getInstance(HistoReadLength.class);
        pl = new ProgressLogger(log);
        this.histo = new HashMap<Integer, Integer>();
    }

    protected int doWork()
    {
        DataOutputStream os = null;
        
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        HashSet<String> readNames = new HashSet<String>();
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        int records = 0;
        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                records++;
                
                String readName = r.getReadName();
                String cigar = r.getCigarString();
                int length = r.getReadLength();
                
                if(!readNames.contains(readName)){
                
                    readNames.add(readName);
                    int c = (int)(length/BIN);

                    //System.out.println(length +"\n" + cigar);

                    if(this.histo.containsKey(c))
                        this.histo.put(c,this.histo.get(c) + 1);
                    else
                        this.histo.put(c,1);
                }
            }
            inputSam.close();
        } catch (Exception e) { e.printStackTrace(); } finally { try {inputSam.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("length\tnumber\n");
            
            int min = MIN/BIN;
            int max = MAX/BIN;
            
            for(int i=min; i<=max; i++){
                int nn = (this.histo.get(i) != null)?this.histo.get(i):0;
                os.writeBytes((i*BIN)+"\t"+ nn + "\n");
            }

            os.close();
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new HistoReadLength().instanceMain(paramArrayOfString));
    }
}
