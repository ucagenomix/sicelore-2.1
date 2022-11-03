package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.util.*;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram; 

@CommandLineProgramProperties(summary = "Select valid cell barcodes ", oneLineSummary = "Select valid cell barcodes", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class SelectValidCellBarcode extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The detected cell barcodes TAB file (BarcodesAssigned.tsv)")
    public File INPUT;
    @Argument(shortName = "O", doc = "The selected cell barcodes file (.csv)")
    public File OUTPUT;
    @Argument(shortName = "MINUMI", doc = "Minimal number of UMI (default=1, no filter)")
    public int MINUMI = 1;
    @Argument(shortName = "ED0ED1RATIO", doc = "ED0 to ED1 ratio (default=1, mean keep cells having more ED0 than ED1 umis (should be the case)", optional=true)
    public double ED0ED1RATIO = 1;
    
    public SelectValidCellBarcode() {
        log = Log.getInstance(SelectValidCellBarcode.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        // Barcode              n Reads with ED<=1 match   ED=0     ED=1
        // CGGACTGTCTTGTACT     140,612                    118,919  21,693
        // AGCCTAAAGGGAAACA     124,745                    88,046   36,699
        // ACCCACTCAGTTCCCT     122,947                    8,499    112,448 --> ED1ED0RATIO > 1 --> potentially filtered by ED1ED0RATIO value
        // ...
        //                      this is the nUMIs
        
        BufferedOutputStream os = null;
        int total_barcodes = 0;
        int kept_barcodes = 0;
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            String line = fichier.readLine();
            line = fichier.readLine();
            while(line != null) {
                line=line.replaceAll(",","");
                String[] tab = line.split("\t");
                
                int total_umi = new Integer(tab[1]).intValue();
                int umi_ed0 = new Integer(tab[2]).intValue();
                int umi_ed1 = new Integer(tab[3]).intValue();
                
                total_barcodes++;
                
                if(total_umi >= MINUMI && new Double(umi_ed0/umi_ed1).doubleValue() >= ED0ED1RATIO){
                    kept_barcodes++;
                    os.write(new String(tab[0] + "\n").getBytes());
                }
                
                line = fichier.readLine();
            }
            fichier.close();
            os.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        log.info(new Object[]{"Total cell barcodes\t\t[" + total_barcodes + "]"});
        log.info(new Object[]{"Valid cell barcodes\t\t[" + kept_barcodes + "]"});
        return 0;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SelectValidCellBarcode().instanceMain(paramArrayOfString));
    }
}
