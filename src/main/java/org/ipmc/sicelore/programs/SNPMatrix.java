package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.*;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SNP/editing events detection/quantification (cellBC/UMI count matrix).", oneLineSummary = "SNP/editing event detection/quantification (cellBC/UMI count matrix).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SNPMatrix extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input molecule SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "SNP", doc = "The SNP/editing events .csv file \n#-----\nchromosome,position,strand,name\n3,80692286,-,Gria2_RG\n3,80706912,-,Gria2_QR\n3,80692286|80706912,-,Gria2_RGQR\n#-----")
    public File SNP;
    //@Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    //public File REFFLAT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=snp)")
    public String PREFIX = "snp";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "RN_min", doc = "Minimum read number to take into account the molecule (default=0, means all)")
    public int RN_min = 0;

    public CellList cellList;
    
    public SNPMatrix() {
        log = Log.getInstance(SNPMatrix.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsReadable(SNP);
        //IOUtil.assertFileIsReadable(REFFLAT);

        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t[" + this.cellList.size() + "]"});
        Matrix matrix = new Matrix(this.cellList);
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary() ;
        
        //editing.csv event list to quantify description file
        //chromosome,position,strand,name
        //11,75300577,-,Rpa1
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(SNP));
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                
                //System.out.println(tok[0] + "\t" + dictionnary.getSequence(tok[0]));
                
                if(dictionnary.getSequence(tok[0]) != null){
                    String chromosome = tok[0];
                    boolean strand = ("-".equals(tok[2]))?true:false;
                    String gene = tok[3];
                    int total=0;
                    nb++;
                    
                    String[] pos = { tok[1] };
                    if(java.util.regex.Pattern.matches(".*\\|.*", tok[1]))
                        pos = tok[1].split("\\|");
                        
                    int[] arr = stringToIntArray(pos);
                    SAMRecordIterator iter = samReader.query(chromosome, arr[0], arr[arr.length - 1], false);
                    while (iter.hasNext()) {
                        SAMRecord r = iter.next();

                        if (strand == r.getReadNegativeStrandFlag()) {
                            pl.record(r);
                            LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, false);

                            if (lrr != null && lrr.getRn() >= RN_min) {
                                Longread lr = new Longread(r.getReadName());
                                lr.addRecord(lrr);
                                String readString = r.getReadString();

                                Integer[] posInt = new Integer[arr.length];
                                for (int i = 0; i < arr.length; i++) {
                                    posInt[i] = r.getReadPositionAtReferencePosition(arr[i]);
                                }

                                int min = Collections.min(Arrays.asList(posInt));
                                int max = Collections.max(Arrays.asList(posInt));

                                if (min > 0 && readString.length() > max) {
                                    String[] nuc = new String[posInt.length];
                                    for (int i = 0; i < arr.length; i++) {
                                        nuc[i] = readString.substring(posInt[i] - 1, posInt[i]);
                                    }

                                    if (r.getReadNegativeStrandFlag()) {
                                        for (int i = 0; i < nuc.length; i++) {
                                            nuc[i] = complementBase(nuc[i].charAt(0));
                                        }
                                    }
                                    
                                    total++;
                                    Molecule molecule = new Molecule(lrr.getBarcode(), lrr.getUmi(), lrr.getRn());
                                    molecule.addLongread(lr);
                                    molecule.setGeneId(gene);
                                    molecule.setTranscriptId(chromosome + ":" + String.join("|", pos) + ".." + String.join("", nuc));
                                    matrix.addMolecule(molecule);
                                }
                            }
                        }
                    }
                    iter.close();
                    log.info(new Object[]{"processing...\t\t" + line + "\t[" + total + " hits]\t"});
                }
                line = fichier.readLine();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }

        if(matrix.getMatrice().size() > 0){
            File MATRIX = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_matrix.txt");
            File METRICS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_metrics.txt");
            File molinfos = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_molinfos.txt");
            matrix.writeIsoformMatrix(MATRIX, METRICS, molinfos, null);
        }
        else
            log.info(new Object[]{"end of processing...\t\tnothing has been detected, check your input parameters, no output files generated\t"});
        
        return 0;
    }
    
    public int[] stringToIntArray(String[] s)
    {
        int size = s.length;
        int [] arr = new int [size];
        for(int i=0; i<size; i++) 
            arr[i] = Integer.parseInt(s[i]);
        
        Arrays.sort(arr);
        return arr;
    }
    
    protected static String complementBase(char base)
    {
        String cc = "";
        if (base == 'A') cc = "T";
        if (base == 'C') cc = "G";
        if (base == 'G') cc = "C";
        if (base == 'T') cc = "A";
	
        return cc;
    }
    
    public String[] complementBaseArray(String[] nuc)
    {
        int size = nuc.length;
        String[] comp = new String[size];
        for(int i=0; i<size; i++) 
            comp[i] = complementBase(nuc[i].charAt(0));
        
        return comp;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SNPMatrix().instanceMain(paramArrayOfString));
    }
}
