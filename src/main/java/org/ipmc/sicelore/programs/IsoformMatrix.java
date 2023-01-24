package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*; 
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

@CommandLineProgramProperties(summary = "Isoform and gene level expression matrices production.", oneLineSummary = "Isoform and gene level expression matrices production.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class IsoformMatrix extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=sicelore)")
    public String PREFIX = "sicelore";
    @Argument(shortName = "ISOBAM", doc = "Wether or not to produce a bam file having geneId (IG) and TranscriptId (IT) SAM flags (default=false)", optional=true)
    public boolean ISOBAM = false;
    @Argument(shortName = "METHOD", doc = "Isoform assignment method (default=STRICT)")
    public String METHOD = "STRICT";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=GE)", optional=true)
    public String GENETAG = "GE";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "POLYASTARTTAG", doc = "PolyA start tag (default=PS)", optional=true)
    public String POLYASTARTTAG = "PS";
    @Argument(shortName = "CDNATAG", doc = "cDNA sequence tag (default=CS)", optional=true)
    public String CDNATAG = "CS";
    @Argument(shortName = "USTAG", doc = "read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "AMBIGUOUS_ASSIGN", doc = "Wether or not to assign the UMI to an isoform if ambiguous (default=false)", optional=true)
    public boolean AMBIGUOUS_ASSIGN = false;
    @Argument(shortName = "MAPQV0", doc = "Wether or not to keep mapqv=0 SAM records (default=false)", optional=true)
    public boolean MAPQV0 = false;
    @Argument(shortName = "TOBULK", doc = "Wether or not to provide a bulk quantification (default=false)", optional=true)
    public boolean TOBULK = false;
    
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;

    public IsoformMatrix() {
        log = Log.getInstance(IsoformMatrix.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        LongreadRecord lrr = new LongreadRecord();
        lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,POLYASTARTTAG,USTAG,CDNATAG,MAXCLIP,RNTAG);
        
        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"\tCells detected\t\t[" + this.cellList.size() + "]"});
        
        if(!"STRICT".equals(METHOD)){
            log.info(new Object[]{"\tIsoform method: [" + METHOD + "] not allowed, only STRICT method allowed (SCORE disabled)"});
            return 0;
        }

        process();

        return 0;
    }

    protected void process()
    {
        LongreadParser bam = new LongreadParser(INPUT, MAPQV0, false, true, true);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        dataset.initModel(REFFLAT);
        dataset.setIsoforms(DELTA, METHOD, AMBIGUOUS_ASSIGN);
        
        File ISOMATRIX   = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isomatrix.txt");
        File ISOMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isometrics.txt");
        File JUNCMATRIX  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_juncmatrix.txt");
        File JUNCMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_juncmetrics.txt");
        File GENEMATRIX  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genematrix.txt");
        File GENEMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genemetrics.txt");
        File CELLMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_cellmetrics.txt");
        File BULKGENE    = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_bulkgene.txt");
        File BULKISO     = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_bulkiso.txt");
        File MOLINFOS    = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_molinfos.txt");
        File LOGS        = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".log");
        File HTML        = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".html");
        File outISOBAM   = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isobam.bam");

        Matrix matrix = dataset.produceMatrix(this.cellList);
        
        log.info(new Object[]{"\twriteIsoformMatrix\t[start]"});
        matrix.writeIsoformMatrix(ISOMATRIX, ISOMETRICS, MOLINFOS, dataset.getModel());
        log.info(new Object[]{"\twriteGeneMatrix\t\t[start]"});
        matrix.writeGeneMatrix(GENEMATRIX, GENEMETRICS);
        log.info(new Object[]{"\twriteCellMetrics\t[start]"});
        matrix.writeCellMetrics(CELLMETRICS);
        log.info(new Object[]{"\twriteJunctionMatrix\t[start]"});
        matrix.writeJunctionMatrix(JUNCMATRIX, JUNCMETRICS);
        
        if(TOBULK)
            matrix.writeBulk(BULKGENE, BULKISO, dataset.getModel());
        
        log.info(new Object[]{"\tMatrix cells size\t[" + matrix.getCellMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix genes size\t[" + matrix.getGeneMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix junctions size\t[" + matrix.getMatriceJunction().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms size\t[" + matrix.getMatrice().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms counts\t[" + matrix.getTotal_count() + "]"});
        log.info(new Object[]{"\tMatrix isoforms define\t[" + matrix.getTotal_isoform_def() + "]"});
        log.info(new Object[]{"\tMatrix isoforms undefine[" + matrix.getTotal_isoform_undef() + "]"});
        
        writeLOGS(LOGS, bam, dataset, matrix);
        writeHTML(HTML, bam, dataset, matrix);
        
        if(ISOBAM){
            log.info(new Object[]{"\tProducing ISOBAM\t[true]"});
            SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
            htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
            samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outISOBAM);
            try{
                for(SAMRecord r : samReader){
                    pl.record(r);
                    String isokey = (String)r.getAttribute(CELLTAG)+":"+(String)r.getAttribute(UMITAG);
                    
                    Molecule molecule = dataset.getMolecule(isokey);
                    if(molecule!=null){
                        r.setAttribute("IG", (molecule.getGeneId()!=null)?molecule.getGeneId():"undef");
                        r.setAttribute("IT", (molecule.getTranscriptId()!=null)?molecule.getTranscriptId():"undef");
                    }
                    else{
                        r.setAttribute("IG", "undef");
                        r.setAttribute("IT", "undef");
                    }
                    samFileWriter.addAlignment(r);
                }
                samReader.close();
                samFileWriter.close();
            } catch (Exception e) { e.printStackTrace(); }
        }
    }
    
    public void writeLOGS(File LOGS, LongreadParser bam, MoleculeDataset dataset, Matrix matrix) {
        BufferedOutputStream os = null;
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(LOGS));
        
            os.write(new String("IsoformMatrix INPUT," + INPUT + "\n").getBytes());
            os.write(new String("IsoformMatrix REFFLAT," + REFFLAT + "\n").getBytes());
            os.write(new String("IsoformMatrix CSV," + CSV + "\n").getBytes());
            os.write(new String("IsoformMatrix DELTA," + DELTA + "\n").getBytes());
            os.write(new String("IsoformMatrix METHOD," + METHOD + "\n").getBytes());
            os.write(new String("IsoformMatrix CELLTAG," + CELLTAG + "\n").getBytes());
            os.write(new String("IsoformMatrix UMITAG," + UMITAG + "\n").getBytes());
            os.write(new String("IsoformMatrix GENETAG," + GENETAG + "\n").getBytes());
            os.write(new String("IsoformMatrix MAXCLIP," + MAXCLIP + "\n").getBytes());
            os.write(new String("IsoformMatrix AMBIGUOUS_ASSIGN," + AMBIGUOUS_ASSIGN + "\n").getBytes());
            os.write(new String("IsoformMatrix MAPQV0," + MAPQV0 + "\n").getBytes());
            
            os.write(new String("Total SAMrecords," + bam.getTotal_records() + "\n").getBytes());
            os.write(new String("SAMrecords valid," + bam.getValid_records() + "\n").getBytes());
            os.write(new String("SAMrecords unvalid," + bam.getUnvalid_records() + "\n").getBytes());
            os.write(new String("SAMrecords mapqv=0," + bam.getMapqv0_records() + "\n").getBytes());
            os.write(new String("SAMrecords no gene," + bam.getGene_unset() + "\n").getBytes());
            os.write(new String("SAMrecords no UMI," + bam.getUmi_unset() + "\n").getBytes());
            os.write(new String("SAMrecords chimeria," + bam.getChimeria_records() + "\n").getBytes());
            os.write(new String("Total reads," + bam.getMapLongreads().size() + "\n").getBytes());
            os.write(new String("Total reads multiSAM," + bam.getMultiRec().size() + "\n").getBytes());

            os.write(new String("Total molecules," + dataset.getMapMolecules().size() + "\n").getBytes());
            os.write(new String("Total molecule reads," + dataset.getTotalReads() + "\n").getBytes());
            os.write(new String("Total molecule multiIG," + dataset.getMultiIG() + "\n").getBytes());

            os.write(new String("UCSCRefFlatParser genes," + dataset.getModel().getMapGenesTranscripts().size() + "\n").getBytes());
            os.write(new String("UCSCRefFlatParser transcripts," + dataset.getModel().getTranscriptsNumber() + "\n").getBytes());
            
            os.write(new String("SetIsoforms monoexon," + dataset.getMonoexon() + "\n").getBytes());
            os.write(new String("SetIsoforms no match," + dataset.getNomatch()  + "\n").getBytes());
            os.write(new String("SetIsoforms one match," + dataset.getOnematch() + "\n").getBytes());
            os.write(new String("SetIsoforms ambiguous," + dataset.getAmbiguous() + "\n").getBytes());
        
            os.write(new String("Matrix cells size," + matrix.getCellMetrics().size() + "\n").getBytes());
            os.write(new String("Matrix genes size," + matrix.getGeneMetrics().size() + "\n").getBytes());
            os.write(new String("Matrix junctions size," + matrix.getMatriceJunction().size() + "\n").getBytes());
            os.write(new String("Matrix isoforms size," + matrix.getMatrice().size() + "\n").getBytes());
            os.write(new String("Matrix isoforms counts," + matrix.getTotal_count() + "\n").getBytes());
            os.write(new String("Matrix isoforms define," + matrix.getTotal_isoform_def() + "\n").getBytes());
            os.write(new String("Matrix isoforms undefined," + matrix.getTotal_isoform_undef() + "\n").getBytes());
           
            
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close();  } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public void writeHTML(File HTML, LongreadParser bam, MoleculeDataset dataset, Matrix matrix)
    {
        BufferedOutputStream os = null;
        
        DefaultCategoryDataset result = new DefaultCategoryDataset();

        result.addValue(bam.getValid_records(), "Valid", "SAM records");
        result.addValue(bam.getMapqv0_records(), "mapqv=0", "SAM records");
        result.addValue(bam.getGene_unset(), "no gene", "SAM records");
        result.addValue(bam.getUmi_unset(), "no UMI", "SAM records");
        result.addValue(bam.getChimeria_records(), "chimeria", "SAM records");
        
        result.addValue(dataset.getOnematch(), "Match", "Molecules");
        result.addValue(dataset.getMonoexon(), "Mono exonic", "Molecules");
        result.addValue(dataset.getAmbiguous(), "Ambiguous", "Molecules");
        result.addValue(dataset.getNomatch(), "No match", "Molecules");
        
        result.addValue(matrix.getTotal_isoform_def(), "defined", "Isoforms");
        result.addValue(matrix.getTotal_isoform_undef(), "undefined", "Isoforms");
        
        final JFreeChart chart = ChartFactory.createStackedBarChart(
            "IsoformMatrix statistics",  // chart title
            "Category",                  // domain axis label
            "#",                     // range axis label
            result,                     // data
            PlotOrientation.HORIZONTAL,    // the plot orientation
            true,                        // legend
            true,                        // tooltips
            false                        // urls
        );
        
        chart.setBackgroundPaint(java.awt.Color.white);
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(HTML));

            ChartUtils.saveChartAsPNG(new File(OUTDIR.getAbsolutePath() + "/histogram.png"), chart, 1000, 1000);
            
            os = new BufferedOutputStream(new FileOutputStream(HTML));
            final PrintWriter writer = new PrintWriter(os);
            writer.println("<HTML>");
            writer.println("<HEAD><TITLE>JFreeChart Image Map Demo 2</TITLE></HEAD>");
            writer.println("<BODY>");
            writer.println("<IMG SRC=\"histogram.png\" WIDTH=\"600\" HEIGHT=\"400\" BORDER=\"0\" USEMAP=\"#chart\">");
            writer.println("</BODY>");
            writer.println("</HTML>");
            writer.close();
            
               os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close();  } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public static void main(String[] args) {
        System.exit(new IsoformMatrix().instanceMain(args));
    }
}
