package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import java.util.regex.Pattern;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Remove duplicate molecule from Fasta/Fastq file.", oneLineSummary = "Remove duplicate molecule from Fasta/Fastq file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class DeduplicateMolecule extends CommandLineProgram {

    private final Log log;

    @Argument(shortName = "I", doc = "The input .fasta or .fastq file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output deduplicate .fasta or .fastq file")
    public File OUTPUT;
    @Argument(shortName = "TSO", doc = "TSO sequence (default=AACGCAGAGTACATGG)")
    public String TSO = "AACGCAGAGTACATGG";
    @Argument(shortName = "MAXPOS", doc = "Number of nucleotides to search for TSO sequence if still in consensus (default=100 first nt.)")
    public int MAXPOS = 100;
 
    public DeduplicateMolecule() {
        log = Log.getInstance(DeduplicateMolecule.class);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        if(Pattern.matches(".*\\.fq", INPUT.getName().toLowerCase()) || Pattern.matches(".*\\.fastq", INPUT.getName().toLowerCase()))
            processFastq();
        else
            processFasta();
        
        return 0;
    }
    
    protected void processFasta()
    {
        int count=0;
        int tsocut=0;
        String line = null;
        BufferedOutputStream os = null;
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"loadFasta\tSTART..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^>.*", line)){
                    String seq = fichier.readLine();
                
                    if(!"null".equals(seq)){
                        count++;
                        
                        int region = (seq.length() < MAXPOS)?seq.length():MAXPOS;
                        
                        // remove remaining TSO if found in the first 100nt.
                        // TSO = AAGCAGTGGTATCAACGCAGAGTACATGG
                        int xx = 0;
                        if((xx = seq.substring(0,region).indexOf(TSO)) > 0){
                            seq = seq.substring(xx+TSO.length());
                            tsocut++;
                        }
                        else if((xx = seq.substring(0,region).indexOf(TSO.substring(0,10))) > 0){
                            seq = seq.substring(xx+16);
                            tsocut++;
                        }
                        else if((xx = seq.substring(0,region).indexOf(TSO.substring(6,16))) > 0){
                            seq = seq.substring(xx+10);
                            tsocut++;
                        }
                        
                        line = line.replace(">", "");
                        line = line.replace("\\|", "-");
                        String[] ids = line.split("-");
                        String key = ids[0]+ids[1];
                        int rn = new Integer(ids[2]).intValue();

                        if(! map.containsKey(key))
                            map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                        else{
                            
                            //log.info(new Object[]{"already seen\t" + ids[0]+"\t"+ids[1]+"\t"+seq.length()+"\t"+rn+"\t"+map.get(key).getConsensus().length()+"\t"+map.get(key).getRn()});
                            
                            if(map.get(key).getRn() < rn)
                                map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                            else if(map.get(key).getRn() == rn){
                                if(map.get(key).getConsensus().length() < seq.length())
                                    map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                            }
                            //log.info(new Object[]{"kept\t\t" +map.get(key).getConsensus().length()+"\t"+map.get(key).getRn()}); 
                        }
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }

        log.info(new Object[]{"loadFasta\tEND..."});
        log.info(new Object[]{"loadFasta\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFasta tso\t"+tsocut});
        log.info(new Object[]{"loadFasta\t" + map.size() + " molecules"});
        log.info(new Object[]{"loadFasta\tEND..."});
        
        log.info(new Object[]{"writeFasta\tSTART..."});
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                Molecule m = (Molecule) map.get(key);
                os.write(new String(">" + m.getBarcode() + "-" + m.getUmi() + "-" + m.getRn() + "\n" + m.getConsensus() + "\n").getBytes());
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFasta\tEND..."});
    }
    
    protected void processFastq()
    {
        int count=0;
        int tsocut=0;
        String line = null;
        BufferedOutputStream os = null;
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"loadFastQ\tSTART..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^@.*", line)){
                    String seq = fichier.readLine();
                    String tmp = fichier.readLine();
                    String qual = fichier.readLine();
                    
                    if(!"null".equals(seq)){
                        count++;
                        
                        int region = (seq.length() < MAXPOS)?seq.length():MAXPOS;
                        
                        // remove remaining TSO if found in the first 100nt.
                        // TSO = AAGCAGTGGTATCAACGCAGAGTACATGG
                        int xx = 0;
                        if((xx = seq.substring(0,region).indexOf(TSO)) > 0){
                            seq = seq.substring(xx+TSO.length());
                            qual = qual.substring(xx+TSO.length());
                            tsocut++;
                        }
                        else if((xx = seq.substring(0,region).indexOf(TSO.substring(0,10))) > 0){
                            seq = seq.substring(xx+16);
                            qual = qual.substring(xx+16);
                            tsocut++;
                        }
                        else if((xx = seq.substring(0,region).indexOf(TSO.substring(6,16))) > 0){
                            seq = seq.substring(xx+10);
                            qual = qual.substring(xx+10);
                            tsocut++;
                        }
                        
                        line = line.replace("@", "");
                        line = line.replace("\\|", "-");
                        String[] ids = line.split("-");
                        String key = ids[0]+ids[1];
                        int rn = new Integer(ids[2]).intValue();

                        if(! map.containsKey(key))
                            map.put(key, new Molecule(ids[0], ids[1], seq, qual, rn));
                        else{
                            if(map.get(key).getRn() < rn)
                                map.put(key, new Molecule(ids[0], ids[1], seq, qual, rn));
                            else if(map.get(key).getRn() == rn){
                                if(map.get(key).getConsensus().length() < seq.length())
                                    map.put(key, new Molecule(ids[0], ids[1], seq, qual, rn));
                            }
                        }
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }

        log.info(new Object[]{"loadFastQ\tEND..."});
        log.info(new Object[]{"loadFastQ\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFastQ tso\t"+tsocut});
        log.info(new Object[]{"loadFastQ\t" + map.size() + " molecules"});
        log.info(new Object[]{"loadFastQ\tEND..."});
        
        log.info(new Object[]{"writeFastQ\tSTART..."});
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                Molecule m = (Molecule) map.get(key);
                os.write(new String("@" + m.getBarcode() + "-" + m.getUmi() + "-" + m.getRn() + "\n" + m.getConsensus() + "\n+\n" + m.getConsensusQV() + "\n").getBytes());
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFastQ\tEND..."});
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new DeduplicateMolecule().instanceMain(paramArrayOfString));
    }
}
