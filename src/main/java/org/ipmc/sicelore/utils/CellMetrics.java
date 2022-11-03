package org.ipmc.sicelore.utils;

import java.util.HashSet;
import java.util.Set;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class CellMetrics
{
    private int isoform_known_count=0;
    private int isoform_undef_count=0;
    private int nb_reads=0; 
    private Set<String> genes = new HashSet<String>();
    private int nb_umis=0;
    
    public CellMetrics(){
        this.genes = new HashSet<String>();
    }
    
    public void addCount(String geneId, String transcriptId, int nb_reads){
        
        this.nb_umis++;
        this.nb_reads+=nb_reads;
        this.genes.add(geneId);
        
        if("undef".equals(transcriptId))
            this.isoform_undef_count++;
        else
            this.isoform_known_count++;
    }
    public int getIsoform_known_count() {
        return isoform_known_count;
    }
    public int getIsoform_undef_count() {
        return isoform_undef_count;
    }
    public int getNb_reads() {
        return nb_reads;
    }
    public int getNb_genes() {
        return this.genes.size();
    }
    public int getNb_umis() {
        return nb_umis;
    }

}
