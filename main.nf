#!/usr/bin/env nextflow

/*
================================================================================
                                  deNovo-circRNA
================================================================================
Started August 2020
--------------------------------------------------------------------------------
Description:
  (To my knowledge) the first circRNA tool to scan RNA-Seq data for circRNAs,
  perform differential gene expression AND conduct circRNA-miRNA predictions
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/BarryDigby/circRNA
 --------------------------------------------------------------------------------
 @Documentation
 Work in progress
--------------------------------------------------------------------------------
*/


/*
================================================================================
                                  Testing Environment
================================================================================
*/

process test{
            echo true
            
            output:
            stdout to out
            
            script:
            """
            python --version
            STAR
            bwa
            python mapsplice.py
            miranda
            bbduk.sh
            CIRCexplorer2
            """
}
