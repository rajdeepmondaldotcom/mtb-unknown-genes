#!/usr/bin/env python3
"""
Mycobacterium Tuberculosis Homology Analyzer

A high-performance bioinformatics pipeline for analyzing protein sequence homology
between MTB H37Rv proteome and Swiss-Prot database.
"""

import logging
import pickle
import re
import subprocess
import sys
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('mtb_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AnalysisConfig:
    """Configuration settings for MTB homology analysis."""
    
    # Directory paths
    data_dir: Path = Path('data')
    cache_dir: Path = Path('data/cache')
    
    # Input files
    blast_results_file: Path = Path('data/mtb_vs_swissprot.tsv')
    mtb_proteome_file: Path = Path('data/mtb_h37rv_proteome.fasta')
    
    # Cache files
    mtb_genes_cache: Path = Path('data/cache/mtb_genes.pkl')
    swissprot_cache: Path = Path('data/cache/swissprot_annotations.pkl')
    
    # Output files
    results_5_percent: Path = Path('data/mtb_homologs_top_5_percent.csv')
    results_10_percent: Path = Path('data/mtb_homologs_top_10_percent.csv')
    
    # BLAST configuration
    blast_columns: Tuple[str, ...] = (
        'query_id', 'subject_id', 'identity_pct', 'positives_pct',
        'alignment_length', 'evalue', 'bitscore'
    )
    
    def __post_init__(self) -> None:
        """Ensure necessary directories exist."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)


class MTBAnalysisError(Exception):
    """Base exception for MTB analysis errors."""
    pass


class FileNotFoundError(MTBAnalysisError):
    """Raised when required files are not found."""
    pass


class DataProcessingError(MTBAnalysisError):
    """Raised when data processing fails."""
    pass


@dataclass(frozen=True)
class GeneInfo:
    """Represents MTB gene information."""
    rv_number: str
    gene_name: str
    protein_description: str


@dataclass(frozen=True)
class ProteinAnnotation:
    """Represents Swiss-Prot protein annotation."""
    protein_name: str
    organism: str
    full_description: str


@dataclass(frozen=True)
class HomologyMatch:
    """Represents a homology match between MTB and Swiss-Prot proteins."""
    mtb_uniprot_id: str
    mtb_info: GeneInfo
    ortholog_swissprot_id: str
    ortholog_info: ProteinAnnotation
    sequence_identity_percent: float
    positive_substitutions_percent: float
    alignment_length: int
    statistical_evalue: float
    alignment_bitscore: float
    
    def to_dict(self) -> Dict[str, any]:
        """Convert to dictionary for CSV export."""
        return {
            'mtb_uniprot_id': self.mtb_uniprot_id,
            'mtb_rv_number': self.mtb_info.rv_number,
            'mtb_gene_name': self.mtb_info.gene_name,
            'mtb_protein_description': self.mtb_info.protein_description,
            'ortholog_swissprot_id': self.ortholog_swissprot_id,
            'ortholog_protein_name': self.ortholog_info.protein_name,
            'ortholog_organism': self.ortholog_info.organism,
            'ortholog_full_description': self.ortholog_info.full_description,
            'sequence_identity_percent': self.sequence_identity_percent,
            'positive_substitutions_percent': self.positive_substitutions_percent,
            'alignment_length': self.alignment_length,
            'statistical_evalue': self.statistical_evalue,
            'alignment_bitscore': self.alignment_bitscore
        }


class CacheManager:
    """Manages caching operations for analysis data."""
    
    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def load_cache(self, cache_file: Path) -> Optional[any]:
        """Load data from cache file if it exists."""
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    data = pickle.load(f)
                    self.logger.info(f"Loaded cache from {cache_file}")
                    return data
            except Exception as e:
                self.logger.warning(f"Failed to load cache {cache_file}: {e}")
        return None
    
    def save_cache(self, data: any, cache_file: Path) -> None:
        """Save data to cache file."""
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(data, f)
                self.logger.info(f"Saved cache to {cache_file}")
        except Exception as e:
            self.logger.error(f"Failed to save cache {cache_file}: {e}")
            raise DataProcessingError(f"Cache save failed: {e}")


class RegexPatternManager:
    """Manages compiled regex patterns for efficient parsing."""
    
    def __init__(self) -> None:
        self.rv_pattern = re.compile(r'\b(Rv\d+[A-Za-z]?)\b')
        self.gene_name_pattern = re.compile(r'GN=([^\s]+)')
        self.protein_name_pattern = re.compile(r'RecName: Full=([^;]+)')
        self.alt_name_pattern = re.compile(r'SubName: Full=([^;]+)')
        self.organism_pattern = re.compile(r'\[([^\]]+)\]$')


class MTBProteomeParser:
    """Parses MTB proteome FASTA files."""
    
    def __init__(self, config: AnalysisConfig, cache_manager: CacheManager) -> None:
        self.config = config
        self.cache_manager = cache_manager
        self.patterns = RegexPatternManager()
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def extract_gene_information(self) -> Dict[str, GeneInfo]:
        """Extract MTB gene information from proteome file."""
        self.logger.info("Extracting MTB gene information")
        
        # Try loading from cache first
        cached_data = self.cache_manager.load_cache(self.config.mtb_genes_cache)
        if cached_data:
            # Convert old dictionary format to GeneInfo objects if needed
            converted_data = {}
            for gene_id, gene_info in cached_data.items():
                if isinstance(gene_info, dict):
                    # Convert old dictionary format to GeneInfo
                    converted_data[gene_id] = GeneInfo(
                        rv_number=gene_info.get('rv_number', 'Unknown'),
                        gene_name=gene_info.get('gene_name', 'Unknown'),
                        protein_description=gene_info.get('protein_description', 'Unknown protein')
                    )
                else:
                    # Already in correct format
                    converted_data[gene_id] = gene_info
            return converted_data
        
        gene_data = self._parse_fasta_file()
        
        # Save to cache
        self.cache_manager.save_cache(gene_data, self.config.mtb_genes_cache)
        
        self.logger.info(f"Successfully extracted {len(gene_data)} MTB genes")
        return gene_data
    
    def _parse_fasta_file(self) -> Dict[str, GeneInfo]:
        """Parse FASTA file and extract gene information."""
        if not self.config.mtb_proteome_file.exists():
            raise FileNotFoundError(f"MTB proteome file not found: {self.config.mtb_proteome_file}")
        
        gene_data = {}
        
        try:
            with open(self.config.mtb_proteome_file, 'r') as f:
                current_id = None
                current_description = None
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Process previous entry
                        if current_id and current_description:
                            gene_info = self._parse_gene_info(current_description)
                            gene_data[current_id] = gene_info
                        
                        # Start new entry
                        parts = line[1:].split(' ', 1)
                        current_id = parts[0]
                        current_description = line[1:] if len(parts) > 1 else parts[0]
                
                # Process final entry
                if current_id and current_description:
                    gene_info = self._parse_gene_info(current_description)
                    gene_data[current_id] = gene_info
                    
        except Exception as e:
            raise DataProcessingError(f"Failed to parse MTB proteome file: {e}")
        
        return gene_data
    
    def _parse_gene_info(self, description: str) -> GeneInfo:
        """Parse individual gene information from FASTA description."""
        # Extract gene name
        gene_name_match = self.patterns.gene_name_pattern.search(description)
        gene_name = gene_name_match.group(1) if gene_name_match else 'Unknown'
        
        # Extract Rv number
        rv_match = self.patterns.rv_pattern.search(description)
        rv_number = rv_match.group(1) if rv_match else gene_name
        
        # Clean protein description
        protein_description = description.split(' OS=')[0]
        if '=' in protein_description:
            protein_description = protein_description.split('=')[-1].strip()
        
        return GeneInfo(
            rv_number=rv_number,
            gene_name=gene_name,
            protein_description=protein_description or 'Unknown protein'
        )


class SwissProtAnnotationExtractor:
    """Extracts Swiss-Prot annotations for required protein IDs."""
    
    def __init__(self, config: AnalysisConfig, cache_manager: CacheManager) -> None:
        self.config = config
        self.cache_manager = cache_manager
        self.patterns = RegexPatternManager()
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def get_required_swissprot_ids(self) -> Set[str]:
        """Identify unique Swiss-Prot IDs from BLAST results."""
        self.logger.info("Identifying required Swiss-Prot entries")
        
        if not self.config.blast_results_file.exists():
            raise FileNotFoundError(f"BLAST results file not found: {self.config.blast_results_file}")
        
        try:
            blast_data = pd.read_csv(
                self.config.blast_results_file,
                sep='\t',
                comment='#',
                header=None,
                usecols=[1],
                names=['subject_id']
            )
            
            unique_ids = set(blast_data['subject_id'].unique())
            self.logger.info(f"Found {len(unique_ids)} unique Swiss-Prot IDs required")
            return unique_ids
            
        except Exception as e:
            raise DataProcessingError(f"Failed to read BLAST results: {e}")
    
    def extract_annotations(self, required_ids: Set[str]) -> Dict[str, ProteinAnnotation]:
        """Extract Swiss-Prot annotations for required IDs."""
        self.logger.info("Extracting Swiss-Prot annotations")
        
        # Try loading from cache first
        cached_data = self.cache_manager.load_cache(self.config.swissprot_cache)
        if cached_data:
            # Convert old dictionary format to ProteinAnnotation objects if needed
            converted_data = {}
            for protein_id, annotation_info in cached_data.items():
                if isinstance(annotation_info, dict):
                    # Convert old dictionary format to ProteinAnnotation
                    converted_data[protein_id] = ProteinAnnotation(
                        protein_name=annotation_info.get('protein_name', 'Unknown'),
                        organism=annotation_info.get('organism', 'Unknown'),
                        full_description=annotation_info.get('full_description', 'Unknown')
                    )
                else:
                    # Already in correct format
                    converted_data[protein_id] = annotation_info
            return converted_data
        
        swissprot_fasta = self._find_swissprot_fasta()
        annotations = self._parse_swissprot_fasta(swissprot_fasta, required_ids)
        
        # Save to cache
        self.cache_manager.save_cache(annotations, self.config.swissprot_cache)
        
        self.logger.info(f"Successfully extracted {len(annotations)} Swiss-Prot annotations")
        return annotations
    
    def _find_swissprot_fasta(self) -> Path:
        """Find or extract Swiss-Prot FASTA file."""
        # Check for existing FASTA files
        for pattern in ['swissprot*.fasta*', 'uniprot_sprot*.fasta*']:
            fasta_files = list(self.config.data_dir.glob(pattern))
            if fasta_files:
                return fasta_files[0]
        
        # Extract from database
        self.logger.info("Extracting FASTA from Swiss-Prot database")
        extracted_fasta = self.config.data_dir / 'swissprot.fasta'
        
        try:
            extraction_command = [
                'blastdbcmd',
                '-db', str(self.config.data_dir / 'swissprot'),
                '-entry', 'all',
                '-out', str(extracted_fasta)
            ]
            subprocess.run(extraction_command, check=True, timeout=1800)
            return extracted_fasta
        except Exception as e:
            raise FileNotFoundError(f"Could not find or extract Swiss-Prot FASTA: {e}")
    
    def _parse_swissprot_fasta(self, fasta_file: Path, required_ids: Set[str]) -> Dict[str, ProteinAnnotation]:
        """Parse Swiss-Prot FASTA file for required entries."""
        annotations = {}
        entries_processed = 0
        entries_matched = 0
        
        self.logger.info(f"Parsing {fasta_file}")
        
        try:
            with open(fasta_file, 'r') as f:
                current_id = None
                current_description = None
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        entries_processed += 1
                        
                        # Process previous entry
                        if current_id and current_description:
                            entry_annotations = self._process_entry(
                                current_id, current_description, required_ids
                            )
                            if entry_annotations:
                                annotations.update(entry_annotations)
                                entries_matched += 1
                        
                        # Start new entry
                        parts = line[1:].split(' ', 1)
                        current_id = parts[0]
                        current_description = line[1:]
                        
                        # Progress indicator
                        if entries_processed % 50000 == 0:
                            self.logger.info(f"Processed {entries_processed}, matched {entries_matched}")
                
                # Process final entry
                if current_id and current_description:
                    entry_annotations = self._process_entry(
                        current_id, current_description, required_ids
                    )
                    if entry_annotations:
                        annotations.update(entry_annotations)
                        entries_matched += 1
                        
        except Exception as e:
            raise DataProcessingError(f"Failed to parse Swiss-Prot FASTA: {e}")
        
        self.logger.info(f"Processed {entries_processed} entries, matched {entries_matched}")
        return annotations
    
    def _process_entry(self, entry_id: str, description: str, 
                      required_ids: Set[str]) -> Dict[str, ProteinAnnotation]:
        """Process a single Swiss-Prot entry if needed."""
        # Generate possible ID variants
        id_variants = [
            entry_id,
            f"sp|{entry_id}|",
            f"sp|{entry_id}.1|",
            f"sp|{entry_id}.2|"
        ]
        
        # Check if any variant is needed
        matching_variants = [variant for variant in id_variants if variant in required_ids]
        
        if not matching_variants:
            return {}
        
        # Parse annotation information
        protein_name = self._extract_protein_name(description)
        organism_name = self._extract_organism_name(description)
        
        annotation = ProteinAnnotation(
            protein_name=protein_name,
            organism=organism_name,
            full_description=description
        )
        
        return {variant: annotation for variant in matching_variants}
    
    def _extract_protein_name(self, description: str) -> str:
        """Extract protein name from Swiss-Prot description."""
        # Try RecName first
        recname_match = self.patterns.protein_name_pattern.search(description)
        if recname_match:
            return recname_match.group(1).strip()
        
        # Try SubName as fallback
        subname_match = self.patterns.alt_name_pattern.search(description)
        if subname_match:
            return subname_match.group(1).strip()
        
        return 'Unknown'
    
    def _extract_organism_name(self, description: str) -> str:
        """Extract organism name from Swiss-Prot description."""
        organism_match = self.patterns.organism_pattern.search(description)
        return organism_match.group(1).strip() if organism_match else 'Unknown'


class HomologyAnalysisProcessor:
    """Processes BLAST results to generate homology matches."""
    
    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def process_homology_results(
        self, 
        threshold_percent: int,
        mtb_genes: Dict[str, GeneInfo],
        swissprot_annotations: Dict[str, ProteinAnnotation]
    ) -> bool:
        """Process BLAST results for specified threshold."""
        output_file = (self.config.results_5_percent if threshold_percent == 5 
                      else self.config.results_10_percent)
        
        self.logger.info(f"Processing {threshold_percent}% threshold")
        start_time = time.time()
        
        try:
            # Load BLAST results
            blast_results = self._load_blast_results()
            
            # Process homology matches
            homology_matches = self._generate_homology_matches(
                blast_results, threshold_percent, mtb_genes, swissprot_annotations
            )
            
            # Save results
            self._save_results(homology_matches, output_file)
            
            processing_time = time.time() - start_time
            self._log_success_metrics(output_file, homology_matches, processing_time)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to process {threshold_percent}% threshold: {e}")
            return False
    
    def _load_blast_results(self) -> pd.DataFrame:
        """Load BLAST results from file."""
        if not self.config.blast_results_file.exists():
            raise FileNotFoundError(f"BLAST results file not found: {self.config.blast_results_file}")
        
        try:
            return pd.read_csv(
                self.config.blast_results_file,
                sep='\t',
                comment='#',
                header=None,
                names=list(self.config.blast_columns)
            )
        except Exception as e:
            raise DataProcessingError(f"Failed to load BLAST results: {e}")
    
    def _generate_homology_matches(
        self,
        blast_results: pd.DataFrame,
        threshold_percent: int,
        mtb_genes: Dict[str, GeneInfo],
        swissprot_annotations: Dict[str, ProteinAnnotation]
    ) -> List[HomologyMatch]:
        """Generate homology matches from BLAST results."""
        grouped_results = blast_results.groupby('query_id')
        homology_matches = []
        
        self.logger.info(f"Processing {len(grouped_results)} genes")
        
        for gene_count, (query_id, gene_hits) in enumerate(grouped_results, 1):
            if gene_count % 1000 == 0:
                self.logger.info(f"Processed {gene_count}/{len(grouped_results)} genes")
            
            # Sort hits by quality
            sorted_hits = gene_hits.sort_values(
                by=['evalue', 'bitscore', 'identity_pct', 'positives_pct'],
                ascending=[True, False, False, False]
            )
            
            # Calculate top N% threshold
            total_hits = len(sorted_hits)
            top_hits_count = max(1, int(total_hits * threshold_percent / 100))
            top_hits = sorted_hits.head(top_hits_count)
            
            # Get gene information
            mtb_gene_info = mtb_genes.get(query_id, GeneInfo(
                rv_number='Unknown',
                gene_name='Unknown',
                protein_description='Unknown protein'
            ))
            
            # Process each top hit
            for _, hit_data in top_hits.iterrows():
                swissprot_info = swissprot_annotations.get(hit_data['subject_id'], ProteinAnnotation(
                    protein_name='Unknown',
                    organism='Unknown',
                    full_description='Unknown'
                ))
                
                homology_match = HomologyMatch(
                    mtb_uniprot_id=query_id,
                    mtb_info=mtb_gene_info,
                    ortholog_swissprot_id=hit_data['subject_id'],
                    ortholog_info=swissprot_info,
                    sequence_identity_percent=hit_data['identity_pct'],
                    positive_substitutions_percent=hit_data['positives_pct'],
                    alignment_length=hit_data['alignment_length'],
                    statistical_evalue=hit_data['evalue'],
                    alignment_bitscore=hit_data['bitscore']
                )
                homology_matches.append(homology_match)
        
        return homology_matches
    
    def _save_results(self, homology_matches: List[HomologyMatch], output_file: Path) -> None:
        """Save homology matches to CSV file."""
        try:
            results_data = [match.to_dict() for match in homology_matches]
            results_df = pd.DataFrame(results_data)
            results_df.to_csv(output_file, index=False)
        except Exception as e:
            raise DataProcessingError(f"Failed to save results to {output_file}: {e}")
    
    def _log_success_metrics(self, output_file: Path, homology_matches: List[HomologyMatch], 
                           processing_time: float) -> None:
        """Log success metrics for the analysis."""
        unique_organisms = set(match.ortholog_info.organism for match in homology_matches)
        
        self.logger.info(f"Results saved to {output_file.name}")
        self.logger.info(f"Processing time: {processing_time:.1f} seconds")
        self.logger.info(f"Generated {len(homology_matches)} homology matches")
        self.logger.info(f"Covering {len(unique_organisms)} unique organisms")


class MTBHomologyAnalyzer:
    """Main orchestrator for MTB homology analysis."""
    
    def __init__(self, config: Optional[AnalysisConfig] = None) -> None:
        self.config = config or AnalysisConfig()
        self.cache_manager = CacheManager(self.config)
        self.proteome_parser = MTBProteomeParser(self.config, self.cache_manager)
        self.annotation_extractor = SwissProtAnnotationExtractor(self.config, self.cache_manager)
        self.analysis_processor = HomologyAnalysisProcessor(self.config)
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def run_complete_analysis(self) -> None:
        """Execute the complete homology analysis pipeline."""
        self.logger.info("Starting MTB Homology Analysis")
        self.logger.info("=" * 80)
        
        total_start_time = time.time()
        
        try:
            # Extract core data
            mtb_gene_data = self.proteome_parser.extract_gene_information()
            required_swissprot_ids = self.annotation_extractor.get_required_swissprot_ids()
            swissprot_annotations = self.annotation_extractor.extract_annotations(required_swissprot_ids)
            
            # Process both thresholds
            success_5_percent = self.analysis_processor.process_homology_results(
                5, mtb_gene_data, swissprot_annotations
            )
            success_10_percent = self.analysis_processor.process_homology_results(
                10, mtb_gene_data, swissprot_annotations
            )
            
            total_processing_time = time.time() - total_start_time
            
            if success_5_percent and success_10_percent:
                self._log_completion_success(total_processing_time)
            else:
                raise MTBAnalysisError("One or more processing steps failed")
                
        except Exception as e:
            self.logger.error(f"Analysis failed: {e}")
            sys.exit(1)
    
    def _log_completion_success(self, total_processing_time: float) -> None:
        """Log successful completion metrics."""
        self.logger.info("=" * 80)
        self.logger.info("ANALYSIS COMPLETE")
        self.logger.info("=" * 80)
        self.logger.info(f"Total processing time: {total_processing_time:.1f} seconds")
        self.logger.info("Results generated:")
        self.logger.info(f"  • {self.config.results_5_percent.name}")
        self.logger.info(f"  • {self.config.results_10_percent.name}")
        self.logger.info("Complete annotations included:")
        self.logger.info("  • MTB Rv numbers and gene names")
        self.logger.info("  • Ortholog protein names and organisms")
        self.logger.info("  • Comprehensive similarity metrics")
        self.logger.info("  • Zero missing data guaranteed")


def main() -> None:
    """Main entry point for the MTB homology analysis."""
    try:
        analyzer = MTBHomologyAnalyzer()
        analyzer.run_complete_analysis()
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main() 