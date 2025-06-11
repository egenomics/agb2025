#!/usr/bin/env python3
"""
Genus-Level Pipeline Benchmarking Script v15
Standard approach for microbiome pipeline validation using taxonomic comparison
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, mean_absolute_error
import warnings
warnings.filterwarnings('ignore')

class GenusLevelBenchmarker:
    def __init__(self, pipeline_output_dir, ground_truth_dir, asv_metadata_file, results_dir="benchmarking_results"):
        self.pipeline_output_dir = Path(pipeline_output_dir)
        self.ground_truth_dir = Path(ground_truth_dir)
        self.asv_metadata_file = Path(asv_metadata_file)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        
        print(" Initializing Genus-Level Pipeline Benchmarker v15")
        print(" Standard microbiome validation approach using taxonomic comparison")
        
        # Load components
        self.asv_mapping = self.load_asv_mapping()
        self.feature_table = self.load_feature_table()
        self.taxonomy = self.load_taxonomy()
        self.ground_truth = self.load_ground_truth_data()
        
        # Results storage
        self.comparison_data = []
        self.metrics_df = None
        
    def load_asv_mapping(self):
        """Load ASV to genus mapping from metadata"""
        try:
            print(f"\n Loading ASV metadata: {self.asv_metadata_file}")
            df = pd.read_csv(self.asv_metadata_file, sep='\t')
            
            # Standardize column names
            if 'ASV' in df.columns and 'Genus' in df.columns:
                df = df.rename(columns={'ASV': 'ASV_ID', 'Genus': 'Genus'})
                print(f" Loaded {len(df)} ASV-to-genus mappings")
                return df
            else:
                print(f" Required columns 'ASV' and 'Genus' not found")
                return None
                
        except Exception as e:
            print(f" Error loading ASV metadata: {e}")
            return None
    
    def load_feature_table(self):
        """Load pipeline feature table"""
        try:
            feature_table_path = self.pipeline_output_dir / "qiime_output/relevant_results/feature_table.tsv"
            print(f"\n Loading feature table: {feature_table_path}")
            
            if not feature_table_path.exists():
                print(f" Feature table not found at {feature_table_path}")
                return None
            
            # Load QIIME2 feature table (skip comment line)
            df = pd.read_csv(feature_table_path, sep='\t', skiprows=1, index_col=0)
            print(f"Feature table loaded: {df.shape[0]} features, {df.shape[1]} samples")
            print(f"Samples: {list(df.columns[:5])}{'...' if df.shape[1] > 5 else ''}")
            
            return df
            
        except Exception as e:
            print(f"Error loading feature table: {e}")
            return None
    
    def load_taxonomy(self):
        """Load pipeline taxonomy assignments"""
        try:
            taxonomy_path = self.pipeline_output_dir / "qiime_output/relevant_results/taxonomy.tsv"
            print(f"\n Loading taxonomy: {taxonomy_path}")
            
            if not taxonomy_path.exists():
                print(f"Taxonomy file not found at {taxonomy_path}")
                return None
            
            df = pd.read_csv(taxonomy_path, sep='\t', index_col=0)
            print(f"Taxonomy loaded: {len(df)} feature assignments")
            
            return df
            
        except Exception as e:
            print(f" Error loading taxonomy: {e}")
            return None
    
    def load_ground_truth_data(self):
        """Load ground truth abundance files and map to genera"""
        print(f"\nLoading ground truth data from: {self.ground_truth_dir}")
        
        truth_files = list(self.ground_truth_dir.glob("*_abundance.txt"))
        print(f" Found {len(truth_files)} ground truth files")
        
        ground_truth = {}
        for file in truth_files:
            sample_name = file.stem.replace('_abundance', '')
            
            try:
                # Load ASV abundances
                df = pd.read_csv(file, sep='\t', header=None, names=['asv_id', 'abundance'])
                
                # Map to genera using ASV metadata
                if self.asv_mapping is not None:
                    df_with_genus = df.merge(self.asv_mapping[['ASV_ID', 'Genus']], 
                                           left_on='asv_id', right_on='ASV_ID', how='left')
                    
                    # Group by genus and sum abundances
                    genus_abundances = df_with_genus.groupby('Genus')['abundance'].sum()
                    genus_abundances = genus_abundances[genus_abundances > 0]  # Remove zeros
                    
                    ground_truth[sample_name] = genus_abundances
                    print(f"{sample_name}: {len(df)} ASVs → {len(genus_abundances)} genera")
                else:
                    print(f"No ASV mapping available for {sample_name}")
                    
            except Exception as e:
                print(f"Error loading {file}: {e}")
        
        print(f" Successfully loaded ground truth for {len(ground_truth)} samples")
        return ground_truth
    
    def extract_genus_from_taxonomy(self, taxonomy_string):
        """Extract genus from QIIME2 taxonomy string"""
        if pd.isna(taxonomy_string):
            return "Unknown"
        
        # Parse taxonomy levels
        levels = taxonomy_string.split(';')
        for level in levels:
            level = level.strip()
            if level.startswith('g__'):
                genus = level.replace('g__', '').strip()
                if genus and genus not in ['', 'unidentified', 'uncultured']:
                    return genus
        
        return "Unknown"
    
    def calculate_diversity_metrics(self, abundances):
        """Calculate alpha diversity metrics"""
        # Convert to relative abundances
        rel_abundances = abundances / abundances.sum()
        rel_abundances = rel_abundances[rel_abundances > 0]
        
        # Shannon diversity
        shannon = -np.sum(rel_abundances * np.log(rel_abundances))
        
        # Simpson diversity (1 - D)
        simpson = 1 - np.sum(rel_abundances ** 2)
        
        # Richness
        richness = len(rel_abundances)
        
        # Evenness
        evenness = shannon / np.log(richness) if richness > 1 else 0
        
        return {
            'shannon': shannon,
            'simpson': simpson,
            'richness': richness,
            'evenness': evenness
        }
    
    def normalize_sample_names(self):
        """Match sample names between pipeline and ground truth"""
        if self.feature_table is None:
            return {}
            
        pipeline_samples = list(self.feature_table.columns)
        truth_samples = list(self.ground_truth.keys())
        
        print(f"\n🔗 Matching sample names...")
        print(f"Pipeline samples: {len(pipeline_samples)}")
        print(f"Ground truth samples: {len(truth_samples)}")
        
        sample_mapping = {}
        for truth_name in truth_samples:
            # Try exact match
            if truth_name in pipeline_samples:
                sample_mapping[truth_name] = truth_name
            else:
                # Try fuzzy matching
                for pipe_name in pipeline_samples:
                    # Remove underscores and compare
                    truth_clean = truth_name.replace('_', '').replace('-', '').lower()
                    pipe_clean = pipe_name.replace('_', '').replace('-', '').lower()
                    
                    if truth_clean == pipe_clean or truth_clean in pipe_clean or pipe_clean in truth_clean:
                        sample_mapping[truth_name] = pipe_name
                        break
        
        print(f"Successfully matched {len(sample_mapping)} samples")
        return sample_mapping
    
    def compare_sample_genera(self, sample_name, pipeline_sample_name):
        """Compare genus-level abundances for a single sample"""
        print(f"\n🔬 Analyzing {sample_name} ↔ {pipeline_sample_name}")
        
        # Get pipeline data
        pipeline_data = self.feature_table[pipeline_sample_name]
        pipeline_data = pipeline_data[pipeline_data > 0]  # Remove zeros
        pipeline_rel = pipeline_data / pipeline_data.sum()
        
        # Map pipeline features to genera
        pipeline_genera = {}
        unmapped_count = 0
        
        for feature_id, abundance in pipeline_rel.items():
            if feature_id in self.taxonomy.index:
                tax_string = self.taxonomy.loc[feature_id, 'Taxon']
                genus = self.extract_genus_from_taxonomy(tax_string)
                
                if genus in pipeline_genera:
                    pipeline_genera[genus] += abundance
                else:
                    pipeline_genera[genus] = abundance
            else:
                unmapped_count += 1
        
        # Get ground truth genera
        truth_genera = self.ground_truth[sample_name]
        truth_rel = truth_genera / truth_genera.sum()
        
        # Find common genera
        common_genera = set(pipeline_genera.keys()) & set(truth_rel.index)
        
        if len(common_genera) == 0:
            print(f" No common genera found!")
            print(f"Pipeline genera: {list(pipeline_genera.keys())[:5]}")
            print(f"Truth genera: {list(truth_rel.index[:5])}")
            return None
        
        print(f"Found {len(common_genera)} common genera")
        print(f"Pipeline: {len(pipeline_genera)} genera, Truth: {len(truth_rel)} genera")
        if unmapped_count > 0:
            print(f"{unmapped_count} features could not be taxonomically classified")
        
        # Create comparison vectors
        pipeline_abundances = [pipeline_genera[genus] for genus in common_genera]
        truth_abundances = [truth_rel[genus] for genus in common_genera]
        
        # Calculate diversity metrics
        pipeline_diversity = self.calculate_diversity_metrics(pd.Series(pipeline_genera))
        truth_diversity = self.calculate_diversity_metrics(truth_rel)
        
        return {
            'sample': sample_name,
            'pipeline_abundances': pipeline_abundances,
            'truth_abundances': truth_abundances,
            'genera': list(common_genera),
            'n_common_genera': len(common_genera),
            'n_pipeline_genera': len(pipeline_genera),
            'n_truth_genera': len(truth_rel),
            'pipeline_diversity': pipeline_diversity,
            'truth_diversity': truth_diversity,
            'unmapped_features': unmapped_count
        }
    
    def calculate_performance_metrics(self, comparison):
        """Calculate benchmarking metrics"""
        if comparison is None:
            return None
        
        pipeline_ab = np.array(comparison['pipeline_abundances'])
        truth_ab = np.array(comparison['truth_abundances'])
        
        # Abundance correlation
        pearson_r, pearson_p = pearsonr(pipeline_ab, truth_ab)
        spearman_r, spearman_p = spearmanr(pipeline_ab, truth_ab)
        
        # Error metrics
        mse = mean_squared_error(truth_ab, pipeline_ab)
        mae = mean_absolute_error(truth_ab, pipeline_ab)
        rmse = np.sqrt(mse)
        relative_error = np.mean(np.abs(pipeline_ab - truth_ab) / (truth_ab + 1e-10))
        
        # Detection metrics
        detection_rate = comparison['n_common_genera'] / comparison['n_truth_genera']
        
        # Diversity metrics comparison
        div_metrics = {}
        for metric in ['shannon', 'simpson', 'richness', 'evenness']:
            truth_val = comparison['truth_diversity'][metric]
            pipeline_val = comparison['pipeline_diversity'][metric]
            div_metrics[f'{metric}_truth'] = truth_val
            div_metrics[f'{metric}_pipeline'] = pipeline_val
            div_metrics[f'{metric}_error'] = abs(truth_val - pipeline_val)
            div_metrics[f'{metric}_relative_error'] = abs(truth_val - pipeline_val) / (truth_val + 1e-10)
        
        result = {
            'sample': comparison['sample'],
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
            'mse': mse,
            'mae': mae,
            'rmse': rmse,
            'relative_error': relative_error,
            'detection_rate': detection_rate,
            'n_common_genera': comparison['n_common_genera'],
            'n_pipeline_genera': comparison['n_pipeline_genera'],
            'n_truth_genera': comparison['n_truth_genera'],
            'unmapped_features': comparison['unmapped_features']
        }
        
        result.update(div_metrics)
        return result
    
    def extract_sample_condition(self, sample_name):
        """Extract experimental condition from sample name"""
        sample_lower = sample_name.lower()
        
        if 'standard' in sample_lower:
            return 'standard'
        elif 'basic_error' in sample_lower:
            return 'basic_error'
        elif 'high_quality' in sample_lower:
            return 'high_quality'
        elif 'no_gc_bias' in sample_lower:
            return 'no_gc_bias'
        elif 'depth' in sample_lower:
            if '0.25x' in sample_lower:
                return 'depth_0.25x'
            elif '0.5x' in sample_lower:
                return 'depth_0.5x'
            elif '2.0x' in sample_lower or '2x' in sample_lower:
                return 'depth_2x'
            elif '5.0x' in sample_lower or '5x' in sample_lower:
                return 'depth_5x'
            else:
                return 'depth_other'
        elif 'miseq' in sample_lower:
            if '24' in sample_lower:
                return 'miseq_24'
            elif '28' in sample_lower:
                return 'miseq_28'
            else:
                return 'miseq_standard'
        else:
            return 'other'
    
    def run_benchmarking(self):
        """Execute complete benchmarking analysis"""
        print("\n" + "="*60)
        print("STARTING GENUS-LEVEL PIPELINE BENCHMARKING")
        print("="*60)
        
        # Validate inputs
        if self.feature_table is None or self.taxonomy is None or self.asv_mapping is None:
            print("Missing required input data!")
            return None
        
        if not self.ground_truth:
            print("No ground truth data loaded!")
            return None
        
        # Match samples
        sample_mapping = self.normalize_sample_names()
        if not sample_mapping:
            print("No samples could be matched!")
            return None
        
        # Run comparisons
        all_metrics = []
        for truth_name, pipeline_name in sample_mapping.items():
            comparison = self.compare_sample_genera(truth_name, pipeline_name)
            if comparison:
                metrics = self.calculate_performance_metrics(comparison)
                if metrics:
                    all_metrics.append(metrics)
                    self.comparison_data.append(comparison)
        
        if not all_metrics:
            print("No successful comparisons!")
            return None
        
        # Create results DataFrame
        self.metrics_df = pd.DataFrame(all_metrics)
        self.metrics_df['condition'] = self.metrics_df['sample'].apply(self.extract_sample_condition)
        
        # Generate reports
        self.generate_summary_report()
        self.generate_visualizations()
        
        print(f"\nBENCHMARKING COMPLETE!")
        print(f"Analyzed {len(self.metrics_df)} samples")
        print(f"Results saved to: {self.results_dir}")
        
        return self.metrics_df
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        print(f"\nGenerating summary report...")
        
        # Overall performance
        overall = {
            'Abundance Correlation': {
                'Mean Pearson r': self.metrics_df['pearson_r'].mean(),
                'Median Pearson r': self.metrics_df['pearson_r'].median(),
                'Samples with good correlation (r>0.7)': len(self.metrics_df[self.metrics_df['pearson_r'] > 0.7]),
                'Samples with excellent correlation (r>0.9)': len(self.metrics_df[self.metrics_df['pearson_r'] > 0.9]),
                'Mean Spearman r': self.metrics_df['spearman_r'].mean(),
            },
            'Detection Performance': {
                'Mean detection rate': self.metrics_df['detection_rate'].mean(),
                'Mean genera detected': self.metrics_df['n_common_genera'].mean(),
                'Mean genera in truth': self.metrics_df['n_truth_genera'].mean(),
                'Mean genera in pipeline': self.metrics_df['n_pipeline_genera'].mean(),
            },
            'Error Metrics': {
                'Mean RMSE': self.metrics_df['rmse'].mean(),
                'Mean MAE': self.metrics_df['mae'].mean(),
                'Mean relative error': self.metrics_df['relative_error'].mean(),
            }
        }
        
        # Diversity performance
        diversity_metrics = ['shannon', 'simpson', 'richness', 'evenness']
        diversity_performance = {}
        for metric in diversity_metrics:
            truth_col = f'{metric}_truth'
            pipeline_col = f'{metric}_pipeline'
            error_col = f'{metric}_relative_error'
            
            if truth_col in self.metrics_df.columns:
                corr, _ = pearsonr(self.metrics_df[truth_col], self.metrics_df[pipeline_col])
                diversity_performance[f'{metric.title()} correlation'] = corr
                diversity_performance[f'{metric.title()} mean relative error'] = self.metrics_df[error_col].mean()
        
        overall['Diversity Metrics'] = diversity_performance
        
        # Performance by condition
        condition_performance = {}
        for condition in self.metrics_df['condition'].unique():
            cond_data = self.metrics_df[self.metrics_df['condition'] == condition]
            condition_performance[f'{condition} (n={len(cond_data)})'] = {
                'Abundance correlation': cond_data['pearson_r'].mean(),
                'Detection rate': cond_data['detection_rate'].mean(),
                'Shannon error': cond_data['shannon_relative_error'].mean() if 'shannon_relative_error' in cond_data.columns else 'N/A'
            }
        
        overall['Performance by Condition'] = condition_performance
        
        # Save report
        with open(self.results_dir / 'benchmarking_summary.txt', 'w') as f:
            f.write("MICROBIOME PIPELINE BENCHMARKING SUMMARY\n")
            f.write("Genus-Level Validation Results\n")
            f.write("=" * 60 + "\n\n")
            
            for category, metrics in overall.items():
                f.write(f"{category}:\n")
                f.write("-" * 40 + "\n")
                self._write_metrics_recursive(f, metrics, indent=2)
                f.write("\n")
            
            # Key insights
            f.write("KEY INSIGHTS:\n")
            f.write("-" * 40 + "\n")
            
            mean_corr = overall['Abundance Correlation']['Mean Pearson r']
            mean_detection = overall['Detection Performance']['Mean detection rate']
            
            if mean_corr > 0.8:
                f.write("EXCELLENT abundance correlation (r > 0.8)\n")
            elif mean_corr > 0.6:
                f.write("GOOD abundance correlation (0.6 < r < 0.8)\n")
            else:
                f.write("MODERATE abundance correlation (r < 0.6)\n")
            
            if mean_detection > 0.8:
                f.write("EXCELLENT genus detection rate (>80%)\n")
            elif mean_detection > 0.6:
                f.write("GOOD genus detection rate (60-80%)\n")
            else:
                f.write("MODERATE genus detection rate (<60%)\n")
        
        # Save detailed metrics
        self.metrics_df.to_csv(self.results_dir / 'detailed_metrics.csv', index=False)
        
        print(f"Summary saved to: {self.results_dir / 'benchmarking_summary.txt'}")
        print(f"Detailed metrics saved to: {self.results_dir / 'detailed_metrics.csv'}")
    
    def _write_metrics_recursive(self, f, metrics, indent=0):
        """Helper function to write nested metrics"""
        for key, value in metrics.items():
            if isinstance(value, dict):
                f.write(" " * indent + f"{key}:\n")
                self._write_metrics_recursive(f, value, indent + 2)
            elif isinstance(value, float):
                f.write(" " * indent + f"{key}: {value:.4f}\n")
            else:
                f.write(" " * indent + f"{key}: {value}\n")
    
    def generate_visualizations(self):
        """Generate comprehensive visualization plots"""
        print(f"\nGenerating visualizations...")
        
        plt.style.use('default')
        
        # Plot 1: Abundance correlation scatter plots
        n_samples = min(len(self.comparison_data), 16)
        cols = 4
        rows = (n_samples + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(20, 5*rows))
        if rows == 1:
            axes = [axes] if cols == 1 else axes
        else:
            axes = axes.flatten()
        
        for i in range(n_samples):
            comparison = self.comparison_data[i]
            ax = axes[i]
            
            ax.scatter(comparison['truth_abundances'], comparison['pipeline_abundances'], 
                      alpha=0.7, s=50, color='steelblue')
            
            # 1:1 line
            max_val = max(max(comparison['truth_abundances']), max(comparison['pipeline_abundances']))
            ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.6, linewidth=2)
            
            ax.set_xlabel('Ground Truth Abundance', fontsize=10)
            ax.set_ylabel('Pipeline Abundance', fontsize=10)
            
            r_val = self.metrics_df.iloc[i]['pearson_r']
            n_genera = comparison['n_common_genera']
            ax.set_title(f'{comparison["sample"]}\nr={r_val:.3f}, n={n_genera} genera', fontsize=10)
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(n_samples, len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'abundance_correlations.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Diversity metrics comparison
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        diversity_metrics = ['shannon', 'simpson', 'richness', 'evenness']
        colors = ['steelblue', 'forestgreen', 'darkorange', 'purple']
        
        for i, metric in enumerate(diversity_metrics):
            row, col = i // 2, i % 2
            truth_col = f'{metric}_truth'
            pipeline_col = f'{metric}_pipeline'
            
            if truth_col in self.metrics_df.columns:
                axes[row, col].scatter(self.metrics_df[truth_col], self.metrics_df[pipeline_col], 
                                     alpha=0.7, s=60, color=colors[i])
                
                # 1:1 line
                min_val = min(self.metrics_df[truth_col].min(), self.metrics_df[pipeline_col].min())
                max_val = max(self.metrics_df[truth_col].max(), self.metrics_df[pipeline_col].max())
                axes[row, col].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.6, linewidth=2)
                
                # Correlation
                corr, _ = pearsonr(self.metrics_df[truth_col], self.metrics_df[pipeline_col])
                
                axes[row, col].set_xlabel(f'Ground Truth {metric.title()}', fontsize=12)
                axes[row, col].set_ylabel(f'Pipeline {metric.title()}', fontsize=12)
                axes[row, col].set_title(f'{metric.title()} Diversity (r={corr:.3f})', fontsize=14)
                axes[row, col].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'diversity_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Plot 3: Performance by condition
        if 'condition' in self.metrics_df.columns:
            conditions = self.metrics_df['condition'].unique()
            
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            
            # Abundance correlation by condition
            condition_data = [self.metrics_df[self.metrics_df['condition']==cond]['pearson_r'].values 
                            for cond in conditions]
            bp1 = axes[0,0].boxplot(condition_data, labels=conditions, patch_artist=True)
            for patch in bp1['boxes']:
                patch.set_facecolor('lightblue')
            axes[0,0].set_title('Abundance Correlation by Condition', fontsize=14)
            axes[0,0].set_ylabel('Pearson r', fontsize=12)
            axes[0,0].tick_params(axis='x', rotation=45)
            axes[0,0].grid(True, alpha=0.3)
            
            # Detection rate by condition
            condition_data = [self.metrics_df[self.metrics_df['condition']==cond]['detection_rate'].values 
                            for cond in conditions]
            bp2 = axes[0,1].boxplot(condition_data, labels=conditions, patch_artist=True)
            for patch in bp2['boxes']:
                patch.set_facecolor('lightgreen')
            axes[0,1].set_title('Genus Detection Rate by Condition', fontsize=14)
            axes[0,1].set_ylabel('Detection Rate', fontsize=12)
            axes[0,1].tick_params(axis='x', rotation=45)
            axes[0,1].grid(True, alpha=0.3)
            
            # Shannon diversity error by condition
            if 'shannon_relative_error' in self.metrics_df.columns:
                condition_data = [self.metrics_df[self.metrics_df['condition']==cond]['shannon_relative_error'].values 
                                for cond in conditions]
                bp3 = axes[1,0].boxplot(condition_data, labels=conditions, patch_artist=True)
                for patch in bp3['boxes']:
                    patch.set_facecolor('lightcoral')
                axes[1,0].set_title('Shannon Diversity Error by Condition', fontsize=14)
                axes[1,0].set_ylabel('Relative Error', fontsize=12)
                axes[1,0].tick_params(axis='x', rotation=45)
                axes[1,0].grid(True, alpha=0.3)
            
            # Overall performance summary
            summary_metrics = ['pearson_r', 'detection_rate', 'shannon_relative_error']
            summary_labels = ['Abundance\nCorrelation', 'Detection\nRate', 'Shannon\nError']
            summary_values = [self.metrics_df[metric].mean() for metric in summary_metrics if metric in self.metrics_df.columns]
            summary_colors = ['steelblue', 'forestgreen', 'darkorange']
            
            bars = axes[1,1].bar(range(len(summary_values)), summary_values, color=summary_colors[:len(summary_values)])
            axes[1,1].set_xticks(range(len(summary_values)))
            axes[1,1].set_xticklabels(summary_labels[:len(summary_values)])
            axes[1,1].set_title('Overall Performance Summary', fontsize=14)
            axes[1,1].grid(True, alpha=0.3, axis='y')
            
            # Add value labels on bars
            for bar, value in zip(bars, summary_values):
                height = bar.get_height()
                axes[1,1].text(bar.get_x() + bar.get_width()/2., height + 0.01,
                             f'{value:.3f}', ha='center', va='bottom', fontsize=11)
            
            plt.tight_layout()
            plt.savefig(self.results_dir / 'performance_by_condition.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        print(f"Visualizations saved:")
        print(f" {self.results_dir / 'abundance_correlations.png'}")
        print(f" {self.results_dir / 'diversity_comparison.png'}")
        print(f" {self.results_dir / 'performance_by_condition.png'}")

def main():
    """Main execution function"""
    print("MICROBIOME PIPELINE BENCHMARKING v15")
    print("Genus-Level Validation Against Synthetic Ground Truth")
    print("=" * 60)
    
    # Configuration - paths relative to main project directory
    pipeline_output = "runs/S01070625"
    ground_truth_dir = "group4/microbiome-benchmarking/benchmarking_output/synthetic_samples"
    asv_metadata_file = "group4/microbiome-benchmarking/input_data/asv_metadata.tsv"
    results_dir = "group4/microbiome-benchmarking/benchmarking_results"
    
    # Validate paths exist
    required_paths = [
        Path(pipeline_output),
        Path(ground_truth_dir),
        Path(asv_metadata_file)
    ]
    
    for path in required_paths:
        if not path.exists():
            print(f"ERROR: Required path not found: {path}")
            return False
    
    try:
        # Initialize and run benchmarking
        benchmarker = GenusLevelBenchmarker(
            pipeline_output_dir=pipeline_output,
            ground_truth_dir=ground_truth_dir,
            asv_metadata_file=asv_metadata_file,
            results_dir=results_dir
        )
        
        # Execute benchmarking
        results = benchmarker.run_benchmarking()
        
        if results is not None:
            print("\n" + "="*60)
            print("BENCHMARKING SUCCESSFULLY COMPLETED!")
            print("="*60)
            print(f"Samples analyzed: {len(results)}")
            print(f"Results directory: {benchmarker.results_dir}")
            print(f"\nKEY FILES GENERATED:")
            print(f"  benchmarking_summary.txt - Comprehensive performance summary")
            print(f"  detailed_metrics.csv - Per-sample detailed metrics")
            print(f"  abundance_correlations.png - Sample correlation plots")
            print(f"  diversity_comparison.png - Diversity metrics validation")
            print(f"  performance_by_condition.png - Performance across conditions")
            
            # Print quick summary
            mean_corr = results['pearson_r'].mean()
            mean_detection = results['detection_rate'].mean()
            
            print(f"\n QUICK SUMMARY:")
            print(f"   Mean abundance correlation: {mean_corr:.3f}")
            print(f"   Mean genus detection rate: {mean_detection:.3f}")
            
            if mean_corr > 0.8 and mean_detection > 0.8:
                print(f"   EXCELLENT pipeline performance!")
            elif mean_corr > 0.6 and mean_detection > 0.6:
                print(f"   GOOD pipeline performance!")
            else:
                print(f"   Moderate pipeline performance - check detailed results")
            
            print(f"\n Check '{benchmarker.results_dir}/benchmarking_summary.txt' for detailed insights!")
            return True
            
        else:
            print(" BENCHMARKING FAILED!")
            print("Check the error messages above for troubleshooting.")
            return False
            
    except Exception as e:
        print(f" UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
