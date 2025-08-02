# Single-Cell RNA-seq Analysis Pipeline

A comprehensive, automated pipeline for single-cell RNA sequencing data analysis using Seurat in R. This pipeline provides end-to-end analysis from raw count matrices to cell type annotation and visualization.

## 🚀 Features

- **Automated Seurat Object Creation**: Build Seurat objects from 10X Genomics outputs with metadata integration
- **Quality Control & Preprocessing**: Automated filtering, normalization, and feature selection
- **Dimensionality Reduction**: PCA, UMAP, and clustering with optimized parameters
- **Data Integration**: RPCA-based integration for multi-sample datasets
- **Cell Type Annotation**: Reference-based cell type prediction using existing atlases
- **Differential Abundance Analysis**: Milo-based detection of cell population changes
- **Cell Communication Analysis**: CellChat-based intercellular signaling inference
- **Comprehensive Visualization**: Automated generation of publication-ready plots
- **Parallel Processing**: Multi-core support for faster computation
- **Organized Output Structure**: Clean directory organization for all results

## 📁 Pipeline Components

### 1. `create_seurat.r`
Creates Seurat objects from 10X Genomics output and integrates sample metadata.

**Usage:**
```bash
Rscript create_seurat.r <input_dir> <metadata.csv> <output_rds>
```

### 2. `automated_seurat_analysis.R`
Comprehensive single-cell analysis including QC, normalization, clustering, and marker discovery.

**Usage:**
```bash
Rscript automated_seurat_analysis.R <seurat_object.rds>
```

### 3. `automated_integration_rpca.R`
RPCA-based integration for multi-sample datasets to remove batch effects.

**Usage:**
```bash
Rscript automated_integration_rpca.R <seurat_object.rds>
```

### 4. `annotation_pipeline.R`
Reference-based cell type annotation using h5ad reference atlases.

**Usage:**
```bash
# Increase system limits first
ulimit -s unlimited
Rscript annotation_pipeline.R <ref_h5ad> <query_rds> <marker_csv>
```

### 5. `automate_milo_analysis.R`
Differential abundance analysis using Milo to identify changes in cell population compositions between conditions.

**Usage:**
```bash
Rscript automate_milo_analysis.R <rds_directory> [group_column] [annotation_column]
```

### 6. `cellchat.R`
Cell-cell communication analysis using CellChat to infer intercellular signaling networks.

**Usage:**
```bash
Rscript cellchat.R <rds_directory> [group_column] [annotation_column]
```

## 🛠️ Installation

### Prerequisites
- R (≥ 4.0.0)
- Required R packages (auto-installed by scripts):
  - `Seurat`
  - `ggplot2`
  - `patchwork`
  - `dplyr`
  - `Matrix`
  - `future`
  - `future.apply`
  - `zellkonverter`
  - `SingleCellExperiment`
  - `BiocManager`
  - `miloR` (for differential abundance)
  - `CellChat` (for cell communication)
  - `ComplexHeatmap`

### Quick Setup
1. Clone this repository:
```bash
git clone https://github.com/yourusername/scrnaseq-pipeline.git
cd scrnaseq-pipeline
```

2. Make scripts executable:
```bash
chmod +x *.R
```

3. Install R dependencies (optional - scripts auto-install):
```r
# In R console
install.packages("BiocManager")
BiocManager::install(c("Seurat", "ggplot2", "patchwork", "dplyr", 
                       "Matrix", "future", "future.apply", 
                       "zellkonverter", "SingleCellExperiment",
                       "miloR", "ComplexHeatmap"))
# CellChat from GitHub
remotes::install_github("sqjin/CellChat")
```

## 📊 Usage Examples

### Basic Workflow

#### Step 1: Create Seurat Object
```bash
# Prepare your input directory with 10X outputs:
# input_dir/
# ├── sample1/
# │   ├── matrix.mtx
# │   ├── genes.tsv
# │   └── barcodes.tsv
# └── sample2/
#     ├── matrix.mtx
#     ├── genes.tsv
#     └── barcodes.tsv

# Create metadata.csv:
# sample_id,condition,batch,patient_id
# sample1,control,batch1,patient_A
# sample2,treatment,batch1,patient_B

Rscript create_seurat.r /path/to/input_dir metadata.csv seurat_object.rds
```

#### Step 2: Quality Control and Analysis
```bash
Rscript automated_seurat_analysis.R seurat_object.rds
```

#### Step 3: Integration (for multi-sample data)
```bash
Rscript automated_integration_rpca.R seurat_object.rds
```

#### Step 4: Cell Type Annotation
```bash
# Download reference atlas (e.g., from cellxgene)
ulimit -s unlimited
Rscript annotation_pipeline.R reference_atlas.h5ad integrated_seurat.rds markers.csv
```

#### Step 5: Differential Abundance Analysis
```bash
# Analyze changes in cell population compositions
# Requires directory with multiple condition-specific RDS files
Rscript automate_milo_analysis.R /path/to/rds_directory condition_column cell_type_column
```

#### Step 6: Cell Communication Analysis
```bash
# Infer intercellular signaling networks
# Uses same directory structure as Milo analysis
Rscript cellchat.R /path/to/rds_directory condition_column cell_type_column
```

### Advanced Usage

#### Custom Parameters
The scripts use sensible defaults but can be modified for specific needs:

- **Quality Control Thresholds**: Edit filtering criteria in `automated_seurat_analysis.R`
- **Clustering Resolution**: Modify resolution parameters for different granularity
- **Integration Features**: Adjust number of variable features for integration
- **Parallel Processing**: Scripts automatically detect cores; modify `workers` variable if needed

## 📂 Output Structure

Each analysis creates an organized output directory:

```
analysis_name_analysis/
├── Visualization/
│   ├── QC plots (before/after filtering)
│   ├── PCA analysis plots
│   ├── UMAP visualizations
│   ├── Clustering results
│   ├── Feature plots
│   ├── Milo differential abundance plots
│   └── CellChat communication networks
└── Data/
    ├── Processed Seurat objects
    ├── Marker gene tables
    ├── Integration results
    ├── Milo analysis results
    └── CellChat objects and networks
```

## 📋 Input Requirements

### For `create_seurat.r`:
- **Input Directory**: Contains subdirectories with 10X Genomics outputs
  - `matrix.mtx`: Count matrix
  - `genes.tsv`: Gene annotations
  - `barcodes.tsv`: Cell barcodes
- **Metadata CSV**: Sample information with `sample_id` column

### For `annotation_pipeline.R`:
- **Reference H5AD**: Single-cell reference atlas in h5ad format
- **Query RDS**: Seurat object to be annotated
- **Marker CSV**: Optional marker gene list

### For `automate_milo_analysis.R`:
- **RDS Directory**: Folder containing condition-specific Seurat RDS files
- **Group Column**: Metadata column defining experimental conditions (default: "stage")
- **Annotation Column**: Metadata column with cell type labels (default: "predicted.id")

### For `cellchat.R`:
- **RDS Directory**: Same structure as Milo analysis
- **Group Column**: Condition metadata column (default: "stage")  
- **Annotation Column**: Cell type annotation column (default: "predicted.id")

## ⚡ Performance Optimization

### System Requirements
- **RAM**: Minimum 16GB, recommended 32GB+ for large datasets
- **CPU**: Multi-core processor recommended
- **Storage**: SSD recommended for faster I/O

### Memory Management
```bash
# Increase system limits before running annotation
ulimit -s unlimited
export R_CStackLimit=10000000
```

### Parallel Processing
Scripts automatically use `detectCores() - 1` workers. Modify if needed:
```r
workers <- 8  # Set specific number of cores
```

## 🔧 Troubleshooting

### Common Issues

1. **Memory Errors**:
   - Increase system memory limits
   - Reduce number of features for integration
   - Use `future.globals.maxSize` parameter

2. **Package Installation**:
   - Scripts auto-install missing packages
   - For manual installation issues, use BiocManager

3. **H5AD Reading Errors**:
   - Ensure zellkonverter and SingleCellExperiment are installed
   - Check h5ad file integrity

4. **Integration Failures**:
   - Verify sample overlap in feature space
   - Adjust integration parameters

### Performance Tips
- Use SSD storage for better I/O performance
- Increase R memory limits for large datasets
- Consider filtering rare genes/cells before integration

## 📈 Example Results

The pipeline generates:
- Quality control metrics and visualizations
- Dimensionality reduction plots (PCA, UMAP)
- Clustering results with different resolutions
- Marker gene identification
- Cell type annotations with confidence scores
- Integration results for batch correction
- Differential abundance analysis with statistical testing
- Cell-cell communication networks and pathway analysis

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{scrnaseq_pipeline,
  title = {Single-Cell RNA-seq Analysis Pipeline},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/scrnaseq-pipeline}
}
```

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/scrnaseq-pipeline/issues)
- **Documentation**: [Wiki](https://github.com/yourusername/scrnaseq-pipeline/wiki)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/scrnaseq-pipeline/discussions)

## 🔄 Version History

- **v1.0.0**: Initial release with core pipeline functionality
- **v1.1.0**: Added RPCA integration support
- **v1.2.0**: Enhanced annotation pipeline with h5ad support
- **v1.3.0**: Added Milo differential abundance analysis
- **v1.4.0**: Integrated CellChat cell communication analysis

---

**Keywords**: single-cell RNA-seq, Seurat, R, pipeline, automation, cell type annotation, data integration, quality control, differential abundance, cell communication, Milo, CellChat