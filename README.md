# Parental-genotype-partitioning (parentalgp)

This R package provides a robust method to partition and infer parental genotypes using VCF files generated from hybrid seed endosperm. It leverages a Beta-binomial model to account for sequencing overdispersion when evaluating allele depth ratios (1:3 and 2:3), which are characteristic of the triploid endosperm (2 maternal : 1 paternal).

## 📦 Installation

You can install the development version of this package directly from GitHub.

```R
# If you haven't installed devtools yet
# install.packages("devtools")

# Install the package
devtools::install_github("zhoutaohao/Parental-genotype-partitioning")
```

## 🚀 Usage

The core function of this package is infer_parental_genotype(). It takes a VCF file containing multiple F1 hybrid endosperm samples and infers the maternal and paternal genotypes for each marker.

```R
library(parentalgp)

# Define the path to your VCF file
# Note: Ensure the VCF contains multiple F1 hybrid samples and the FORMAT column includes "AD" (Allele Depth).
vcf_path <- "path/to/your/hybrid_endosperm.vcf"

# Run the inference function
result <- infer_parental_genotype(
  vcf_file = vcf_path, 
  output_file = "inferred_parents.tsv", 
  min_depth = 150, 
  rho = 0.001, 
  alpha = 0.05
)

# View the first few rows of the result
head(result)
```

## 📊 Output
The function returns a data frame (and saves a TSV file) containing:

markerID: The genomic position of the SNP (e.g., Chr_Pos).

[SampleName]_maternal: The inferred maternal genotype (0 = REF, 2 = ALT).

[SampleName]_paternal: The inferred paternal genotype (0 = REF, 2 = ALT).

## ⚙️ Parameters
vcf_file: Character. Path to the input VCF file.

output_file: Character. Path to save the output table. Default is "parental_genotype.csv".

min_depth: Integer. Minimum total read depth (REF + ALT) required to consider a site. Default is 150.

rho: Numeric. The overdispersion parameter for the Beta-binomial model. Default is 0.001.

alpha: Numeric. The cutoff for the statistical test. Default is 0.05.

## 📄 License
This project is licensed under the MIT License.



