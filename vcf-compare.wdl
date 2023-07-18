version 1.0

# Copyright (c) 2023 Leiden University Medical Center

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/bcftools.wdl" as bcftools
import "tasks/snpeff.wdl" as snpeff

struct CompareUnit {
    File vcf1
    File vcf1index
    File vcf2
    File vcf2index
    String? name1
    String? name2
}

workflow VcfCompare {
    input {
        Array[CompareUnit] units
        String outputDir = "."
        String snpEffGenomeVersion
        File snpEffDatadirZip
        Array[String] snpEffConfigOptions = []
        File? regionsFile
    }

    scatter (unit in units) {
        String name1 = select_first([unit.name1, "1"])
        String name2 = select_first([unit.name2, "2"])

        call bcftools.Norm as normalizeVcf1 {
            input:
                inputVcf = unit.vcf1,
                inputVcfIndex = unit.vcf1index,
                outputPath = outputDir + "/" + name1 + ".normalized.vcf.gz",
                regionsFile=regionsFile
        }
        call bcftools.Norm as normalizeVcf2 {
            input:
                inputVcf = unit.vcf2,
                inputVcfIndex = unit.vcf2index,
                outputPath = outputDir + "/" + name2 + ".normalized.vcf.gz",
                regionsFile=regionsFile
        }

        call bcftools.View as filterHomRefVcf1 {
            input:
                inputFile=normalizeVcf1.outputVcf,
                outputPath = outputDir + "/" + name1 + ".normalized.nohomref.vcf.gz",
                include="'GT[*]=\"alt\"'"
        }

        call bcftools.View as filterHomRefVcf2 {
            input:
                inputFile=normalizeVcf2.outputVcf,
                outputPath = outputDir + "/" + name2 + ".normalized.nohomref.vcf.gz",
                include="'GT[*]=\"alt\"'"
        }

        call bcftools.Isec as intersect {
            input: 
                aVcf = filterHomRefVcf1.outputVcf,
                aVcfIndex = select_first([filterHomRefVcf1.outputVcfIndex]),
                bVcf = filterHomRefVcf2.outputVcf,
                bVcfIndex = select_first([filterHomRefVcf2.outputVcfIndex]),
                prefix = outputDir + "/" + name1
        }
        
        call vennStats as vennStats {
            input:
                uniqueAVcf = intersect.privateAVcf,
                uniqueBVcf = intersect.privateBVcf,
                sharedVcf = intersect.sharedAVcf,
                nameA = name1,
                nameB = name2,
                outputPath = outputDir + "/" + name1 + "." + "venn.txt"
        }

        call snpeff.SnpEff as annotateUnique1 {
            input:
                vcf = intersect.privateAVcf,
                vcfIndex = intersect.privateAVcfIndex,
                genomeVersion = snpEffGenomeVersion,
                datadirZip = snpEffDatadirZip,
                outputPath = outputDir + "/" + name1 + ".annotated.vcf",
                configOptions = snpEffConfigOptions,
        }

        call snpeff.SnpEff as annotateUnique2 {
            input:
                vcf = intersect.privateBVcf,
                vcfIndex = intersect.privateBVcfIndex,
                genomeVersion = snpEffGenomeVersion,
                datadirZip = snpEffDatadirZip,
                outputPath = outputDir + "/" + name2 + ".annotated.vcf",
                configOptions = snpEffConfigOptions,
        }
    }

    output {
        Array[File] privateVcf1 = intersect.privateAVcf
        Array[File] privateVcf1Index = intersect.privateAVcfIndex
        Array[File] privateVcf2 = intersect.privateBVcf
        Array[File] privateVcf2Index = intersect.privateBVcfIndex
        Array[File] sharedVcf1 = intersect.sharedAVcf
        Array[File] sharedVcf1Index = intersect.sharedAVcfIndex
        Array[File] sharedVcf2 = intersect.sharedBVcf
        Array[File] sharedVcf2Index = intersect.sharedBVcfIndex
        Array[File] bcftoolsIsecReadmes = intersect.readme
        Array[File] vennFiles = vennStats.out
        Array[File] annotatedPrivateVcf1 = annotateUnique1.outputVcf
        Array[File] annotatedPrivateVcf2 = annotateUnique2.outputVcf
    }
}


task vennStats {
    input {
        File uniqueAVcf
        File uniqueBVcf
        File sharedVcf
        String nameA = basename(uniqueAVcf)
        String nameB = basename(uniqueBVcf)
        String outputPath = "venn.txt"
    }

    command <<<
        mkdir -p $(dirname '~{outputPath}')
        ONE_UNIQUES=$(bcftools view --no-header ~{uniqueAVcf} | wc -l)
        TWO_UNIQUES=$(bcftools view --no-header ~{uniqueBVcf} | wc -l)
        SHARED_UNIQUES=$(bcftools view --no-header ~{sharedVcf} | wc -l)
        echo "$ONE_UNIQUES" > one.txt
        echo "$TWO_UNIQUES" > two.txt
        echo "$SHARED_UNIQUES" > shared.txt
        echo "~{nameA}: $ONE_UNIQUES" > ~{outputPath}
        echo "~{nameB}: $TWO_UNIQUES" >> ~{outputPath}
        echo "shared: $SHARED_UNIQUES" >> ~{outputPath} 
    >>>
    output {
        Int uniqueA = read_int("one.txt")
        Int uniqueB = read_int("two.txt")
        Int shared = read_int("shared.txt")
        File out = outputPath
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        # 10 minutes should be ample as only the header and the first line are processed.
        time_minutes: 30
        memory: "512MiB"
    }
}
