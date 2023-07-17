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
    }

    scatter (unit in units) {
        call retrieveSamplesFromVcf as retrieveSample1Id {
            input:
                vcf = unit.vcf1 
        }
        String prefix = outputDir + "/" + retrieveSample1Id.samples[0]
        String name1 = select_first([unit.name1, "1"])
        String name2 = select_first([unit.name2, "2"])

        call retrieveSamplesFromVcf as retrieveSample2Id {
            input:
                vcf = unit.vcf2
        }

        call bcftools.Norm as normalizeVcf1 {
            input:
                inputVcf = unit.vcf1,
                inputVcfIndex = unit.vcf1index,
                outputPath = prefix + "." + name1 + ".normalized.vcf.gz"
        }
        call bcftools.Norm as normalizeVcf2 {
            input:
                inputVcf = unit.vcf2,
                inputVcfIndex = unit.vcf2index,
                outputPath = prefix + "." + name2 + ".normalized.vcf.gz"
        }

        call bcftools.Isec as intersect {
            input: 
                aVcf = normalizeVcf1.outputVcf,
                aVcfIndex = normalizeVcf1.outputVcfIndex,
                bVcf = normalizeVcf2.outputVcf,
                bVcfIndex = normalizeVcf2.outputVcfIndex,
                prefix = prefix 
        }
        
        call vennStats as vennStats {
            input:
                uniqueAVcf = intersect.privateAVcf,
                uniqueBVcf = intersect.privateBVcf,
                sharedVcf = intersect.sharedAVcf,
                nameA = name1,
                nameB = name2,
                outputPath = prefix + "/venn.txt"
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
    }
}

task retrieveSamplesFromVcf {
    input {
        File vcf
    }

    command <<<
        bcftools query --list-samples ~{vcf}
    >>>

    output {
        Array[String] samples = read_lines(stdout())
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        # 10 minutes should be ample as only the header and the first line are processed.
        time_minutes: 10
        memory: "512MiB"
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
