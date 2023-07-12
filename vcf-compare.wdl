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

}

workflow VcfCompare {
    input {
        Array[CompareUnit] units
    }

    scatter (unit in units) {
        call bcftools.Norm as normalizeVcf1 {
            input:
                inputVcf = unit.vcf1,
                inputVcfIndex = unit.vcf1index,
        }
        call bcftools.Norm as normalizeVcf2 {
            input:
                inputVcf = unit.vcf2,
                inputVcfIndex = unit.vcf2index,
        }

        call bcftools.Isec as intersect {
            input: 
                aVcf = normalizeVcf1.outputVcf,
                aVcfIndex = normalizeVcf1.outputVcfIndex,
                bVcf = normalizeVcf2.outputVcf,
                bVcfIndex = normalizeVcf2.outputVcfIndex,
        }
    }

    output {
        Array[File] privateVcf1 = intersect.privateAVcf
        Array[File] privateVcf1Index = intersect.privateAVcfIndex
        Array[File] privateVcf2 = intersect.privateBVcf
        Array[File] privateVcf2Index = intersect.privateBVcfIndex
    }
}