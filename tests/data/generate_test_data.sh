#!/bin/bash
# Generate test data for the integration tests

GENERATOR="$(dirname "$0")/generate_mutations.py"


function generate {
    $GENERATOR $@ > "${1}.fasta"
}

rm -f -- *.fasta

#K417N
generate K417N 22813:C
#K417T
generate K417T 22812:C
#N439K
generate N439K 22879:A
#N440K
generate N440K 22882:A
#G446S
generate G446S 22898:A
#Y449H
generate Y449H 22907:C
#Y449N
generate Y449N 22907:A
#L452R
generate L452R 22917:G
#Y453F
generate Y453F 22920:T
#S477N
generate S477N 22992:A
#T478K
generate T478K 22995:A
#T478R
generate T478R 22995:G
#E484A
generate E484A 23013:C
#E484K
generate E484K 23012:A
#E484Q
generate E484Q 23012:C
#F490R
generate F490R 23030:C 23031:G
#Q493K
generate Q493K 23039:A
#N501Y
generate N501Y 23063:T
#A570D
generate A570D 23271:A
#Q613H
generate Q613H 23401:C
#D614G
generate D614G 23403:G
#A626S
generate A626S 23438:T
#H655Y
generate H655Y 23525:T
#Q677H
generate Q677H 23593:C
#P681H
generate P681H 23604:A
#P681R
generate P681R 23604:G
#I692V
generate I692V 23636:G
#A701V
generate A701V 23664:T
#T716I
generate T716I 23709:T
#T732A
generate T732A 23756:G

