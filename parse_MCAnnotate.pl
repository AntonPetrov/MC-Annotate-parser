#!/usr/bin/perl -w

    # http://rna.bgsu.edu/FR3D/AnalyzedStructures/1J5E/1J5E_stacking.html
    # interactions are listed only in one direction - fixed.
    # insertion codes are present, but alternate ids are missing
    # for stacking only positions in chains are given - A1040
    # check that all stacking correspondences match up - done. Used 1J5E_A.

    # Alternative implementation: java parser written by Jose (Dumontier's lab);
    # http://code.google.com/p/semanticscience/source/browse/branches/jose/java/MC-Annotator/trunk/MC-Annotator/src/main/java/com/dumontierlab/mcannotator/bin/ParseMCAnnotate.java?r=677

    # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
    # does not recognize alternate ids - tested on 1DK1 (residue 27). Verify that it's always A. How?

use Switch;

if ( scalar(@ARGV) < 1 ) {
	die( "Input: text file created by MC-Annotate\n" );
}
chomp(@ARGV);
open( IN, '<', $ARGV[0] ) or die("Could not open $ARGV[0].\n");
open( BPS, '>', 'MC-Annotate_basePairs.csv' ) or die('Could not open MC-Annotate_basePairs.csv');
open( NBPS, '>', 'MC-Annotate_nearBasePairs.csv' ) or die('Could not open MC-Annotate_nearBasePairs.csv');
open( BST, '>', 'MC-Annotate_baseStacks.csv' ) or die('Could not open MC-Annotate_baseStacks.csv');

if ( $ARGV[0] !~ /(?<pdbId>[A-z0-9]{4})_(?<type>A|\d)/ ) {
    die('Check the input filename');
}
$pdbId = $+{pdbId};
if ( $+{type} eq 'A' ) {
    $pdbType = 'AU';
} else {
    $pdbType = "BA$+{type}";
}
print $pdbId , '_' , $pdbType , "\n";

$model = -1;
$/ = "Residue conformations";
while ( $record = <IN> ) {

    $model++;
    # skip the first block because it's always empty
    if ( $model == 0 ) {
        next;
    }
    print "Model $model\n";

    @lines = split("\n", $record);

    %ntMap = ();
    foreach $line (@lines) {

        # Parse residue conformations and store residue names
        # A1003 : G C3p_endo anti OR A1003.A : G C3p_endo anti
        if ( $line =~ /^(?<mcId>(?<ch>\w)(?<num>\d+)(?<insCode>\.\w)?)\s+:\s+(?<res>\w)\s+(?<puckerAtom>\w+)_(?<puckerQual>\w+)\s+(?<conf>\w+)/ ) {

            # store the conformations if necessary
            $ntMap{$+{mcId}} = $+{res};

        }

        # Parse base pairs
        # A419-A424 : C-G Ww/Ww pairing antiparallel cis XIX OR A1030.A-A1030.C : G-G O2'/Hh Ss/O2P pairing
        elsif ( $line =~ /^(?<ch1>\w)(?<res1>\d+)(?<insCode1>\.\w)?-(?<ch2>\w)(?<res2>\d+)(?<insCode2>\.\w)?\s+:\s+(?<base1>\w)-(?<base2>\w)\s+(?<MCpair>(w|h|s){2}\/(w|h|s){2})(.*?)pairing\s+(?<orientation>antiparallel|parallel)\s+(?<dir>cis|trans)/i ) {

            $insertionCode1 = ( defined($+{insCode1}) ) ? substr($+{insCode1},1,1) : '';
            $insertionCode2 = ( defined($+{insCode2}) ) ? substr($+{insCode2},1,1) : '';

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$+{ch1},$+{res1},$+{base1},$insertionCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$+{ch2},$+{res2},$+{base2},$insertionCode2);

            $LWpair  = substr($+{dir},0,1) . substr($+{MCpair},0,1) . substr($+{MCpair},3,1);
            $rLWpair = substr($+{dir},0,1) . substr($+{MCpair},3,1) . substr($+{MCpair},0,1);
            $rMCpair = substr($+{MCpair},3,2) . '/' . substr($+{MCpair},0,2);

            print BPS join("\t",$ntId1,$LWpair,$+{MCpair},$ntId2) , "\n";
            print BPS join("\t",$ntId2,$rLWpair,$rMCpair,$ntId1) , "\n";

        }

        # Other MC-Annotate pairs
        # A1266-A1268 : G-A O2'/Hh Hh/O2P pairing OR A1281-A1282 : U-C O2P/Bh adjacent_5p pairing
        elsif ( $line =~ /^(?<ch1>\w)(?<res1>\d+)(?<insCode1>\.\w)?-(?<ch2>\w)(?<res2>\d+)(?<insCode2>\.\w)?\s+:\s+(?<base1>\w)-(?<base2>\w)\s+(?<nearMCpair>(.*?))(\s+adjacent_.p){0,1}(\s+.{2,4}ward){0,1}(\s+pairing){0,1}\s+/ ) {

            $nearMCpair = $+{nearMCpair};
            $insertionCode1 = ( defined($+{insCode1}) ) ? substr($+{insCode1},1,1) : '';
            $insertionCode2 = ( defined($+{insCode2}) ) ? substr($+{insCode2},1,1) : '';

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$+{ch1},$+{res1},$+{base1},$insertionCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$+{ch2},$+{res2},$+{base2},$insertionCode2);

#            $rNearMCpair = '';
#            @pairs = split(/ /,$nearMCpair);
#            foreach $pair (@pairs) {
#                $pair =~ m/(.+)\/(.+)/;
#                $rNearMCpair .= $2 . '/' . $1  . ' ';
#            }
#            $rNearMCpair =~ s/\s+$//;
#            
#            print $nearMCpair , ' ' , $rNearMCpair , "\n";

            print NBPS join("\t",$ntId1,$nearMCpair,$ntId2) , "\n";
#                print NBPS join("\t",$ntId2,$rNearMCpair,$ntId1) , "\n";

        }

        # Parse base stacking
        # adjacent base stacking:
        # A1028-A1029 : adjacent_5p upward OR A1030.B-A1030.C : adjacent_5p upward
        # non-adjacent stacking:
        # A1346-A1348 : outward OR A1347-A1373 : inward pairing
        elsif ( $line =~ /^(?<mcId1>(?<ch1>\w)(?<res1>\d+)(?<insCode1>\.\w)?)-(?<mcId2>(?<ch2>\w)(?<res2>\d+)(?<insCode2>\.\w)?)\s+:\s+(\w+\s+)?(?<Xward>outward|inward|downward|upward)/ ) {

            $insertionCode1 = ( defined($+{insCode1}) ) ? substr($+{insCode1},1,1) : '';
            $insertionCode2 = ( defined($+{insCode2}) ) ? substr($+{insCode2},1,1) : '';

            # upward s35
            # downward s53
            # inward s33
            # outward s55
            switch ($+{Xward}) {
            	case "upward"	{ $LWstack = 's35'; $rLWstack = 's53'; $rXward = 'downward'; }
            	case "downward"	{ $LWstack = 's53'; $rLWstack = 's35'; $rXward = 'upward'; }
            	case "inward"	{ $LWstack = 's33'; $rLWstack = 's33'; $rXward = 'inward'; }
            	case "outward"	{ $LWstack = 's55'; $rLWstack = 's55'; $rXward = 'outward'; }
            }

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_',$pdbId,$pdbType,$model,$+{ch1},$+{res1},$ntMap{$+{mcId1}},$insertionCode1);
            $ntId2 = join('_',$pdbId,$pdbType,$model,$+{ch2},$+{res2},$ntMap{$+{mcId2}},$insertionCode2);
            print BST "$ntId1\t$LWstack\t$+{Xward}\t$ntId2\n$ntId2\t$rLWstack\t$rXward\t $ntId1\n";

        }

        else {
            # inspect the lines that were not parsed
            if ( length($line) > 15 and $line !~ /-+/ ) {
                print $line , "\n";
            }
            next;
        }

    }

    print "++++++++++++++++++\n";

}

close(IN);
close(BPS);
close(BST);
close(NBPS);

#            if ( defined($+{insCode1}) ) {
#                $insertionCode1 = substr($+{insCode1},1,1);
#            } else {
#                $insertionCode1 = '';
#            }
#            if ( defined($+{insCode2}) ) {
#                $insertionCode2 = substr($+{insCode2},1,1);
#            } else {
#                $insertionCode2 = '';
#            }
