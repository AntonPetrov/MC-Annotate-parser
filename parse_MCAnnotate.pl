#!/usr/bin/perl -w

    # http://rna.bgsu.edu/FR3D/AnalyzedStructures/1J5E/1J5E_stacking.html
    # interactions are listed only in one direction - fixed.
    # insertion codes are present, but alternate ids are missing
    # for stacking only positions in chains are given - A1040
    # check that all stacking correspondences match up - done. Used 1J5E_A.
    # Had to replace named capture groups with $1 etc because perl on the server is at v5.8

    # Alternative implementation: java parser written by Jose (Dumontier's lab);
    # http://code.google.com/p/semanticscience/source/browse/branches/jose/java/MC-Annotator/trunk/MC-Annotator/src/main/java/com/dumontierlab/mcannotator/bin/ParseMCAnnotate.java?r=677

    # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
    # does not recognize alternate ids - tested on 1DK1 (residue 27). Verify that it's always A. How?

use Switch;

if ( scalar(@ARGV) < 1 ) {
	die( "Input: text file created by MC-Annotate\n" );
}
chomp(@ARGV);
open( IN,  '<', $ARGV[0] ) or die("Could not open $ARGV[0].\n");
open( BPS, '>', 'MC-Annotate_basePairs.csv' ) or die('Could not open MC-Annotate_basePairs.csv');
open( NBPS,'>', 'MC-Annotate_nearBasePairs.csv' ) or die('Could not open MC-Annotate_nearBasePairs.csv');
open( BST, '>', 'MC-Annotate_baseStacks.csv' ) or die('Could not open MC-Annotate_baseStacks.csv');

if ( $ARGV[0] !~ /([A-z0-9]{4}) # $1, pdbId
                  _
                  (A|\d+)       # $2, A for Asymmetric unit, numbers for biological assemblies
                  /x ) {
    die('Check the input filename');
}
$pdbId = $1;
if ( $2 eq 'A' ) {
    $pdbType = 'AU';
} else {
    $pdbType = "BA$2";
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
        if ( $line =~ /^((\w)(\d+)(\.\w)?) # $1-mcID,$2-chain,$3-number,$4-insCode
                       \s+:\s+
                       (\w)\s+             # $5-residue
                       (\w+)               # $6-puckerAtom
                       _
                       (\w+)\s+            # $7-puckerQual
                       (\w+)               # $8-conf
                       /x ) {

            # store the conformations if necessary
            $ntMap{$1} = $5;
        }

        # Parse base pairs
        # A419-A424 : C-G Ww/Ww pairing antiparallel cis XIX OR A1030.A-A1030.C : G-G O2'/Hh Ss/O2P pairing
        elsif ( $line =~ /^(\w)      # $1-ch1
                           (\d+)     # $2-res1
                           (\.\w)?-  # $3-insCode1
                           (\w)      # $4-ch2
                           (\d+)     # $5-res2
                           (\.\w)?   # $6-insCode2
                           \s:\s
                           (\w)-     # $7-base1
                           (\w)\s    # $8-base2
                           ((?:w|h|s){2}\/(?:w|h|s){2}) # $9-MCpair
                           (?:.*?)pairing\s             # additional descriptions
                           (antiparallel|parallel)\s    # $10-orientation
                           (cis|trans)                  # $11-cis,trans
                           /ix ) {

            $insertionCode1 = ( defined($3) ) ? substr($3,1,1) : '';
            $insertionCode2 = ( defined($6) ) ? substr($6,1,1) : '';

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$1,$2,$7,$insertionCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$4,$5,$8,$insertionCode2);

            $LWpair  = substr($11,0,1) . substr($9,0,1) . substr($9,3,1);
            $rLWpair = substr($11,0,1) . substr($9,3,1) . substr($9,0,1);
            $rMCpair = substr($9,3,2) . '/' . substr($9,0,2);

            print BPS join("\t",$ntId1,$LWpair,$9,$ntId2) , "\n";
            print BPS join("\t",$ntId2,$rLWpair,$rMCpair,$ntId1) , "\n";

        }

        # Other MC-Annotate pairs
        # A1266-A1268 : G-A O2'/Hh Hh/O2P pairing OR A1281-A1282 : U-C O2P/Bh adjacent_5p pairing
        elsif ( $line =~ /^(\w)      # $1-ch1
                           (\d+)     # $2-res1
                           (\.\w)?-  # $3-insCode1
                           (\w)      # $4-ch2
                           (\d+)     # $5-res2
                           (\.\w)?   # $6-insCode2
                           \s:\s
                           (\w)-     # $7-base1
                           (\w)\s    # $8-base2
                           ((?:(?:\w|'){2,3}\/(?:\w|'){2,3}\s)+)  # $9-nearMCpair
                           (?:adjacent_.p){0,1}   # don't include adjacency
                           (?:\s.{2,4}ward){0,1}  # don't include stacking
                           (?:\spairing){0,1}
                           /x ) {

            $nearMCpair = $9;
            $insertionCode1 = ( defined($3) ) ? substr($3,1,1) : '';
            $insertionCode2 = ( defined($6) ) ? substr($6,1,1) : '';

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$1,$2,$7,$insertionCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$4,$5,$8,$insertionCode2);

            $rNearMCpair = '';
            @pairs = split(/\s/,$nearMCpair);
            foreach $pair (@pairs) {
                $pair =~ m/(.+)\/(.+)/;
                $rNearMCpair .= $2 . '/' . $1  . ' ';
            }
            $rNearMCpair =~ s/\s+$//;
            
#            print $nearMCpair , ' ' , $rNearMCpair , "\n";

            print NBPS join("\t",$ntId1,$nearMCpair,$ntId2) , "\n";
            print NBPS join("\t",$ntId2,$rNearMCpair,$ntId1) , "\n";

        }

        # Parse base stacking
        # adjacent base stacking:
        # A1028-A1029 : adjacent_5p upward OR A1030.B-A1030.C : adjacent_5p upward
        # non-adjacent stacking:
        # A1346-A1348 : outward OR A1347-A1373 : inward pairing
        elsif ( $line =~ /^((\w)(\d+)(\.\w)?)- # $1-mcId1,$2-ch1,$3-res1,$4-insCode1
                           ((\w)(\d+)(\.\w)?)  # $5-mcId2,$6-ch2,$7-res2,$8-insCode2
                           \s:\s
                           (?:\w+\s)?
                           (outward|inward|downward|upward) # $9-Xward
                           /x ) {

            $insertionCode1 = ( defined($4) ) ? substr($4,1,1) : '';
            $insertionCode2 = ( defined($8) ) ? substr($8,1,1) : '';

            # upward=s35, downward=s53, inward=s33, outward=s55
            switch ($9) {
            	case "upward"	{ $LWstack = 's35'; $rLWstack = 's53'; $rXward = 'downward'; }
            	case "downward"	{ $LWstack = 's53'; $rLWstack = 's35'; $rXward = 'upward'; }
            	case "inward"	{ $LWstack = 's33'; $rLWstack = 's33'; $rXward = 'inward'; }
            	case "outward"	{ $LWstack = 's55'; $rLWstack = 's55'; $rXward = 'outward'; }
            }

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_',$pdbId,$pdbType,$model,$2,$3,$ntMap{$1},$insertionCode1);
            $ntId2 = join('_',$pdbId,$pdbType,$model,$6,$7,$ntMap{$5},$insertionCode2);
            print BST join("\t",$ntId1,$LWstack,$9,$ntId2), "\n";
            print BST join("\t",$ntId2,$rLWstack,$rXward,$ntId1), "\n";

        }

        else {
            # inspect the lines that were not parsed
            if ( length($line) > 15 and $line !~ /-{2,}/ ) {
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
