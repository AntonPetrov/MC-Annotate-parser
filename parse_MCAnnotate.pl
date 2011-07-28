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

$bps = 'MC-Annotate_basePairs_' . $pdbId . '_' . $pdbType . '.csv';
$bst = 'MC-Annotate_baseStacks_' . $pdbId . '_' . $pdbType . '.csv';
open( BPS, '>', $bps ) or die("Could not open $bps");
open( BST, '>', $bst ) or die("Could not open $bst");
#open( NBPS,'>', 'MC-Annotate_nearBasePairs.csv' ) or die('Could not open MC-Annotate_nearBasePairs.csv');

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
        if ( $line =~ /^((\w|'\d')(\d+)(\.\w)?) # $1-mcId,$2-chain,$3-number,$4-insCode
                       \s:\s
                       (\w)\s        # $5-residue
                       (?:\w+)_      # $6-puckerAtom
                       (?:\w+)\s     # $7-puckerQual
                       (?:\w+)       # $8-conf
                       /x ) {
                        
            $mcId = $1;
            $res  = $5;
            $mcId =~ s/('|\.)//g; #'      
            $ntMap{$mcId} = $res;
            # store the conformations if necessary      
        }
        

        # Parse base pairs
        # A419-A424 : C-G Ww/Ww pairing antiparallel cis XIX OR A1030.A-A1030.C : G-G O2'/Hh Ss/O2P pairing
        elsif ( $line =~ /^(\w|'\d') # $1-ch1
                           (\d+)     # $2-res1
                           (\.\w)?-  # $3-insCode1
                           (\w|'\d') # $4-ch2
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

            $ch1  = $1; $ch2 = $4; 
            $res1 = $2; $res2 = $5;
            $insCode1 = ( defined($3) ) ? substr($3,1,1) : '';
            $insCode2 = ( defined($6) ) ? substr($6,1,1) : '';
            $base1 = $7; $base2= $8;
            $MCpair = $9;
            $cistrans = $11;
            $ch1 =~ s/'//g; 
            $ch2 =~ s/'//g;

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$ch1,$res1,$base1,$insCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$ch2,$res2,$base2,$insCode2);

            $LWpair  = substr($cistrans,0,1) . substr($MCpair,0,1) . substr($MCpair,3,1);
            $rLWpair = substr($cistrans,0,1) . substr($MCpair,3,1) . substr($MCpair,0,1);
            $rMCpair = substr($MCpair,3,2) . '/' . substr($MCpair,0,2);

            print BPS join("\t",$ntId1,$LWpair,$MCpair,$ntId2) , "\n";
            print BPS join("\t",$ntId2,$rLWpair,$rMCpair,$ntId1) , "\n";

        }
        

        # Pares other MC-Annotate pairs
        # A1266-A1268 : G-A O2'/Hh Hh/O2P pairing OR A1281-A1282 : U-C O2P/Bh adjacent_5p pairing
        elsif ( $line =~ /^(\w|'\d') # $1-ch1
                           (\d+)     # $2-res1
                           (\.\w)?-  # $3-insCode1
                           (\w|'\d') # $4-ch2
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
                            
            $ch1  = $1; $ch2 = $4; 
            $res1 = $2; $res2 = $5;
            $insCode1 = ( defined($3) ) ? substr($3,1,1) : '';
            $insCode2 = ( defined($6) ) ? substr($6,1,1) : '';
            $base1 = $7; $base2= $8;
            $nearMCpair = $9;
            $ch1 =~ s/'//g; 
            $ch2 =~ s/'//g;

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_', $pdbId,$pdbType,$model,$ch1,$res1,$base1,$insCode1);
            $ntId2 = join('_', $pdbId,$pdbType,$model,$ch2,$res2,$base2,$insCode2);

            $rNearMCpair = '';
            @pairs = split(/\s/,$nearMCpair);
            foreach $pair (@pairs) {
                $pair =~ m/(.+)\/(.+)/;
                $rNearMCpair .= $2 . '/' . $1  . ' ';
            }
            $rNearMCpair =~ s/\s+$//;            

            $LWcompatible = '';

            print BPS join("\t",$ntId1,$LWcompatible,$nearMCpair,$ntId2) , "\n";
            print BPS join("\t",$ntId2,$LWcompatible,$rNearMCpair,$ntId1) , "\n";

        }


        # Parse base stacking
        # adjacent base stacking:
        # A1028-A1029 : adjacent_5p upward OR A1030.B-A1030.C : adjacent_5p upward
        # non-adjacent stacking:
        # A1346-A1348 : outward OR A1347-A1373 : inward pairing
        elsif ( $line =~ /^(\w|'\d')(\d+)(\.\w)?- # $1-ch1,$2-res1,$3-insCode1
                           (\w|'\d')(\d+)(\.\w)?  # $4-ch2,$5-res2,$6-insCode2
                           \s:\s
                           (?:\w+\s)?
                           (outward|inward|downward|upward) # $7-Xward
                           /x ) {
        
            $ch1 = $1;  $ch2 = $4;
            $res1 = $2; $res2 = $5;
            $insCode1 = ( defined($3) ) ? substr($3,1,1) : '';
            $insCode2 = ( defined($6) ) ? substr($6,1,1) : '';            
            $Xward = $7;                                   
            $ch1 =~ s/'//g;
            $ch2 =~ s/'//g;

            $mcId1 = $ch1 . $res1 . $insCode1;
            $mcId2 = $ch2 . $res2 . $insCode2; 

            # upward=s35, downward=s53, inward=s33, outward=s55
            switch ($Xward) {
            	case "upward"	{ $LWstack = 's35'; $rLWstack = 's53'; $rXward = 'downward'; }
            	case "downward"	{ $LWstack = 's53'; $rLWstack = 's35'; $rXward = 'upward'; }
            	case "inward"	{ $LWstack = 's33'; $rLWstack = 's33'; $rXward = 'inward'; }
            	case "outward"	{ $LWstack = 's55'; $rLWstack = 's55'; $rXward = 'outward'; }
            }

            # pdb_id,pdb_type,NT.ModelNum,NT.Chain,NT.Number,NT.Base,insCode,alternateId
            $ntId1 = join('_',$pdbId,$pdbType,$model,$ch1,$res1,$ntMap{$mcId1},$insCode1);
            $ntId2 = join('_',$pdbId,$pdbType,$model,$ch2,$res2,$ntMap{$mcId2},$insCode2);
            print BST join("\t",$ntId1,$LWstack,$Xward,$ntId2), "\n";
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
#close(NBPS);
