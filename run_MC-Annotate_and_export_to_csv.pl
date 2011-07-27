#!/usr/bin/perl -w

$MCAnnotate = '/net/users/apetrov/software/mcannotate/MC-Annotate';
$Parser = '/net/users/apetrov/software/mcannotate/MC-Annotate/parser/parse_MCAnnotate.pl';

# folder with pdb files
$sourceFolder = shift(@ARGV);
# folder with MC-Annotate files and csv files
$destinationFolder = shift(@ARGV);
unless ( -d $destinationFolder ) {
    mkdir($destinationFolder);    
}

opendir(DIR, $sourceFolder) or die $!;

while (my $file = readdir(DIR)) {

    # Use a regular expression to ignore files beginning with a period
    next if ($file =~ m/^\./ or $file !~ m/\.pdb/);

    $command = $MCAnnotate . ' ' . $file . ' > ' . $destinationFolder . '/' . $file . '.txt';
    print $command , "\n";
#    system($command);

    $command = 'perl ' . $Parser . ' ' . $destinationFolder . '/' . $file . '.txt';
    print $command , "\n";
    
}

closedir(DIR);




