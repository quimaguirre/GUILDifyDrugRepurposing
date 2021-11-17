#!/usr/bin/perl -w
#
# Script to read an input file and create R.script to be submitted for guildify network expansion
# it submit batches of '$limit' jobs and then sleeps for 30 minutes before to submitting new ones
#
# PARAMETERS
my $counter = 0;   
my $limit   = 20; # if combo --> the less than 10;
my $sleep   = 700; # 1800 sec --> 30 minutes

# PATHS

my $inputfile    = $ARGV[0]; # Input file with the info
my $Rscriptspath = "/home/quim/PHD/Projects/GUILDify/drug_repurposing/RScripts";


#
open(IN, "<$inputfile") || die "Cannot open $inputfile\n";

# Loop through file, create R.scripts and submit

while(<IN>){
#C0024305	lymphoma, non-hodgkin	20	20	22	chorioretinal abnormality	20	20	21
#C0018916	hemangioma	20	20	22
#C0024305	lymphoma, non-hodgkin	20	20	22

	next if ($_ =~ /^#/);
	chomp($_);
	my @l = split(/\t/);
	my $rscriptfile = $Rscriptspath . "/". $l[0] .".R";
	if (-e $rscriptfile){ # Skip if present
		#print $l[0] . " skipped as R.script is present\n";
		next;
	}
	# Create Rscript
	#create_Rscript_combo($rscriptfile, $_);
	create_Rscript_netscore($rscriptfile, $l[0]);
	#create_Rscript_zscore($rscriptfile, $_);
	#create_Rscript_short($rscriptfile, $_);
	# Execute Rscript
	#die;
	print "Submitting Rscript for BIANA id: ". $l[0] ."\n"; 
#die;
	system "Rscript '$rscriptfile'";
	#die;
	$counter++;
	if ($counter == $limit){
		$counter = 0;
		sleep($sleep);
	}
}


close(IN);

exit;



########### SUBRUTINES


sub create_Rscript_combo {

my ($lfile, $lkw) = @_; 

#library(guildifyR)
#a = data.frame(id=c("9385"), in.network=1)
#job.id =submit.job(a, species="9606", tissue="All", scoring.options=list(netscore=T, repetitionSelector=3, iterationSelector=2))

open(TMP, ">$lfile") || die "Cannot open $lfile\n";

print TMP qq[library(guildifyR)\n];
print TMP qq[result.table = query("$lkw", species="9606", tissue="all", network.source = "BIANA")\n];
print TMP qq[submit.job(result.table, species="9606", tissue="all", network.source="BIANA", scoring.options = list(netcombo = T))\n];

close(TMP);


}


################################################################################################################################################

sub create_Rscript_netscore {

my ($lfile, $lkw) = @_;

#library(guildifyR)
#a = data.frame(id=c("9385"), in.network=1)
#job.id =submit.job(a, species="9606", tissue="All", scoring.options=list(netscore=T, repetitionSelector=3, iterationSelector=2))
open(TMP, ">$lfile") || die "Cannot open $lfile\n";
print TMP qq[library(guildifyR)\n];
print TMP qq[result.table = query("$lkw", species="9606", tissue="all", network.source = "BIANA")\n];
print TMP qq[submit.job(result.table, species="9606", tissue="all", network.source="BIANA", scoring.options = list(netscore=T, repetitionSelector=3, iterationSelector=2))\n];
close(TMP);

}

################################################################################################################################################

sub create_Rscript_zscore {

my ($lfile, $lkw) = @_;

#library(guildifyR)
#a = data.frame(id=c("9385"), in.network=1)
#job.id =submit.job(a, species="9606", tissue="All", scoring.options=list(netscore=T, repetitionSelector=3, iterationSelector=2))
open(TMP, ">$lfile") || die "Cannot open $lfile\n";
print TMP qq[library(guildifyR)\n];
print TMP qq[result.table = query("$lkw", species="9606", tissue="all", network.source = "BIANA")\n];
print TMP qq[submit.job(result.table, species="9606", tissue="all", network.source="BIANA", scoring.options = list(netzcore=T, repetitionZelector=5))\n];
close(TMP);

}

#######################################################################################################################################

sub create_Rscript_short {

my ($lfile, $lkw) = @_;

#library(guildifyR)
#a = data.frame(id=c("9385"), in.network=1)
#job.id =submit.job(a, species="9606", tissue="All", scoring.options=list(netscore=T, repetitionSelector=3, iterationSelector=2))
open(TMP, ">$lfile") || die "Cannot open $lfile\n";
print TMP qq[library(guildifyR)\n];
print TMP qq[result.table = query("$lkw", species="9606", tissue="all", network.source = "BIANA")\n];
print TMP qq[submit.job(result.table, species="9606", tissue="all", network.source="BIANA", scoring.options = list(netshort=T))\n];
close(TMP);

}

