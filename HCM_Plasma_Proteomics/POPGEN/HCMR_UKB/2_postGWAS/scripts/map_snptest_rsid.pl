use strict;
#Script takes time to load master db. If many files are to be processed, contact Anuj. Can loop over the files to avoid loading hash multiple times!
#29 Jan 2021.

#Usage:
#Go to the folder where you have your input file that needs appending with rsid
#perl /well/PROCARDIS/datasets/database/map_snptest_rsid.pl <inputfilename> <outputfilename> <colnum containing snptestid -1>
#column number where snptestid is stored minus 1 (perl speak)

#If Memory::Usage is not working after the 2 commands, then just comment out lines which has "####### this line"
#script takes 32GB RAM
#make sure you run these 2 commands before 
#module load Perl/5.28.1-GCCcore-8.2.0
#export PERL5LIB=/well/PROCARDIS/bin/perl_modules/lib/perl5/site_perl/5.28.1/:$PERL5LIB


use Memory::Usage;                       ####### this line

my $mu = Memory::Usage->new();           ####### this line
my $in=shift;#with header
my $out=shift;
my $snptest_col=shift;#column number where snptestid is stored minus 1
#read master db
my %rs=();#snptestid as key, rsid as value (the first rsid if multiple rsids encountered)
#1	10177	A	AC	1:10177_AC_A	1:10177_A_AC	+-----------	NA
#1	10178	A	AC	1:10178_AC_A	1:10178_A_AC	+-----------	rs99
#1	10235	T	TA	1:10235_T_TA	1:10235_TA_T	+-----------	rs100,rs101
open (M,"/well/PROCARDIS/datasets/database/master1.db.b37.txt");
open (IF, $in) or die "Cannot open input file";
open (OF, ">$out") or die "Cannot open output file";
my $h=<IF>;
chomp $h;
print OF "rsid\t$h\n";
$mu->record('start');                    ####### this line
while (my $l=<M>){
  chomp $l;
  my @cells=split(/\t/, $l);
  my @rsids=split(",",$cells[7]);
  $rs{$cells[4]}=$rsids[0];
  $rs{$cells[5]}=$rsids[0];
}
close M;
print "Master database read\n";
print "Hash size=",scalar keys %rs,"\n";
$mu->record('after');                    ####### this line
$mu->dump();                             ####### this line
#read input dataset
my $lno=0;
my $f=0;
while (my $l=<IF>){
  $lno++;
  chomp $l;
  my @cells=split(/\s+/,$l);
  my $st=$cells[$snptest_col];
  if (exists $rs{$st} && $rs{$st} ne "NA"){
    $f++;
    print OF "$rs{$st}\t$l\n";
  }else{
    print OF "$st\t$l\n";
  }
}
print "Input file read\nNumber of lines=$lno\nNumber of snptest ids found=$f\n";
close IF;
close OF;
exit;
