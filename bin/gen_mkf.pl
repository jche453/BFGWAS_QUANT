#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

Generate bfGWAS Makefile

=head1 SYNOPSIS

./gen_mkf.pl [options]

Options:
  -h      help messages
  -d      degub
  --mf    output make file

=head1 OPTIONS

=over 8

=item B<-help>

  -w         work directory : location for all output files
  -t         directory for the C++ executable file for Estep (running MCMC)
  --geno     specify genotype format: (gzipped) vcf, (gzipped) genotxt, or bed
  --gd       genotype file directory
  --ad       annotation file directory
  --ac       annotation classification code file
  --pheno    phenotype file
  --hyp      initial hyper parameter value file
  -f         file with a list of fileheads for genotype files
  -G         genotype format: GT(genotype data 0/0, 0/1, 1/1) or EC (dosage data)
  --maf      maf threshold: default 0.5%
  --pp       specify prior for the causal probability: default 1e-6
  --abgamma  specify inverse gamma prior for the effect size variance: default 0.1
  --win      window size of the neighborhood: default 100
  --em       number of EM iterations: default 5
  -b         number of burn ins: default 50,000
  -N         number of MCMC iterations: default 50,000
  --NL       number of MCMC iterations for the last EM iteration: default 50,000
  -c         compress genotype data (1) or not (0): default 0 (do not compress)
  --initype  specify initial model (1:start with top signal), (2: start with genome-wide significant signals), or default (3: stepwise selected signals)
  --rv       fixed residual variance value : default 1 (recomended)
  --smin     minimum number of variates per block in the model: default 0
  --smax     maximum number of variates per block in the model: default 5
  --mem      specify maximum memory usage: default 3000MB
  --time     specify time of runing MCMC per block per MCMC iteration: default 24hr
  -l         specify how job will be submited (options: slurm, mosix, local): default slurm
  --mf       output make file
  --nice     SLURM option for scheduling priority
  --xnode    compuation nodes to be excluded
  -j         specify SLURM job names
  --wnode    specify compuation nodes to be used
  --part     specify SLURM partition
  --AnnoNumber specify number of Annotation files

=back

=head1 DESCRIPTION

B<gen_mkf.pl> will generate a makefile for conducting Bayesian Functional GWASs by bfGWAS.

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "qsub";
my $wkDir=getcwd();
my $makeFile = "bfGWAS.mk";
my $genofile = "vcf";

my $toolE="/home/jchen/bfGWAS/BFGWAS_QUAN_ver2/bin/Estep_mcmc";
my $rs="/home/jchen/bfGWAS/BFGWAS_QUAN_ver2/bin/Mstep.r";
my $annoDir="/home/jyang/Collaborations/IrwinSAGE/BU_GWASs_CDSymptomDimensions/Anno_Files/PriceAnno/AA";
my $genoDir = "";
my $pheno="";
my $hyppar="/home/jyang/Collaborations/IrwinSAGE/BFGWAS_Analysis/hypval.current";
my $filelist = "/home/jyang/Collaborations/IrwinSAGE/SummaryData/FileHead.txt";


my $EM=5;
my $GTfield="GT";
my $maf="0.005";
my $rho="1";
my $smin="0";
my $smax="5";
my $win="100";
my $burnin="50000";
my $Nmcmc="50000";
my $NmcmcLast="50000";
my $compress=0;
my $initype="3";
my $rv="1";
my $pp="1e-6";
my $abgamma="0.1";

my $maxmem = "8000";
my $time = "24:00:00";
my $nice = "0";
my $jobid="";
my $xnode="";
my $wnode="";
my $part="nomosix";

my $seed="2017";
my $refld=0;
my $usextxld=0;
my $r2_level="0.001";
my $LDdir="/home/jyang/Collaborations/IrwinSAGE/BFGWAS_Analysis/RefLD/AA";
my $Scoredir="/mnt/YangFSS/data/BU_GWASs_CDSymptomDimensions/BFGWAS/AA/ScoreStat";

my $N="3272";
my $pv="6.598844";
my $EMdir="";

my $AnnoNumber="0";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
                'w:s'=>\$wkDir, 'Estep:s' =>\$toolE, 'ad:s'=>\$annoDir,
                'geno:s'=>\$genofile,
                'gd:s'=>\$genoDir, 'pheno:s'=>\$pheno,
                'hyp:s'=>\$hyppar, 'G:s'=>\$GTfield, 'maf:s'=>\$maf,
                'smin:s'=>\$smin, 'smax:s'=>\$smax, 'win:s'=>\$win,
                'b:s'=>\$burnin, 'N:s'=>\$Nmcmc, 'NL:s'=>\$NmcmcLast,
                'c:i'=>\$compress, 'n:i'=>\$N,
                'initype:s'=>\$initype, 'rv:s'=>\$rv,
                'pp:s'=>\$pp,'abgamma:s'=>\$abgamma,
                'mem:s'=>\$maxmem, 'time:s'=>\$time, 'f:s'=>\$filelist, 'em:i'=>\$EM, 'rs:s'=>\$rs,
                'l:s'=>\$launchMethod, 'mf:s'=>\$makeFile, 'nice:s'=>\$nice,
                'j:s' =>\$jobid, 'xnode:s'=>\$xnode, 'wnode:s'=>\$wnode, 'part:s'=>\$part,
                'refLD:s'=>\$refld, 'seed:s'=>\$seed, 'r2:s'=>\$r2_level,
                'LDdir:s'=>\$LDdir, 'Scoredir:s'=>\$Scoredir, 'EMdir:s'=>\$EMdir,
                'pv:s'=>\$pv, 'usextxLD:i'=>\$usextxld, 'AnnoNumber:s'=>\$AnnoNumber )
  || !defined($wkDir) || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


if ($launchMethod ne "local" && $launchMethod ne "slurm" && $launchMethod ne "mosix" && $launchMethod ne "qsub")
{
    print STDERR "Launch method has to be local, slurm, mosix, or qsub\! \n";
    exit(1);
}

##############
#print options
##############
printf("Options\n");
printf("\n");
printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkDir);
print "Estep: ", $toolE, "\n";
print "genoDir: ", $genoDir,
      "\nannoDir: ", $annoDir,
        "\npheno: ", $pheno,
      "\n LDdir: ", $LDdir,
       "\n Scoredir: ", $Scoredir,
        "hyppar: ", $hyppar,
        "\nfileheads: ", $filelist, "\n",
        "Rscript: ", $rs, "\n",
        "run_Estep.sh and run_Mstep.sh directory: ", $EMdir, "\n",
        "AnnoNumber: ", $AnnoNumber, "\n",
      #  "GTfield ", $GTfield, "; maf ", $maf, "; smin ", $smin, "\n",
      #  "smax ", $smax, "; win ", $win, "; burnin ", $burnin, "; Nmcmc ", $Nmcmc, "\n",
      #  "NmcmcLast ", $NmcmcLast, "; compress ", $compress, "; initype ", $initype, "\n",
       "rv ", $rv, "; pp ", $pp, "; abgamma ", $abgamma, "\n";
printf("\n");

my $comp = "";
if ($compress != 0){
  $comp = "-comp"
}

my $refLD = "";
if ($refld != 0){
  $refLD = "-refLD"
}

my $usextxLD = "";
if ($usextxld != 0){
  $usextxLD = "-usextxLD"
}

#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkDir);

########################################
# Initial Set up before EM iterations
########################################

### prepare files before EM_MCMC
my $hypcurrent="$wkDir/hypval.current";
$tgt = "$wkDir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkDir/output $wkDir/Eoutput $wkDir/OUT";
push(@cmd, "mkdir -p $wkDir/output $wkDir/Eoutput $wkDir/OUT");
push(@cmd, "cp -f $hyppar $hypcurrent");
push(@cmd, "> $wkDir/Eoutput/EM\_result.txt");
push(@cmd, "> $wkDir/Rout.txt");
makeJob("local", $tgt, $dep, $wkDir, @cmd);


###### EM step 0 without dependencies ###########
my $i;
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filelist)
    or die "Can not open $filelist \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
}
else{ print "Total \# of fileheads: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $i=0;
        $line = $filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/pre_em.OK";
        if($genofile eq "vcf"){
           @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -p $pheno -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp  -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "genotxt"){
          @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "bed") {
          @cmd = "$toolE -bfile $genoDir/$line -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
        }elsif($genofile eq "sumstat"){
          @cmd = "$EMdir/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $annoDir $AnnoNumber $initype ";
          # @cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
        }
        else{
          die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
        }

        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

my $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
my $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
my $logfile="$wkDir/Eoutput/log$i.txt";
my $Annofile="$wkDir/Eoutput/AnnoScore$i.txt";

$tgt = "$wkDir/Eoutput/cp_param$i.OK";
$dep = "$premcmcOK";
@cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort\` > $paramfile";
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep AnnoScore | sort\` > $Annofile");
makeJob("local", $tgt, $dep, $wkDir, @cmd);

$tgt = "$wkDir/R$i.OK";
$dep = "$wkDir/Eoutput/cp_param$i.OK $wkDir/pre_em.OK";
@cmd = "$EMdir/run_Mstep.sh $pv $i $wkDir $hypcurrent $Annofile $abgamma";
# @cmd = "Rscript --vanilla $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent >> $wkDir/Rout.txt";
makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);


####### With dependencies of previous output
my $ipre="";

for $i (1..$EM){

    $ipre=$i-1; $premcmcOK="";

    for(my $j=0; $j< @filehead; ++$j){
        $line=$filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.$i.OK ";
        $tgt = "$wkDir/OUT/$line.$i.OK";
        $dep = "$wkDir/R$ipre.OK";
        if($i < $EM){
          if($genofile eq "vcf"){
            @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -p $pheno -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt ";
          }elsif ($genofile eq "genotxt"){
            @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
          }elsif ($genofile eq "bed") {
            @cmd = "$toolE -bfile $genoDir/$line -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
          }elsif($genofile eq "sumstat"){
            @cmd = "$EMdir/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $annoDir $AnnoNumber $initype";
           # @cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
          }else{
            die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
          }
        } elsif ($i == $EM){
            if($genofile eq "vcf"){
              @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -a $annoDir/Anno\_$line.gz -p $pheno -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $NmcmcLast $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt ";
            }elsif ($genofile eq "genotxt"){
              @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $NmcmcLast $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
            }elsif ($genofile eq "bed") {
              @cmd = "$toolE -bfile $genoDir/$line -a $annoDir/Anno\_$line.gz -hfile $hypcurrent -GTfield $GTfield -maf $maf -bvsrm -rmin $rho -rmax $rho -smin $smin -smax $smax -win $win -o $line -w $burnin -s $NmcmcLast $comp -initype $initype -rv $rv -AnnoNumber $AnnoNumber > $wkDir/OUT/$line.output.txt";
            }elsif($genofile eq "sumstat"){
              @cmd = "$EMdir/run_Estep.sh $wkDir $line $N $pv $burnin $Nmcmc $LDdir $Scoredir $hypcurrent $annoDir $AnnoNumber $initype";
           #@cmd = "$toolE -score $genoDir/$line.SS.score.txt.gz -inputSS $refLD -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a $annoDir/Anno\_$line.gz -fcode $annoCode -hfile $hypcurrent -n $N -pv $pv -maf $maf -r2 $r2_level -bvsrm -smin $smin -smax $smax -win $win -o $line -w $burnin -s $Nmcmc $comp -initype $initype -seed ${seed} > $wkDir/OUT/$line.output.txt";
          }else{
            die "genoDir need to be one of the vcf, genotxt, bed, or sumstat file types\!\n"
          }
      }
        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

    $paramfile="$wkDir/Eoutput/paramtemp$i.txt";
    $hypfile="$wkDir/Eoutput/hyptemp$i.txt";
    $logfile="$wkDir/Eoutput/log$i.txt";
    $Annofile="$wkDir/Eoutput/AnnoScore$i.txt";

  $tgt = "$wkDir/Eoutput/cp_param$i.OK";
  $dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkDir/output/** | grep paramtemp | sort \` > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep hyptemp | sort\` > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
  push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep AnnoScore | sort\` > $Annofile");
  makeJob("local", $tgt, $dep, $wkDir, @cmd);

  $tgt = "$wkDir/R$i.OK";
  $dep = "$wkDir/Eoutput/cp_param$i.OK";
  @cmd = "$EMdir/run_Mstep.sh $pv $i $wkDir $hypcurrent $Annofile $abgamma";
  #@cmd = "Rscript --vanilla $rs $hypfile $i $pp $abgamma $wkDir/Eoutput/EM_result.txt $hypcurrent >> $wkDir/Rout.txt";
  makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);

}


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
mkpath("$wkDir/tmp");
push(@tgts, "clean_err");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/tmp/*.err");


push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/*.OK $wkDir/Eoutput/*.OK $wkDir/OUT/*.OK $wkDir/tmp/*.err");


for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, $wkd, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, $wkd, @cmd);
    }
    elsif($method eq "mosix")
    {
      makeMosix($tgt, $dep, $wkd, @cmd);
    }
    elsif($method eq "qsub")
    {
      makeQsub($tgt, $dep, $wkd, @cmd);
    }
}

#run qsub jobs
sub makeQsub
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tqsub -q b.q -j y -wd $wkd/ -N $jobid -o ${wkd}/tmp/\\\$\$JOB_NAME.o\\\$\$JOB_ID.err -sync y $c \n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}


#run mosix jobs
sub makeMosix
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
      if($wnode){
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -j$wnode " . $c . "\n";
      }else{
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -b " . $c . "\n";
      }

    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tsrun --exclude=$xnode --partition=$part --mem-per-cpu\=$maxmem --time\=$time --nice\=$nice --error\=$wkDir/tmp/\%N.\%j.err -J $jobid -D $wkd $c \n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
