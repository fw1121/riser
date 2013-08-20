#!/usr/bin/perl 

use strict;
use Bio::SeqIO;
use List::Util qw[min max];


if ($#ARGV != 2) {
    print "Error: no files provided\n";
    exit(1);
}

my $gbkfile=$ARGV[0];
my $outdir=$ARGV[1];
my $name=$ARGV[2];

open(OUTTRANSCRIPTABLE, ">$outdir/transcripts.txt"); #open for write, overwrite
open(OUTSEQ, ">$outdir/sequence.fa"); #open for write, overwrite

my $in  = Bio::SeqIO->new(-file => $gbkfile,
         '-format' => 'genbank');

# loop over diff sequences 
while (my $seq_object = $in->next_seq) {
    # print("$refgenome_id\n");
    my @gene_names = ();
    my @cds_start_jctn = ();
    my @cds_end_jctn = ();
    my @cds_start_g = ();
    my @cds_end_g = ();
    my @cds_start = ();
    my @cds_end = ();
    my @cds_strand=();
    my @gene_names2 = ();
    my @gene_start = ();
    my @gene_end = ();
    my @gene_start_jctn = ();
    my @gene_end_jctn = ();
    my @gene_strand=();
    my @gene_GI=();
    my @gene_ID=();
    my @cds_GENEID=();
    my @cds_GI=();
    my $cds;
    my @tag_values;
    my @grepGeneID;
    my @GeneIDnum;
    my $gene;
    my @GInum;
    my @grepGI;
    my @grepGENEID;
    my @GENIDnum;
    my $gene_name;
    my $transcript_strand;
    my $gene_id;
    my $cdsStart;
    my $cdcsEnd;
    my $exonCount;
    my $exonStarts;
    my $cdcsEnd;
    my $exonEnds;
    my $gene_index;

# get the accession nmb

    my $acc_nmb = $seq_object->accession_number;
    my $seqout = $seq_object->seq();


    print OUTSEQ ">"."$acc_nmb\n";
    print OUTSEQ "$seqout";


# get the features 
my @feats = $seq_object->get_all_SeqFeatures(); 

F1:foreach $cds (@feats) {
    next F1 unless ($cds->primary_tag() eq 'CDS');
if ($cds->location->isa('Bio::Location::SplitLocationI'))
{
for my $location ($cds->location->sub_Location ) {
    push @cds_start_jctn, $location->start;    
    push @cds_end_jctn, $location->end;
    push @cds_start_g,'NA';
    push @cds_end_g,'NA';
    push @cds_start, $location->start;    
    push @cds_end, $location->end;
    push @cds_strand, $location->strand;
    if ($cds->has_tag('gene'))
    {
	push @gene_names, $cds->get_tag_values('gene');
    } elsif ($cds->has_tag('product'))
    {
	push @gene_names, $cds->get_tag_values('product');
    } else {push(@gene_names, 'NA');}


    if ($cds->has_tag('db_xref'))
    {
	@tag_values = $cds->get_tag_values('db_xref');
   
	if(grep(/^GeneID:/, @tag_values))
	{
	    @grepGeneID = grep(/^GeneID:/, @tag_values);
	    @GeneIDnum = split(/:/,$grepGeneID[0]);
	    push @cds_GENEID,$GeneIDnum[1];
	    
	} else {push(@cds_GENEID, 'NA');}
    
	if(grep(/^GI:/, @tag_values))
	{
	    @grepGI = grep(/^GI:/, @tag_values);
	    @GInum = split(/:/,$grepGI[0]);
	    push @cds_GI,$GInum[1];} else {push(@cds_GI, 'NA');} 
    } else {push(@cds_GENEID, 'NA'); push(@cds_GI, 'NA');}

}
     
} else {

    push @cds_start, $cds->location->start;
    push @cds_end, $cds->location->end;
    push @cds_start_g, $cds->location->start;
    push @cds_end_g, $cds->location->end;
    push @cds_start_jctn, 'NA';
    push @cds_end_jctn, 'NA';	
    push @cds_strand, $cds->location->strand;	
    if ($cds->has_tag('gene'))
    {
	push @gene_names, $cds->get_tag_values('gene');
    } elsif ($cds->has_tag('product'))
    {
	push @gene_names, $cds->get_tag_values('product');
    } else {push(@gene_names, 'NA');}


   if ($cds->has_tag('db_xref'))
    {
	@tag_values = $cds->get_tag_values('db_xref');
   
	if(grep(/^GeneID:/, @tag_values))
	{
	    @grepGeneID = grep(/^GeneID:/, @tag_values);
	    @GeneIDnum = split(/:/,$grepGeneID[0]);
	    push @cds_GENEID,$GeneIDnum[1];
	    
	} else {push(@cds_GENEID, 'NA');}
    
	if(grep(/^GI:/, @tag_values))
	{
	    @grepGI = grep(/^GI:/, @tag_values);
	    @GInum = split(/:/,$grepGI[0]);
	    push @cds_GI,$GInum[1];} else {push(@cds_GI, 'NA');} 
    } else {push(@cds_GENEID, 'NA'); push(@cds_GI, 'NA');}


}

}


F2:foreach $gene (@feats) {
    next F2 unless ($gene->primary_tag() eq 'gene');

if ($gene->location->isa('Bio::Location::SplitLocationI'))
{
for my $location ($gene->location->sub_Location ) {
    push @gene_start, $location->start;    
    push @gene_end, $location->end;
    push @gene_strand, $location->strand;
    if ($gene->has_tag('gene')) {push @gene_names2, $gene->get_tag_values('gene');} else {push(@gene_names2, 'NA');} 
    if ($gene->has_tag('db_xref')){
	@tag_values = $gene->get_tag_values('db_xref');
        @grepGI = grep(/^GI:/, @tag_values);
        
    if(grep(/^GI:/, @tag_values))
     {
     @grepGI = grep(/^GI:/, @tag_values);
     @GInum = split(/:/,$grepGI[0]);
     push @gene_GI,@GInum[1];} else {push(@gene_GI, 'NA');} 
     
    if(grep(/^GeneID:/, @tag_values))
     {
     @grepGENEID = grep(/^GeneID:/, @tag_values);
     @GENIDnum = split(/:/,$grepGENEID[0]);
     push @gene_ID,@GENIDnum[1];} else {push(@gene_ID, 'NA');} 
    } 
   
    else {push(@gene_GI, 'NA');
          push(@gene_ID, 'NA');
    }
 
  
}

} else {
    push @gene_start, $gene->location->start;
    push @gene_end, $gene->location->end;
    push @gene_strand, $gene->location->strand;
    if ($gene->has_tag('gene')) {push @gene_names2, $gene->get_tag_values('gene');} else {push(@gene_names2, 'NA');} 
    if ($gene->has_tag('db_xref')){
	@tag_values = $gene->get_tag_values('db_xref');

    if(grep(/^GI:/, @tag_values))
     {
     @grepGI = grep(/^GI:/, @tag_values);
     @GInum = split(/:/,$grepGI[0]);
     push @gene_GI,@GInum[1];} else {push(@gene_GI, 'NA');} 
    
     if(grep(/^GeneID:/, @tag_values))
     {
     @grepGENEID = grep(/^GeneID:/, @tag_values);
     @GENIDnum = split(/:/,$grepGENEID[0]);
     push @gene_ID,@GENIDnum[1];} else {push(@gene_ID, 'NA');} 
    } 
 
     else {
           push(@gene_GI, 'NA');
           push(@gene_ID, 'NA'); 
     }

   }

}


# print "@cds_GENEID\n";
#print "size: ". @gene_names ."\n";
#print "size: ". @gene_names2 ."\n";
#print "@gene_names\n";
#print "@cds_start\n";

if(@gene_names2 == 0)
{
    print "Warning: no genes detected in the GenBank file. Will only consider CDSs.";
    my $gene_index = 0; 

    foreach $gene_name (@gene_names) {

	if($gene_name ne 'NA')
	{
        
	    # get the info for each CDS  
	    my $search_for =   $gene_name;
	    my @index = grep { $gene_names[$_] eq $search_for } 0..$#gene_names;   
	    my $cdsStart = min(@cds_start[@index]);    
	    my $cdcsEnd = max(@cds_end[@index]);       
	    my $exonCount = @index;
	    my $exonStarts = join(',', @cds_start[@index]);
	    my $exonEnds = join(',', @cds_end[@index]);
   
	    #print "@cds_start[@index]\n";
	    #print "$exonStarts\n";

	    if (@cds_strand[$gene_index] == 1) {$transcript_strand="+"}    
	    if (@cds_strand[$gene_index] == -1) {$transcript_strand="+"}   
   
	    print OUTTRANSCRIPTABLE "@cds_GI[$gene_index]"."\t"."$gene_name"."\t"."$acc_nmb"."\t"."$transcript_strand"."\t"."$cdsStart"."\t"."$cdcsEnd"."\t"."$cdsStart"."\t"."$cdcsEnd"."\t"."$exonCount"."\t"."$exonStarts"."\t"."$exonEnds"."\n";
    

	}

	$gene_index++;

    }


} elsif (@gene_names2 > 0) {

    foreach $gene_id (@gene_ID) {
	if($gene_id ne 'NA')
	{

	    my $search_for =   $gene_id;
    	    my @index = grep { $cds_GENEID[$_] eq $search_for } 0..$#cds_GENEID;     
	    
	    if (@cds_start_g[@index[0]] ne 'NA'){
		$cdsStart = @cds_start_g[@index];  
		$cdcsEnd = @cds_end_g[@index];       
		$exonCount = @index;
		$exonStarts = $cdsStart;
		$exonEnds = $cdcsEnd;
	    } else {
		$exonStarts = join(',', @cds_start_jctn[@index]);
		$exonEnds = join(',', @cds_end_jctn[@index]);
		$exonCount = @index;
		$cdsStart = min(@cds_start_jctn[@index]);
		$cdcsEnd = max(@cds_end_jctn[@index]);
	    }
	
	    $gene_name=@gene_names[@index];

 
	    if (@gene_strand[$gene_index] == 1) {$transcript_strand="+"}    
	    if (@gene_strand[$gene_index] == -1) {$transcript_strand="-"}   

	    print OUTTRANSCRIPTABLE "@gene_ID[$gene_index]"."\t"."$gene_name"."\t"."$acc_nmb"."\t"."$transcript_strand"."\t"."@gene_start[$gene_index]"."\t"."@gene_end[$gene_index]"."\t"."$cdsStart"."\t"."$cdcsEnd"."\t"."$exonCount"."\t"."$exonStarts"."\t"."$exonEnds"."\n";

	}
 

	$gene_index++;

     }

  }


}

close(OUTTRANSCRIPTABLE)
